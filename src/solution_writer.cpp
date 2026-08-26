#include "gravityx/solution_writer.hpp"

#include <charconv>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <system_error>

namespace gravityx {
namespace {

void require_size(
    std::size_t actual,
    std::size_t expected,
    const char* name) {
    if (actual != expected) {
        throw std::runtime_error(
            std::string("submission state ") + name + " size " +
            std::to_string(actual) + " does not match case size " +
            std::to_string(expected));
    }
}

void reject_onedrive(const std::string& path) {
    std::string normalized = path;
    for (char& character : normalized) {
        character = static_cast<char>(std::tolower(
            static_cast<unsigned char>(character)));
    }
    if (normalized.find("onedrive") != std::string::npos) {
        throw std::runtime_error("refusing OneDrive path: " + path);
    }
}

template <typename Integer>
void append_integer(std::string& output, Integer value) {
    char buffer[32];
    const auto conversion = std::to_chars(
        buffer, buffer + sizeof(buffer), value);
    if (conversion.ec != std::errc{}) {
        throw std::runtime_error("failed to format submission integer");
    }
    output.append(buffer, conversion.ptr);
}

void append_double(std::string& output, double value) {
    char buffer[64];
    const auto conversion = std::to_chars(
        buffer, buffer + sizeof(buffer), value,
        std::chars_format::general,
        std::numeric_limits<double>::max_digits10);
    if (conversion.ec != std::errc{}) {
        throw std::runtime_error("failed to format submission floating-point value");
    }
    output.append(buffer, conversion.ptr);
}

void append_separator(std::string& output) {
    output.append(", ");
}

void append_static_section(
    std::string& output,
    const std::string& section,
    const std::unordered_map<int, std::pair<std::size_t, std::size_t>>&
        outage_rows,
    int outage_index) {
    const auto omission = outage_rows.find(outage_index);
    if (omission == outage_rows.end()) {
        output.append(section);
        return;
    }
    const auto [offset, length] = omission->second;
    output.append(section.data(), offset);
    const auto suffix_offset = offset + length;
    output.append(
        section.data() + suffix_offset,
        section.size() - suffix_offset);
}

}  // namespace

GoSolutionWriter::GoSolutionWriter(const CaseData& data)
    : data_(data),
      bus_prefixes_(data.buses.size()),
      load_prefixes_(data.loads.size()),
      generator_prefixes_(data.generators.size()),
      shunt_prefixes_(data.shunts.size()) {
    for (std::size_t position = 0; position < data_.buses.size(); ++position) {
        const auto& bus = data_.buses[position];
        if (!bus.present) {
            continue;
        }
        append_integer(bus_prefixes_[position], bus.index);
        append_separator(bus_prefixes_[position]);
    }
    for (std::size_t position = 0; position < data_.loads.size(); ++position) {
        const auto& load = data_.loads[position];
        if (!load.present) {
            continue;
        }
        append_integer(load_prefixes_[position], load.source_bus);
        append_separator(load_prefixes_[position]);
        load_prefixes_[position].append(load.source_id);
        append_separator(load_prefixes_[position]);
    }
    for (std::size_t position = 0;
         position < data_.generators.size(); ++position) {
        const auto& generator = data_.generators[position];
        if (!generator.present) {
            continue;
        }
        append_integer(
            generator_prefixes_[position], generator.source_bus);
        append_separator(generator_prefixes_[position]);
        generator_prefixes_[position].append(generator.source_id);
        append_separator(generator_prefixes_[position]);
    }
    for (std::size_t position = 0; position < data_.shunts.size(); ++position) {
        const auto& shunt = data_.shunts[position];
        if (!shunt.present || !shunt.dispatchable) {
            continue;
        }
        append_integer(shunt_prefixes_[position], shunt.source_bus);
    }

    line_section_.append(
        "--line section\n"
        "iorig, idest, id, x\n");
    transformer_section_.append(
        "--transformer section\n"
        "iorig, idest, id, x, xst\n");
    for (const auto& branch : data_.branches) {
        if (!branch.present) {
            continue;
        }
        auto& section = branch.transformer
            ? transformer_section_ : line_section_;
        const auto offset = section.size();
        append_integer(section, branch.source_from);
        append_separator(section);
        append_integer(section, branch.source_to);
        append_separator(section);
        section.append(branch.source_id);
        append_separator(section);
        append_integer(section, branch.status != 0 ? 1 : 0);
        if (branch.transformer) {
            int step = 0;
            if (branch.control_mode == 1 || branch.control_mode == -1) {
                step = branch.tm_step;
            } else if (branch.control_mode == 3 ||
                       branch.control_mode == -3) {
                step = branch.ta_step;
            }
            append_separator(section);
            append_integer(section, step);
        }
        section.push_back('\n');
        auto& outage_rows = branch.transformer
            ? transformer_outage_rows_ : line_outage_rows_;
        const auto inserted = outage_rows.emplace(
            branch.index, RowSpan{offset, section.size() - offset});
        if (!inserted.second) {
            throw std::runtime_error(
                "duplicate branch index in solution writer: " +
                std::to_string(branch.index));
        }
    }
}

std::string GoSolutionWriter::text(
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency) const {
    require_size(state.vm.size(), data_.buses.size(), "vm");
    require_size(state.va.size(), data_.buses.size(), "va");
    require_size(state.demand_factor.size(), data_.loads.size(),
                 "demand_factor");
    require_size(state.pg.size(), data_.generators.size(), "pg");
    require_size(state.qg.size(), data_.generators.size(), "qg");
    require_size(commitment.size(), data_.generators.size(), "commitment");

    const bool generator_outage = contingency != nullptr &&
        contingency->type == ContingencyType::Generator;
    const bool branch_outage = contingency != nullptr &&
        contingency->type == ContingencyType::Branch;
    const int outage_index = contingency != nullptr
        ? contingency->source_index : -1;

    std::string output;
    const std::size_t estimated_rows =
        data_.buses.size() + data_.loads.size() + data_.generators.size() +
        data_.shunts.size();
    output.reserve(
        line_section_.size() + transformer_section_.size() +
        512 + 48 * estimated_rows);
    output.append("--bus section\n"
                  "i, v, theta\n");
    for (std::size_t position = 0; position < data_.buses.size(); ++position) {
        const auto& bus = data_.buses[position];
        if (bus.present) {
            output.append(bus_prefixes_[position]);
            append_double(output, state.vm[position]);
            append_separator(output);
            append_double(output, state.va[position]);
            output.push_back('\n');
        }
    }

    output.append("--load section\n"
                  "i, id, t\n");
    for (std::size_t position = 0; position < data_.loads.size(); ++position) {
        const auto& load = data_.loads[position];
        if (load.present) {
            output.append(load_prefixes_[position]);
            append_double(output, state.demand_factor[position]);
            output.push_back('\n');
        }
    }

    output.append("--generator section\n"
                  "i, id, p, q, x\n");
    for (std::size_t position = 0;
         position < data_.generators.size(); ++position) {
        const auto& generator = data_.generators[position];
        if (!generator.present ||
            (generator_outage && generator.index == outage_index)) {
            continue;
        }
        output.append(generator_prefixes_[position]);
        append_double(output, state.pg[position]);
        append_separator(output);
        append_double(output, state.qg[position]);
        append_separator(output);
        append_integer(output, commitment[position]);
        output.push_back('\n');
    }

    append_static_section(
        output, line_section_, line_outage_rows_,
        branch_outage ? outage_index : -1);
    append_static_section(
        output, transformer_section_, transformer_outage_rows_,
        branch_outage ? outage_index : -1);

    output.append(
        "--switched shunt section\n"
        "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8\n");
    for (std::size_t position = 0; position < data_.shunts.size(); ++position) {
        const auto& shunt = data_.shunts[position];
        if (!shunt.present || !shunt.dispatchable) {
            continue;
        }
        const std::vector<int>* steps = &shunt.steps;
        if (position < state.shunt_steps.size()) {
            steps = &state.shunt_steps[position];
        }
        output.append(shunt_prefixes_[position]);
        for (const int step : *steps) {
            append_separator(output);
            append_integer(output, step);
        }
        output.push_back('\n');
    }
    return output;
}

void GoSolutionWriter::write(
    const std::string& path,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency) const {
    reject_onedrive(path);
    const std::filesystem::path output(path);
    if (output.has_parent_path()) {
        std::filesystem::create_directories(output.parent_path());
    }
    const auto temporary = output.string() + ".tmp";
    {
        std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
        if (!stream) {
            throw std::runtime_error(
                "cannot open submission output: " + temporary);
        }
        const auto solution_text = text(state, commitment, contingency);
        stream.write(
            solution_text.data(),
            static_cast<std::streamsize>(solution_text.size()));
        if (!stream) {
            throw std::runtime_error(
                "failed while writing submission output: " + temporary);
        }
    }
    std::filesystem::rename(temporary, output);
}

std::string go_solution_text(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency) {
    return GoSolutionWriter(data).text(state, commitment, contingency);
}

void write_go_solution(
    const std::string& path,
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency) {
    GoSolutionWriter(data).write(path, state, commitment, contingency);
}

}  // namespace gravityx
