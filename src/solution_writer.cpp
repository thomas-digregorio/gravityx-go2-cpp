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

}  // namespace

std::string go_solution_text(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency) {
    require_size(state.vm.size(), data.buses.size(), "vm");
    require_size(state.va.size(), data.buses.size(), "va");
    require_size(state.demand_factor.size(), data.loads.size(),
                 "demand_factor");
    require_size(state.pg.size(), data.generators.size(), "pg");
    require_size(state.qg.size(), data.generators.size(), "qg");
    require_size(commitment.size(), data.generators.size(), "commitment");

    const bool generator_outage = contingency != nullptr &&
        contingency->type == ContingencyType::Generator;
    const bool branch_outage = contingency != nullptr &&
        contingency->type == ContingencyType::Branch;
    const int outage_index = contingency != nullptr
        ? contingency->source_index : -1;

    std::string output;
    const std::size_t estimated_rows =
        data.buses.size() + data.loads.size() + data.generators.size() +
        data.branches.size() + data.shunts.size();
    output.reserve(512 + 48 * estimated_rows);
    output.append("--bus section\n"
                  "i, v, theta\n");
    for (std::size_t position = 0; position < data.buses.size(); ++position) {
        const auto& bus = data.buses[position];
        if (bus.present) {
            append_integer(output, bus.index);
            append_separator(output);
            append_double(output, state.vm[position]);
            append_separator(output);
            append_double(output, state.va[position]);
            output.push_back('\n');
        }
    }

    output.append("--load section\n"
                  "i, id, t\n");
    for (std::size_t position = 0; position < data.loads.size(); ++position) {
        const auto& load = data.loads[position];
        if (load.present) {
            append_integer(output, load.source_bus);
            append_separator(output);
            output.append(load.source_id);
            append_separator(output);
            append_double(output, state.demand_factor[position]);
            output.push_back('\n');
        }
    }

    output.append("--generator section\n"
                  "i, id, p, q, x\n");
    for (std::size_t position = 0;
         position < data.generators.size(); ++position) {
        const auto& generator = data.generators[position];
        if (!generator.present ||
            (generator_outage && generator.index == outage_index)) {
            continue;
        }
        append_integer(output, generator.source_bus);
        append_separator(output);
        output.append(generator.source_id);
        append_separator(output);
        append_double(output, state.pg[position]);
        append_separator(output);
        append_double(output, state.qg[position]);
        append_separator(output);
        append_integer(output, commitment[position]);
        output.push_back('\n');
    }

    output.append("--line section\n"
                  "iorig, idest, id, x\n");
    for (const auto& branch : data.branches) {
        if (!branch.present || branch.transformer ||
            (branch_outage && branch.index == outage_index)) {
            continue;
        }
        append_integer(output, branch.source_from);
        append_separator(output);
        append_integer(output, branch.source_to);
        append_separator(output);
        output.append(branch.source_id);
        append_separator(output);
        append_integer(output, branch.status != 0 ? 1 : 0);
        output.push_back('\n');
    }

    output.append("--transformer section\n"
                  "iorig, idest, id, x, xst\n");
    for (const auto& branch : data.branches) {
        if (!branch.present || !branch.transformer ||
            (branch_outage && branch.index == outage_index)) {
            continue;
        }
        int step = 0;
        if (branch.control_mode == 1 || branch.control_mode == -1) {
            step = branch.tm_step;
        } else if (branch.control_mode == 3 ||
                   branch.control_mode == -3) {
            step = branch.ta_step;
        }
        append_integer(output, branch.source_from);
        append_separator(output);
        append_integer(output, branch.source_to);
        append_separator(output);
        output.append(branch.source_id);
        append_separator(output);
        append_integer(output, branch.status != 0 ? 1 : 0);
        append_separator(output);
        append_integer(output, step);
        output.push_back('\n');
    }

    output.append(
        "--switched shunt section\n"
        "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8\n");
    for (std::size_t position = 0; position < data.shunts.size(); ++position) {
        const auto& shunt = data.shunts[position];
        if (!shunt.present || !shunt.dispatchable) {
            continue;
        }
        const std::vector<int>* steps = &shunt.steps;
        if (position < state.shunt_steps.size()) {
            steps = &state.shunt_steps[position];
        }
        append_integer(output, shunt.source_bus);
        for (const int step : *steps) {
            append_separator(output);
            append_integer(output, step);
        }
        output.push_back('\n');
    }
    return output;
}

void write_go_solution(
    const std::string& path,
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency) {
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
        const auto text = go_solution_text(
            data, state, commitment, contingency);
        stream.write(text.data(), static_cast<std::streamsize>(text.size()));
        if (!stream) {
            throw std::runtime_error(
                "failed while writing submission output: " + temporary);
        }
    }
    std::filesystem::rename(temporary, output);
}

}  // namespace gravityx
