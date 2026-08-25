#include "gravityx/solution_writer.hpp"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

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

    std::ostringstream stream;
    stream << std::setprecision(17) << std::defaultfloat;
    stream << "--bus section\n"
           << "i, v, theta\n";
    for (std::size_t position = 0; position < data.buses.size(); ++position) {
        const auto& bus = data.buses[position];
        if (bus.present) {
            stream << bus.index << ", " << state.vm[position] << ", "
                   << state.va[position] << '\n';
        }
    }

    stream << "--load section\n"
           << "i, id, t\n";
    for (std::size_t position = 0; position < data.loads.size(); ++position) {
        const auto& load = data.loads[position];
        if (load.present) {
            stream << load.source_bus << ", " << load.source_id
                   << ", " << state.demand_factor[position] << '\n';
        }
    }

    stream << "--generator section\n"
           << "i, id, p, q, x\n";
    for (std::size_t position = 0;
         position < data.generators.size(); ++position) {
        const auto& generator = data.generators[position];
        if (!generator.present ||
            (generator_outage && generator.index == outage_index)) {
            continue;
        }
        stream << generator.source_bus << ", "
               << generator.source_id << ", " << state.pg[position] << ", "
               << state.qg[position] << ", " << commitment[position] << '\n';
    }

    stream << "--line section\n"
           << "iorig, idest, id, x\n";
    for (const auto& branch : data.branches) {
        if (!branch.present || branch.transformer ||
            (branch_outage && branch.index == outage_index)) {
            continue;
        }
        stream << branch.source_from << ", "
               << branch.source_to << ", " << branch.source_id
               << ", " << (branch.status != 0 ? 1 : 0) << '\n';
    }

    stream << "--transformer section\n"
           << "iorig, idest, id, x, xst\n";
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
        stream << branch.source_from << ", "
               << branch.source_to << ", " << branch.source_id
               << ", " << (branch.status != 0 ? 1 : 0) << ", " << step
               << '\n';
    }

    stream << "--switched shunt section\n"
           << "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8\n";
    for (std::size_t position = 0; position < data.shunts.size(); ++position) {
        const auto& shunt = data.shunts[position];
        if (!shunt.present || !shunt.dispatchable) {
            continue;
        }
        const std::vector<int>* steps = &shunt.steps;
        if (position < state.shunt_steps.size()) {
            steps = &state.shunt_steps[position];
        }
        stream << shunt.source_bus;
        for (const int step : *steps) {
            stream << ", " << step;
        }
        stream << '\n';
    }
    return stream.str();
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
