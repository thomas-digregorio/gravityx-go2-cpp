#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gravityx {

class GoSolutionWriter {
public:
    explicit GoSolutionWriter(const CaseData& data);

    std::string text(
        const AcState& state,
        const std::vector<int>& commitment,
        const Contingency* contingency = nullptr) const;

    void write(
        const std::string& path,
        const AcState& state,
        const std::vector<int>& commitment,
        const Contingency* contingency = nullptr) const;

private:
    using RowSpan = std::pair<std::size_t, std::size_t>;

    const CaseData& data_;
    std::vector<std::string> bus_prefixes_;
    std::vector<std::string> load_prefixes_;
    std::vector<std::string> generator_prefixes_;
    std::vector<std::string> shunt_prefixes_;
    std::string line_section_;
    std::string transformer_section_;
    std::unordered_map<int, RowSpan> line_outage_rows_;
    std::unordered_map<int, RowSpan> transformer_outage_rows_;
};

std::string go_solution_text(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency = nullptr);

void write_go_solution(
    const std::string& path,
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& commitment,
    const Contingency* contingency = nullptr);

}  // namespace gravityx
