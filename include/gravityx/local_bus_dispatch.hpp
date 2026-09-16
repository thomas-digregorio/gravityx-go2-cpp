#pragma once

#include "gravityx/ac_model.hpp"

namespace gravityx {

struct LocalBusDispatchResult {
    int active_generation_changes{};
    int reactive_generation_changes{};
    int load_changes{};
    double predicted_native_gain{};
    bool changed() const {
        return active_generation_changes + reactive_generation_changes + load_changes > 0;
    }
};

// Bound to one immutable case, original within-run base, and commitment.
// Only source-derived intervals/curves are cached, never candidate solutions.
class LocalBusDispatchCache {
public:
    LocalBusDispatchCache(const CaseData& data, const AcState& original_base,
                          const std::vector<int>& commitment);
    LocalBusDispatchResult improve(
        const Contingency& contingency, const std::vector<double>& network_p,
        const std::vector<double>& network_q, AcState& candidate, int passes = 2) const;
private:
    const CaseData& data_;
    const AcState& original_base_;
    std::vector<int> commitment_;
    std::vector<double> lower_g_, upper_g_, lower_t_, upper_t_;
    std::vector<std::vector<PwlPoint>> gen_curves_, load_curves_;
};

// Candidate generation only. Keep all network controls/flows fixed, and use
// ORIGINAL-base corrective bounds. The caller must rebuild and independently
// validate the complete state and retain its incumbent on any failed check.
LocalBusDispatchResult improve_fixed_network_bus_dispatch(
    const CaseData& data, const AcState& original_base,
    const std::vector<int>& commitment, const Contingency& contingency,
    const std::vector<double>& network_p, const std::vector<double>& network_q,
    AcState& candidate, int passes = 2);

void run_local_bus_dispatch_regression();

}  // namespace gravityx
