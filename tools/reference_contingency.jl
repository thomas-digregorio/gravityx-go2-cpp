#!/usr/bin/env julia

# Development-only semantic oracle.  This script solves one corrective case
# with the pinned official Julia formulation so the Gravity C++ formulation can
# be compared component-by-component.  Production optimization remains C++.

using Ipopt
using JSON
using JuMP
using PowerModels
using PowerModelsSecurityConstrained

const PM = PowerModels
const PMSC = PowerModelsSecurityConstrained

length(ARGS) == 4 || error(
    "usage: reference_contingency.jl <case-dir> <cpp-base.json> <label> <output.json>")

case_dir = abspath(ARGS[1])
base_path = abspath(ARGS[2])
label = ARGS[3]
output_path = abspath(ARGS[4])
for path in (case_dir, base_path, output_path)
    occursin("onedrive", lowercase(path)) && error("refusing OneDrive path: $(path)")
end

function ordered_ids(items::Dict)
    return sort(collect(keys(items)), by=id -> Int(items[id]["index"]))
end

raw_path = joinpath(case_dir, "case.raw")
con_path = joinpath(case_dir, "case.con")
json_path = joinpath(case_dir, "case.json")
goc_data = PMSC.parse_c2_files(
    raw_path, con_path, json_path; case_id="semantic_tiny", scenario_id="01")
data = PMSC.build_c2_pm_model(goc_data)
PMSC.correct_c2_solution!(data)
PMSC.normalize_tm_values!(data)

base = JSON.parsefile(base_path)
base["success"] == true || error("C++ base result was not successful")
state = base["selected_state"]
commitment = Int.(base["commitment"])

bus_ids = ordered_ids(data["bus"])
gen_ids = ordered_ids(data["gen"])
load_ids = ordered_ids(data["load"])
branch_ids = ordered_ids(data["branch"])

for (position, id) in enumerate(bus_ids)
    bus = data["bus"][id]
    bus["vm"] = state["vm"][position]
    bus["va"] = state["va"][position]
    bus["vm_start"] = state["vm"][position]
    bus["va_start"] = state["va"][position]
end
for (position, id) in enumerate(gen_ids)
    gen = data["gen"][id]
    gen["gen_status"] = commitment[position]
    gen["status_prev"] = commitment[position]
    gen["pg"] = state["pg"][position]
    gen["qg"] = state["qg"][position]
    gen["pg_prev"] = state["pg"][position]
    gen["qg_prev"] = state["qg"][position]
    gen["pg_start"] = state["pg"][position]
    gen["qg_start"] = state["qg"][position]
end
for (position, id) in enumerate(load_ids)
    load = data["load"][id]
    factor = state["demand_factor"][position]
    load["pd"] = load["pd_nominal"] * factor
    load["qd"] = load["qd_nominal"] * factor
    load["pd_prev"] = load["pd"]
    load["qd_prev"] = load["qd"]
    load["z_demand_start"] = factor
end
for (position, id) in enumerate(branch_ids)
    branch = data["branch"][id]
    for key in ("pf", "qf", "pt", "qt")
        branch["$(key)_start"] = state[key][position]
    end
end

contingencies = vcat(data["gen_contingencies"], data["branch_contingencies"])
matches = filter(item -> item.label == label, contingencies)
length(matches) == 1 || error("expected one contingency named $(label)")
contingency = only(matches)
component_group = contingency.type == "gen" ? "gen" : "branch"
component_id = string(contingency.idx)
data[component_group][component_id]["present"] = false
if contingency.type == "gen"
    data[component_group][component_id]["gen_status"] = 0
else
    data[component_group][component_id]["br_status"] = 0
end

optimizer = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "tol" => 1e-6,
    "print_level" => 0,
    "linear_solver" => "mumps",
)
result = PMSC.run_c2_opf_soft_ctg(data, PM.ACPPowerModel, optimizer)

solution = result["solution"]
state_out = Dict(
    "vm" => [solution["bus"][id]["vm"] for id in bus_ids],
    "va" => [solution["bus"][id]["va"] for id in bus_ids],
    "pg" => [haskey(solution["gen"], id) ? solution["gen"][id]["pg"] : 0.0 for id in gen_ids],
    "qg" => [haskey(solution["gen"], id) ? solution["gen"][id]["qg"] : 0.0 for id in gen_ids],
)
output = Dict(
    "label" => label,
    "termination_status" => string(result["termination_status"]),
    "objective" => result["objective"],
    "state" => state_out,
)
mkpath(dirname(output_path))
open(output_path, "w") do io
    JSON.print(io, output, 2)
    println(io)
end
println(JSON.json(Dict("output" => output_path, "objective" => result["objective"])))
