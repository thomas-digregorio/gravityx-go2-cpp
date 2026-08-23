#!/usr/bin/env julia

using JSON
using PowerModelsSecurityConstrained

length(ARGS) == 4 || error("usage: export_case.jl <case-dir> <output.json> <case-id> <scenario-id>")

case_dir = abspath(ARGS[1])
output_path = abspath(ARGS[2])
case_id = ARGS[3]
scenario_id = ARGS[4]

for path in (case_dir, output_path)
    occursin("onedrive", lowercase(path)) && error("refusing OneDrive path: $(path)")
end

raw_path = joinpath(case_dir, "case.raw")
con_path = joinpath(case_dir, "case.con")
json_path = joinpath(case_dir, "case.json")
all(isfile, (raw_path, con_path, json_path)) || error("case directory must contain case.raw, case.con, and case.json")

goc_data = parse_c2_files(raw_path, con_path, json_path; case_id=case_id, scenario_id=scenario_id)
pm_data = build_c2_pm_model(goc_data)
correct_c2_solution!(pm_data)
normalize_tm_values!(pm_data)
set_va_start_values!(pm_data)

mkpath(dirname(output_path))
open(output_path, "w") do io
    JSON.print(io, pm_data)
    println(io)
end

println(JSON.json(Dict(
    "case_id" => case_id,
    "scenario_id" => scenario_id,
    "buses" => length(pm_data["bus"]),
    "generators" => length(pm_data["gen"]),
    "loads" => length(pm_data["load"]),
    "branches" => length(pm_data["branch"]),
    "output" => output_path,
)))
