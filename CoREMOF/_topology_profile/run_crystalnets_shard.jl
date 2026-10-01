#!/usr/bin/env julia

using CrystalNets
using JSON3
using SHA

const SCHEMA_VERSION = "crystalnets-topology-result/1.0"

CrystalNets.toggle_export(false)
CrystalNets.toggle_warning(false)

const OPTIONS = CrystalNets.Options(
    structure=CrystalNets.StructureType.MOF,
    clusterings=[
        CrystalNets.Clustering.SingleNodes,
        CrystalNets.Clustering.AllNodes,
    ],
)

function clean_error(error, cif_path)
    message = sprint(showerror, error)
    return replace(message, cif_path => "<CIF>")
end

function topology_payload(topology)
    if topology === missing
        return Dict("status" => "NOT_AVAILABLE")
    end
    label = string(topology)
    return Dict(
        "status" => "SUCCESS",
        "dimension" => ndims(topology.genome),
        "topology" => label,
        "is_named_net" => !startswith(label, "UNKNOWN"),
    )
end

function consensus(values)
    isempty(values) && return nothing
    first_value = first(values)
    return all(value == first_value for value in values) ? first_value : nothing
end

function process_one(structure_id, cif_path, cif_sha256)
    started_ns = time_ns()
    common = Dict(
        "schema_version" => SCHEMA_VERSION,
        "structure_id" => structure_id,
        "cif_sha256" => cif_sha256,
        "software" => Dict(
            "julia_version" => string(VERSION),
            "crystalnets_version" => string(Base.pkgversion(CrystalNets)),
        ),
        "method" => Dict(
            "structure_type" => "MOF",
            "clusterings" => ["SingleNodes", "AllNodes"],
            "exports_enabled" => false,
            "warnings_enabled" => false,
            "interpenetration_definition" => (
                "number of subnets in InterpenetratedTopologyResult"
            ),
        ),
    )
    try
        result = CrystalNets.determine_topology(cif_path, OPTIONS)
        subnets = []
        for (index, subnet) in enumerate(result)
            single = topology_payload(
                subnet[CrystalNets.Clustering.SingleNodes]
            )
            allnode = topology_payload(
                subnet[CrystalNets.Clustering.AllNodes]
            )
            push!(
                subnets,
                Dict(
                    "subnet_index" => index,
                    "single_node" => single,
                    "all_node" => allnode,
                    "single_all_agree" => (
                        single == allnode
                        && get(single, "status", "") == "SUCCESS"
                    ),
                ),
            )
        end
        single_dimensions = [
            subnet["single_node"]["dimension"]
            for subnet in subnets
            if get(subnet["single_node"], "status", "") == "SUCCESS"
        ]
        all_dimensions = [
            subnet["all_node"]["dimension"]
            for subnet in subnets
            if get(subnet["all_node"], "status", "") == "SUCCESS"
        ]
        single_topologies = [
            subnet["single_node"]["topology"]
            for subnet in subnets
            if get(subnet["single_node"], "status", "") == "SUCCESS"
        ]
        all_topologies = [
            subnet["all_node"]["topology"]
            for subnet in subnets
            if get(subnet["all_node"], "status", "") == "SUCCESS"
        ]
        single_dimension = consensus(single_dimensions)
        all_dimension = consensus(all_dimensions)
        network_dimension = (
            single_dimension !== nothing
            && single_dimension == all_dimension
        ) ? single_dimension : nothing
        merge!(
            common,
            Dict(
                "execution_status" => "SUCCESS",
                "runtime_seconds" => (time_ns() - started_ns) / 1.0e9,
                "interpenetrated_subnet_count" => length(result),
                "catenation_degree" => length(result),
                "network_dimension" => network_dimension,
                "single_node_net" => consensus(single_topologies),
                "all_node_net" => consensus(all_topologies),
                "all_subnets_single_all_agree" => all(
                    get(subnet, "single_all_agree", false)
                    for subnet in subnets
                ),
                "subnets" => subnets,
            ),
        )
    catch error
        merge!(
            common,
            Dict(
                "execution_status" => "ERROR",
                "runtime_seconds" => (time_ns() - started_ns) / 1.0e9,
                "error_type" => string(typeof(error)),
                "error_message" => clean_error(error, cif_path),
                "interpenetrated_subnet_count" => nothing,
                "catenation_degree" => nothing,
                "network_dimension" => nothing,
                "single_node_net" => nothing,
                "all_node_net" => nothing,
                "all_subnets_single_all_agree" => nothing,
                "subnets" => [],
            ),
        )
    end
    return common
end

function main()
    length(ARGS) == 4 || error(
        "usage: run_crystalnets_shard.jl MANIFEST SHARD_INDEX SHARD_COUNT OUTPUT"
    )
    manifest_path = ARGS[1]
    shard_index = parse(Int, ARGS[2])
    shard_count = parse(Int, ARGS[3])
    output_path = ARGS[4]
    0 <= shard_index < shard_count || error("invalid shard index/count")
    lines = readlines(manifest_path)
    isempty(lines) && error("empty manifest")
    header = split(first(lines), '\t')
    header == ["structure_id", "cif_path", "cif_sha256"] || error(
        "unexpected manifest header"
    )
    temporary = output_path * ".tmp." * string(getpid())
    processed = 0
    open(temporary, "w") do output
        for (row_index, line) in enumerate(Iterators.drop(lines, 1))
            (row_index - 1) % shard_count == shard_index || continue
            fields = split(line, '\t')
            length(fields) == 3 || error("malformed manifest row $(row_index)")
            record = process_one(fields[1], fields[2], fields[3])
            JSON3.write(output, record)
            write(output, '\n')
            flush(output)
            processed += 1
            processed % 10 == 0 && GC.gc(false)
        end
    end
    mv(temporary, output_path; force=true)
    println("shard=$(shard_index) processed=$(processed) output=$(output_path)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
