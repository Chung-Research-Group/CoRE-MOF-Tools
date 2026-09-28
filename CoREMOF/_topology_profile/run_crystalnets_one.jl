#!/usr/bin/env julia

include(joinpath(@__DIR__, "run_crystalnets_shard.jl"))

function main_one()
    length(ARGS) == 4 || error(
        "usage: run_crystalnets_one.jl STRUCTURE_ID CIF_PATH CIF_SHA256 OUTPUT"
    )
    output_path = ARGS[4]
    temporary = output_path * ".tmp." * string(getpid())
    record = process_one(ARGS[1], ARGS[2], ARGS[3])
    open(temporary, "w") do output
        JSON3.write(output, record)
        write(output, '\n')
    end
    mv(temporary, output_path; force=true)
end

main_one()
