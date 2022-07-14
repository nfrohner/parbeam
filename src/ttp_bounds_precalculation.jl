#!/usr/bin/env julia

@everywhere begin
    using Distributed
    using SharedArrays
    include("TTP/ttp_cvrp.jl")
end

using DelimitedFiles
using Random
using Printf: @printf
using JLSO

function main()
    if(length(ARGS) != 4)
        println("Usage: $PROGRAM_FILE <instance-file> <streak-limit> <bounds-file> <cvrp-h-bounds>")
        exit(1)
    end

    instance_file = ARGS[1]
    streak_limit = parse(Int, ARGS[2])
    bounds_file = ARGS[3]
    d = readdlm(instance_file, Int)
    n = size(d)[1]
    cvrp_h_bounds = parse(Bool, ARGS[4])

    if cvrp_h_bounds
        bounds_by_state = SharedArray(zeros(UInt32, (n, 2^(n-1), n, streak_limit, n)))
        @sync @distributed for team in 1:n
            TTPCVRP.calculate_bounds_for_team(n, team, d, streak_limit, bounds_by_state)
        end
    else
        bounds_by_state = SharedArray(zeros(UInt32, (n, 2^(n-1), n, streak_limit)))
        @time @sync @distributed for team in 1:n
            TTPCVRP.calculate_bounds_for_team(n, team, d, streak_limit, bounds_by_state)
        end
    end

    @printf("bounds sum %d\n", sum(bounds_by_state))
    @time JLSO.save(bounds_file, :bounds => Array(bounds_by_state))
end

@time main()
