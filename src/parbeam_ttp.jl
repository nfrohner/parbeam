include("lib/parbeam.jl")
include("TTP/ttp_cvrp.jl")
include("TTP/ttp_guidance.jl")

using DelimitedFiles
using Random
using Printf: @printf
using JLD2
using CodecZlib
using JLSO

# specialized types
mutable struct AdditionalData <: AbstractAdditionalData
    games::Array{UInt16, 2}
    forbidden_opponents::Array{UInt8, 2}
    rounds::Array{UInt8, 2}
    positions::Array{UInt8, 2}
    possible_away_streaks::Array{UInt8, 2}
    possible_home_stands::Array{UInt8, 2}
    heuristic_estimates::Array{Int, 2}
    number_of_away_games_left::Array{UInt8, 2}
    number_of_home_games_left::Array{UInt8, 2}
    away_teams_left_masks::Array{UInt32, 2}
end

mutable struct AdditionalSuccessorData <: AbstractAdditionalSuccessorData
    heuristic_estimates_delta_away_team::Array{Int, 1}
    heuristic_estimates_delta_home_team::Array{Int, 1}
end

TTPNodeData = NodeData{UInt32, Float64, Int16}
TTPSuccessorData = SuccessorData{UInt32, Float64, UInt16, AdditionalSuccessorData, Int16}
TTPBeam = Beam{TTPNodeData, AdditionalData}

mutable struct TTPStatistics
    eval_successors_prep_time::Float64
    eval_successors_loop_time::Float64
    TTPStatistics() = new(0.0, 0.0)
end

struct AuxiliaryData
    number_of_away_games_left::Array{Int, 1}
    number_of_home_games_left::Array{Int, 1}
    away_teams_left_masks::Array{Int, 1}
    number_of_permitted_games_left::Array{Int, 1}
    statistics::TTPStatistics
    bounds_by_state::Union{Array{UInt32, 4}, Array{UInt32, 5}, Nothing}
end

struct Instance <: AbstractInstance
    n::Int
    d::Array{Int, 2}
    streak_limit::Int
    no_repeat::Bool
    games::Vector{Tuple{Int,Int}}
    aux_by_thread::Vector{AuxiliaryData}
    teams_permutation::Array{Int, 1}
end

function maximum_layer_number(instance::Instance)
    instance.n*(instance.n-1)
end

function maximum_number_of_successors_per_node(instance::Instance)
    2*(instance.n-1)
end

function number_of_games(instance::Instance)
    instance.n*(instance.n-1)
end

function initialize_beam(configuration::Configuration, instance::Instance)
    layer = Array{Int16, 1}(undef, configuration.beam_width)
    costs_so_far = Array{UInt32, 1}(undef, configuration.beam_width)
    priority = Array{Float64, 1}(undef, configuration.beam_width)
    nodes = TTPNodeData(layer, costs_so_far, priority)

    games = Array{UInt16, 2}(undef, number_of_games(instance), configuration.beam_width)
    forbidden_opponents = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    rounds = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    positions = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    possible_away_streaks = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    possible_home_stands = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    heuristic_estimates = Array{Int, 2}(undef, instance.n, configuration.beam_width)
    number_of_away_games_left = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    number_of_home_games_left = Array{UInt8, 2}(undef, instance.n, configuration.beam_width)
    away_teams_left_masks = Array{UInt32, 2}(undef, instance.n, configuration.beam_width)

    additional_data = AdditionalData(games, forbidden_opponents, rounds, positions, possible_away_streaks, possible_home_stands, heuristic_estimates, number_of_away_games_left, number_of_home_games_left, away_teams_left_masks)
    TTPBeam(nodes, additional_data)
end

function initialize_successor_specs(configuration::Configuration, instance::Instance)::TTPSuccessorData
    layer = Array{Int16, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    costs_after_move = Array{UInt32, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    priority = Array{Float64, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    Q_idx = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    successor_move = Array{UInt16, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    heuristic_estimates_delta_away_team = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    heuristic_estimates_delta_home_team = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)

    successor_data = TTPSuccessorData(layer, costs_after_move, priority, Q_idx, successor_move, AdditionalSuccessorData(heuristic_estimates_delta_away_team, heuristic_estimates_delta_home_team))
    successor_data
end

function initialize_root(instance::Instance, Q::TTPBeam)
    Q.nodes.layer[1] = 0
    Q.nodes.costs_so_far[1] = 0
    Q.nodes.priority[1] = 0

    for i in 1:number_of_games(instance)
        Q.additional_data.games[i,1] = i
    end

    #prepare_auxiliary_data(instance, Q, 1, 0)

    for i in 1:instance.n
        Q.additional_data.forbidden_opponents[i,1] = 0
        Q.additional_data.rounds[i,1] = 0
        Q.additional_data.positions[i,1] = i
        Q.additional_data.possible_away_streaks[i,1] = instance.streak_limit
        Q.additional_data.possible_home_stands[i,1] = instance.streak_limit
        Q.additional_data.heuristic_estimates[i,1] = TTPGuidance.heuristic_estimate(instance.streak_limit, instance.d, i, 2^(instance.n-1) - 1, instance.n-1, instance.n-1, i, 0, instance.aux_by_thread[Threads.threadid()].bounds_by_state)
        Q.additional_data.number_of_away_games_left[i,1] = instance.n-1
        Q.additional_data.number_of_home_games_left[i,1] = instance.n-1
        Q.additional_data.away_teams_left_masks[i,1] = 2^(instance.n - 1) -1
    end
end

function initial_estimate(configuration::Configuration, instance::Instance, beam::TTPBeam)
    if configuration.guidance_function == "cvrph" || configuration.guidance_function == "cvrp"
        heuristics_sum = 0
        for i in 1:instance.n
            heuristics_sum += TTPGuidance.heuristic_estimate(instance.streak_limit, instance.d, i, 2^(instance.n-1) - 1, instance.n-1, instance.n-1, i, 0, instance.aux_by_thread[Threads.threadid()].bounds_by_state)
        end
        heuristics_sum
    else
        # TODO: use solution of approximation algorithm here 
        0
    end
end

function initial_successor_moves(instance::Instance, depth::Int)
    moves_list = Vector{Vector{Int}}()
    push!(moves_list, [])
    if depth > 0
        throw(ArgumentError("subtree splitting not yet implemented for the TTP"))
    end
    moves_list
end

function create_objectives_and_solutions()
    Vector{Tuple{Float64,Array{Int}}}()
end

function create_all_objectives(size::Int)
    all_objectives = Array{Array{Float64}}(undef, size)
    for i in 1:size
        all_objectives[i] = zeros(Int, 1)
    end
    all_objectives
end

function create_all_solutions(instance::Instance, size::Int)
    all_solutions = Vector{Array{Int, 2}}(undef, size)
    for i in 1:size
        all_solutions[i] = zeros(Int, 2*(instance.n - 1), instance.n)
    end
    all_solutions
end

function get_solution(instance::Instance, beam::TTPBeam, best_terminal::Int)
    best_solution = beam.additional_data.games[:,best_terminal]

    schedule = zeros(Int, 2*(instance.n - 1), instance.n)
    round = 0
    half_number_of_teams = convert(Int, instance.n / 2)
    for (i, j) in enumerate(best_solution)
        if (i - 1) % half_number_of_teams == 0
            round += 1
        end
        game = instance.games[j]
        schedule[round, game[1]] = -game[2]
        schedule[round, game[2]] = game[1]
    end

    schedule
end

function read_in_TTP(guidance_function::String, variable_ordering::String, instance_name::String)
    file_name = string("insts/TTP/", instance_name, ".txt")
    d = readdlm(file_name, Int)
    n = size(d)[1]::Int
    streak_limit = 3
    no_repeat = true

    games = Vector{Tuple{Int,Int}}()
    sizehint!(games, n*(n-1))
    for i in 1:n
        for j in 1:n
            if i != j
                push!(games, (i, j))
            end
        end
    end

    if guidance_function == "cvrph"
        bounds_file = string("data/TTP/", instance_name, "_cvrph.jld2")
        if isfile(bounds_file)
            bounds_by_state = load(bounds_file)["bounds"]::Array{UInt32, 5}
        else
            bounds_file = string("data/TTP/", instance_name, "_cvrph.jlso")
            if isfile(bounds_file)
                bounds_by_state = JLSO.load(bounds_file)[:bounds]::Array{UInt32, 5}
            else
                throw(ArgumentError("no cvrph bounds file found"))
            end
        end
    elseif guidance_function == "cvrp"
        bounds_file = string("data/TTP/", instance_name, "_cvrp.jld2")
        if isfile(bounds_file)
            bounds_by_state = load(bounds_file)["bounds"]::Array{UInt32, 4}
        else
            bounds_file = string("data/TTP/", instance_name, "_cvrp.jlso")
            if isfile(bounds_file)
                bounds_by_state = convert(Array{UInt32, 4}, JLSO.load(bounds_file)[:bounds])::Array{UInt32, 4}
            else
                throw(ArgumentError("no cvrp bounds file found"))
            end
        end
    else
        bounds_by_state = nothing
    end

    aux_by_thread = Vector{AuxiliaryData}()
    for _ in 1:Threads.nthreads()
        push!(aux_by_thread, AuxiliaryData(zeros(Int, n), zeros(Int, n), zeros(Int, n), zeros(Int, n), TTPStatistics(), bounds_by_state))
    end

    if variable_ordering == "random"
        teams_permutation = convert(Array{Int}, randperm(n))
    else
        teams_permutation = convert(Array{Int}, 1:n)
    end

    Instance(n, d, streak_limit, no_repeat, games, aux_by_thread, teams_permutation)
end

function read_instance(guidance_function::String, variable_ordering::String, instance_description::String)
    read_in_TTP(guidance_function, variable_ordering, instance_description)
end

function instance_csv_string(configuration::Configuration, instance::Instance)
    @sprintf("%d", instance.n)
end

function next_team_dynamic(instance::Instance, beam::TTPBeam, Q_idx::Int, auxiliary_data::AuxiliaryData, min_round::Int)
    min_bound_by_team = ones(Int, instance.n)*typemax(Int)
    bound_sum_by_team = zeros(Int, instance.n)
    count_by_team = zeros(Int, instance.n)
    #min_bound_by_team = zeros(Int, instance.n)
    layer = beam.nodes.layer[Q_idx]
    heuristic_estimate = sum(@view beam.additional_data.heuristic_estimates[:,Q_idx])
    for game_idx in (layer+1):number_of_games(instance)
        j = beam.additional_data.games[game_idx, Q_idx]
        game = instance.games[j]
        heuristic_estimates_delta_away_team = 0
        heuristic_estimates_delta_home_team = 0
        if beam.additional_data.rounds[game[1],Q_idx] == min_round && game_allowed(instance, beam.additional_data, Q_idx, game[1], game[2]) && delta_feasibility_check(instance, beam.additional_data, Q_idx, game[1], game[2])
            heuristic_estimates_delta_away_team, heuristic_estimates_delta_home_team = calc_heuristic_estimates_delta(instance, beam.additional_data, Q_idx, game[1], game[2], auxiliary_data)
            priority = travel_distance_after_next_game(instance, beam, Q_idx, game[1], game[2]) + heuristic_estimate + heuristic_estimates_delta_away_team + heuristic_estimates_delta_home_team
            min_bound_by_team[game[1]] = min(priority, min_bound_by_team[game[1]])
            bound_sum_by_team[game[1]] += priority
            count_by_team[game[1]] += 1
            min_bound_by_team[game[2]] = min(priority, min_bound_by_team[game[2]])
            count_by_team[game[2]] += 1
            bound_sum_by_team[game[2]] += priority
        end
    end
    max_team = 1
    best_bound = 0
    for i in 1:instance.n
        if bound_sum_by_team[i] > 0
            bound_avg = bound_sum_by_team[i]/count_by_team[i]
            if bound_avg > best_bound
                best_bound = bound_avg
                max_team = i
            end
        end
    end
    #println(min_bound_by_team, max_team)
    max_team
end

function eval_successors(configuration::Configuration, instance::Instance, beam::TTPBeam, successor_specs::TTPSuccessorData, start_index::Int, Q_idx::Int, sigma::Float64)
    successor_count = 0
    layer = beam.nodes.layer[Q_idx]
    pivot_team, min_round, max_round = next_team(instance, beam.additional_data, Q_idx)

    #instance.aux_by_thread[Threads.threadid()].statistics.eval_successors_prep_time += @elapsed prepare_auxiliary_data(instance, beam, Q_idx, layer)
    auxiliary_data = instance.aux_by_thread[Threads.threadid()]
    
    if !dead_team_check_without_aux(instance, beam, Q_idx, layer, min_round, max_round)
        return successor_count
    end

    #pivot_team = next_team_dynamic(instance, beam, Q_idx, auxiliary_data, min_round)

    heuristic_estimate = sum(@view beam.additional_data.heuristic_estimates[:,Q_idx])

    cached_heuristic_estimates_delta_home_team = nothing
    for game_idx in (layer+1):number_of_games(instance)
        j = beam.additional_data.games[game_idx, Q_idx]
        game = instance.games[j]
        heuristic_estimates_delta_away_team = 0
        heuristic_estimates_delta_home_team = 0
        if (game[1] == pivot_team || game[2] == pivot_team) && game_allowed(instance, beam.additional_data, Q_idx, game[1], game[2]) && delta_feasibility_check(instance, beam.additional_data, Q_idx, game[1], game[2])
            i = start_index + successor_count
            @assert i <= configuration.maximum_number_of_successors_by_node*configuration.beam_width
            successor_specs.costs_after_move[i] = travel_distance_after_next_game(instance, beam, Q_idx, game[1], game[2])
            if configuration.guidance_function == "short"
                successor_specs.priority[i] = successor_specs.costs_after_move[i]
                heuristic_estimates_delta_away_team = 0
                heuristic_estimates_delta_home_team = 0
            elseif configuration.guidance_function == "cvrp" || configuration.guidance_function == "cvrph"
                if game[2] != pivot_team
                    heuristic_estimates_delta_away_team, heuristic_estimates_delta_home_team = calc_heuristic_estimates_delta(instance, beam.additional_data, Q_idx, game[1], game[2], auxiliary_data)
                elseif cached_heuristic_estimates_delta_home_team === nothing && game[2] == pivot_team
                    heuristic_estimates_delta_away_team, heuristic_estimates_delta_home_team = calc_heuristic_estimates_delta(instance, beam.additional_data, Q_idx, game[1], game[2], auxiliary_data)
                    cached_heuristic_estimates_delta_home_team = heuristic_estimates_delta_home_team
                elseif cached_heuristic_estimates_delta_home_team !== nothing && game[2] == pivot_team
                    heuristic_estimates_delta_away_team = calc_heuristic_estimates_delta_for_away_team(instance, beam.additional_data, Q_idx, game[1], game[2], auxiliary_data)
                    heuristic_estimates_delta_home_team = cached_heuristic_estimates_delta_home_team
                else
                    @assert false
                end
                successor_specs.priority[i] = successor_specs.costs_after_move[i] + heuristic_estimate + heuristic_estimates_delta_away_team + heuristic_estimates_delta_home_team
            else
                throw(ArgumentError("unknown guidance"))
            end
            @assert successor_specs.priority[i] >= 0
            successor_specs.layer[i] = layer + 1
            if sigma > 0.0
                successor_specs.priority[i] += randn()*sigma
            end
            successor_specs.Q_idx[i] = Q_idx
            successor_specs.successor_move[i] = j
            successor_specs.additional_data.heuristic_estimates_delta_away_team[i] = heuristic_estimates_delta_away_team
            successor_specs.additional_data.heuristic_estimates_delta_home_team[i] = heuristic_estimates_delta_home_team
            successor_count += 1
        end
    end
    successor_count
end

function calc_heuristic_estimates_delta(instance::Instance, additional_data::AdditionalData, Q_idx::Int, away_team::Int, home_team::Int, auxiliary_data::AuxiliaryData)
    heuristic_estimate_delta_away_team = -additional_data.heuristic_estimates[away_team,Q_idx]
    streak = instance.streak_limit - (additional_data.possible_away_streaks[away_team,Q_idx] - 1)
    if home_team < away_team
        teams_left_mask = additional_data.away_teams_left_masks[away_team,Q_idx] - 2 ^ (home_team-1)
    else
        teams_left_mask = additional_data.away_teams_left_masks[away_team,Q_idx] - 2 ^ (home_team-2)
    end

    heuristic_estimate_delta_away_team += TTPGuidance.heuristic_estimate(instance.streak_limit, instance.d, away_team, teams_left_mask*1, additional_data.number_of_away_games_left[away_team,Q_idx]-1, additional_data.number_of_home_games_left[away_team,Q_idx]+0, home_team, streak, auxiliary_data.bounds_by_state)
    
    heuristic_estimate_delta_home_team = -additional_data.heuristic_estimates[home_team,Q_idx]
    streak = instance.streak_limit - (additional_data.possible_home_stands[home_team,Q_idx] - 1)
    teams_left_mask = additional_data.away_teams_left_masks[home_team,Q_idx]
    heuristic_estimate_delta_home_team += TTPGuidance.heuristic_estimate(instance.streak_limit, instance.d, home_team, teams_left_mask*1, additional_data.number_of_away_games_left[home_team,Q_idx]+0, additional_data.number_of_home_games_left[home_team,Q_idx]-1, home_team, streak, auxiliary_data.bounds_by_state)

    heuristic_estimate_delta_away_team, heuristic_estimate_delta_home_team
end

function calc_heuristic_estimates_delta_for_away_team(instance::Instance, additional_data::AdditionalData, Q_idx::Int, away_team::Int, home_team::Int, auxiliary_data::AuxiliaryData)
    heuristic_estimate_delta_away_team = -additional_data.heuristic_estimates[away_team,Q_idx]
    streak = instance.streak_limit - (additional_data.possible_away_streaks[away_team,Q_idx] - 1)
    if home_team < away_team
        teams_left_mask = additional_data.away_teams_left_masks[away_team,Q_idx] - 2 ^ (home_team-1)
    else
        teams_left_mask = additional_data.away_teams_left_masks[away_team,Q_idx] - 2 ^ (home_team-2)
    end

    heuristic_estimate_delta_away_team += TTPGuidance.heuristic_estimate(instance.streak_limit, instance.d, away_team, teams_left_mask*1, additional_data.number_of_away_games_left[away_team,Q_idx]-1, additional_data.number_of_home_games_left[away_team,Q_idx]+0, home_team, streak, auxiliary_data.bounds_by_state)
    heuristic_estimate_delta_away_team
end

function prepare_auxiliary_data(instance::Instance, beam::TTPBeam, Q_idx::Int, layer::Int)
    for i in 1:instance.n
        instance.aux_by_thread[Threads.threadid()].number_of_away_games_left[i] = 0
        instance.aux_by_thread[Threads.threadid()].number_of_home_games_left[i] = 0
        instance.aux_by_thread[Threads.threadid()].number_of_permitted_games_left[i] = 0
        instance.aux_by_thread[Threads.threadid()].away_teams_left_masks[i] = 0
    end

    for game_idx in (layer+1):number_of_games(instance)
        j = beam.additional_data.games[game_idx,Q_idx]
        game = instance.games[j]
        instance.aux_by_thread[Threads.threadid()].number_of_away_games_left[game[1]] += 1
        instance.aux_by_thread[Threads.threadid()].number_of_home_games_left[game[2]] += 1
        if game[2] < game[1]
            instance.aux_by_thread[Threads.threadid()].away_teams_left_masks[game[1]] += 2 ^ (game[2]-1)
        else
            instance.aux_by_thread[Threads.threadid()].away_teams_left_masks[game[1]] += 2 ^ (game[2]-2)
        end
        if game_allowed(instance, beam.additional_data, Q_idx, game[1], game[2])
            instance.aux_by_thread[Threads.threadid()].number_of_permitted_games_left[game[1]] += 1
            instance.aux_by_thread[Threads.threadid()].number_of_permitted_games_left[game[2]] += 1
        end
    end
end

function next_team(instance::Instance, additional_data::AdditionalData, Q_idx::Int)
    #argmin(map(x -> (node.state.rounds[x], teams_permutation[x]), 1:ttp_instance.n))
    min_team = instance.teams_permutation[1]
    min_round = additional_data.rounds[min_team,Q_idx]
    max_round = additional_data.rounds[min_team,Q_idx]
    for i in 2:instance.n
        team = instance.teams_permutation[i]
        if additional_data.rounds[team,Q_idx] < min_round
            min_team = team
            min_round -= 1
            break
        end
    end
    min_team, min_round, max_round
end

function game_allowed(instance::Instance, additional_data::AdditionalData, Q_idx::Int, away_team::Int, home_team::Int)
    additional_data.rounds[away_team,Q_idx] == additional_data.rounds[home_team,Q_idx] && additional_data.possible_away_streaks[away_team,Q_idx] > 0 && additional_data.possible_home_stands[home_team,Q_idx] > 0 && (!instance.no_repeat || (additional_data.forbidden_opponents[away_team,Q_idx] != home_team && additional_data.forbidden_opponents[home_team,Q_idx] != away_team))
end

function travel_distance_after_next_game(instance::Instance, beam::TTPBeam, Q_idx::Int, away_team::Int, home_team::Int)
    beam.nodes.costs_so_far[Q_idx] + instance.d[beam.additional_data.positions[away_team,Q_idx],home_team] + instance.d[beam.additional_data.positions[home_team,Q_idx],home_team]
end

function is_terminal(instance::Instance, beam::TTPBeam, Q_idx::Int)
    beam.nodes.layer[Q_idx] == number_of_games(instance)
end

function terminal_costs(instance::Instance, beam::TTPBeam, Q_idx::Int)
    remaining_terminal_costs = 0
    for i in 1:instance.n
        remaining_terminal_costs += instance.d[beam.additional_data.positions[i,Q_idx], i]
    end
    remaining_terminal_costs
end

function copy(instance::Instance, from_Q::TTPBeam, to_Q::TTPBeam, from_idx::Int, to_idx::Int)
    to_Q.nodes.layer[to_idx] = from_Q.nodes.layer[from_idx]
    to_Q.nodes.costs_so_far[to_idx] = from_Q.nodes.costs_so_far[from_idx]
    to_Q.nodes.priority[to_idx] = from_Q.nodes.priority[from_idx]
    for i in 1:number_of_games(instance)
        to_Q.additional_data.games[i,to_idx] = from_Q.additional_data.games[i,from_idx]
    end
    for i in 1:instance.n
        to_Q.additional_data.rounds[i,to_idx] = from_Q.additional_data.rounds[i,from_idx]
        to_Q.additional_data.forbidden_opponents[i,to_idx] = from_Q.additional_data.forbidden_opponents[i,from_idx]
        to_Q.additional_data.positions[i,to_idx] = from_Q.additional_data.positions[i,from_idx]
        to_Q.additional_data.possible_away_streaks[i,to_idx] = from_Q.additional_data.possible_away_streaks[i,from_idx]
        to_Q.additional_data.possible_home_stands[i,to_idx] = from_Q.additional_data.possible_home_stands[i,from_idx]
        to_Q.additional_data.heuristic_estimates[i,to_idx] = from_Q.additional_data.heuristic_estimates[i,from_idx]
        to_Q.additional_data.number_of_away_games_left[i,to_idx] = from_Q.additional_data.number_of_away_games_left[i,from_idx]
        to_Q.additional_data.number_of_home_games_left[i,to_idx] = from_Q.additional_data.number_of_home_games_left[i,from_idx]
        to_Q.additional_data.away_teams_left_masks[i,to_idx] = from_Q.additional_data.away_teams_left_masks[i,from_idx]
    end
end

function copy(from_successor::TTPSuccessorData, to_successors::TTPSuccessorData, from_idx::Int, to_idx::Int)
    to_successors.layer[to_idx] = from_successor.layer[from_idx]
    to_successors.costs_after_move[to_idx] = from_successor.costs_after_move[from_idx]
    to_successors.priority[to_idx] = from_successor.priority[from_idx]
    to_successors.Q_idx[to_idx] = from_successor.Q_idx[from_idx]
    to_successors.successor_move[to_idx] = from_successor.successor_move[from_idx]
    to_successors.additional_data.heuristic_estimates_delta_away_team[to_idx] = from_successor.additional_data.heuristic_estimates_delta_away_team[from_idx]
    to_successors.additional_data.heuristic_estimates_delta_home_team[to_idx] = from_successor.additional_data.heuristic_estimates_delta_home_team[from_idx]
end

function make_transition(configuration::Configuration, instance::Instance, successor_specs::TTPSuccessorData, Q::TTPBeam, Q_idx::Int, spec_idx::Int)
    Q.nodes.layer[Q_idx] += 1
    Q.nodes.costs_so_far[Q_idx] = successor_specs.costs_after_move[spec_idx]
    Q.nodes.priority[Q_idx] = successor_specs.priority[spec_idx]
    j = successor_specs.successor_move[spec_idx][1]
    heuristic_estimates_delta_away_team = successor_specs.additional_data.heuristic_estimates_delta_away_team[spec_idx]
    heuristic_estimates_delta_home_team = successor_specs.additional_data.heuristic_estimates_delta_home_team[spec_idx]
    make_move(instance, Q, Q_idx, j, heuristic_estimates_delta_away_team, heuristic_estimates_delta_home_team)
end

function make_move(instance::Instance, Q::TTPBeam, Q_idx::Int, next_j::UInt16, heuristic_estimates_delta_away_team::Int, heuristic_estimates_delta_home_team::Int)
    next_game = instance.games[next_j]
    additonal_data_after_next_game!(instance, Q, Q_idx, next_game[1], next_game[2], heuristic_estimates_delta_away_team, heuristic_estimates_delta_home_team)
    layer = Q.nodes.layer[Q_idx]
    j_to_swap = Q.additional_data.games[layer,Q_idx]
    @assert 1 <= next_j <= number_of_games(instance)
    @assert 1 <= j_to_swap <= number_of_games(instance)
    if next_j != j_to_swap
        target_idx = -1
        for j_idx in (layer+1):number_of_games(instance)
            if Q.additional_data.games[j_idx,Q_idx] == next_j
                target_idx = j_idx
                break 
            end
        end
        @assert target_idx != -1
        Q.additional_data.games[layer,Q_idx], Q.additional_data.games[target_idx,Q_idx] = Q.additional_data.games[target_idx,Q_idx], Q.additional_data.games[layer,Q_idx]
    end
end

function additonal_data_after_next_game!(instance::Instance, Q::TTPBeam, idx::Int, away_team::Int, home_team::Int, heuristic_estimates_delta_away_team::Int, heuristic_estimates_delta_home_team::Int)
    #print_state_and_game(instance, Q, idx, away_team, home_team)
    # heuristic estimates delta
    Q.additional_data.heuristic_estimates[away_team,idx] += heuristic_estimates_delta_away_team
    Q.additional_data.heuristic_estimates[home_team,idx] += heuristic_estimates_delta_home_team

    # rounds
    Q.additional_data.rounds[away_team,idx] += 1
    Q.additional_data.rounds[home_team,idx] += 1
        
    # streaks
    max_away_streak_away_team = min(Q.additional_data.number_of_away_games_left[away_team,idx], instance.streak_limit)
    max_home_stand_away_team = min(Q.additional_data.number_of_home_games_left[away_team,idx], instance.streak_limit)
    max_away_streak_home_team = min(Q.additional_data.number_of_away_games_left[home_team,idx], instance.streak_limit)
    max_home_stand_home_team = min(Q.additional_data.number_of_home_games_left[home_team,idx], instance.streak_limit)

    if Q.additional_data.positions[away_team,idx] == away_team
        Q.additional_data.possible_away_streaks[away_team,idx] = max_away_streak_away_team - 1
        Q.additional_data.possible_home_stands[away_team,idx] = max_home_stand_away_team
    else
        Q.additional_data.possible_away_streaks[away_team,idx] -= 1
        Q.additional_data.possible_home_stands[away_team,idx] = max_home_stand_away_team
    end

    if Q.additional_data.positions[home_team,idx] == home_team
        Q.additional_data.possible_away_streaks[home_team,idx] = max_away_streak_home_team
        Q.additional_data.possible_home_stands[home_team,idx] -= 1
    else
        Q.additional_data.possible_away_streaks[home_team,idx] = max_away_streak_home_team
        Q.additional_data.possible_home_stands[home_team,idx] = max_home_stand_home_team - 1
    end

    # number of games left
    Q.additional_data.number_of_away_games_left[away_team,idx] -= 1
    Q.additional_data.number_of_home_games_left[home_team,idx] -= 1

    repeater_exists = false
    if home_team < away_team
        if Q.additional_data.away_teams_left_masks[home_team,idx] & 2 ^ (away_team-2) > 0
            repeater_exists = true
        end
        Q.additional_data.away_teams_left_masks[away_team,idx] -= 2 ^ (home_team-1)
    else
        if Q.additional_data.away_teams_left_masks[home_team,idx] & 2 ^ (away_team-1) > 0
            repeater_exists = true
        end
        Q.additional_data.away_teams_left_masks[away_team,idx] -= 2 ^ (home_team-2)
    end

    # positions
    Q.additional_data.positions[away_team,idx] = home_team
    Q.additional_data.positions[home_team,idx] = home_team

    # forbidden opponents
    if repeater_exists
        Q.additional_data.forbidden_opponents[away_team,idx] = home_team
        Q.additional_data.forbidden_opponents[home_team,idx] = away_team
    else
        Q.additional_data.forbidden_opponents[away_team,idx] = 0
        Q.additional_data.forbidden_opponents[home_team,idx] = 0
    end

    for team in 1:instance.n
        if team != away_team && team != home_team && (Q.additional_data.forbidden_opponents[team,idx] == away_team || Q.additional_data.forbidden_opponents[team,idx] == home_team)
            Q.additional_data.forbidden_opponents[team,idx] = 0
        end
    end

    @assert 0 <= Q.additional_data.possible_away_streaks[away_team,idx] <= instance.streak_limit
    @assert 0 <= Q.additional_data.possible_away_streaks[home_team,idx] <= instance.streak_limit
    @assert 0 <= Q.additional_data.possible_home_stands[away_team,idx] <= instance.streak_limit
    @assert 0 <= Q.additional_data.possible_home_stands[home_team,idx] <= instance.streak_limit
end

function print_state_and_game(instance, beam, Q_idx, away_team, home_team)
    @printf("game (%d, %d) to be played\n", away_team, home_team)
    @printf("rounds: %s\n", beam.additional_data.rounds[:,Q_idx])
    @printf("positions: %s\n", beam.additional_data.positions[:,Q_idx])
    @printf("forbidden opponents: %s\n", beam.additional_data.forbidden_opponents[:,Q_idx])
    @printf("possible away streaks: %s\n", beam.additional_data.possible_away_streaks[:,Q_idx])
    @printf("possible home stands: %s\n", beam.additional_data.possible_home_stands[:,Q_idx])
end

function print_custom_stats(instance::Instance)
    for i in 1:Threads.nthreads()
        @printf("thread %d, eval succ prep time: %.02f\n", i, instance.aux_by_thread[i].statistics.eval_successors_prep_time)
        @printf("thread %d, eval succ loop time: %.02f\n", i, instance.aux_by_thread[i].statistics.eval_successors_loop_time)
    end
end

function delta_feasibility_check(instance::Instance, additional_data::AdditionalData, Q_idx::Int, away_team::Int, home_team::Int)
    away_team_home_games_left = additional_data.number_of_home_games_left[away_team,Q_idx]*1
    away_team_away_games_left = additional_data.number_of_away_games_left[away_team,Q_idx] - 1
    home_team_home_games_left = additional_data.number_of_home_games_left[home_team,Q_idx] - 1
    home_team_away_games_left = additional_data.number_of_away_games_left[home_team,Q_idx]*1

    if !delta_at_most_check(instance, away_team_home_games_left, away_team_away_games_left, home_team_home_games_left, home_team_away_games_left)
        return false
    end

    return true
end

function delta_at_most_check(instance::Instance, away_team_home_games_left::Int, away_team_away_games_left::Int, home_team_home_games_left::Int, home_team_away_games_left::Int)
    if home_team_home_games_left * instance.streak_limit + instance.streak_limit < home_team_away_games_left
        return false
    end

    if away_team_away_games_left * instance.streak_limit + instance.streak_limit < away_team_home_games_left
        return false
    end

    return true
end

function dead_team_check(instance::Instance, beam::TTPBeam, Q_idx::Int, min_round::Number, max_round::Number)
    for i in 1:instance.n
        if (min_round == max_round || beam.additional_data.rounds[i,Q_idx] < max_round) && instance.aux_by_thread[Threads.threadid()].number_of_permitted_games_left[i] == 0 && instance.aux_by_thread[Threads.threadid()].number_of_home_games_left[i] + instance.aux_by_thread[Threads.threadid()].number_of_away_games_left[i] > 0
            return false
        end
    end
    
    return true
end

function dead_team_check_without_aux(instance::Instance, beam::TTPBeam, Q_idx::Int, layer::Number, min_round::Number, max_round::Number)
    teams_mask = 0
    for i in 1:instance.n
        if beam.additional_data.number_of_home_games_left[i,Q_idx] + beam.additional_data.number_of_away_games_left[i,Q_idx] == 0 || (min_round != max_round && beam.additional_data.rounds[i,Q_idx] == max_round)
            teams_mask |= 2^(i-1)
        end
    end
    for game_idx in (layer+1):number_of_games(instance)
        j = beam.additional_data.games[game_idx,Q_idx]
        game = instance.games[j]
        if game_allowed(instance, beam.additional_data, Q_idx, game[1], game[2])
            teams_mask |= 2^(game[1]-1)
            teams_mask |= 2^(game[2]-1)
        end
    end
    return teams_mask == 2^instance.n - 1
end

function state_hash(configuration::Configuration, instance::Instance, Q::TTPBeam, Q_idx::Int, layer::Number)
    away_teams_left_mask = @view Q.additional_data.away_teams_left_masks[:,Q_idx]
    forbidden_opponents = @view Q.additional_data.forbidden_opponents[:,Q_idx]
    positions = @view Q.additional_data.positions[:,Q_idx]
    possible_away_streaks = @view Q.additional_data.possible_away_streaks[:,Q_idx]
    possible_home_stands = @view Q.additional_data.possible_home_stands[:,Q_idx]

    hash(away_teams_left_mask, hash(forbidden_opponents, hash(positions, hash(possible_away_streaks, hash(possible_home_stands)))))
end

function isequal_by_idx(Q_a::TTPBeam, Q_b::TTPBeam, Q_idx_a::Int, Q_idx_b::Int, layer::Number)
    away_teams_left_mask_a = @view Q_a.additional_data.away_teams_left_masks[:,Q_idx_a]
    away_teams_left_mask_b = @view Q_b.additional_data.away_teams_left_masks[:,Q_idx_b]
    forbidden_opponents_a = @view Q_a.additional_data.forbidden_opponents[:,Q_idx_a]
    forbidden_opponents_b = @view Q_b.additional_data.forbidden_opponents[:,Q_idx_b]
    positions_a = @view Q_a.additional_data.positions[:,Q_idx_a]
    positions_b = @view Q_b.additional_data.positions[:,Q_idx_b]
    possible_away_streaks_a = @view Q_a.additional_data.possible_away_streaks[:,Q_idx_a]
    possible_away_streaks_b = @view Q_b.additional_data.possible_away_streaks[:,Q_idx_b]
    possible_home_stands_a = @view Q_a.additional_data.possible_home_stands[:,Q_idx_a]
    possible_home_stands_b = @view Q_b.additional_data.possible_home_stands[:,Q_idx_b]

    isequal(away_teams_left_mask_a, away_teams_left_mask_b) && isequal(forbidden_opponents_a, forbidden_opponents_b) && isequal(positions_a, positions_b) && isequal(possible_away_streaks_a, possible_away_streaks_b) && isequal(possible_home_stands_a, possible_home_stands_b)
end

function isdominated_by_idx(Q_a::TTPBeam, Q_b::TTPBeam, Q_idx_a::Int, Q_idx_b::Int)
    away_teams_left_mask_a = @view Q_a.additional_data.away_teams_left_masks[:,Q_idx_a]
    away_teams_left_mask_b = @view Q_b.additional_data.away_teams_left_masks[:,Q_idx_b]
    forbidden_opponents_a = @view Q_a.additional_data.forbidden_opponents[:,Q_idx_a]
    forbidden_opponents_b = @view Q_b.additional_data.forbidden_opponents[:,Q_idx_b]
    positions_a = @view Q_a.additional_data.positions[:,Q_idx_a]
    positions_b = @view Q_b.additional_data.positions[:,Q_idx_b]
    possible_away_streaks_a = @view Q_a.additional_data.possible_away_streaks[:,Q_idx_a]
    possible_away_streaks_b = @view Q_b.additional_data.possible_away_streaks[:,Q_idx_b]
    possible_home_stands_a = @view Q_a.additional_data.possible_home_stands[:,Q_idx_a]
    possible_home_stands_b = @view Q_b.additional_data.possible_home_stands[:,Q_idx_b]

    Q_a.nodes.priority[Q_idx_b] <= Q_a.nodes.priority[Q_idx_a] && isequal(away_teams_left_mask_a, away_teams_left_mask_b) && isequal(forbidden_opponents_a, forbidden_opponents_b) && isequal(positions_a, positions_b) && all(possible_away_streaks_a .<= possible_away_streaks_b) && all(possible_home_stands_a .<= possible_home_stands_b)
end

function successors_size(successor_specs::TTPSuccessorData)
    (sizeof(successor_specs.layer)+sizeof(successor_specs.costs_after_move)+sizeof(successor_specs.priority)+sizeof(successor_specs.Q_idx)+sizeof(successor_specs.successor_move)+sizeof(successor_specs.additional_data.heuristic_estimates_delta_away_team)+sizeof(successor_specs.additional_data.heuristic_estimates_delta_home_team))/(1024^3)
end

function beam_size(Q::TTPBeam)
    (sizeof(Q.nodes.layer)+sizeof(Q.nodes.costs_so_far)+sizeof(Q.nodes.priority)+sizeof(Q.additional_data.games)+sizeof(Q.additional_data.forbidden_opponents)+sizeof(Q.additional_data.rounds)+sizeof(Q.additional_data.positions)+sizeof(Q.additional_data.possible_away_streaks)+sizeof(Q.additional_data.possible_home_stands)+sizeof(Q.additional_data.number_of_away_games_left)+sizeof(Q.additional_data.number_of_home_games_left)+sizeof(Q.additional_data.heuristic_estimates)+sizeof(Q.additional_data.away_teams_left_masks))/(1024^3)
end

function instance_size(instance::Instance)
    sizeof(instance.aux_by_thread[1].bounds_by_state)/1024^3
end

function pre_run_instance_name(instance_name::String)::String
    "NL/nl6"
end
   
@time main()