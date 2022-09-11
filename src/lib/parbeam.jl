# parbeam.jl: a parallel beam search framework for combinatorial optimization problems

import Base.isequal
import Base.isless
using Printf
using MPI
using Logging
using Random
using Printf: @printf, @sprintf

# abstract types
abstract type AbstractInstance end
abstract type AbstractNodeData end
abstract type AbstractSuccessorData end
abstract type AbstractBeam end
abstract type AbstractAdditionalData end
abstract type AbstractAdditionalSuccessorData end
abstract type AbstractAdditionalConfiguration end

# structs
@enum ProblemType begin
    mini = 1
    maxi = 2
end

mutable struct Configuration
    problem_type::ProblemType
    beam_width::Int
    maximum_number_of_successors_by_node::Int
    maximum_layer_number::Int
    variable_ordering::String
    guidance_function::String
    tie_breaking::String
    shrink_successor_specs::Bool
    cacheline_length_int::Int
    successor_creation_load_balancing::String
    states_filtering::String
    additional_configuration::AbstractAdditionalConfiguration
end

struct EmptyAdditionalConfiguration <: AbstractAdditionalConfiguration end

mutable struct NodeData{T<:Number,U<:Number,V<:Number} <: AbstractNodeData
    layer::Array{V, 1}
    costs_so_far::Array{T, 1}
    priority::Array{U, 1}
end

mutable struct SuccessorData{T<:Number,U<:Number,V<:Number,W<:AbstractAdditionalSuccessorData,X<:Number} <: AbstractSuccessorData
    layer::Array{X, 1}
    costs_after_move::Array{T, 1}
    priority::Array{U, 1}
    Q_idx::Array{Int, 1}
    successor_move::Array{V, 1}
    additional_data::W
end

mutable struct Beam{T<:AbstractNodeData,U<:AbstractAdditionalData} <: AbstractBeam
    nodes::T
    additional_data::U
end

mutable struct Statistics
    main_loop_time::Float64
    instance_preparation_time::Float64
    memory_init_time::Float64
    successor_creation_time::Float64
    compress_successors_time::Float64
    surviving_successors_time::Float64
    transitions_time::Float64
    terminal_check_time::Float64
    sequential_counting_sort_time::Float64
    parallel_counting_sort_time_min_max::Float64
    parallel_counting_sort_time_binning::Float64
    parallel_counting_sort_time_cutting::Float64
    tie_breaking_time::Float64
    successor_specs_shrinking_time::Float64
    number_of_non_tie_successors::Int
    number_of_tie_successors::Int
    successors_creation_time_loop::Float64
    successors_creation_time_prep::Float64
    number_of_successors::Int
    number_of_surviving_successors::Int
    states_filtering_time::Float64
    number_of_states_filtered::Int
    Statistics() = new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0.0, 0, 0, 0.0, 0)
end

# abstract functions
function additional_configuration_type()::DataType
    EmptyAdditionalConfiguration
end

function problem_type(instance::AbstractInstance)
    mini
end

function read_instance(guidance_function::String, variable_ordering::String, instance_description::String)::AbstractInstance
end

function prepare_beam_and_successor_specs(configuration::Configuration, instance::AbstractInstance)
end

function initial_estimate(configuration::Configuration, instance::AbstractInstance, beam::AbstractBeam)::Number
end

function is_terminal(instance::AbstractInstance, beam::AbstractBeam, Q_idx::Int)::Bool
end

function maximum_number_of_successors_per_node(instance::AbstractInstance)
end

function maximum_layer_number(instance::AbstractInstance)
    throw("not implemented")
end

function initialize_beam(configuration::Configuration, instance::AbstractInstance)
end

function initialize_successor_specs(configuration::Configuration, instance::AbstractInstance)::AbstractSuccessorData
end

function eval_successors(configuration::Configuration, instance::AbstractInstance, beam::AbstractBeam, successor_specs::AbstractSuccessorData, start_index::Int, Q_idx::Int, sigma::Float64)
end

function copy(instance::AbstractInstance, from_Q::AbstractBeam, to_Q::AbstractBeam, from_idx::Int, to_idx::Int)
end

function copy(from_successor::AbstractSuccessorData, to_successors::AbstractSuccessorData, from_idx::Int, to_idx::Int)
end

function make_transition(configuration::Configuration, instance::AbstractInstance, successor_specs::AbstractSuccessorData, Q::AbstractBeam, Q_idx::Int, spec_idx::Int)
end

function initial_move_down(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, successor_moves)
end

function initial_successor_moves(instance::AbstractInstance, depth::Int)
end

function create_objectives_and_solutions()::Vector{Tuple{Float64,Array{Int}}}
end

function create_all_objectives(size::Number)::Array{Array{Float64}}
end

function create_all_solutions(instance::AbstractInstance, size::Number)::Vector{Array{Int}}
end

function get_solution(instance::AbstractInstance, beam::AbstractBeam, Q_idx::Int)::Array{Int}
end

function terminal_costs(instance::AbstractInstance, beam::AbstractBeam, Q_idx::Int)
end

function print_custom_stats(instance::AbstractInstance)
end

function state_hash(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_idx::Int, layer::Number)
    throw("not implemented")
end

function isequal_by_idx(Q_a::AbstractBeam, Q_b::AbstractBeam, Q_idx_a::Int, Q_idx_b::Int, layer::Number)
    throw("not implemented")
end

function isdominated_by_idx(Q_a::AbstractBeam, Q_b::AbstractBeam, Q_idx_a::Int, Q_idx_b::Int)
    throw("not implemented")
end

function beam_size(Q::AbstractBeam)
    0
end

function successors_size(successor_specs::AbstractSuccessorData)
    0
end

function instance_size(instance::AbstractInstance)
    0
end

# concrete functions
function initialize_successor_spec_layers(successor_data::AbstractSuccessorData)
    Threads.@threads for i in 1:length(successor_data.layer)
        successor_data.layer[i] = 0
    end
end

function create_successors(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, successor_specs::AbstractSuccessorData, sigma::Float64, layer::Number)
    number_of_successors_arr = zeros(Int, Threads.nthreads())
    #time_per_thread = zeros(Float64, Threads.nthreads())
    total_time = @elapsed if Threads.nthreads() > 1 && configuration.successor_creation_load_balancing == "shuffle"
        partitioned_ranges = collect(Iterators.partition(1:Q_len, 100*configuration.cacheline_length_int))
        random_ordering = randperm(length(partitioned_ranges))

        Threads.@threads for j::Int in 1:length(partitioned_ranges::Vector{UnitRange{Int64}})
            for i::Int in partitioned_ranges[random_ordering[j]::Int]::UnitRange{Int64}
                if Q.nodes.layer[i] == layer
                    successors_count = eval_successors(configuration, instance, Q, successor_specs, (i-1)*configuration.maximum_number_of_successors_by_node+1, i, sigma)
                    number_of_successors_arr[Threads.threadid()] += successors_count
                end
            end
        end
    else
        if configuration.states_filtering == "none"
            Threads.@threads for i in 1:Q_len
                #time_per_thread[Threads.threadid()] += @elapsed successors_count = eval_successors(configuration, instance, Q, successor_specs, (i-1)*configuration.maximum_number_of_successors_by_node+1, i, sigma)
                successors_count = eval_successors(configuration, instance, Q, successor_specs, (i-1)*configuration.maximum_number_of_successors_by_node+1, i, sigma)
                number_of_successors_arr[Threads.threadid()] += successors_count
            end
        else
            Threads.@threads for i in 1:Q_len
                if Q.nodes.layer[i] == layer
                    successors_count = eval_successors(configuration, instance, Q, successor_specs, (i-1)*configuration.maximum_number_of_successors_by_node+1, i, sigma)
                    number_of_successors_arr[Threads.threadid()] += successors_count
                end
            end
        end

    end
    #@printf("[CSV-threads-successors] %d;%d;%.02f;%s;%s\n", layer, Threads.nthreads(), total_time, join(number_of_successors_arr, ';'), join(time_per_thread, ';'))
    sum(number_of_successors_arr)
end

function compress_successors(successor_specs::AbstractSuccessorData, successor_count::Int, layer::Number)
    dst = 1
    src = length(successor_specs.layer)
    next_layer = layer + 1

    while dst < src
        # find next empty destination (from front)
        while successor_specs.layer[dst] == next_layer && dst <= successor_count && dst < src
            dst += 1
        end
        if dst >= src || dst > successor_count
            break
        end
        @assert successor_specs.layer[dst] != next_layer
        # find next non-empty source (from back)
        while successor_specs.layer[src] != next_layer && src > dst
            src -= 1
        end
        if src <= dst
            break
        end
        @assert successor_specs.layer[src] == next_layer
        #@printf("copy from %d to %d\n", src, dst)
        copy(successor_specs, successor_specs, src, dst)
        src -= 1
        @assert successor_specs.layer[dst] == next_layer
    end
end

function surviving_successors(configuration::Configuration, instance::AbstractInstance, layer::Number, successor_specs::AbstractSuccessorData, current_successor_spec_count::Int, Q::AbstractBeam, Q_next::AbstractBeam, stats::Statistics)
    min_prios = ones(Int, Threads.nthreads())*typemax(Int)
    max_prios = -ones(Int, Threads.nthreads())*typemax(Int)

    stats.parallel_counting_sort_time_min_max += @elapsed Threads.@threads for i in 1:length(successor_specs.layer)
        if successor_specs.layer[i] == layer + 1
            converted_prio = convert(Int, floor(successor_specs.priority[i]))
            if converted_prio < min_prios[Threads.threadid()] 
                min_prios[Threads.threadid()] = converted_prio
            end
            if converted_prio > max_prios[Threads.threadid()] 
                max_prios[Threads.threadid()] = converted_prio
            end
        end
    end

    min_prio = minimum(min_prios)
    max_prio = maximum(max_prios)

    #println(min_prio)
    #println(max_prio)

    #println(current_successor_spec_count)
    min_len = min(current_successor_spec_count, configuration.beam_width)

    if min_len > 0
        cut_off_prio = surviving_successor_specs_from_counting_sort!(configuration, instance, successor_specs, current_successor_spec_count, Q, Q_next, layer, min_len, min_prio, max_prio, stats)
    end

    min_len
end

function surviving_successor_specs_from_counting_sort!(configuration::Configuration, instance::AbstractInstance, successor_specs::AbstractSuccessorData, current_successor_spec_count::Int, Q::AbstractBeam, Q_next::AbstractBeam, layer::Number, k::Int, min_priority::Int, max_priority::Int, stats::Statistics)::Number
    @assert max_priority >= min_priority
    span = max_priority - min_priority + 1
    # avoid false sharing, put safety cacheline length between counting array regions, watch also out for column major ordering
    stats.sequential_counting_sort_time += @elapsed counting_arrs = zeros(Int, (span+configuration.cacheline_length_int, Threads.nthreads()))
    stats.parallel_counting_sort_time_binning += @elapsed Threads.@threads for i in 1:length(successor_specs.layer)
        if successor_specs.layer[i] == layer + 1
            converted_prio = convert(Int, floor(successor_specs.priority[i]))
            counting_arrs[converted_prio - min_priority + 1, Threads.threadid()] += 1
        end
    end
    stats.sequential_counting_sort_time += @elapsed counting_arr = sum(counting_arrs, dims=2)
    #prefix_sum = cumsum(counting_arr, dims=2)
    #@printf("[CSV3] %d;%d;%d;%s\n", layer, min_priority, max_priority, join(prefix_sum, ","))

    counter = 0
    number_of_ties = 0
    maximum_number_of_ties = -1
    cut_off_prio = nothing
    number_of_non_ties = -1

    stats.sequential_counting_sort_time += @elapsed for i in 1:span
        current_prio = min_priority + i - 1
        if counter + counting_arr[i] >= k
            cut_off_prio = current_prio
            number_of_ties = counting_arr[i]
            maximum_number_of_ties = k - counter
            number_of_non_ties = counter
            #@printf("cut off prio: %d with ties: %d can only select ties: %d\n", current_prio, number_of_ties, maximum_number_of_ties)
            break
        end
        counter += counting_arr[i]
    end
    @assert cut_off_prio !== nothing
    @assert number_of_non_ties != -1
    @assert maximum_number_of_ties > 0
    @assert number_of_ties > 0

    stats.number_of_non_tie_successors += number_of_non_ties
    stats.number_of_tie_successors += number_of_ties

    make_transitions_for_surviving_successor_specs_from_counting_sort!(configuration, instance, successor_specs, current_successor_spec_count, Q, Q_next, layer, k, cut_off_prio, number_of_ties, number_of_non_ties, maximum_number_of_ties, stats)
    cut_off_prio
end

function make_transitions_for_surviving_successor_specs_from_counting_sort!(configuration::Configuration, instance::AbstractInstance, successor_specs::AbstractSuccessorData, current_successor_spec_count::Int, Q::AbstractBeam, Q_next::AbstractBeam, layer::Number, k::Int, cut_off_priority::Int, number_of_ties::Int, number_of_non_ties::Int, maximum_number_of_ties::Int, stats::Statistics)
    #non_tie_successor_spec_indices = Array{Int, 1}(undef, number_of_non_ties)
    non_tie_successor_spec_indices_by_thread = Array{Vector{Int}}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        non_tie_successor_spec_indices_by_thread[i] = Vector{Int}()
        sizehint!(non_tie_successor_spec_indices_by_thread, convert(Int, ceil(number_of_non_ties/Threads.nthreads())))
    end
    #tie_breaking_successor_spec_indices = Array{Int, 1}(undef, number_of_ties)
    tie_breaking_successor_spec_indices_by_thread = Array{Vector{Int}}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        tie_breaking_successor_spec_indices_by_thread[i] = Vector{Int}()
        sizehint!(tie_breaking_successor_spec_indices_by_thread, convert(Int, ceil(number_of_ties/Threads.nthreads())))
    end
    #consumed_specs = Threads.Atomic{Int}(0)
    #consumed_ties = Threads.Atomic{Int}(0)
    # work_by_threads_count = zeros(Float64, Threads.nthreads())
    # work_by_threads_count[Threads.threadid()] += @elapsed 

    stats.parallel_counting_sort_time_cutting += @elapsed Threads.@threads for i in 1:length(successor_specs.layer)
        if successor_specs.layer[i] == layer + 1
            converted_prio = convert(Int, floor(successor_specs.priority[i]))
            if converted_prio < cut_off_priority
                #j = Threads.atomic_add!(consumed_specs, 1) + 1
                #non_tie_successor_spec_indices[j] = i
                push!(non_tie_successor_spec_indices_by_thread[Threads.threadid()], i)
            elseif converted_prio == cut_off_priority
                #j = Threads.atomic_add!(consumed_ties, 1) + 1
                #tie_breaking_successor_spec_indices[j] = i
                push!(tie_breaking_successor_spec_indices_by_thread[Threads.threadid()], i)
            end
        end
    end

    tie_breaking_successor_spec_indices = vcat(tie_breaking_successor_spec_indices_by_thread...)
    non_tie_successor_spec_indices = vcat(non_tie_successor_spec_indices_by_thread...)

    #consumed_specs = sum(map(length, non_tie_successor_spec_indices_arr))
    #non_tie_successor_spec_indices = vcat(non_tie_successor_spec_indices_arr...)
    #println(work_by_threads_count)

    #@assert consumed_ties[] == number_of_ties
    stats.tie_breaking_time += @elapsed tie_breaking_successor_spec_priorities = @view successor_specs.priority[tie_breaking_successor_spec_indices]
    #@assert consumed_specs[] + maximum_number_of_ties == k

    # non ties state transitions
    stats.transitions_time += @elapsed Threads.@threads for j in 1:number_of_non_ties
        i = non_tie_successor_spec_indices[j]
        copy(instance, Q, Q_next, successor_specs.Q_idx[i], j)
        make_transition(configuration, instance, successor_specs, Q_next, j, i)
    end

    # tie breaking
    stats.tie_breaking_time += @elapsed if configuration.tie_breaking == "lex"
    elseif configuration.tie_breaking == "random"
        sorted_tie_breaking_indices = @view randperm(number_of_ties)[1:maximum_number_of_ties]
    elseif configuration.tie_breaking == "quicksort"
        sorted_tie_breaking_indices = partialsortperm(tie_breaking_successor_spec_priorities, 1:maximum_number_of_ties)
    else
        throw(@sprintf("tie breaking %s not implemented", configuration.tie_breaking))
    end
    if configuration.tie_breaking == "lex"
        stats.transitions_time += @elapsed Threads.@threads for i in 1:maximum_number_of_ties
            successor_spec_idx = tie_breaking_successor_spec_indices[i]
            copy(instance, Q, Q_next, successor_specs.Q_idx[successor_spec_idx], number_of_non_ties + i)
            make_transition(configuration, instance, successor_specs, Q_next, number_of_non_ties + i, successor_spec_idx)
        end
    else
        stats.transitions_time += @elapsed Threads.@threads for i in 1:maximum_number_of_ties
            successor_spec_idx = tie_breaking_successor_spec_indices[sorted_tie_breaking_indices[i]]
            copy(instance, Q, Q_next, successor_specs.Q_idx[successor_spec_idx], number_of_non_ties + i)
            make_transition(configuration, instance, successor_specs, Q_next, number_of_non_ties + i, successor_spec_idx)
        end
    end
end

#function make_transitions(instance::AbstractInstance, surviving_successor_spec_indices::Array{Int}, successor_specs, Q, Q_next, Q_len::Int)
#    Threads.@threads for i in 1:Q_len
#        successor_spec_index = surviving_successor_spec_indices[i]
#        successor_spec = successor_specs[successor_spec_index]
#        copy(Q[successor_spec.Q_idx], Q_next[i])
#        make_transition(instance, successor_spec, Q_next[i])
#    end
#end

function is_new_objective_strictly_better(configuration::Configuration, new_objective, old_objective)
    if configuration.problem_type == mini
        new_objective < old_objective
    elseif configuration.problem_type == maxi
        new_objective > old_objective
    end
end

function check_for_terminals(configuration::Configuration, instance::AbstractInstance, beam::AbstractBeam, Q_len::Int, best_solution, best_objective)
    for i in 1:Q_len
        if is_terminal(instance, beam, i)
            remaining_terminal_costs = terminal_costs(instance, beam, i)
            if is_new_objective_strictly_better(configuration, beam.nodes.costs_so_far[i] + remaining_terminal_costs, best_objective)
                best_solution = get_solution(instance, beam, i)
                best_objective = beam.nodes.costs_so_far[i] + remaining_terminal_costs
            end
        end
    end

    best_solution, best_objective
end

function parallel_check_for_terminals(configuration::Configuration, instance::AbstractInstance, beam::AbstractBeam, Q_len::Int)
    best_objectives = ones(Int, configuration.cacheline_length_int*Threads.nthreads())*typemax(Int)
    best_terminals = -ones(Int, configuration.cacheline_length_int*Threads.nthreads())

    Threads.@threads for i in 1:Q_len
        if is_terminal(instance, beam, i)
            remaining_terminal_costs = terminal_costs(instance, beam, i)
            if beam.nodes.costs_so_far[i] + remaining_terminal_costs < best_objectives[configuration.cacheline_length_int*Threads.threadid()]
                best_terminals[configuration.cacheline_length_int*Threads.threadid()] = i
                best_objectives[configuration.cacheline_length_int*Threads.threadid()] = beam.nodes.costs_so_far[i] + remaining_terminal_costs
            end
        end
    end

    min_index = argmin(best_objectives)
    if best_terminals[min_index] == -1
        nothing, Inf
    else
        best_terminals[min_index], best_objectives[min_index]
    end
end

function flip_beams(Q, Q_up, Q_down)
    (Q == Q_up) ? (Q_down, Q_up) : (Q_up, Q_down)
end

function flip_beams(is_up::Bool, Q_up, Q_down)
    (is_up) ? (!is_up, Q_down, Q_up) : (!is_up, Q_up, Q_down)
end

function shrink_successor_specs(configuration::Configuration, instance::AbstractInstance, layer::Number, successor_specs::AbstractSuccessorData)
    configuration.maximum_number_of_successors_by_node = instance.n - layer
    @assert configuration.maximum_number_of_successors_by_node >= 0
    resize!(successor_specs.layer, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    resize!(successor_specs.costs_after_move, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    resize!(successor_specs.priority, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    resize!(successor_specs.Q_idx, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    resize!(successor_specs.successor_move, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
end

function print_solution(solution)
    @info solution
end

function instance_csv_string(configuration::Configuration, instance::AbstractInstance)::String
end

function print_stats(stats::Statistics)
    @printf("instance prep time: %.02f\n", stats.instance_preparation_time)
    @printf("main loop time: %.02f\n", stats.main_loop_time)
    @printf("memory init time: %.02f\n", stats.memory_init_time)
    @printf("successor creation time: %.02f\n", stats.successor_creation_time)
    @printf("successor compression time: %.02f\n", stats.compress_successors_time)
    @printf("surviving successor time: %.02f\n", stats.surviving_successors_time)
    @printf("transition time: %.02f\n", stats.transitions_time)
    @printf("terminal check time: %.02f\n", stats.terminal_check_time)
    @printf("sequential counting sort time: %.02f\n", stats.sequential_counting_sort_time)
    @printf("parallel counting sort time min max: %.02f\n", stats.parallel_counting_sort_time_min_max)
    @printf("parallel counting sort time binning: %.02f\n", stats.parallel_counting_sort_time_binning)
    @printf("parallel counting sort time cutting: %.02f\n", stats.parallel_counting_sort_time_cutting)
    @printf("tie breaking time: %.02f\n", stats.tie_breaking_time)
    @printf("successor specs shrinking time: %.02f\n", stats.successor_specs_shrinking_time)
    @printf("states filtering time: %.02f\n", stats.states_filtering_time)
    if stats.number_of_non_tie_successors > 0
        @printf("number of tie successors: %d\n", stats.number_of_tie_successors)
        @printf("number of non-tie successors: %d\n", stats.number_of_non_tie_successors)
        @printf("tie rate: %.04f\n", stats.number_of_tie_successors/(stats.number_of_non_tie_successors + stats.number_of_tie_successors))
    end
    @printf("number of successors: %d\n", stats.number_of_successors)
    @printf("number of surviving successors: %d\n", stats.number_of_surviving_successors)
    @printf("number of states filtered: %d\n", stats.number_of_states_filtered)
end

function print_additional_output(configuration::Configuration, instance::AbstractInstance, best_solution)
end

include("parbeam_states_filtering.jl")

function filter_states(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, layer::Number, lks::Array{Threads.SpinLock}, hash_table::Array{Vector{Int}}, hash_locks::Array{Threads.SpinLock}, m::Int)
    if configuration.states_filtering == "duplicates_spinlocked"
        number_of_states_filtered = filter_duplicate_states_spinlocked(configuration, instance, Q, Q_len, lks, layer)
    elseif configuration.states_filtering == "duplicates_spinlocked_collect"
        number_of_states_filtered = filter_duplicate_states_spinlocked_collect(configuration, instance, Q, Q_len, lks, layer)
    elseif configuration.states_filtering == "duplicates_dict_merge"
        number_of_states_filtered = filter_duplicate_states_dict_merge(configuration, instance, Q, Q_len, layer)
    elseif configuration.states_filtering == "duplicates_separate_dicts"
        number_of_states_filtered = filter_duplicate_states_separate_dicts(configuration, instance, Q, Q_len, layer)
    elseif configuration.states_filtering == "duplicates_separate_dicts_final_merge"
        number_of_states_filtered = filter_duplicate_states_separate_dicts_final_merge(configuration, instance, Q, Q_len, layer)
    elseif configuration.states_filtering == "duplicates_concurrent_queues"
        number_of_states_filtered = filter_duplicate_states_concurrent_queues(configuration, instance, Q, Q_len, hash_table, hash_locks, m, layer)
    elseif configuration.states_filtering == "duplicates_spinlocked_non_dominated_lists"
        number_of_states_filtered = filter_duplicate_states_spinlocked_non_dominated_lists(configuration, instance, Q, Q_len, hash_table, hash_locks, m, layer)
    elseif configuration.states_filtering == "dominated"
        number_of_states_filtered = filter_dominated_states(configuration, instance, Q, Q_len)
    elseif configuration.states_filtering == "none"
    else
        throw(@sprintf("state filtering %s not implemented", configuration.states_filtering))
    end

    number_of_states_filtered
end

function construct(configuration::Configuration, instance::AbstractInstance, Q_up::AbstractBeam, Q_down::AbstractBeam, successor_specs::AbstractSuccessorData, sigma::Float64, stats::Statistics)
    #root_node = root(instance)

    best_solution = nothing
    if configuration.problem_type == mini
        best_objective = Inf
    elseif configuration.problem_type == maxi
        best_objective = -Inf
    end

    #Q_up.nodes[1] = root_node
    Q_len = 1
    Q = Q_up
    Q_next = Q_down
    layer = Q.nodes.layer[1]
    is_up = true
    stats.number_of_surviving_successors += Q_len

    stats.states_filtering_time += @elapsed begin
        if configuration.states_filtering != "none"
            lks = Array{Threads.SpinLock}(undef, 97)
            for i in 1:length(lks)
                lks[i] = Threads.SpinLock()
            end
            hash_table, hash_locks, m = initialize_hash_table(configuration)
        end
    end

    @info "Starting beam search"
    stats.main_loop_time = @elapsed while layer < configuration.maximum_layer_number
        #@debug @sprintf("current layer: %d, Q_len: %d\n", layer, Q_len)
        if configuration.states_filtering != "none"
            stats.states_filtering_time += @elapsed stats.number_of_states_filtered += filter_states(configuration, instance, Q, Q_len, layer, lks, hash_table, hash_locks, m)
        end
        stats.successor_creation_time += @elapsed successors_count = create_successors(configuration, instance, Q, Q_len, successor_specs, sigma, layer)
        #@time stats.compress_successors_time += @elapsed compress_successors(successor_specs, successors_count, layer)
        stats.surviving_successors_time += @elapsed Q_len = surviving_successors(configuration, instance, layer, successor_specs, successors_count, Q, Q_next, stats)

        #@debug @sprintf("new successors %d, survived have: %d\n", successors_count, Q_len)
        stats.number_of_successors += successors_count
        stats.number_of_surviving_successors += Q_len
        if Q_len == 0
            #@warn "no surving successors found!"
            @info @sprintf("No more successors at layer %d\n", layer)
            break
        end
        #stats.transitions_time += @elapsed make_transitions(instance, surviving_successor_spec_indices, successor_specs, Q, Q_next, Q_len)
        stats.terminal_check_time += @elapsed best_solution, best_objective = check_for_terminals(configuration, instance, Q_next, Q_len, best_solution, best_objective)
        if configuration.shrink_successor_specs
            stats.successor_specs_shrinking_time += @elapsed shrink_successor_specs(configuration, instance, layer+1, successor_specs)
        end
        is_up, Q, Q_next = flip_beams(is_up, Q_up, Q_down)

        layer += 1
    end
    #@info @sprintf("Finished beam search after %.02f seconds\n", stats.main_loop_time)

    @debug begin
        print_stats(stats)
        print_custom_stats(instance)
    end

    if best_solution !== nothing
        best_solution, best_objective, stats
    else
        [], best_objective, stats
    end
end

# from the tutorial https://juliafolds.github.io/data-parallelism/tutorials/mutations/
# by Takafumi Arakaki, licensed under https://creativecommons.org/licenses/by-sa/4.0/
function get_cacheline_length()
    try
        parse(Int, read("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size", String))
    catch err
        #@warn "cannot read cache line size" exception = (err, catch_backtrace())
        64
    end
end

function main_mpi(pre_run::Bool=false)
    start_time = time()

    AdditionalConfigurationType = additional_configuration_type()
    additional_configuration_help_string = join(map(x -> "<" * string(x) * ">", fieldnames(AdditionalConfigurationType)), " ")

    if(length(ARGS) != 10 + length(fieldnames(AdditionalConfigurationType)))
        @printf("Usage: %s <instance> <beam_width> <guidance_function> <tie_breaking> <depth_for_subtree_splitting> <number-of-runs> <sigma_rel> <variable_ordering> <successor-creation-load-balancing> <duplicate-states-filtering> %s\n", "$PROGRAM_FILE", additional_configuration_help_string)
        exit(1)
    end

    stats = Statistics()

    instance_description = ARGS[1]
    beam_width = parse(Int, ARGS[2])
    guidance_function = ARGS[3]
    tie_breaking = ARGS[4]
    depth_for_subtree_splitting = parse(Int, ARGS[5])
    number_of_runs = parse(Int, ARGS[6])
    sigma_rel = parse(Float64, ARGS[7])
    variable_ordering = ARGS[8]
    successor_creation_load_balancing = ARGS[9]
    states_filtering = ARGS[10]
    arg_idx = 11

    additional_parameters = []
    for fieldtype in fieldtypes(AdditionalConfigurationType)
        push!(additional_parameters, parse(fieldtype, ARGS[arg_idx]))
        arg_idx += 1
    end

    @info "Initializing MPI"
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)

    # to make runs deterministic
    Random.seed!(7 + rank)

    # run with small beam width for time measurement without compilation time
    if pre_run
        beam_width = 10
        number_of_runs = size
        instance_description = pre_run_instance_name(instance_description)
    end

    @info "Preparing instance"
    stats.instance_preparation_time = @elapsed instance = read_instance(guidance_function, variable_ordering, instance_description)
    @info @sprintf("Finished preparing instance after %.02f s", stats.instance_preparation_time)
    configuration = Configuration(problem_type(instance), beam_width, maximum_number_of_successors_per_node(instance), maximum_layer_number(instance), variable_ordering, guidance_function, tie_breaking, false, cld(get_cacheline_length(), sizeof(Int)), successor_creation_load_balancing, states_filtering, AdditionalConfigurationType(additional_parameters...))
    if !pre_run
        @info @sprintf("configuration: %s\n", configuration)
    else
        @sprintf("configuration: %s\n", configuration)
    end

    objectives_and_solutions = create_objectives_and_solutions()

    @info "Initializing beams and successor specs"
    stats.memory_init_time += @elapsed Q_up = initialize_beam(configuration, instance)
    stats.memory_init_time += @elapsed Q_down = initialize_beam(configuration, instance)
    stats.memory_init_time += @elapsed successor_specs = initialize_successor_specs(configuration, instance)
    
    if !pre_run
        @info @sprintf("one beam size: %f GB\n", beam_size(Q_up))
        @info @sprintf("successors size: %f GB\n", successors_size(successor_specs))
        @info @sprintf("sorting indices size: %f GB\n", sizeof(Int)*configuration.beam_width/1024^3)
        @info @sprintf("instance size: %f GB\n", instance_size(instance))
        @info @sprintf("total estimated memory demand: %f GB", 2*beam_size(Q_up)+successors_size(successor_specs)+instance_size(instance)+sizeof(Int)*configuration.beam_width/1024^3)
    else
        @sprintf("one beam size: %f GB\n", beam_size(Q_up))
        @sprintf("successors size: %f GB\n", successors_size(successor_specs))
        @sprintf("sorting indices size: %f GB\n", sizeof(Int)*configuration.beam_width/1024^3)
        @sprintf("instance size: %f GB\n", instance_size(instance))
        @sprintf("total estimated memory demand: %f GB", 2*beam_size(Q_up)+successors_size(successor_specs)+instance_size(instance)+sizeof(Int)*configuration.beam_width/1024^3)
    end

    stats.memory_init_time += @elapsed initialize_root(instance, Q_up)
    root_estimate = initial_estimate(configuration, instance, Q_up)
    sigma = root_estimate*sigma_rel

    moves_list = initial_successor_moves(instance, depth_for_subtree_splitting)
    for (i, initial_moves) in enumerate(moves_list)
        for k in 1:number_of_runs
            if (number_of_runs*(i - 1) + k - 1) % size == rank
                stats.memory_init_time += @elapsed initialize_root(instance, Q_up)
                stats.memory_init_time += @elapsed initialize_successor_spec_layers(successor_specs)
                initial_priority = initial_move_down(configuration, instance, Q_up, initial_moves)
            
                if !pre_run
                    @info @sprintf("solving instance %s using beam search with beam width %d\n", instance_description, beam_width)
                else
                    @sprintf("solving instance %s using beam search with beam width %d\n", instance_description, beam_width)
                end
                root_estimate_run = initial_estimate(configuration, instance, Q_up)
                #@assert initial_lower_bound == initial_priority
                @info @sprintf("initial estimate of run: %d\n", root_estimate_run)
            
                construction_time = @elapsed best_solution, best_objective, stats = construct(configuration, instance, Q_up, Q_down, successor_specs, sigma, stats)

                @info @sprintf("best solution: %s\n", best_solution)
                @info @sprintf("best objective: %f\n", best_objective)

                push!(objectives_and_solutions, (best_objective, best_solution))
            
                #@printf("[CSV] %s;%d;%s;%f;%f;%f;%f;%f;%f;%f\n", instance_name, beam_width, guidance_function, construction_time, stats.main_loop_time, stats.memory_init_time, stats.successor_creation_time, stats.surviving_successors_time, stats.transitions_time, stats.terminal_check_time)
                #instance_description, beam_width, guidance_function, best_objective_over_subtrees, construction_time, stats
                if !pre_run
                    println(@sprintf("[CSV] %d;%d;%s;%s;%d;%s;%d;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%d;%d", rank, Threads.nthreads(), instance_description, instance_csv_string(configuration, instance), beam_width, guidance_function, root_estimate, best_objective, construction_time, stats.main_loop_time, stats.memory_init_time, stats.instance_preparation_time, stats.successor_creation_time, stats.surviving_successors_time, stats.transitions_time, stats.terminal_check_time, stats.sequential_counting_sort_time, stats.tie_breaking_time, stats.parallel_counting_sort_time_min_max, stats.parallel_counting_sort_time_binning, stats.parallel_counting_sort_time_cutting, stats.states_filtering_time, stats.number_of_surviving_successors, stats.number_of_states_filtered))
                else
                    @sprintf("[CSV] %d;%d;%s;%s;%d;%s;%d;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%d;%d", rank, Threads.nthreads(), instance_description, instance_csv_string(configuration, instance), beam_width, guidance_function, root_estimate, best_objective, construction_time, stats.main_loop_time, stats.memory_init_time, stats.instance_preparation_time, stats.successor_creation_time, stats.surviving_successors_time, stats.transitions_time, stats.terminal_check_time, stats.sequential_counting_sort_time, stats.tie_breaking_time, stats.parallel_counting_sort_time_min_max, stats.parallel_counting_sort_time_binning, stats.parallel_counting_sort_time_cutting, stats.states_filtering_time, stats.number_of_surviving_successors, stats.number_of_states_filtered)
                end
            end
        end
    end

    if rank == 0
        all_objectives = create_all_objectives(size)
        all_solutions = create_all_solutions(instance, size)
        rreqs = Vector{MPI.Request}()
        for src in 1:size
            push!(rreqs, MPI.Irecv!(all_objectives[src], src-1, src-1, comm))
            push!(rreqs, MPI.Irecv!(all_solutions[src], src-1, src-1, comm))
        end
    end

    #min_obj_and_sol = minimum(objectives_and_solutions)
    if configuration.problem_type == mini
        min_obj_and_sol = sort(objectives_and_solutions, by=first)[1]
    elseif configuration.problem_type == maxi
        min_obj_and_sol = sort(objectives_and_solutions, by=first)[end]
    end
    if !pre_run
        @info @sprintf("minimum objective and solution: %s", min_obj_and_sol)
    else
        @sprintf("minimum objective and solution: %s", min_obj_and_sol)
    end
    sreq_obj = MPI.Isend(min_obj_and_sol[1], 0, rank, comm)
    sreq_sol = MPI.Isend(min_obj_and_sol[2], 0, rank, comm)

    MPI.Waitall!([sreq_obj, sreq_sol])
    if rank == 0
        MPI.Waitall!(rreqs)
    end

    MPI.Barrier(comm)

    if rank == 0
        if configuration.problem_type == mini
            best_obj_and_sol_idx = argmin(all_objectives)
        elseif configuration.problem_type == maxi
            best_obj_and_sol_idx = argmax(all_objectives)
        end
        best_objective = all_objectives[best_obj_and_sol_idx][1]
        best_solution = all_solutions[best_obj_and_sol_idx]
        end_time = time()
        runtime_without_startup = end_time - start_time
        if !pre_run
            @info @sprintf("best solution over all runs: %s\n", best_solution)
            println(@sprintf("[CSV2] %s;%s;%d;%s;%s;%s;%s;%d;%f;%d;%d;%f;%f;%d;%d", instance_description, instance_csv_string(configuration, instance), beam_width, guidance_function, tie_breaking, variable_ordering, states_filtering, depth_for_subtree_splitting, sigma_rel, root_estimate, best_objective, stats.instance_preparation_time, runtime_without_startup, size, Threads.nthreads()))
        else
            @sprintf("best solution over all runs: %s\n", best_solution)
            @sprintf("[CSV2] %s;%s;%d;%s;%s;%s;%s;%d;%f;%d;%d;%f;%f;%d;%d", instance_description, instance_csv_string(configuration, instance), beam_width, guidance_function, tie_breaking, variable_ordering, states_filtering, depth_for_subtree_splitting, sigma_rel, root_estimate, best_objective, stats.instance_preparation_time, runtime_without_startup, size, Threads.nthreads())
        end
        #@time print_solution(best_solution)
        print_additional_output(configuration, instance, best_solution)
    end
end

function main()
    @info "Performing pre-run"
    old_logger = global_logger(Logging.NullLogger())
    main_mpi(true)
    global_logger(old_logger)
    @time main_mpi(false)
end