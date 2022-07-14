include("lib/parbeam.jl")

using Printf: @printf, @sprintf

# specialized types
mutable struct AdditionalData <: AbstractAdditionalData
    independent_set::Array{Bool, 2}
    remaining_vertices::Array{Bool, 2}
    remaining_n::Array{Int, 1}
    remaining_m::Array{Int, 1}
end

mutable struct AdditionalSuccessorData <: AbstractAdditionalSuccessorData
end

MISPNodeData = NodeData{Int, Float64, Int}
MISPSuccessorData = SuccessorData{Int, Float64, Bool, AdditionalSuccessorData, Int}
MISPBeam = Beam{MISPNodeData, AdditionalData}

# instance
struct Instance <: AbstractInstance
    n::Int
    m::Int
    adjacency_lists::Array{Vector{Int}, 1}
    perm::Array{Int, 1}
end

function problem_type(instance::Instance)
    maxi
end

function read_in_MISP(guidance_function::String, variable_ordering::String, instance_name::String)
    file_name = string("insts/MISP/", instance_name, ".clq")
    lines = readlines(file_name)

    count_edges = 0
    n = 0
    m = 0
    adjacency_lists = nothing
    for line in lines
        if line[1] == 'c'
            continue
        elseif line[1] == 'p'
            split_p_line = split(line, " ")
            n = parse(Int, split_p_line[3])
            m = parse(Int, split_p_line[4])
            adjacency_lists = Array{Vector{Int}, 1}(undef, n)
            for i in 1:n
                adjacency_lists[i] = Vector{Int}()
            end
        elseif line[1] == 'e'
            split_e_line = split(line, " ")
            v_1 = parse(Int, split_e_line[2])
            v_2 = parse(Int, split_e_line[3])
            @assert v_1 != v_2
            push!(adjacency_lists[v_1], v_2)
            push!(adjacency_lists[v_2], v_1)
            count_edges += 1
        else
            throw("unknown clq file line type")
        end
    end

    @assert count_edges == m

    if variable_ordering == "lex"
        perm = convert(Array{Int, 1}, 1:n)
    else
        throw(@sprintf("variable ordering %s does not exist", variable_ordering))
    end

    Instance(n, m, adjacency_lists, perm)
end

function read_instance(guidance_function::String, variable_ordering::String, instance_description::String)
    read_in_MISP(guidance_function, variable_ordering, instance_description)
end

function instance_csv_string(instance::Instance)
    @sprintf("%d;%d", instance.n, instance.m)
end

# initialization
function initialize_beam(configuration::Configuration, instance::Instance)
    layer = Array{Int, 1}(undef, configuration.beam_width)
    costs_so_far = Array{Int, 1}(undef, configuration.beam_width)
    priority = Array{Float64, 1}(undef, configuration.beam_width)
    nodes = MISPNodeData(layer, costs_so_far, priority)

    independent_set = Array{Bool, 2}(undef, configuration.beam_width, instance.n)
    remaining_vertices = Array{Bool, 2}(undef, configuration.beam_width, instance.n)
    remaining_n = Array{Int, 1}(undef, configuration.beam_width)
    remaining_m = Array{Int, 1}(undef, configuration.beam_width)

    additional_data = AdditionalData(independent_set, remaining_vertices, remaining_n, remaining_m)
    MISPBeam(nodes, additional_data)
end

function initialize_successor_specs(configuration::Configuration, instance::Instance)::MISPSuccessorData
    layer = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    costs_after_move = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    priority = Array{Float64, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    Q_idx = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    successor_move = Array{Bool, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)

    successor_data = MISPSuccessorData(layer, costs_after_move, priority, Q_idx, successor_move, AdditionalSuccessorData())
    successor_data
end

function initialize_root(instance::Instance, Q::MISPBeam)
    Q.nodes.layer[1] = 0
    Q.nodes.costs_so_far[1] = 0
    Q.nodes.priority[1] = 0

    for i in 1:instance.n
        Q.additional_data.remaining_vertices[1,i] = true
        Q.additional_data.independent_set[1,i] = false
    end

    Q.additional_data.remaining_n[1] = instance.n
    Q.additional_data.remaining_m[1] = instance.m
end

function initial_lb(configuration::Configuration, instance::Instance, beam::MISPBeam)
    convert(Int, floor(1//2 + sqrt(1//4 + instance.n^2 - instance.n - 2*instance.m)))
end

function initial_successor_moves(instance::Instance, depth::Int)
    moves_list = Vector{Vector{Int}}()
    push!(moves_list, [])
    if depth > 0
        throw(ArgumentError("subtree splitting not yet implemented for the MISP"))
    end
    moves_list
end

# objectives and solution
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
    all_solutions = Vector{Array{Int}}(undef, size)
    for i in 1:size
        all_solutions[i] = zeros(Int, instance.n)
    end
    all_solutions
end

function get_solution(instance::Instance, beam::MISPBeam, best_terminal::Int)
    best_solution = beam.additional_data.independent_set[best_terminal,:]

    selected_vertices = Vector{Int}()
    for i in 1:instance.n
        if best_solution[i]
            push!(selected_vertices, i)
        end
    end

    selected_vertices
end

# successors and transitions
function maximum_number_of_successors_per_node(instance::Instance)
    2
end

function maximum_layer_number(instance::Instance)
    instance.n
end

function eval_successors(configuration::Configuration, instance::Instance, beam::MISPBeam, successor_specs::MISPSuccessorData, start_index::Int, Q_idx::Int, sigma::Float64)
    successor_count = 0
    layer = beam.nodes.layer[Q_idx]

    if is_terminal(instance, beam, Q_idx)
        return successor_count
    end

    next_vertex = instance.perm[layer+1]

    i = start_index + successor_count
    @assert i <= configuration.maximum_number_of_successors_by_node*configuration.beam_width

    # select vertex if possible
    if beam.additional_data.remaining_vertices[Q_idx,next_vertex]
        eval_successor(configuration, instance, layer, beam, Q_idx, successor_specs, i, true, next_vertex, sigma)
        successor_count += 1
    end
    
    i = start_index + successor_count
    @assert i <= configuration.maximum_number_of_successors_by_node*configuration.beam_width

    # do not select vertex
    eval_successor(configuration, instance, layer, beam, Q_idx, successor_specs, i, false, next_vertex, sigma)
    successor_count += 1

    successor_count
end

function eval_successor(configuration::Configuration, instance::Instance, layer::Int, beam::MISPBeam, Q_idx::Int, successor_specs::MISPSuccessorData, i::Int, choice::Bool, next_vertex::Int, sigma::Float64)
    if choice
        costs_delta = 1
    else
        costs_delta = 0
    end

    successor_specs.costs_after_move[i] = beam.nodes.costs_so_far[Q_idx] + costs_delta
    if configuration.guidance_function == "g"
        successor_specs.priority[i] = -successor_specs.costs_after_move[i]
    elseif configuration.guidance_function == "lex_g_trivial_ub"
        remaining_n = beam.additional_data.remaining_n[Q_idx]
        if beam.additional_data.remaining_vertices[Q_idx,next_vertex]
            remaining_n -= 1
        end
        if choice
            for neighbor in instance.adjacency_lists[next_vertex]
                if beam.additional_data.remaining_vertices[Q_idx,neighbor]
                    remaining_n -= 1
                end
            end
        end
        successor_specs.priority[i] = -(successor_specs.costs_after_move[i] + 10^-6*remaining_n)
    else
        throw(ArgumentError("unknown guidance"))
    end
    successor_specs.layer[i] = layer + 1
    if sigma > 0.0
        successor_specs.priority[i] += randn()*sigma
    end
    successor_specs.Q_idx[i] = Q_idx
    successor_specs.successor_move[i] = choice
end

function is_terminal(instance::Instance, beam::MISPBeam, Q_idx::Int)
    beam.additional_data.remaining_n[Q_idx] == 0
end

function terminal_costs(instance::Instance, beam::MISPBeam, Q_idx::Int)
    0
end

function copy(instance::Instance, from_Q::MISPBeam, to_Q::MISPBeam, from_idx::Int, to_idx::Int)
    to_Q.nodes.layer[to_idx] = from_Q.nodes.layer[from_idx]
    to_Q.nodes.costs_so_far[to_idx] = from_Q.nodes.costs_so_far[from_idx]
    to_Q.nodes.priority[to_idx] = from_Q.nodes.priority[from_idx]
    to_Q.additional_data.remaining_n[to_idx] = from_Q.additional_data.remaining_n[from_idx]
    to_Q.additional_data.remaining_m[to_idx] = from_Q.additional_data.remaining_m[from_idx]
    for i in 1:instance.n
        to_Q.additional_data.independent_set[to_idx,i] = from_Q.additional_data.independent_set[from_idx,i]
        to_Q.additional_data.remaining_vertices[to_idx,i] = from_Q.additional_data.remaining_vertices[from_idx,i]
    end
end

function copy(from_successor::MISPSuccessorData, to_successors::MISPSuccessorData, from_idx::Int, to_idx::Int)
    to_successors.layer[to_idx] = from_successor.layer[from_idx]
    to_successors.costs_after_move[to_idx] = from_successor.costs_after_move[from_idx]
    to_successors.priority[to_idx] = from_successor.priority[from_idx]
    to_successors.Q_idx[to_idx] = from_successor.Q_idx[from_idx]
    to_successors.successor_move[to_idx] = from_successor.successor_move[from_idx]
end

function make_transition(instance::Instance, successor_specs::MISPSuccessorData, Q::MISPBeam, Q_idx::Int, spec_idx::Int)
    Q.nodes.layer[Q_idx] += 1
    Q.nodes.costs_so_far[Q_idx] = successor_specs.costs_after_move[spec_idx]
    Q.nodes.priority[Q_idx] = successor_specs.priority[spec_idx]
    j = successor_specs.successor_move[spec_idx][1]
    make_move(instance, Q, Q_idx, j)
end

function make_move(instance::Instance, Q::MISPBeam, Q_idx::Int, next_j::Bool)
    layer = Q.nodes.layer[Q_idx]
    vertex = instance.perm[layer]
    if Q.additional_data.remaining_vertices[Q_idx,vertex]
        Q.additional_data.remaining_n[Q_idx] -= 1
    end
    Q.additional_data.remaining_vertices[Q_idx,vertex] = false
    Q.additional_data.independent_set[Q_idx,vertex] = next_j
    if next_j
        for neighbor in instance.adjacency_lists[vertex]
            if Q.additional_data.remaining_vertices[Q_idx,neighbor]
                Q.additional_data.remaining_n[Q_idx] -= 1
            end
            Q.additional_data.remaining_vertices[Q_idx,neighbor] = false
        end
    end

    @assert Q.additional_data.remaining_n[Q_idx] >= 0
end

# duplicate states/dominance
function state_hash(configuration::Configuration, instance::Instance, Q::MISPBeam, Q_idx::Int, layer::Number)
    remaining_vertices = @view Q.additional_data.remaining_vertices[Q_idx,:]

    hash(remaining_vertices)
end

function isequal_by_idx(Q_a::MISPBeam, Q_b::MISPBeam, Q_idx_a::Int, Q_idx_b::Int, layer::Number)
    remaining_vertices_a = @view Q_a.additional_data.remaining_vertices[Q_idx_a,:]
    remaining_vertices_b = @view Q_b.additional_data.remaining_vertices[Q_idx_b,:]

    isequal(remaining_vertices_a, remaining_vertices_b)
end

# memory footprint

# pre-run
function pre_run_instance_name(instance_name::String)::String
    "bergman2013/random_graph_200_10_1"
end
   
@time main()