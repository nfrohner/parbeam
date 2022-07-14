# TTP CVRP(H) bounds precalculation module by recursive state space traversal and subsequent bottom up sweep

module TTPCVRP
    using DataStructures: Queue, enqueue!, dequeue!, Deque
    using SharedArrays
    import Base.hash
    import Base.isequal

    using Printf: @printf

    struct State
        teams_left::BitSet
        position::Int
        streak::Int
    end

    struct Arc{T}
        destination::T
        weight::Int
    end

    mutable struct Node
        layer::Int
        shortest_path_length::Int
        lower_bound::Int
        constrained_lower_bounds::Union{Vector{Int}, Nothing}
        parent::Union{Node, Nothing}
        forward_arcs::Vector{Arc{Node}}
        state::Union{State, Nothing}
        Node() = new(0, 0, typemax(UInt32), nothing, nothing, Vector{Arc{Node}}(), nothing)
    end

    function mask_teams_left(team, teams_left)
        set_mask = 0
        for away_team in teams_left
            if away_team < team
                set_mask += 2 ^ (away_team-1)
            else
                set_mask += 2 ^ (away_team-2)
            end
        end
        set_mask+1
    end

    function hash(state::State, h::UInt)
        hash(state.teams_left, hash(state.position, hash(state.streak)))
    end

    function isequal(a::State, b::State)
        isequal(a.teams_left, b.teams_left) && isequal(a.position, b.position) && isequal(a.streak, b.streak)
    end

    function move_to_team(d::Array{Int, 2}, node::Node, to_team::Int)
        new_node = Node()
        weight = d[node.state.position, to_team]
        new_node.shortest_path_length = node.shortest_path_length + weight
        new_node.state = State(delete!(copy(node.state.teams_left), to_team), to_team, node.state.streak+1)
        new_node, weight
    end

    function move_to_team_and_home(d::Array{Int, 2}, node::Node, to_team::Int, home::Int)
        new_node = Node()
        weight = d[node.state.position, to_team] + d[to_team, home]
        new_node.shortest_path_length = node.shortest_path_length + weight
        new_node.state = State(delete!(copy(node.state.teams_left), to_team), home, 0)
        new_node, weight
    end

    function incorporate(parent::Node, new_node::Node, weight::Int, nodes::Dict{State, Node}, nodes_by_layers::Dict{Int, Deque{Node}}, Q::Queue{Node})
        new_node.layer = parent.layer+1
        if haskey(nodes, new_node.state)
            existing_node = nodes[new_node.state]
            if new_node.shortest_path_length < existing_node.shortest_path_length
                existing_node.shortest_path_length = new_node.shortest_path_length
                existing_node.parent = parent
            end
            push!(parent.forward_arcs, Arc{Node}(existing_node, weight))
        else
            nodes[new_node.state] = new_node
            push!(nodes_by_layers[new_node.layer], new_node)
            push!(parent.forward_arcs, Arc{Node}(new_node, weight))
            enqueue!(Q, new_node)
            #push!(Q, new_node)
        end
    end

    function construct(team::Int, n::Int, d::Array{Int, 2}, streak_limit::Int)
        root = Node()
        root.shortest_path_length = 0
        root.state = State(delete!(BitSet(1:n), team), team, 0)

        terminal = Node()
        terminal.shortest_path_length = typemax(UInt32)
        terminal.lower_bound = 0
        terminal.constrained_lower_bounds = ones(Int, n)*typemax(UInt32)
        terminal.constrained_lower_bounds[1] = 0
        terminal.state = State(BitSet(), team, 0)

        Q = Queue{Node}()
        enqueue!(Q, root)
        #push!(Q, root)

        nodes = Dict{State, Node}()
        nodes[root.state] = root
        nodes[terminal.state] = terminal

        nodes_by_layers = Dict{Int, Deque{Node}}(i => Deque{Node}() for i = 0:n-1)
        push!(nodes_by_layers[0], root)
        push!(nodes_by_layers[n-1], terminal)

        transitions = 0

        while length(Q) > 0
            node = dequeue!(Q)
            #node = pop!(Q)
            for to_team in node.state.teams_left
                if length(node.state.teams_left) > 1 && node.state.streak < streak_limit - 1
                    new_node, weight = move_to_team(d, node, to_team)
                    transitions += 1
                    incorporate(node, new_node, weight, nodes, nodes_by_layers, Q)
                end

                new_node, weight = move_to_team_and_home(d, node, to_team, team)
                transitions += 1
                incorporate(node, new_node, weight, nodes, nodes_by_layers, Q)
            end
        end

        @printf("%d transitions\n", transitions)
        @printf("%d nodes\n", length(nodes))

        nodes, nodes_by_layers, terminal.shortest_path_length
    end

    function calculate_bounds_for_team(n::Int, team::Int, d::Array{Int, 2}, streak_limit::Int, bounds_by_state::SharedArray{UInt32,4})
        @printf("calculating team %d\n", team)
        @time nodes, nodes_by_layers, shortest_path = construct(team, n, d, streak_limit)

        @time for i = n-2:-1:0
            for node in nodes_by_layers[i]
                node.lower_bound = minimum(map(x -> x.destination.lower_bound + x.weight, node.forward_arcs))
            end
        end

        @time for (node_state, node) in nodes
            bounds_by_state[team, mask_teams_left(team, node.state.teams_left), node.state.position, node.state.streak+1] = node.lower_bound
        end

        flush(stdout)
    end

    function calculate_bounds_for_team_cvrp_sparse(n::Int, team::Int, d::Array{Int, 2}, streak_limit::Int)
        @time nodes, nodes_by_layers, shortest_path = construct(team, n, d, streak_limit)

        @time for i = n-2:-1:0
            for node in nodes_by_layers[i]
                node.lower_bound = minimum(map(x -> x.destination.lower_bound + x.weight, node.forward_arcs))
            end
        end

        sparse_matrix_update = Vector{Tuple{NTuple{4, Int}, Int}}()
        sizehint!(sparse_matrix_update, length(nodes))
        @time for (node_state, node) in nodes
            push!(sparse_matrix_update, ((team, mask_teams_left(team, node.state.teams_left), node.state.position, node.state.streak+1), node.lower_bound))
        end

        flush(stdout)

        sparse_matrix_update
    end

    function sum_with_potential_infinity(a::Int, b::Int)::UInt32
        if a == typemax(UInt32) || b == typemax(UInt32)
            typemax(UInt32)
        else
            a + b
        end
    end

    function calculate_bounds_for_team(n::Int, team::Int, d::Array{Int, 2}, streak_limit::Int, bounds_by_state::SharedArray{UInt32,5})
        @printf("calculating team %d\n", team)
        @time nodes, nodes_by_layers, shortest_path = construct(team, n, d, streak_limit)

        @time for i = n-2:-1:0
            for node in nodes_by_layers[i]
                node.constrained_lower_bounds = ones(Int, n)*typemax(UInt32)
                if node.state.position == team
                    #node.constrained_lower_bounds[2:end] = min.(node.constrained_lower_bounds[2:end], map(x -> [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n-1], node.forward_arcs)...)
                    for x in node.forward_arcs
                        node.constrained_lower_bounds[2:end] = min.(node.constrained_lower_bounds[2:end], [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n-1])
                    end
                else
                    #node.constrained_lower_bounds = min.(node.constrained_lower_bounds, map(x -> [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n], node.forward_arcs)...)
                    for x in node.forward_arcs
                        node.constrained_lower_bounds = min.(node.constrained_lower_bounds, [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n])
                    end
                end
            end
        end

        @time for (node_state, node) in nodes
            bounds_by_state[team, mask_teams_left(team, node.state.teams_left), node.state.position, node.state.streak+1, :] = node.constrained_lower_bounds
        end

        flush(stdout)
    end

    function calculate_bounds_for_team_sparse(n::Int, team::Int, d::Array{Int, 2}, streak_limit::Int)
        @printf("calculating team %d\n", team)
        @time nodes, nodes_by_layers, shortest_path = construct(team, n, d, streak_limit)

        @time for i = n-2:-1:0
            for node in nodes_by_layers[i]
                node.constrained_lower_bounds = ones(Int, n)*typemax(UInt32)
                if node.state.position == team
                    #node.constrained_lower_bounds[2:end] = min.(node.constrained_lower_bounds[2:end], map(x -> [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n-1], node.forward_arcs)...)
                    for x in node.forward_arcs
                        node.constrained_lower_bounds[2:end] = min.(node.constrained_lower_bounds[2:end], [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n-1])
                    end
                else
                    #node.constrained_lower_bounds = min.(node.constrained_lower_bounds, map(x -> [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n], node.forward_arcs)...)
                    for x in node.forward_arcs
                        node.constrained_lower_bounds = min.(node.constrained_lower_bounds, [sum_with_potential_infinity(x.destination.constrained_lower_bounds[i], x.weight) for i in 1:n])
                    end
                end
            end
        end

        sparse_matrix_update = Vector{Tuple{NTuple{4, Int}, UInt32}}()
        sizehint!(sparse_matrix_update, length(nodes))
        @time for (node_state, node) in nodes
            push!(sparse_matrix_update, ((team, mask_teams_left(team, node.state.teams_left), node.state.position, node.state.streak+1), node.lower_bound))
        end

        flush(stdout)

        sparse_matrix_update
    end
end
