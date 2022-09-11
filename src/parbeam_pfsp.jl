import Base.copy

include("lib/parbeam.jl")

using DelimitedFiles
using Random
using Printf: @printf

mutable struct AdditionalData <: AbstractAdditionalData
    perms::Array{Int, 2}
    idle_times::Array{Int, 2}
    time_markers::Array{Int, 2}
    remaining_processing_times::Array{Int, 2}
end

mutable struct AdditionalSuccessorData <: AbstractAdditionalSuccessorData
end

PFSPNodeData = NodeData{Int, Float64, Int}
PFSPSuccessorData = SuccessorData{Int, Float64, Int, AdditionalSuccessorData, Int}
PFSPBeam = Beam{PFSPNodeData, AdditionalData}

struct AuxiliaryData
    remaining_jobs::Array{Int, 1}
    current_johnson_schedules_with_machine_m::Array{Int, 2}
end

struct Instance <: AbstractInstance
    n::Int
    m::Int
    PTM::Array{Int, 2}
    cPTM::Array{Int, 2}
    rcPTM::Array{Int, 2}
    time_lags::Array{Int, 3}
    time_lags_with_last_machine::Array{Int, 2}
    johnson_schedules::Array{Int, 3}
    johnson_schedules_with_last_machine::Array{Int, 2}
    time_markers_by_thread::Array{Int, 2}
    remaining_jobs_by_thread::Array{Int, 2}
    aux_by_thread::Vector{AuxiliaryData}
end

function root(instance::Instance)
    layer = 0
    costs_so_far = 0
    priority = 0
    PFSPNode(layer, costs_so_far, priority)
end

function initial_estimate(configuration::Configuration, instance::Instance, beam::PFSPBeam)
    #lb_2(instance, beam, 1, beam.additional_data)
    if configuration.guidance_function == "flow_time_bound_and_idle_time"
        0
    else
        lb_taillard(instance, beam, 1, beam.additional_data)
    end
end

function is_terminal(instance::Instance, beam::PFSPBeam, Q_idx::Int)
    beam.nodes.layer[Q_idx] == instance.n
end

function terminal_costs(instance::Instance, beam::PFSPBeam, Q_idx::Int)
    0
end

function maximum_number_of_successors_per_node(instance::Instance)
    instance.n
end

function maximum_layer_number(instance::Instance)
    instance.n
end

function initialize_beam(configuration::Configuration, instance::Instance)
    #@assert isbitstype(PFSPNode)
    layer = Array{Int, 1}(undef, configuration.beam_width)
    costs_so_far = Array{Int, 1}(undef, configuration.beam_width)
    priority = Array{Float64, 1}(undef, configuration.beam_width)
    nodes = PFSPNodeData(layer, costs_so_far, priority)

    perms = Array{Int, 2}(undef, configuration.beam_width, instance.n)
    time_markers = Array{Int, 2}(undef, configuration.beam_width, instance.m)
    idle_times = Array{Int, 2}(undef, configuration.beam_width, instance.m)
    remaining_processing_times = Array{Int, 2}(undef, configuration.beam_width, instance.m)
    additional_data = AdditionalData(perms, time_markers, idle_times, remaining_processing_times)

    PFSPBeam(nodes, additional_data)
end

function initialize_root(instance::Instance, Q::PFSPBeam)
    Q.nodes.layer[1] = 0
    Q.nodes.costs_so_far[1] = 0
    Q.nodes.priority[1] = 0

    for i in 1:instance.n
        Q.additional_data.perms[1,i] = i
    end

    for i in 1:instance.m
        Q.additional_data.time_markers[1,i] = 0
        Q.additional_data.idle_times[1,i] = 0
        Q.additional_data.remaining_processing_times[1,i] = sum(@view instance.PTM[i,:])
    end
end

function initialize_successor_specs(configuration::Configuration, instance::Instance)::PFSPSuccessorData
    #@assert isbitstype(PFSPSuccessorSpecification)
    layer = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    costs_after_move = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    priority = Array{Float64, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    Q_idx = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    successor_move = Array{Int, 1}(undef, configuration.beam_width*configuration.maximum_number_of_successors_by_node)
    successor_data = PFSPSuccessorData(layer, costs_after_move, priority, Q_idx, successor_move, AdditionalSuccessorData())
    successor_data
end

function eval_successors(configuration::Configuration, instance::Instance, beam::PFSPBeam, successor_specs::PFSPSuccessorData, start_index::Int, Q_idx::Int, sigma::Float64)
    successor_count = 0
    layer = beam.nodes.layer[Q_idx]
    
    if configuration.guidance_function == "LB2_and_idle_time" || configuration.guidance_function == "LB2"
        #remaining_jobs = zeros(Int, instance.n)
        #for j_idx in (layer+1):instance.n
        #    j = beam.additional_data.perms[Q_idx,j_idx]
        #    remaining_jobs[j] = 1
        #end
        for j_idx in 1:layer
            j = beam.additional_data.perms[Q_idx,j_idx]
            instance.aux_by_thread[Threads.threadid()].remaining_jobs[j] = 0
        end
        for j_idx in (layer+1):instance.n
            j = beam.additional_data.perms[Q_idx,j_idx]
            instance.aux_by_thread[Threads.threadid()].remaining_jobs[j] = 1
        end

        #current_johnson_schedules_with_machine_m = Array{Int, 2}(undef, (instance.m-1, instance.n - layer))
        for m_1 in 1:(instance.m - 1)
            dst = 1
            for j_idx in 1:instance.n
                j = instance.johnson_schedules_with_last_machine[m_1,j_idx]
                @assert j >= 1 && j <= instance.n
                if instance.aux_by_thread[Threads.threadid()].remaining_jobs[j] == 1
                    instance.aux_by_thread[Threads.threadid()].current_johnson_schedules_with_machine_m[m_1,dst] = j
                    dst += 1
                end
            end
            @assert dst == instance.n - layer + 1
        end
    end

    for j_idx in (layer+1):instance.n
        j = beam.additional_data.perms[Q_idx, j_idx]
        @assert j >= 1 && j <= instance.n
        i = start_index + successor_count
        #successor_spec = successor_specs[start_index+j-1]
        successor_specs.layer[i] = layer + 1
        if configuration.guidance_function == "flow_time_bound_and_idle_time"
            successor_specs.costs_after_move[i] = flow_time_after_next_job(instance, beam, Q_idx, j)
        else
            successor_specs.costs_after_move[i] = makespan_after_next_job(instance, beam, Q_idx, j)
        end
        if configuration.guidance_function == "cmax"
            successor_specs.priority[i] = successor_specs.costs_after_move[i]
        #elseif configuration.guidance_function == "bound"
        #    successor_specs.priority[i] = incremental_bound(instance, node, j, successor_spec.costs_after_move)
        elseif configuration.guidance_function == "LB2"
            #successor_specs.priority[i] = lb_2(instance, beam, Q_idx, beam.additional_data, j)
            lb_2_bound, _ = incremental_lb_2_and_idle_time(instance, beam, Q_idx, j, layer)
            successor_specs.priority[i] = lb_2_bound
        #elseif configuration.guidance_function == "idle_time"
        #    successor_specs.priority[i] = incremental_idle_time(instance, node, j)
        elseif configuration.guidance_function == "bound_and_idle_time"
            successor_specs.priority[i] = incremental_bound_and_idle_time(instance, beam, Q_idx, j)
        elseif configuration.guidance_function == "LB2_and_idle_time"
            alpha = (beam.nodes.layer[Q_idx] + 1) / instance.n
            C = layer / instance.m
            lb_2_bound, total_idle_time = incremental_lb_2_and_idle_time(instance, beam, Q_idx, j, layer)
            successor_specs.priority[i] = alpha * lb_2_bound + (1 - alpha) * C * total_idle_time
        #elseif configuration.guidance_function == "lb_taillard"
        #    successor_specs.priority[i] = incremental_lb_taillard(instance, node, j)
        elseif configuration.guidance_function == "flow_time_bound_and_idle_time"
            successor_specs.priority[i] = incremental_flow_time_bound_and_idle_time(instance, beam, Q_idx, j)
        else
            throw(ArgumentError("unknown guidance"))
        end
        if sigma > 0.0
            successor_specs.priority[i] += randn()*sigma
        end
        successor_specs.Q_idx[i] = Q_idx
        successor_specs.successor_move[i] = j
        successor_count += 1
    end
    successor_count
end

function makespan_after_next_job(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int)
    last_marker = beam.additional_data.time_markers[Q_idx, 1] + instance.PTM[1, next_job]
    for machine in 2:instance.m
        last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx, machine]) + instance.PTM[machine, next_job]
    end
    last_marker
end

function flow_time_after_next_job(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int) 
    last_marker = beam.additional_data.time_markers[Q_idx, 1] + instance.PTM[1, next_job]
    for machine in 2:instance.m
        last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx, machine]) + instance.PTM[machine, next_job]
    end
    beam.nodes.costs_so_far[Q_idx] + last_marker
end

function incremental_bound_and_idle_time(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int)
    layer = beam.nodes.layer[Q_idx] + 1
    alpha = layer / instance.n
    @assert alpha >= 0
    @assert alpha <= 1
    #C = alpha * instance.n / instance.m
    #C = 1
    C = layer / instance.m
    total_idle_time = beam.additional_data.idle_times[Q_idx,1]
    last_marker = beam.additional_data.time_markers[Q_idx,1] + instance.PTM[1, next_job]
    max_marker = last_marker + incremental_remaining_jobs_processing_time_on_machine(instance, beam, Q_idx, next_job, 1)
    for machine in 2:instance.m
        total_idle_time += beam.additional_data.idle_times[Q_idx,machine] + max(last_marker - beam.additional_data.time_markers[Q_idx,machine], 0)
        last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx,machine]) + instance.PTM[machine, next_job]
        max_marker = max(max_marker, last_marker + incremental_remaining_jobs_processing_time_on_machine(instance, beam, Q_idx, next_job, machine))
    end
    alpha * max_marker + (1 - alpha) * C * total_idle_time
end

function incremental_flow_time_bound_and_idle_time(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int)
    layer = beam.nodes.layer[Q_idx] + 1
    alpha = layer / instance.n
    C = layer / instance.m
    total_idle_time = beam.additional_data.idle_times[Q_idx,1]
    last_marker = beam.additional_data.time_markers[Q_idx,1] + instance.PTM[1, next_job]
    max_marker = last_marker + incremental_remaining_jobs_processing_time_on_machine(instance, beam, Q_idx, next_job, 1)
    for machine in 2:instance.m
        total_idle_time += beam.additional_data.idle_times[Q_idx,machine] + max(last_marker - beam.additional_data.time_markers[Q_idx,machine], 0)
        last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx,machine]) + instance.PTM[machine, next_job]
        max_marker = max(max_marker, last_marker + incremental_remaining_jobs_processing_time_on_machine(instance, beam, Q_idx, next_job, machine))
    end
    alpha * (beam.nodes.costs_so_far[Q_idx] + last_marker) + (1 - alpha) * C * total_idle_time
end

function incremental_lb_2_and_idle_time(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int, layer::Int)
    #alpha = (layer + 1) / instance.n
    lb_2_bound = 0
    #@assert alpha >= 0
    #@assert alpha <= 1
    
    total_idle_time = beam.additional_data.idle_times[Q_idx,1]
    last_marker = beam.additional_data.time_markers[Q_idx,1] + instance.PTM[1, next_job]
    lb_2_bound = max(lb_2_bound, last_marker + incremental_johnson_makespan(instance, 1, next_job, layer))
    for machine in 2:(instance.m-1)
        total_idle_time += beam.additional_data.idle_times[Q_idx,machine] + max(last_marker - beam.additional_data.time_markers[Q_idx,machine], 0)
        last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx,machine]) + instance.PTM[machine, next_job]
        lb_2_bound = max(lb_2_bound, last_marker + incremental_johnson_makespan(instance, machine, next_job, layer))
    end
    total_idle_time += beam.additional_data.idle_times[Q_idx,instance.m] + max(last_marker - beam.additional_data.time_markers[Q_idx,instance.m], 0)
    last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx,instance.m]) + instance.PTM[instance.m, next_job]
    #println(lb_2_bound)
    #println(lb_2(instance, beam, Q_idx, beam.additional_data, next_job))
    #@assert lb_2_bound == lb_2(instance, beam, Q_idx, beam.additional_data, next_job)
    #alpha * lb_2_bound + (1 - alpha) * total_idle_time
    lb_2_bound, total_idle_time
end

function incremental_johnson_makespan(instance::Instance, m_1::Int, next_job::Int, layer::Int)
    last_marker_m_1 = 0
    last_marker_m_2 = 0
    for j_idx in 1:(instance.n - layer)
        j = instance.aux_by_thread[Threads.threadid()].current_johnson_schedules_with_machine_m[m_1,j_idx]
        if j != next_job
            last_marker_m_1 += instance.PTM[m_1,j]
            last_marker_m_2 = max(last_marker_m_1 + instance.time_lags_with_last_machine[m_1,j], last_marker_m_2) + instance.PTM[instance.m,j]
        end
    end
    last_marker_m_2
end

function lb_2_and_idle_time(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int)
    alpha = (beam.nodes.layer[Q_idx] + 1) / instance.n
    @assert alpha >= 0
    @assert alpha <= 1
    
    total_idle_time = beam.additional_data.idle_times[Q_idx,1]
    last_marker = beam.additional_data.time_markers[Q_idx,1] + instance.PTM[1, next_job]
    for machine in 2:instance.m
        total_idle_time += max(last_marker - beam.additional_data.time_markers[Q_idx,machine], 0)
        last_marker = max(last_marker, beam.additional_data.time_markers[Q_idx,machine]) + instance.PTM[machine, next_job]
    end
    lb_2_bound = lb_2(instance, beam, Q_idx, beam.additional_data, next_job)
    alpha * lb_2_bound + (1 - alpha) * total_idle_time
end

function incremental_remaining_jobs_processing_time_on_machine(instance::Instance, beam::PFSPBeam, Q_idx::Int, next_job::Int, machine::Int)
    beam.additional_data.remaining_processing_times[Q_idx,machine] - instance.PTM[machine,next_job]
end

function lb_taillard(instance::Instance, beam::Beam, node_idx::Int, additional_data::AdditionalData)
    remaining_jobs = Vector{Int}()
    for j_idx in (beam.nodes.layer[node_idx]+1):instance.n
        j = additional_data.perms[node_idx,j_idx]
        push!(remaining_jobs, j)
    end   

    println(remaining_jobs)

    if length(remaining_jobs) == 0
        return beam.nodes.costs_so_far[node_idx]
    end

    b_i = zeros(Int, instance.m)
    a_i = zeros(Int, instance.m)
    T_i = sum(instance.PTM[:,remaining_jobs], dims=2)

    b_i[2:instance.m] .= minimum(instance.cPTM[:,remaining_jobs], dims=2)[1:instance.m-1]
    a_i[1:instance.m-1] .= reverse(minimum(instance.rcPTM[:,remaining_jobs], dims=2)[1:instance.m-1])

    println(additional_data.time_markers[node_idx,:])
    max(maximum(additional_data.time_markers[node_idx,:] + b_i + T_i + a_i), additional_data.time_markers[node_idx,1] + maximum(sum(instance.PTM[remaining_jobs], dims=1)))
end

function lb_1(instance::Instance, beam::Beam, node_idx::Int, additional_data::AdditionalData)
    if beam.nodes.layer[node_idx] == 0
        bounds_by_machine = zeros(Int, instance.m)
        bounds_by_machine[1] = sum(instance.PTM[1,:])
        for machine in 2:instance.m
            bounds_by_machine[machine] = minimum(instance.PTM[machine-1,:]) + sum(instance.PTM[machine,:])
        end
        maximum(bounds_by_machine)
    end
end

function flow_time_lb_1(instance::Instance, beam::Beam, node_idx::Int, additional_data::AdditionalData)
    if beam.nodes.layer[node_idx] == 0
        bounds_by_machine = zeros(Int, instance.m)
        bounds_by_machine[1] = sum(instance.PTM[1,:])
        for machine in 2:instance.m
            bounds_by_machine[machine] = minimum(instance.PTM[machine-1,:]) + sum(instance.PTM[machine,:])
        end
        maximum(bounds_by_machine)
    end
end

function lb_2(instance::Instance, beam::Beam, node_idx::Int, additional_data::AdditionalData, next_job::Int=0) 
    max_bound = 0
    for j in 1:instance.m
        instance.time_markers_by_thread[Threads.threadid(), j] = additional_data.time_markers[node_idx,j]
    end
    if next_job != 0
        instance.time_markers_by_thread[Threads.threadid(), 1] += instance.PTM[1,next_job]
        for m in 2:instance.m
            instance.time_markers_by_thread[Threads.threadid(), m] = max(instance.time_markers_by_thread[Threads.threadid(), m-1], instance.time_markers_by_thread[Threads.threadid(), m]) + instance.PTM[m,next_job]
        end
    end

    for j_idx in 1:beam.nodes.layer[node_idx]
        j = additional_data.perms[node_idx,j_idx]
        instance.remaining_jobs_by_thread[Threads.threadid(), j] = 0
    end

    for j_idx in (beam.nodes.layer[node_idx]+1):instance.n
        j = additional_data.perms[node_idx,j_idx]
        if j != next_job
            instance.remaining_jobs_by_thread[Threads.threadid(), j] = 1
        else
            instance.remaining_jobs_by_thread[Threads.threadid(), j] = 0
        end
    end

    for m_1 in 1:(instance.m-1)
        #for m_2 in (m_1 + 1):instance.m
        for m_2 in instance.m:instance.m
            last_marker_m_1 = 0
            last_marker_m_2 = 0
            for j_idx in 1:instance.n
                #j = instance.johnson_schedules[m_1,m_2,j_idx]
                j = instance.johnson_schedules_with_last_machine[m_1,j_idx]
                @assert j >= 1 && j <= instance.n
                if instance.remaining_jobs_by_thread[Threads.threadid(), j] == 1
                    last_marker_m_1 += instance.PTM[m_1,j]
                    #last_marker_m_2 = max(last_marker_m_1 + instance.time_lags[m_1,m_2,j], last_marker_m_2) + instance.PTM[m_2,j]
                    last_marker_m_2 = max(last_marker_m_1 + instance.time_lags_with_last_machine[m_1,j], last_marker_m_2) + instance.PTM[m_2,j]
                end
            end
            
            if beam.nodes.layer[node_idx] == 0 && next_job == 0 && m_1 > 1
                r_by_machine = minimum(instance.PTM[m_1-1,:])
            else
                r_by_machine = instance.time_markers_by_thread[Threads.threadid(), m_1]
            end

            if r_by_machine + last_marker_m_2 > max_bound
                max_bound = last_marker_m_2 + r_by_machine
            end
        end
    end

    max_bound
end

function time_markers_idle_times_and_remaining_processing_times_after_next_job!(instance::Instance, Q::PFSPBeam, idx::Int, next_job::Int)
    Q.additional_data.time_markers[idx,1] = Q.additional_data.time_markers[idx,1] + instance.PTM[1, next_job]
    Q.additional_data.remaining_processing_times[idx,1] -= instance.PTM[1,next_job]
    for machine in 2:instance.m
        Q.additional_data.idle_times[idx,machine] += max(Q.additional_data.time_markers[idx,machine-1] - Q.additional_data.time_markers[idx,machine], 0)
        Q.additional_data.time_markers[idx,machine] = max(Q.additional_data.time_markers[idx,machine-1], Q.additional_data.time_markers[idx,machine]) + instance.PTM[machine, next_job]
        Q.additional_data.remaining_processing_times[idx,machine] -= instance.PTM[machine,next_job]
    end
end

function copy(instance::Instance, from_Q::PFSPBeam, to_Q::PFSPBeam, from_idx::Int, to_idx::Int)
    to_Q.nodes.layer[to_idx] = from_Q.nodes.layer[from_idx]
    to_Q.nodes.costs_so_far[to_idx] = from_Q.nodes.costs_so_far[from_idx]
    to_Q.nodes.priority[to_idx] = from_Q.nodes.priority[from_idx]
    for i in 1:instance.n
        to_Q.additional_data.perms[to_idx,i] = from_Q.additional_data.perms[from_idx,i]
    end
    for i in 1:instance.m
        to_Q.additional_data.time_markers[to_idx,i] = from_Q.additional_data.time_markers[from_idx,i]
        to_Q.additional_data.idle_times[to_idx,i] = from_Q.additional_data.idle_times[from_idx,i]
        to_Q.additional_data.remaining_processing_times[to_idx,i] = from_Q.additional_data.remaining_processing_times[from_idx,i]
    end
end

function copy(from_successor::PFSPSuccessorData, to_successors::PFSPSuccessorData, from_idx::Int, to_idx::Int)
    to_successors.layer[to_idx] = from_successor.layer[from_idx]
    to_successors.costs_after_move[to_idx] = from_successor.costs_after_move[from_idx]
    to_successors.priority[to_idx] = from_successor.priority[from_idx]
    to_successors.Q_idx[to_idx] = from_successor.Q_idx[from_idx]
    to_successors.successor_move[to_idx] = from_successor.successor_move[from_idx]
end

function make_transition(configuration::Configuration, instance::Instance, successor_specs::PFSPSuccessorData, Q::PFSPBeam, Q_idx::Int, spec_idx::Int)
    Q.nodes.layer[Q_idx] += 1
    Q.nodes.costs_so_far[Q_idx] = successor_specs.costs_after_move[spec_idx]
    Q.nodes.priority[Q_idx] = successor_specs.priority[spec_idx]
    next_job = successor_specs.successor_move[spec_idx]
    make_move(instance, Q, Q_idx, next_job)
end

function make_move(instance::Instance, Q::PFSPBeam, Q_idx::Int, next_job::Int)
    @assert next_job >= 1 && next_job <= instance.n
    time_markers_idle_times_and_remaining_processing_times_after_next_job!(instance, Q, Q_idx, next_job)
    #@time time_markers_and_idle_times_after_next_job!(instance, Q, idx, next_job)
    #node.remaining_jobs[next_job] = false
    #push!(node.perm, next_job)
    layer = Q.nodes.layer[Q_idx]
    j_to_swap = Q.additional_data.perms[Q_idx,layer]
    if next_job != j_to_swap
        target_idx = -1
        for j_idx in (layer+1):instance.n
            if Q.additional_data.perms[Q_idx,j_idx] == next_job
                target_idx = j_idx
                break 
            end
        end
        @assert target_idx != -1
        Q.additional_data.perms[Q_idx,layer], Q.additional_data.perms[Q_idx,target_idx] = Q.additional_data.perms[Q_idx,target_idx], Q.additional_data.perms[Q_idx,layer]
    end
    #@assert node.priority == lb_taillard(instance, node)
end

function initial_move_down(configuration::Configuration, instance::Instance, Q::PFSPBeam, successor_moves::Vector{Int})
    Q_idx = 1
    #Q.nodes.priority[Q_idx] = lb_2(instance, Q, Q_idx, Q.additional_data)
    Q.nodes.priority[Q_idx] = initial_estimate(configuration, instance, Q)
    for j in successor_moves
        if configuration.guidance_function == "flow_time_bound_and_idle_time"
            Q.nodes.costs_so_far[Q_idx] = flow_time_after_next_job(instance, Q, Q_idx, j)
        else
            Q.nodes.costs_so_far[Q_idx] = makespan_after_next_job(instance, Q, Q_idx, j)
        end
        if configuration.guidance_function == "cmax"
            Q.nodes.priority[Q_idx] = Q.nodes.costs_so_far[Q_idx]
        elseif configuration.guidance_function == "bound_and_idle_time"
            Q.nodes.priority[Q_idx] = incremental_bound_and_idle_time(instance, Q, Q_idx, j)
        elseif configuration.guidance_function == "LB2_and_idle_time"
            Q.nodes.priority[Q_idx] = lb_2_and_idle_time(instance, Q, Q_idx, j)
        elseif configuration.guidance_function == "LB2"
            Q.nodes.priority[Q_idx] = lb_2(instance, Q, Q_idx, Q.additional_data, j)
        elseif configuration.guidance_function == "flow_time_bound_and_idle_time"
            Q.nodes.priority[Q_idx] = incremental_flow_time_bound_and_idle_time(instance, Q, Q_idx, j)
        else
            throw(ArgumentError("unknown guidance"))
        end
        Q.nodes.layer[Q_idx] += 1
        make_move(instance, Q, Q_idx, j)
    end
    Q.nodes.priority[Q_idx]
end

function initial_successor_moves(instance::Instance, depth::Int)
    moves_list = Vector{Vector{Int}}()
    moves = Vector{Int}(undef, 0)
    initial_successor_moves_enumerate(instance, moves_list, moves, depth)
    moves_list
end

function initial_successor_moves_enumerate(instance::Instance, moves_list::Vector{Vector{Int}}, moves::Vector{Int}, depth::Int)
    if length(moves) == depth
        push!(moves_list, moves)
    else
        for j in 1:instance.n
            if !(j in moves)
                new_moves = copy(moves)
                push!(new_moves, j)
                initial_successor_moves_enumerate(instance, moves_list, new_moves, depth)
            end
        end
    end
end

#function shrink_successor_specs(configuration::Configuration, instance::Instance, layer::Int, successor_specs)
#end

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

function get_solution(instance::Instance, beam::PFSPBeam, best_terminal::Int)
    beam.additional_data.perms[best_terminal,:]
end

#function print_solution(solution::Array{Int})
#    @printf("best solution as permutation of jobs: %s\n", solution)
#end

function read_in_FSP(guidance_function::String, instance_name::String)
    VFR = false
    if occursin(r"^ta", instance_name)
        filename = string("insts/PFSP/taillard_instances/", instance_name, ".dat")
    elseif occursin(r"^VFR", instance_name)
        filename = string("insts/PFSP/VFR/", instance_name, "_Gap.txt")
        VFR = true
    else
        throw(ArgumentError("unknown instance class"))
    end
    file = open(filename)

    str = readline(file)
    n = parse(Int, split(str)[1])
    m = parse(Int, split(str)[2])
    PTM = readdlm(file, Int)
    if VFR
        PTM = permutedims(PTM)
        PTM = PTM[2*(convert(Array{Int}, 1:m)), :]
    else
        # dummy execution
        PTM_dummy = permutedims(PTM)
        PTM_dummy = PTM_dummy[2*(convert(Array{Int}, 1:m)), :]
    end
    cPTM = cumsum(PTM, dims=1)
    rcPTM = cumsum(reverse(PTM, dims=1), dims=1)

    # precalculate time lags
    time_lags = zeros(Int, (m, m, n))
    time_lags_with_last_machine = zeros(Int, (m-1, n))
    if guidance_function == "LB2_and_idle_time" || guidance_function == "LB2"
        for m_1 in 1:(m - 1)
            for m_2 in (m_1 + 1):m
                for j in 1:n
                    for machine_between in (m_1 + 1):(m_2 - 1)
                        time_lags[m_1,m_2,j] += PTM[machine_between,j]
                        if m_2 == m
                            time_lags_with_last_machine[m_1,j] += PTM[machine_between,j]
                        end
                    end
                end
            end
        end
    end

    # precalculate johnson schedules
    johnson_schedules = zeros(Int, (m, m, n))
    johnson_schedules_with_last_machine = zeros(Int, (m-1, n))
    remaining_jobs = convert(Vector{Int}, 1:n)
    if guidance_function == "LB2_and_idle_time" || guidance_function == "LB2"
        for m_1 in 1:(m-1)
            for m_2 in (m_1 + 1):m
                job_processing_times_with_machines = vcat(map(x -> (x[2] + time_lags[m_1,m_2,remaining_jobs[x[1]]], remaining_jobs[x[1]], 0), enumerate(PTM[m_1,remaining_jobs])), map(x -> (x[2] + time_lags[m_1,m_2,remaining_jobs[x[1]]], remaining_jobs[x[1]], 1), enumerate(PTM[m_2,remaining_jobs])))
                #println(job_processing_times_with_machines)
                sort!(job_processing_times_with_machines)
                i_front = 1
                i_back = length(remaining_jobs)
                scheduled_jobs = Set{Int}()
                #print(job_processing_times_with_machines)
                for (_, j, m_val) in job_processing_times_with_machines
                    if j in scheduled_jobs
                        continue
                    end
                    
                    push!(scheduled_jobs, j)
                    if m_val == 0
                        johnson_schedules[m_1,m_2,i_front] = j
                        if m_2 == m
                            johnson_schedules_with_last_machine[m_1,i_front] = j
                        end
                        i_front += 1
                    else
                        johnson_schedules[m_1,m_2,i_back] = j
                        if m_2 == m
                            johnson_schedules_with_last_machine[m_1,i_back] = j
                        end
                        i_back -= 1
                    end
                end
                @assert i_front == i_back + 1
            end
        end
    end

    time_markers_by_thread = zeros(Int, (Threads.nthreads(), m))
    remaining_jobs_by_thread = zeros(Int, (Threads.nthreads(), n))

    aux_by_thread = Vector{AuxiliaryData}()
    for _ in 1:Threads.nthreads()
        push!(aux_by_thread, AuxiliaryData(zeros(Int, n), zeros(Int, (m-1, n))))
    end

    Instance(n, m, PTM, cPTM, rcPTM, time_lags, time_lags_with_last_machine, johnson_schedules, johnson_schedules_with_last_machine, time_markers_by_thread, remaining_jobs_by_thread, aux_by_thread)
end

function read_instance(guidance_function::String, variable_ordering::String, instance_description::String)::Instance
    read_in_FSP(guidance_function, instance_description)
end

function instance_csv_string(configuration::Configuration, instance::Instance)
    @sprintf("%d;%d", instance.n, instance.m)
end

function state_hash(configuration::Configuration, instance::Instance, Q::PFSPBeam, Q_idx::Int, layer::Int)
    selected_jobs = Set(@view Q.additional_data.perms[Q_idx,(layer+1):end])
    time_markers = @view Q.additional_data.time_markers[Q_idx,:]

    hash(selected_jobs, hash(time_markers))
end

function isequal_by_idx(Q_a::PFSPBeam, Q_b::PFSPBeam, Q_idx_a::Int, Q_idx_b::Int, layer::Int)
    time_markers_a = @view Q_a.additional_data.time_markers[Q_idx_a,:]
    time_markers_b = @view Q_b.additional_data.time_markers[Q_idx_b,:]

    if isequal(time_markers_a, time_markers_b)
        selected_jobs_a = Set(@view Q_a.additional_data.perms[Q_idx_a,(layer+1):end])
        selected_jobs_b = Set(@view Q_b.additional_data.perms[Q_idx_b,(layer+1):end])
        isequal(selected_jobs_a, selected_jobs_b)
    else
        false
    end
end

function pre_run_instance_name(instance_name::String)::String
    "ta001"
end

@time main()
#total_main_time = @elapsed @time instance_name, beam_width, guidance_function, best_objective, construction_time, stats = main()
#@printf("[CSV] %s;%d;%s;%d;%d;%f;%f;%f;%f;%f;%f;%f;%f\n", instance_name, beam_width, guidance_function, Threads.nthreads(), best_objective, total_main_time, construction_time, stats.main_loop_time, stats.memory_init_time, stats.successor_creation_time, stats.surviving_successors_time, stats.transitions_time, stats.terminal_check_time)
