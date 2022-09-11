using Primes
using ConcurrentCollections

# for duplicates counter
mutable struct Atomic{T}; @atomic x::T; end

# helper functions
function initialize_hash_table(configuration::Configuration)
    m = prevprime(configuration.beam_width)

    hash_table = Array{Vector{Int}}(undef, m)
    hash_locks = Array{Threads.SpinLock}(undef, m)

    Threads.@threads for i in 1:m
        hash_table[i] = Vector{Int}()
        hash_locks[i] = Threads.SpinLock()
    end

    hash_table, hash_locks, m
end


function initialize_hash_table_old(configuration::Configuration)
    m = prevprime(configuration.beam_width)

    hash_table = Array{ConcurrentQueue{Int}}(undef, m)

    Threads.@threads for i in 1:m
        hash_table[i] = ConcurrentQueue{Int}()
    end

    hash_table, m
end

function hash_dict_get(H::Dict{UInt64,Int}, hash::UInt64, lk::Threads.SpinLock)
    lock(lk) do
        get(H, hash, 0)
    end
end

function hash_dict_set(H::Dict{UInt64,Int}, hash::UInt64, val::Int, lk::Threads.SpinLock)
    lock(lk) do
        H[hash] = val
    end
end

function conditionally_filter_dominated_duplicate!(Q::AbstractBeam, idx_1::Int, idx_2::Int, layer::Number, duplicates_filtered::Atomic{Int})
    if isequal_by_idx(Q, Q, idx_1, idx_2, layer)
        #@atomic duplicates_filtered.x += 1

        if Q.nodes.priority[idx_1] < Q.nodes.priority[idx_2]
            #Q.nodes.layer[idx_2] -= 1
            if Q.nodes.layer[idx_2] == layer
                Q.nodes.layer[idx_2] -= 1
                @atomic duplicates_filtered.x += 1
            end
            #push!(filtered_indices, idx_2)
        elseif Q.nodes.priority[idx_1] == Q.nodes.priority[idx_2]
            if idx_1 < idx_2
                #Q.nodes.layer[idx_2] -= 1
                if Q.nodes.layer[idx_2] == layer
                    Q.nodes.layer[idx_2] -= 1
                    @atomic duplicates_filtered.x += 1
                end
                #push!(filtered_indices, idx_2)
            else
                #Q.nodes.layer[idx_1] -= 1
                if Q.nodes.layer[idx_1] == layer
                    Q.nodes.layer[idx_1] -= 1
                    @atomic duplicates_filtered.x += 1
                end
                #push!(filtered_indices, idx_1)
            end
        else
            #Q.nodes.layer[idx_1] -= 1
            if Q.nodes.layer[idx_1] == layer
                Q.nodes.layer[idx_1] -= 1
                @atomic duplicates_filtered.x += 1
            end
            #push!(filtered_indices, idx_1)
        end
    end
end

##
# concurrent hashing with global spinlocks
##

# use an array of spinlocks to manage access to dict, dict maps hash to index to currently dominating state
function filter_duplicate_states_spinlocked(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, lks::Array{Threads.SpinLock}, layer::Number)
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i, layer)
    end

    H = Dict{UInt64,Int}()
    duplicates_filtered = Atomic{Int}(0)
    lks_len = length(lks)

    Threads.@threads for i in 1:Q_len
        current_hash = hashes[i]
        lk_idx = current_hash % lks_len + 1
        lk = lks[lk_idx]
        existing_i = hash_dict_get(H, current_hash, lk)
        @assert 1 <= existing_i <= configuration.beam_width
        if existing_i != 0
            if isequal_by_idx(Q, Q, i, existing_i, layer)
                @atomic duplicates_filtered.x += 1

                lock(lk) do
                    if Q.nodes.priority[existing_i] <= Q.nodes.priority[i]
                        Q.nodes.layer[i] -= 1
                    else
                        Q.nodes.layer[existing_i] -= 1
                        H[current_hash] = i
                    end
                end
            end
        else
            hash_dict_set(H, current_hash, i, lk)
        end
    end

    duplicates_filtered.x
end

# test function to measure collection of states into lists using spinlocks
function filter_duplicate_states_spinlocked_collect(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, lks::Array{Threads.SpinLock}, layer::Number)
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i, layer)
    end

    H = Dict{UInt64,Vector{Int}}()
    lks_len = length(lks)
    Threads.@threads for i in 1:Q_len
        current_hash = hashes[i]

        lk_idx = current_hash % lks_len + 1
        lk = lks[lk_idx]
        lock(lk) do
            if haskey(H, current_hash)
                push!(H[current_hash], i)
            else
                H[current_hash] = [i]::Vector{Int}
            end
        end
    end

    duplicates_filtered = Atomic{Int}(0)
    duplicates_filtered.x
end

##
# concurrent hashing with separate dicts
##

# own dict for each thread collecting all states with same hash, thread which has found lowest index for a state hash collects and checks potential duplicates
# deterministic/independent of thread number but can grow quadratically in the number of duplicates 
function filter_duplicate_states_dict_merge(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, layer::Number)
    hashes = Array{UInt64}(undef, Q_len)
    first_seen = Array{Int}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i, layer)
    end

    Hs = Array{Dict{UInt64,Vector{Int}}}(undef, Threads.nthreads())

    for i in 1:Threads.nthreads()
        Hs[i] = Dict{UInt64,Vector{Int}}()
    end

    Threads.@threads for i in 1:Q_len
        H = Hs[Threads.threadid()]
        current_hash = hashes[i]

        if haskey(H, current_hash)
            push!(H[current_hash], i)
            first_seen[i] = H[current_hash][1]
        else
            H[current_hash] = [i]::Vector{Int}
            first_seen[i] = i
        end
    end

    duplicates_filtered = Atomic{Int}(0)

    Threads.@threads for i in 1:Q_len
        if first_seen[i] != i
            continue
        end

        current_hash = hashes[i]
        found_count = 0
        assigned_thread = 0
        for j in 1:length(Hs)
            if haskey(Hs[j], current_hash)
                found_count += length(Hs[j][current_hash])
                if assigned_thread == 0
                    assigned_thread = j
                end
            end
        end
        @assert found_count >= 1
        @assert assigned_thread >= 1

        if found_count > 1 && assigned_thread == Threads.threadid()
            possible_duplicates = Vector{Int}()
            for j in 1:length(Hs)
                if haskey(Hs[j], current_hash)
                    append!(possible_duplicates, Hs[j][current_hash])
                end
            end
            @assert length(possible_duplicates) == found_count

            for k in 1:(length(possible_duplicates)-1)
                for l in (k+1):length(possible_duplicates)
                    idx_1 = possible_duplicates[k]
                    idx_2 = possible_duplicates[l]
    
                    @assert l > k
                    @assert idx_1 != idx_2
                    conditionally_filter_dominated_duplicate!(Q, idx_1, idx_2, layer, duplicates_filtered)
                end
            end
        end
    end
#=     merged_H = Hs[1]
    for i in 2:length(Hs)
        merged_H = merge(union, merged_H, Hs[i])
    end

    possible_duplicates = Vector{Vector{Int}}()
    for (current_hash, indices) in merged_H
        if length(indices) > 1
            push!(possible_duplicates, indices)
        end
    end

    Threads.@threads for k in 1:length(possible_duplicates)
        indices = possible_duplicates[k]
        for i in 1:(length(indices)-1)
            for j in (i+1):length(indices)
                idx_1 = indices[i]
                idx_2 = indices[j]

                if isequal_by_idx(Q, Q, idx_1, idx_2)
                    @atomic duplicates_filtered.x += 1
    
                    if Q.nodes.priority[idx_1] <= Q.nodes.priority[idx_2]
                        Q.nodes.layer[idx_2] -= 1
                    else
                        Q.nodes.layer[idx_1] -= 1
                    end
                end
            end
        end
    end =#

    duplicates_filtered.x
end

# duplicate detection only within threads, parallelizes well but beam search thread dependent
function filter_duplicate_states_separate_dicts(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, layer::Number)
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i, layer) 
    end

    Hs = Array{Dict{UInt64,Int}}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        Hs[i] = Dict{UInt64,Int}()
    end

    duplicates_filtered = Atomic{Int}(0)
    Threads.@threads for i in 1:Q_len
        H = Hs[Threads.threadid()]
        current_hash = hashes[i]

        existing_i = get(H, current_hash, 0)
        if existing_i != 0
            if isequal_by_idx(Q, Q, i, existing_i, layer)
                @atomic duplicates_filtered.x += 1

                if Q.nodes.priority[existing_i] <= Q.nodes.priority[i]
                    Q.nodes.layer[i] -= 1
                else
                    Q.nodes.layer[existing_i] -= 1
                    H[current_hash] = i
                end
            end
        else
            H[current_hash] = i
        end
    end    

    duplicates_filtered.x
end

# first duplicate detection within threads, then check between threads iterating over dicts for a given hash
function filter_duplicate_states_separate_dicts_final_merge(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, layer::Number)
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i, layer) 
    end

    Hs = Array{Dict{UInt64,Int}}(undef, Threads.nthreads())
    for i in 1:Threads.nthreads()
        Hs[i] = Dict{UInt64,Int}()
    end

    duplicates_filtered = Atomic{Int}(0)
    Threads.@threads for i in 1:Q_len
        H = Hs[Threads.threadid()]
        current_hash = hashes[i]

        existing_i = get(H, current_hash, 0)
        if existing_i != 0
            @assert i != existing_i
            if isequal_by_idx(Q, Q, i, existing_i, layer)
                @atomic duplicates_filtered.x += 1

                if Q.nodes.priority[existing_i] <= Q.nodes.priority[i]
                    Q.nodes.layer[i] -= 1
                else
                    Q.nodes.layer[existing_i] -= 1
                    H[current_hash] = i
                end
            end
        else
            H[current_hash] = i
        end
    end

    if length(Hs) > 1
        Threads.@threads for i in 1:Q_len
            current_hash = hashes[i]
            previous_i = 0

            for j in 1:length(Hs)
                existing_i = get(Hs[j], current_hash, 0)
                if existing_i != 0
                    @assert existing_i != previous_i
                    # seen it already
                    if previous_i != 0
                        # possible race conditions in counting if two threads access same data, but dominance check is still correct
                        if Q.nodes.layer[existing_i] == layer && Q.nodes.layer[previous_i] == layer && isequal_by_idx(Q, Q, previous_i, existing_i, layer)
                            @atomic duplicates_filtered.x += 1
                        
                            if Q.nodes.priority[previous_i] <= Q.nodes.priority[existing_i]
                                Q.nodes.layer[existing_i] -= 1
                            else
                                Q.nodes.layer[previous_i] -= 1
                                previous_i = existing_i
                            end
                        end
                    # see it for the first time    
                    else
                        previous_i = existing_i
                    end
                end
            end
        end
    end

    duplicates_filtered.x
end

##
# concurrent hashing with thread-safe lists for each slot
##

function filter_duplicate_states_concurrent_queues_old(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, hash_table::Array{ConcurrentQueue{Int}}, m::Int, layer::Number)
    Threads.@threads for i in 1:Q_len
        hash_val = state_hash(configuration, instance, Q, i, layer)
        idx = 1 + hash_val % m
        push!(hash_table[idx], i)
    end

    duplicates_filtered = Atomic{Int}(0)

    Threads.@threads for i in 1:m
        possible_duplicates = Vector{Int}()

        while true
            maybe_first = maybepopfirst!(hash_table[i])
            if maybe_first === nothing
                break
            end
            push!(possible_duplicates, maybe_first.value)
        end

        if length(possible_duplicates) > 1
            for k in 1:(length(possible_duplicates)-1)
                for l in (k+1):length(possible_duplicates)
                    idx_1 = possible_duplicates[k]
                    idx_2 = possible_duplicates[l]

                    @assert l > k
                    @assert idx_1 != idx_2
                    conditionally_filter_dominated_duplicate!(Q, idx_1, idx_2, layer, duplicates_filtered)
                end
            end
        end
    end

    duplicates_filtered.x
end

function filter_duplicate_states_concurrent_queues(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, hash_table::Array{Vector{Int}}, hash_locks::Array{Threads.SpinLock}, m::Int, layer::Number)
    if Threads.nthreads() == 1
        for i in 1:Q_len
            hash_val = state_hash(configuration, instance, Q, i, layer)
            idx = 1 + hash_val % m
            push!(hash_table[idx], i)
        end
    else
        Threads.@threads for i in 1:Q_len
            hash_val = state_hash(configuration, instance, Q, i, layer)
            idx = 1 + hash_val % m
            lock(hash_locks[idx]) do
                push!(hash_table[idx], i)
            end
        end
    end

    duplicates_filtered = Atomic{Int}(0)

    #max_lengths = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:m
        possible_duplicates = hash_table[i]
        #max_lengths[Threads.threadid()] = max(max_lengths[Threads.threadid()], length(possible_duplicates))
        #filtered_indices = Set()

        if length(possible_duplicates) > 1
            for k in 1:(length(possible_duplicates)-1)
                for l in (k+1):length(possible_duplicates)
                    idx_1 = possible_duplicates[k]
                    idx_2 = possible_duplicates[l]

                    @assert l > k
                    @assert idx_1 != idx_2
                    conditionally_filter_dominated_duplicate!(Q, idx_1, idx_2, layer, duplicates_filtered)
                end
            end
        end

        empty!(possible_duplicates)
        #@atomic duplicates_filtered.x += length(filtered_indices)

        
        #@assert length(hash_table[i]) == 0
    end
    #println(layer)
    #println(maximum(max_lengths))

    duplicates_filtered.x
end

# each hash bucket holds list of currently non-dominated states
# on insert protected by spinlock, check whether new state dominates existing (then replace and stop) or is dominated (then deactive state and stop)
#    if not, append to list, since it is currently non-dominated
# by domination we mean here the same state with better priority
# algorithm is deterministic and independent from number of threads
function maintain_non_dominated_list!(Q::AbstractBeam, hash_table::Array{Vector{Int}}, i::Int, idx::UInt, current_priority::Number, layer::Number, duplicates_filtered::Atomic{Int})
    non_dominated = true
    for j in 1:length(hash_table[idx])
        existing_i = hash_table[idx][j]
        @assert existing_i != i
        @assert Q.nodes.layer[existing_i] == layer
        if isequal_by_idx(Q, Q, existing_i, i, layer)
            if Q.nodes.priority[existing_i] < current_priority
                Q.nodes.layer[i] -= 1
                @atomic duplicates_filtered.x += 1
                non_dominated = false
                break
            elseif Q.nodes.priority[existing_i] == current_priority && existing_i < i
                Q.nodes.layer[i] -= 1
                @atomic duplicates_filtered.x += 1
                non_dominated = false
                break
            else
                Q.nodes.layer[existing_i] -= 1
                hash_table[idx][j] = i
                @atomic duplicates_filtered.x += 1
                non_dominated = false
                break
            end
        end
    end
    if non_dominated
        push!(hash_table[idx], i)
    end
end

function filter_duplicate_states_spinlocked_non_dominated_lists(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, hash_table::Array{Vector{Int}}, hash_locks::Array{Threads.SpinLock}, m::Int, layer::Number)
    duplicates_filtered = Atomic{Int}(0)

    if Threads.nthreads() == 1
        for i in 1:Q_len
            hash_val = state_hash(configuration, instance, Q, i, layer)
            idx = 1 + hash_val % m
            current_priority = Q.nodes.priority[i]
            maintain_non_dominated_list!(Q, hash_table, i, idx, current_priority, layer, duplicates_filtered)
        end
    else
        Threads.@threads for i in 1:Q_len
            hash_val = state_hash(configuration, instance, Q, i, layer)
            idx = 1 + hash_val % m
            current_priority = Q.nodes.priority[i]
            lock(hash_locks[idx]) do
                maintain_non_dominated_list!(Q, hash_table, i, idx, current_priority, layer, duplicates_filtered)
            end
        end
    end

    Threads.@threads for i in 1:m
        empty!(hash_table[i])
    end

    duplicates_filtered.x
end

##
# pair-wise domination check
##

function filter_dominated_states(configuration, instance, Q, Q_len)
    duplicates_filtered = Atomic{Int}(0)

    Threads.@threads for i in 1:Q_len
        for j in 1:Q_len
            if i != j && isdominated_by_idx(Q, Q, i, j)
                @atomic duplicates_filtered.x += 1
                Q.nodes.layer[i] -= 1
            end
        end
    end

    duplicates_filtered.x
end