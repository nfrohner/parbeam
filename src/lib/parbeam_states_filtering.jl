using Primes
using ConcurrentCollections

mutable struct Atomic{T}; @atomic x::T; end

function filter_duplicate_states_spinlocked(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, lks::Array{Threads.SpinLock})
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i) 
    end

    H = Dict{UInt64,Int}()
    duplicates_filtered = Atomic{Int}(0)
    lks_len = length(lks)

    Threads.@threads for i in 1:Q_len
        current_hash = hashes[i]
        lk_idx = current_hash % lks_len + 1
        lk = lks[lk_idx]
        existing_i = hash_dict_get(H, current_hash, lk)
        if existing_i != 0
            if isequal_by_idx(Q, Q, i, existing_i)
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

function filter_duplicate_states_spinlocked_collect(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, lks::Array{Threads.SpinLock})
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i) 
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

function filter_duplicate_states_dict_merge(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int)
    hashes = Array{UInt64}(undef, Q_len)
    first_seen = Array{Int}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i) 
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
                    if isequal_by_idx(Q, Q, idx_1, idx_2)
                        @atomic duplicates_filtered.x += 1
        
                        if Q.nodes.priority[idx_1] < Q.nodes.priority[idx_2]
                            Q.nodes.layer[idx_2] -= 1
                        elseif Q.nodes.priority[idx_1] == Q.nodes.priority[idx_2]
                            if idx_1 < idx_2
                                Q.nodes.layer[idx_2] -= 1
                            else
                                Q.nodes.layer[idx_1] -= 1
                            end
                        else
                            Q.nodes.layer[idx_1] -= 1
                        end
                    end
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

function filter_duplicate_states_separate_dicts(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int)
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i) 
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
            if isequal_by_idx(Q, Q, i, existing_i)
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

function filter_duplicate_states_separate_dicts_final_merge(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int)
    hashes = Array{UInt64}(undef, Q_len)
    Threads.@threads for i in 1:Q_len
        hashes[i] = state_hash(configuration, instance, Q, i) 
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
            if isequal_by_idx(Q, Q, i, existing_i)
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
                        if isequal_by_idx(Q, Q, previous_i, existing_i)
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

function initialize_hash_table_old(configuration::Configuration)
    m = prevprime(configuration.beam_width)

    hash_table = Array{ConcurrentQueue{Int}}(undef, m)
    
    Threads.@threads for i in 1:m
        hash_table[i] = ConcurrentQueue{Int}()
    end

    hash_table, m
end

function filter_duplicate_states_concurrent_queues_old(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, hash_table::Array{ConcurrentQueue{Int}}, m::Int)
    Threads.@threads for i in 1:Q_len
        hash_val = state_hash(configuration, instance, Q, i)
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
        end
    end

    duplicates_filtered.x
end

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

function filter_duplicate_states_concurrent_queues(configuration::Configuration, instance::AbstractInstance, Q::AbstractBeam, Q_len::Int, hash_table::Array{Vector{Int}}, hash_locks::Array{Threads.SpinLock}, m::Int, layer::Number)
    if Threads.nthreads() == 1
        Threads.@threads for i in 1:Q_len
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

    Threads.@threads for i in 1:m
        possible_duplicates = hash_table[i]
        #filtered_indices = Set()

        if length(possible_duplicates) > 1
            for k in 1:(length(possible_duplicates)-1)
                for l in (k+1):length(possible_duplicates)
                    idx_1 = possible_duplicates[k]
                    idx_2 = possible_duplicates[l]

                    @assert l > k
                    @assert idx_1 != idx_2
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
            end
        end

        empty!(possible_duplicates)
        #@atomic duplicates_filtered.x += length(filtered_indices)

        #@assert length(hash_table[i]) == 0
    end

    duplicates_filtered.x
end

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