# TTP node estimators based on cvrp and cvrph lower bounds

module TTPGuidance
    using Printf: @printf

    export heuristic_estimate

    # heuristic estimates based on precalculated exact cvrp bounds
    function heuristic_estimate(streak_limit::Int, d::Array{Int, 2}, team::Int, teams_left_mask::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int, bounds_by_state::Array{UInt32,4})
        if teams_left_mask == 0
            return d[position, team]
        end

        if team == position
            streak = 0
        end
        if streak == streak_limit
            d[position, team] + bounds_by_state[team, teams_left_mask+1, team, 1]
        else
            min(d[position, team] + bounds_by_state[team, teams_left_mask+1, team, 1], bounds_by_state[team, teams_left_mask+1, position, streak+1])
        end
    end

    # heuristic estimates based on precalculated exact cvrph bounds
    function heuristic_estimate(streak_limit::Int, d::Array{Int, 2}, team::Int, teams_left_mask::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int, bounds_by_state::Array{UInt32,5})
        if teams_left_mask == 0
            return d[position, team]
        end

        if team == position
            away_streak = 0
        else
            away_streak = streak
        end

        # we hit the away streak limit and have to return home, adding a detour
        if team != position && streak == streak_limit
            detour = d[position, team]
            away_streak = 0
            streak = 0
            position = team
        else
            detour = 0
        end

        minimum_home_stands_needed, maximum_home_stands_allowed = min_max_home_stands(streak_limit, team, number_of_away_games_left, number_of_home_games_left, position, streak)

        best_bound_direct = minimum(@view bounds_by_state[team, teams_left_mask+1, position, away_streak+1, minimum_home_stands_needed:maximum_home_stands_allowed])
        if team != position
            minimum_home_stands_needed_detour, maximum_home_stands_allowed_detour = min_max_home_stands(streak_limit, team, number_of_away_games_left, number_of_home_games_left, team, 0)
            best_bound_via_home = d[position, team] + minimum(@view bounds_by_state[team, teams_left_mask+1, team, 1, minimum_home_stands_needed_detour:maximum_home_stands_allowed_detour])
        else
            best_bound_via_home = typemax(Int)
        end

        return detour + min(best_bound_direct, best_bound_via_home)
    end

    # vehicle bounds

    # the minimum number and maximum number of vehicles needed for solving the CVRPH problem for a team, inferred by the feasibilty and optimality considerations
    #     the current streak counts as a vehicle
    function min_max_vehicles(streak_limit::Int, team::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int)
        minimum_vehicles_needed = max(min_vehicles_by_home_games(streak_limit, team, number_of_away_games_left, number_of_home_games_left, position, streak), min_vehicles_by_away_games(streak_limit, team, number_of_away_games_left, number_of_home_games_left, position, streak))
        maximum_vehicles_allowed = min(max_vehicles_by_home_games(streak_limit, team, number_of_away_games_left, number_of_home_games_left, position, streak), max_vehicles_by_away_games(streak_limit, team, number_of_away_games_left, number_of_home_games_left, position, streak))
        maximum_vehicles_allowed = max(maximum_vehicles_allowed, minimum_vehicles_needed)

        minimum_vehicles_needed, maximum_vehicles_allowed
    end

    # this is similar to min max vehicles, but for the precalculated lower bound values we do not count the vehicles (i.e. streaks) but the home stands, where the return home at the end also counts as home stand
    function min_max_home_stands(streak_limit::Int, team::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int)
        minimum_vehicles_needed, maximum_vehicles_allowed = min_max_vehicles(streak_limit, team, number_of_away_games_left, number_of_home_games_left, position, streak)

        if team == position
            minimum_vehicles_needed + 1, maximum_vehicles_allowed + 1
        else
            minimum_vehicles_needed, maximum_vehicles_allowed
        end
    end

    # the remaining number of home games imply a minimum of vehicles (away streaks) for feasibility...
    function min_vehicles_by_home_games(streak_limit::Int, team::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int)
        if team == position
            max(convert(Int64, ceil((streak + number_of_home_games_left)/streak_limit))-1, 0)
        else
            convert(Int64, ceil(number_of_home_games_left/streak_limit))
        end
    end

    # .. so do the remaining away games
    function min_vehicles_by_away_games(streak_limit::Int, team::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int)
        if team == position
            convert(Int64, ceil(number_of_away_games_left/streak_limit))
        else
            convert(Int64, ceil((streak + number_of_away_games_left)/streak_limit))
        end
    end

    # the number of home games left imply a maximum number of vehicles
    function max_vehicles_by_home_games(streak_limit::Int, team::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int)
        1 + number_of_home_games_left
    end

    # if the streak limit > 1, it is never optimal to have more the one streak with length one, since they can be merged
    function max_vehicles_by_away_games(streak_limit::Int, team::Int, number_of_away_games_left::Int, number_of_home_games_left::Int, position::Int, streak::Int)
        if team == position
            convert(Int64, ceil(number_of_away_games_left/min(streak_limit, 2)))
        else
            if streak == 1 && streak_limit > 1
                1 + convert(Int64, ceil((number_of_away_games_left-1)/2))
            else
                1 + convert(Int64, ceil(number_of_away_games_left/min(streak_limit, 2)))
            end
        end
    end
end
