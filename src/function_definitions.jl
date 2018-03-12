#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
Test if trial solution is better than the incumbent
"""
comparison(::Type{Minimisation}, trial, incumbent) = trial < incumbent
comparison(::Type{Maximisation}, trial, incumbent) = trial > incumbent

"""
Get worst case objective
"""
worstcase(::Type{Minimisation}) = 1e100  #  Inf
worstcase(::Type{Maximisation}) = -1e100 # -Inf

dynamics!(sp, f::Function) = (sp[:dynamics] = f)
dynamics!(f::Function, sp) = dynamics!(sp, f)

terminalobjective!(sp, f::Function) = (sp[:terminalcost] = f)
terminalobjective!(f::Function, sp) = terminalobjective!(sp, f)

constraints!(sp, f::Function) = (sp[:isfeasible] = f)
constraints!(f::Function, sp) = constraints!(sp, f)

"""
Function for adding states in Model definition
"""
function addstates!(sp; kwargs...)
    sp[:statespace] = GenericSpace(Float64;kwargs...)
end

"""
Function for adding controls in Model definition
"""
function addcontrols!(sp; kwargs...)
    sp[:controlspace] = GenericSpace(;kwargs...)
end

"""
Function for adding noises in Model definition
"""
function addnoises!(sp; kwargs...)
    sp[:noisespace] = GenericSpace(WeightedProbability;kwargs...)
end

"""
Interpolate the bellman surface using Interpolations.jl
"""
function interpolatesurface!(stage)
    stage.interpolatedsurface = interpolate(stage.statespace.dimensions, stage.bellmansurface, Gridded(Linear()))
end

"""
Calculate the end of horizon reward as current stage cost + terminal cost
"""
function calculatereward(::Type{TerminalReward}, m, t::Int, stage, newstate, state, control, noise)
    return stage.reward(state, control, noise) + stage.terminalcost(newstate)
end

function calculatereward(::Type{TerminalReward}, m, t, newstate, reward)
    reward + m.stages[t].terminalcost(newstate)
end

"""
Calculate the middle of horizon reward as current stage cost + interpolated value to go
"""
function calculatereward(::Type{InterpolatedReward}, m, t, stage, newstate, state, control, noise)
    stage.reward(state, control, noise) + getinterpolatedvalue(m.stages[t+1], newstate)
end

function calculatereward(::Type{InterpolatedReward}, m, t, newstate, reward)
    reward + getinterpolatedvalue(m.stages[t+1], newstate)
end

"""
Lookup the interpolated value at a state
"""
getinterpolatedvalue(stage, state) = stage.interpolatedsurface[state...]

"""
Create an array of WeightedProbability
"""
function DiscreteDistribution{T}(observations::AbstractVector{T}, probabilities::AbstractVector)
    @assert length(x) == length(weights)
    if !isapprox(sum(weights) , 1.0)
        warn("Weight vector is not normalised. Sum = $(sum(weights)). Normalising.")
        weights ./= sum(weights)
    end
    y = WeightedProbability{T}[]
    for i=1:length(x)
        push!(y, WeightedProbability(x[i], weights[i]))
    end
    y
end

"""
Risk measure for nested CVaR
"""
function NestedCVaR(isminimisation::Bool, beta, lambda)
    if lambda < 0.999
        return (x,p) -> lambda * dot(x,p) + (1 - lambda) * CVaR(x, p, beta, isminimisation)
    else
        return dot
    end
end
NestedCVaR(::Type{Minimisation}, beta, lambda) = NestedCVaR(true, beta, lambda)
NestedCVaR(::Type{Maximisation}, beta, lambda) = NestedCVaR(false, beta, lambda)

"""
Calculate CVaR of array
"""
function CVaR(x::AbstractVector{Float64},  p::AbstractVector{Float64}, beta::Float64, isminimisation::Bool)
    @assert length(x) == length(p)
    q = 0.
    cp = 0.
    cvar = 0.
    for i in sortperm(x, rev=isminimisation)
            if q <= beta
                cp = min(p[i], beta-q)
                q += cp
                cvar += cp * x[i]
            else
                break
            end
    end
    cvar / beta
end

"""
Here and Now models choose an action, then observe the noise.
"""
function _innerloop{T,N}(::Type{HereAndNow}, rewardtype::RewardType, m::SDPModel{T, N}, t::Int, state, newstate, riskmeasure::Function, outcomes::Vector{Float64}, probabilities::Vector{Float64})
    stage = m.stages[t]
    bestobj = worstcase(m.sense)
    reward = 0.0
    for control in product(stage.controlspace)
        n = 1
        for noise in product(stage.noisespace)
            prob = cumulateprobability(noise)
            if stage.isfeasible(state, control, noise)
                reward = stage.dynamics!(newstate, state, control, noise)::Float64
                if withinbounds(newstate, stage.statespace)
                    outcomes[n] = calculatereward(rewardtype, m, t, newstate, reward)::Float64
                else
                    outcomes[n] = worstcase(m.sense)
                end
            else
                outcomes[n] = worstcase(m.sense)
            end
            probabilities[n] = prob
            n += 1
        end
        tmp_obj = riskmeasure(outcomes, probabilities)
        if comparison(m.sense, tmp_obj, bestobj)
            bestobj = tmp_obj
        end
    end
    return bestobj::Float64
end
function cumulateprobability(noise)
    prob = 1.0
    for i in noise
        prob *= i.probability
    end
    prob::Float64
end
"""
Wait and see models observe the noise, then choose the best control.
The value is the expectation of taking those controls
"""
function _innerloop{T,N}(::Type{WaitAndSee}, rewardtype::RewardType, m::SDPModel{T, N}, t::Int, state, newstate, riskmeasure::Function, outcomes::Vector{Float64}, probabilities::Vector{Float64})

    stage = m.stages[t]
    n = 1
    reward = 0.0
    for noise in product(stage.noisespace)
        prob = cumulateprobability(noise)
        bestobj = worstcase(m.sense)

        for control in product(stage.controlspace)
            if stage.isfeasible(state, control, noise)
                reward = stage.dynamics!(newstate, state, control, noise)::Float64
                if withinbounds(newstate, stage.statespace)
                    tmp_obj = calculatereward(rewardtype, m, t, newstate, reward)::Float64
                    if comparison(m.sense, tmp_obj, bestobj)
                        bestobj = tmp_obj::Float64
                    end
                end
            end
        end
        outcomes[n] = bestobj
        probabilities[n] = prob
        n += 1
    end

    riskadjustedobjective = riskmeasure(outcomes, probabilities)
    if isnan(riskadjustedobjective)
        return worstcase(m.sense)
    end
    return riskadjustedobjective::Float64
end

function expectedvalue(x)
    ev = 0.
    for i in x
        ev += i.value * i.probability
    end
    ev
end
"""
ExpectedValue problem
"""
function _innerloop{T,N}(::Type{ExpectedValue}, rewardtype::RewardType, m::SDPModel{T, N}, t::Int, state, newstate, riskmeasure::Function, outcomes::Vector{Float64}, probabilities::Vector{Float64})
    stage = m.stages[t]
    noise = tuple([expectedvalue(dim) for dim in stage.noisespace.dimensions]...)
    bestobj = worstcase(m.sense)
    reward = 0.0
    for control in product(stage.controlspace)
        if stage.isfeasible(state, control, noise)
            reward = stage.dynamics!(newstate, state, control, noise)
            if withinbounds(newstate, stage.statespace)
                tmp_obj = calculatereward(rewardtype, m, t, newstate, reward)::Float64
                if comparison(m.sense, tmp_obj, bestobj)
                    bestobj = tmp_obj
                end
            end
        end
    end
    bestobj::Float64
end

function withinbounds(state, statespace)
    for i=1:length(state)
        if statespace.bounded[i]
            if (state[i] < statespace.minimum[i]) || (state[i] > statespace.maximum[i])
                return false
            end
        end
    end
    return true
end

"""
Solve end of horizon stage
"""
solveterminal!{T, N}(paralleltype, modtype, m::SDPModel{T, N},riskmeasure::Function) = solvestage!(paralleltype, modtype, TerminalReward, m, N, riskmeasure)

"""
Solve stage
"""
function solvestage!(::Type{Serial}, modtype, rewardtype::RewardType, m::SDPModel, t::Int,riskmeasure::Function)
    stage = m.stages[t]
    newstate = zeros(length(stage.statespace))
    outcomes = zeros(length(product(stage.noisespace)))
    probabilities = zeros(length(product(stage.noisespace)))

    map!(state->_innerloop(modtype, rewardtype, m, t, state, newstate, riskmeasure, outcomes, probabilities), stage.bellmansurface, collect(product(stage.statespace)))

    interpolatesurface!(stage)
end

function _innerloop(T, N, realisationtype, rewardtype, t, state, riskmeasure)
    mm = DynamicProgramming.m::SDPModel{T,N}
    stage = mm.stages[t]
    outcomes = DynamicProgramming._outcomes::Vector{Float64}
    probabilities = DynamicProgramming._probabilities::Vector{Float64}
    _innerloop(realisationtype, rewardtype, mm, t, state, zeros(length(stage.statespace)), riskmeasure, outcomes, probabilities)
end

function pmap!(results, f, lst)
    np = nprocs()  # determine the number of processes available
    n = length(lst)
    # results = Vector{Any}(n)
    i = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
                        results[idx] = remotecall_fetch(f, p, lst[idx])
                    end
                end
            end
        end
    end
end

function solvestage!{T, N}(::Type{Parallel}, modtype, rewardtype::RewardType, m::SDPModel{T, N}, t::Int, riskmeasure::Function)
    sendtoall(_probabilities=zeros(length(product(m.stages[t].noisespace))))
    sendtoall(_outcomes=zeros(length(product(m.stages[t].noisespace))))

    pmap!(m.stages[t].bellmansurface, state->_innerloop(T, N, modtype, rewardtype, t, state, riskmeasure), collect(product(m.stages[t].statespace)))

    interpolatesurface!(m.stages[t])

    distribute_work_void!((A)->(
            (DynamicProgramming.m::SDPModel{T,N}).stages[t].interpolatedsurface = A
        ),
        deepcopy(m.stages[t].interpolatedsurface)
    )
end

function distribute_work_void!(f::Function, args...)
    @sync begin
        for (i, procid) in enumerate(workers())
            @async begin
                remotecall_fetch(f, procid, args...)
            end
        end
    end
end
function sendto(p::Int, env; args...)
    for (nm, val) in args
        @spawnat(p, eval(env, Expr(:(=), nm, val)))
    end
end

function sendtoall(env=DynamicProgramming; args...)
    nprocs() == 1 && return false
    for p in procs()
        sendto(p, env; args...)
    end
    return true
end

"""
Solve SDPModel
"""
function solve(m::SDPModel;
    realisation::ModelType=WaitAndSee,
    riskmeasure::NestedCVaRType=Expectation()
    )
    solvetype = Serial

    if nprocs() >= 3
        solvetype = Parallel
        sendtoall(m = deepcopy(m))
    end

    risk_measure_function = NestedCVaR(m.sense, riskmeasure.beta, riskmeasure.lambda)
    totaltime = [0.0]
    # solve final stage
    printheader()
    tic()
    solveterminal!(solvetype, realisation, m, risk_measure_function)
    totaltime[1] += toq()
    printlog(length(m.stages), totaltime[1])
    # backwards recursion
    for t in (length(m.stages)-1):-1:1
        tic()
        solvestage!(solvetype, realisation, InterpolatedReward, m, t, risk_measure_function)
        totaltime[1] += toq()
        printlog(t, totaltime[1])
    end
end

function printheader()
    println("""-------------------------------------------------------------------------------
                       DynamicProgramming.jl © Oscar Dowson, 2018
    -------------------------------------------------------------------------------""")
    println("Stage | Elapsed Time")
    println("-------------------------------------------------------------------------------""")
end
function printlog(t::Int, time::Float64)
    println(humanize(t, "5d"), "| ", humanize(time, "8.3f"))
end
"""
Initialise storage for simulation
"""
function initialiseresult!(results::Dict{Symbol, Any}, replications, stages, key::Symbol, Ty::DataType)
    results[key] = Array{Ty}((stages, replications))
end
initialiseresult!{T}(results, replications, stages, key, x::AbstractVector{T}) = initialiseresult!(results, replications, stages, key, T)
initialiseresult!{T}(results, replications, stages, key, x::AbstractVector{WeightedProbability{T}}) = initialiseresult!(results, replications, stages, key, T)

"""
Simulate policy
"""
function simulate{T,N}(m::SDPModel{T,N},  n::Int; kwargs...)
    results = Dict{Symbol, Any}(
    :n_stages => N,
    :n_replications => n,
    :objective => zeros(n)
    )
    stage0 = m.stages[1]
    for space in [:statespace, :controlspace, :noisespace]
        for (key, idx) in getfield(stage0, space).nameindices
            initialiseresult!(results, n, N, key, getfield(stage0, space).dimensions[idx])
        end
    end
    stagecost = 0.
    futurecostcost = 0.
    beststagecost = 0.
    bestcontrol = Float64[] # TODO: Type stability
    for replication = 1:n
        state = zeros(length(stage0.statespace))
        newstate = similar(state)
        for (key, val) in kwargs
            if !haskey(stage0.statespace.nameindices, key)
                error("State not found")
            else
                state[stage0.statespace.nameindices[key]] = val
            end
        end
        for (t, stage) in enumerate(m.stages)
            noise = map(i->WeightedProbability(rand(i), 0.), stage.noisespace.dimensions)
            bestobj = worstcase(m.sense)
            beststagecost = worstcase(m.sense)
            bestcontrol = zeros(length(stage.controlspace))

            for control in product(stage.controlspace)
                if stage.isfeasible(state, control, noise)
                    stagecost = stage.dynamics!(newstate, state, control, noise)
                    if withinbounds(newstate, stage.statespace)
                        if t < length(m.stages)
                            futurecost = getinterpolatedvalue(m.stages[t+1], newstate)
                        else
                            stagecost += stage.terminalcost(newstate)
                            futurecost = 0.
                        end
                        if comparison(m.sense, stagecost + futurecost, bestobj)
                            bestobj       = stagecost + futurecost
                            beststagecost = stagecost
                            bestcontrol   = identity(control)
                        end
                    end
                end
            end
            stage.dynamics!(newstate, state, bestcontrol, noise)
            results[:objective][replication] += beststagecost
            for (key, index) in stage.statespace.nameindices
                results[key][t, replication] = newstate[index]
            end
            for (key, index) in stage.controlspace.nameindices
                results[key][t, replication] = bestcontrol[index]
            end
            for (key, index) in stage.noisespace.nameindices
                results[key][t, replication] = noise[index].value
            end
            copy!(state, newstate)

        end
    end
    results
end
