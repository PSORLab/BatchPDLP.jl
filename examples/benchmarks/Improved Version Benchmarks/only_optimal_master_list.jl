
# Bring in the packages we need
using SourceCodeMcCormick

# Create a "problem" struct that contains useful information.
mutable struct Problem
    name::String
    vars::Vector{Num}
    nvars::Int
    ncons::Int
    lvbs::Vector{Float64}
    uvbs::Vector{Float64}
    obj::Num
    eq::Vector{Num}
    leq::Vector{Num}
    geq::Vector{Num}
    sense::String
end

function Problem(name::String, vars::Vector{Num}, obj::Num, lvbs, uvbs; eq=Num[], leq=Num[], geq=Num[], sense="min")
    return Problem(name, vars, length(vars), length(eq)+length(leq)+length(geq), Float64.(lvbs), Float64.(uvbs), obj, eq, leq, geq, sense)
end

mutable struct LoadedProblem
    name::String
    vars::Vector{Num}
    nvars::Int
    ncons::Int
    lvbs::Vector{Float64}
    uvbs::Vector{Float64}
    obj_fun::Function
    eq_cons::Vector{Function}
    leq_cons::Vector{Function}
    geq_cons::Vector{Function}
    obj_sp::Vector{Bool}
    eq_sp::Vector{Vector{Bool}}
    leq_sp::Vector{Vector{Bool}}
    geq_sp::Vector{Vector{Bool}}
    sense::String
end

function LoadedProblem(prob::Problem; overwrite::Bool=false)
    return LoadedProblem(prob.name, prob.vars, prob.nvars, prob.ncons, prob.lvbs, prob.uvbs, 
        kgen(prob.obj, prob.vars, overwrite=overwrite),
        kgen.(prob.eq, Ref(prob.vars), overwrite=overwrite),
        kgen.(prob.leq, Ref(prob.vars), overwrite=overwrite),
        kgen.(prob.geq, Ref(prob.vars), overwrite=overwrite),
        [x in string.(pull_vars(prob.obj)) ? true : false for x in string.(prob.vars)],
        [[x in string.(pull_vars(prob.eq[i])) ? true : false for x in string.(prob.vars)] for i in eachindex(prob.eq)],
        [[x in string.(pull_vars(prob.leq[i])) ? true : false for x in string.(prob.vars)] for i in eachindex(prob.leq)],
        [[x in string.(pull_vars(prob.geq[i])) ? true : false for x in string.(prob.vars)] for i in eachindex(prob.geq)],
        prob.sense)
end


# Create symbolics variables for the examples
for i=1:1001 
    var = Symbol("x$i")
    @eval Symbolics.@variables $var
end

# Create oft-used functions
sqr(x) = x^2

# List of included instances, and the number of variables and constraints
included = [
            "ex4_1_1"                  1     0;
            "rbrock"                   2     0;
            "ex8_4_2"                 24    10;
            "prob09"                   3     1;
            "eq6_1"                   16    60;
            "ex4_1_3"                  1     0;
            "ex4_1_4"                  1     0;
            "ex8_1_1"                  1     1;
            "ex8_1_2"                  1     0;
            "ex8_1_6"                  2     0;
            "ex8_4_7"                 62    40;
            "ex8_6_2"                 30     0;
            "kriging_peaks_red010"     2     0;
            "mathopt5_1"               1     0;
            "mathopt5_2"               1     0;
            "mathopt5_3"               1     0;
            "mathopt5_4"               1     0;
            "mathopt5_5"               1     0;
            "mathopt5_7"               1     0;
            "mathopt5_8"               1     0;
            "mathopt6"                 1     0;
            "maxmin"                  27    78; 
]



for i in included[:,1]
    try
        include("./all_examples/$i.jl")
    catch
        println("$i failed")
    end
end