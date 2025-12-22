
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
original_included = [
 "ex4_1_1"                  1     0;
 "ex4_1_2"                  1     0;
 "ex4_1_3"                  1     0;
 "ex4_1_4"                  1     0;
 "ex4_1_6"                  1     0;
 "ex4_1_7"                  1     0;
 "ex8_1_2"                  1     0;
 "mathopt5_1"               1     0;
 "mathopt5_2"               1     0;
 "mathopt5_3"               1     0;
 "mathopt5_4"               1     0;
 "mathopt5_5"               1     0;
 "mathopt5_7"               1     0;
 "mathopt5_8"               1     0;
 "mathopt6"                 1     0;
 "ex8_1_1"                  1     1;
 "trig"                     1     1;
 "ex4_1_5"                  2     0;
 "ex8_1_3"                  2     0; # Some singular errors
 "ex8_1_4"                  2     0;
 "ex8_1_5"                  2     0; # Lots of "large complementary slackness residual" warnings
 "ex8_1_6"                  2     0;
 "kriging_peaks_red010"     2     0;
#  "kriging_peaks_red020"     2     0; # Objective length is 9578
#  "kriging_peaks_red030"     2     0; # Objective length is 14362
#  "kriging_peaks_red050"     2     0; # Objective length is 23762
#  "kriging_peaks_red100"     2     0; # Objective length is 47766
#  "kriging_peaks_red200"     2     0; # Objective length is 95485
#  "kriging_peaks_red500"     2     0; # Objective length is 238721
 "rbrock"                   2     0;
 "ex4_1_8"                  2     1;
 "ex4_1_9"                  2     2;
 "ex14_1_9"                 2     2;
 "mathopt1"                 2     2;
 "mathopt4"                 2     2;
 "st_e19"                   2     2;
 "trigx"                    2     2;
 "mathopt2"                 2     4;
 "least"                    3     0; # This causes one of the other solvers to crash vscode
 "ex6_2_6"                  3     1;
 "ex6_2_8"                  3     1;
 "ex6_2_11"                 3     1;
 "hs62"                     3     1;
 "prob09"                   3     1;
 "st_e11"                   3     2;
 "st_e06"                   3     3;
 "ex14_1_1"                 3     4;
 "ex14_1_3"                 3     4;
 "ex14_1_4"                 3     4;
 "ex14_1_8"                 3     4;
 "st_e37"                   4     1;
 "ex6_2_9"                  4     2;
 "ex6_2_12"                 4     2;
 "ex6_2_14"                 4     2;
 "st_e04"                   4     2;
 "st_e41"                   4     2;
 "chance"                   4     3;
 "ex6_1_2"                  4     3;
 "st_e12"                   4     3;
 "ex14_2_2"                 4     5;
 "ex14_2_5"                 4     5;
 "ex14_2_8"                 4     5;
 "ex14_2_9"                 4     5;
 "ex7_3_1"                  4     7;
 "ex7_3_2"                  4     7; # This causes one of the other solvers to crash vscode
 "mhw4d"                    5     3;
 "ex8_5_3"                  5     4; # This gets NaN constraints
 "ex8_5_4"                  5     4; # This gets NaN constraints
 "ex8_5_5"                  5     4; # This gets NaN constraints
 "ex8_1_7"                  5     5;
 "ex14_2_1"                 5     7;
 "ex14_2_4"                 5     7;
 "ex14_2_6"                 5     7;
 "ex6_2_10"                 6     3;
 "ex6_2_13"                 6     3;
 "ex6_1_4"                  6     4;
 "ex8_5_1"                  6     4; # This gets NaN constraints
 "ex8_5_2"                  6     4; # This gets NaN constraints
 "ex8_5_6"                  6     4; # This gets NaN constraints
 "ex7_2_2"                  6     5;
 "ex14_1_5"                 6     6; # This one crashes Julia
 "st_e21"                   6     6;
 "wallfix"                  6     6; # This gets NaN constraints
 "mathopt3"                 6     7; # Causes Julia to crash
 "ex14_1_2"                 6     9;
 "ex14_2_3"                 6     9;
 "ex14_2_7"                 6     9;
 "ex7_2_1"                  7    14;
 "ex7_2_4"                  8     4;
 "ex6_1_1"                  8     6; # Seems to be taking a VERY long time (1-2 hours?)
 "ex7_2_3"                  8     6;
 "inscribedsquare01"        8     8;
 "inscribedsquare02"        8     8;
 "ex6_2_5"                  9     3;
 "ex6_2_7"                  9     3;
 "like"                     9     3; # This gets NaN constraints
 "process"                 10     7;
 "st_e03"                  10     7;
 "alkylation"              10    11;
 "ex14_1_7"                10    17; # This crashes Julia (maybe in Gurobi or HiGHS PDLP?)
 "shiporig"                10    17; # This gets NaN constraints
 "chem"                    11     4;
 "ann_fermentation_exp"    12     9; # This gets NaN constraints (due to a problem in McCormick)
 "ex6_1_3"                 12     9;
 "st_e16"                  12     9;
 "ex7_3_4"                 12    17; # Seems to be taking a VERY long time (overnight, 10+ hours?)
 "ex7_3_5"                 13    15;
 "alkyl"                   14     7;
 "ex8_4_6"                 14     8;
 "ex8_4_5"                 15    11; # This gets NaN constraints
 "prob07"                  15    36;
 "ex5_4_3"                 16    13;
 "eq6_1"                   16    60;
 "ex7_3_6"                 17    17;
 "ex8_4_2"                 24    10;
 "orth_d3m6"               25    62;
 "kriging_peaks_full010"   26    24;
 "ex5_4_4"                 27    19;
 "maxmin"                  27    78; # This breaks Julia (maybe in Gurobi or HiGHS PDLP?)
 "hhfair"                  29    25; # This gets NaN constraints
 "ex8_6_2"                 30     0;
 "ramsey"                  33    22;
 "launch"                  38    28;
 "ex8_4_8_bnd"             42    30;
 "orth_d4m6_pl"            42    86;
 "orth_d3m6_pl"            42   127;
 "chenery"                 43    38;
 "kriging_peaks_full020"   46    44;
 "pricing050"              50     5;
 "polygon25"               50   324;
 "ex8_4_3"                 52    25;
 "ex8_2_1b"                57    33;
 "powerflow0009p"          60   139;
 "ex8_2_4b"                61    87;
 "ex8_4_7"                 62    40;
 "chakra"                  62    41;
 "minlphi"                 64    79; # This gets NaN constraints
 "kriging_peaks_full030"   66    64;
#  "elec25"                  75    25; # Objective length is 17303
 "ex8_6_1"                 75    45; # This gets NaN constraints
 "gsg_0001"                78   112;
 "lakes"                   90    78;
 "korcns"                  96    78; # This gets NaN constraints
 "ann_compressor_exp"      96    95;
 "etamac"                  97    70;
 "ann_peaks_exp"          100    98;
#  "polygon50"              100  1274; # Over 1000 constraints
 "chain50"                102    51;
 "otpop"                  103    76;
 "kriging_peaks_full050"  106   104;
 "ex8_3_13"               115    72;
 "kall_ellipsoids_tc02b"  124   128;
 "ex8_3_7"                126    92; # This gets NaN constraints
 "btest14"                135    93; # This gets NaN constraints
#  "elec50"                 150    50; # Objective length is 73418
 "hybriddynamic_varcc"    151   110;
 "kall_ellipsoids_tc03c"  193   196;
 "chain100"               202   101;
 "kriging_peaks_full100"  206   204;
#  "infeas1"                272  1614; # Objective length is 37958
 "camcns"                 280   243;
#  "elec100"                300   100; # Objective length is 311512
 "rocket50"               307   252;
#  "cesam2log"              316   165; # Objective length is 7402
 "gancns"                 357   274;
 "chain200"               402   201;
 "kriging_peaks_full200"  406   404;
 "kall_ellipsoids_tc05a"  464   461;
 "steenbrf"               468   108;
#  "arki0019"               510     2; # Objective length is 93706
#  "elec200"                600   200; # Objective length is 1281701
 "rocket100"              607   502;
 "glider50"               665   609;
 "ann_cumene_exp"         794   790;
#  "chain400"               802   401; # Objective length is 9921
]

good_examples_list = ["ex4_1_1", "ex4_1_2", "ex4_1_3", "ex4_1_4", "ex4_1_6", "ex4_1_7", 
               "mathopt5_4", "mathopt5_7", "mathopt5_8", "ex14_1_9", "ex4_1_5", 
               "ex4_1_8", "ex4_1_9", "ex8_1_3", "ex8_1_4", "ex8_1_5", "ex8_1_6", 
               "kriging_peaks_red010", "mathopt1", "mathopt2", "rbrock", "st_e19", 
               "ex14_1_1", "ex14_1_3", "ex14_1_8", "ex6_2_11", "ex6_2_6", "ex6_2_8", 
               "hs62", "prob09", "st_e06", "st_e11", "chance", "ex14_2_2", "ex14_2_5", 
               "ex14_2_8", "ex14_2_9", "ex6_1_2", "ex6_2_12", "ex6_2_14", "ex6_2_9", 
               "ex7_3_1", "st_e04", "st_e12", "st_e37", "st_e41", "ex14_2_1", 
               "ex14_2_4", "ex14_2_6", "ex8_1_7", "mhw4d", "ex14_1_2", 
               #"ex14_1_5", # Breaks Julia
               "ex14_2_3", "ex14_2_7", "ex6_1_4", "ex6_2_10", "ex6_2_13", "ex7_2_2", 
               "st_e21", "ex7_2_1", "ex6_1_1", "ex7_2_3", "ex7_2_4", "ex6_2_5", 
               "ex6_2_7", "alkylation", 
               #"ex14_1_7", # Breaks Julia
               "process", "st_e03", "chem", 
               "ex6_1_3", "ex7_3_4", "st_e16", "ex7_3_5", "alkyl", "ex8_4_6", "prob07", 
               "eq6_1", "ex5_4_3", "ex7_3_6", "ex8_4_2", "orth_d3m6", "kriging_peaks_full010", 
               "ex5_4_4", 
               # "maxmin", # Breaks Julia
               "ex8_6_2", "ramsey", "launch", "ex8_4_8_bnd", 
               "orth_d3m6_pl", "orth_d4m6_pl", "chenery", "kriging_peaks_full020", 
               "pricing050", "ex8_4_3", "ex8_2_1b", "ex8_2_4b", "chakra", "ex8_4_7", 
               "kriging_peaks_full030", "gsg_0001", "lakes", "ann_compressor_exp", 
               "etamac", "ann_peaks_exp", 
            #    "chain50", "otpop", "kriging_peaks_full050", "ex8_3_13", "kall_ellipsoids_tc02b", 
               # New examples with sin/cos:
               "ex8_1_2", "mathopt5_1", "mathopt5_2", 
               "mathopt5_3", "mathopt5_5", "mathopt6", "trig", "ex8_1_1", "mathopt4", "trigx", "ex14_1_4", 
               # "mathopt3", # Breaks Julia
               "inscribedsquare01", "inscribedsquare02", "polygon25", "powerflow0009p", 
               # "polygon50", # Too many constraints (based on limit set in paper)
]

included = original_included[in.(original_included[:,1], Ref(good_examples_list)),:]

for i in included[:,1]
    try
        include("./all_examples/$i.jl")
    catch
        println("$i failed")
    end
end

