
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
 "ex8_1_3"                  2     0; # Note: Solvers disagree on solutions
 "ex8_1_4"                  2     0; # Note: Solvers disagree on solutions
 "ex8_1_5"                  2     0; # Note: Solvers disagree on solutions (lots of "large complementary slackness residual" warnings)
 "ex8_1_6"                  2     0;
 "kriging_peaks_red010"     2     0;
 "kriging_peaks_red020"     2     0; # Exclude: Objective length is 9578
 "kriging_peaks_red030"     2     0; # Exclude: Objective length is 14362
 "kriging_peaks_red050"     2     0; # Exclude: Objective length is 23762
 "kriging_peaks_red100"     2     0; # Exclude: Objective length is 47766
 "kriging_peaks_red200"     2     0; # Exclude: Objective length is 95485
 "kriging_peaks_red500"     2     0; # Exclude: Objective length is 238721
 "rbrock"                   2     0;
 "ex4_1_8"                  2     1;
 "ex4_1_9"                  2     2;
 "ex14_1_9"                 2     2;
 "mathopt1"                 2     2;
 "mathopt4"                 2     2;
 "st_e19"                   2     2;
 "trigx"                    2     2;
 "mathopt2"                 2     4;
 "least"                    3     0; # Run separately: Crashes Julia when run
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
 "ex7_3_2"                  4     7; # Run separately: Crashes Julia when run
 "mhw4d"                    5     3;
 "ex8_5_3"                  5     4; # Exclude: This gets NaN constraints
 "ex8_5_4"                  5     4; # Exclude: This gets NaN constraints
 "ex8_5_5"                  5     4; # Exclude: This gets NaN constraints
 "ex8_1_7"                  5     5;
 "ex14_2_1"                 5     7;
 "ex14_2_4"                 5     7;
 "ex14_2_6"                 5     7;
 "ex6_2_10"                 6     3;
 "ex6_2_13"                 6     3;
 "ex6_1_4"                  6     4;
 "ex8_5_1"                  6     4; # Exclude: This gets NaN constraints
 "ex8_5_2"                  6     4; # Exclude: This gets NaN constraints
 "ex8_5_6"                  6     4; # Exclude: This gets NaN constraints
 "ex7_2_2"                  6     5;
 "ex14_1_5"                 6     6; # Run separately: Crashes Julia when run
 "st_e21"                   6     6;
 "wallfix"                  6     6; # Exclude: This gets NaN constraints
 "mathopt3"                 6     7; # Run separately: Crashes Julia when run
 "ex14_1_2"                 6     9;
 "ex14_2_3"                 6     9;
 "ex14_2_7"                 6     9;
 "ex7_2_1"                  7    14;
 "ex7_2_4"                  8     4;
 "ex6_1_1"                  8     6; # Note: Took a VERY long time in testing (1-2 hours?)
 "ex7_2_3"                  8     6;
 "inscribedsquare01"        8     8;
 "inscribedsquare02"        8     8;
 "ex6_2_5"                  9     3;
 "ex6_2_7"                  9     3;
 "like"                     9     3; # Exclude: This gets NaN constraints
 "process"                 10     7;
 "st_e03"                  10     7;
 "alkylation"              10    11;
 "ex14_1_7"                10    17; # Run separately: Crashes Julia when run
 "shiporig"                10    17; # Exclude: This gets NaN constraints
 "chem"                    11     4;
 "ann_fermentation_exp"    12     9; # Exclude: This gets NaN constraints
 "ex6_1_3"                 12     9;
 "st_e16"                  12     9;
 "ex7_3_4"                 12    17; # Note: Took a VERY long time in testing (overnight, 10+ hours?)
 "ex7_3_5"                 13    15;
 "alkyl"                   14     7;
 "ex8_4_6"                 14     8;
 "ex8_4_5"                 15    11; # Exclude: This gets NaN constraints
 "prob07"                  15    36;
 "ex5_4_3"                 16    13;
 "eq6_1"                   16    60;
 "ex7_3_6"                 17    17;
 "ex8_4_2"                 24    10;
 "orth_d3m6"               25    62;
 "kriging_peaks_full010"   26    24;
 "ex5_4_4"                 27    19;
 "maxmin"                  27    78; # Run separately: Crashes Julia when run
 "hhfair"                  29    25; # Exclude: This gets NaN constraints
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
 "minlphi"                 64    79; # Exclude: This gets NaN constraints
 "kriging_peaks_full030"   66    64;
 "elec25"                  75    25; # Exclude: Objective length is 17303
 "ex8_6_1"                 75    45; # Exclude: This gets NaN constraints
 "gsg_0001"                78   112;
 "lakes"                   90    78;
 "korcns"                  96    78; # Exclude: This gets NaN constraints
 "ann_compressor_exp"      96    95;
 "etamac"                  97    70;
 "ann_peaks_exp"          100    98;
 "polygon50"              100  1274; # Exclude: Over 1000 constraints
 "chain50"                102    51; # Exclude: Over 100 variables
 "otpop"                  103    76; # Exclude: Over 100 variables
 "kriging_peaks_full050"  106   104; # Exclude: Over 100 variables
 "ex8_3_13"               115    72; # Exclude: Over 100 variables
 "kall_ellipsoids_tc02b"  124   128; # Exclude: Over 100 variables
 "ex8_3_7"                126    92; # Exclude: This gets NaN constraints
 "btest14"                135    93; # Exclude: This gets NaN constraints
 "elec50"                 150    50; # Exclude: Objective length is 73418
 "hybriddynamic_varcc"    151   110; # Exclude: Over 100 variables
 "kall_ellipsoids_tc03c"  193   196; # Exclude: Over 100 variables
 "chain100"               202   101; # Exclude: Over 100 variables
 "kriging_peaks_full100"  206   204; # Exclude: Over 100 variables
 "infeas1"                272  1614; # Exclude: Objective length is 37958
 "camcns"                 280   243; # Exclude: Over 100 variables
 "elec100"                300   100; # Exclude: Objective length is 311512
 "rocket50"               307   252; # Exclude: Over 100 variables
 "cesam2log"              316   165; # Exclude: Objective length is 7402
 "gancns"                 357   274; # Exclude: Over 100 variables
 "chain200"               402   201; # Exclude: Over 100 variables
 "kriging_peaks_full200"  406   404; # Exclude: Over 100 variables
 "kall_ellipsoids_tc05a"  464   461; # Exclude: Over 100 variables
 "steenbrf"               468   108; # Exclude: Over 100 variables
 "arki0019"               510     2; # Exclude: Objective length is 93706
 "elec200"                600   200; # Exclude: Objective length is 1281701
 "rocket100"              607   502; # Exclude: Over 100 variables
 "glider50"               665   609; # Exclude: Over 100 variables
 "ann_cumene_exp"         794   790; # Exclude: Over 100 variables
 "chain400"               802   401; # Exclude: Objective length is 9921
]

exclude_list = ["kriging_peaks_red020", "kriging_peaks_red030", "kriging_peaks_red050", 
                "kriging_peaks_red100", "kriging_peaks_red200", "kriging_peaks_red500", 
                "least", "ex7_3_2", "ex8_5_3", "ex8_5_4", "ex8_5_5", "ex8_5_1", "ex8_5_2", 
                "ex8_5_6", "ex14_1_5", "wallfix", "mathopt3", "like", "ex14_1_7", "shiporig", 
                "ann_fermentation_exp", "ex8_4_5", "maxmin", "hhfair", "minlphi", 
                "elec25", "ex8_6_1", "korcns", "polygon50", "chain50", "otpop", 
                "kriging_peaks_full050", "ex8_3_13", "kall_ellipsoids_tc02b", 
                "ex8_3_7", "btest14", "elec50", "hybriddynamic_varcc", 
                "kall_ellipsoids_tc03c", "chain100", "kriging_peaks_full100", "infeas1", 
                "camcns", "elec100", "rocket50", "cesam2log", "gancns", "chain200", 
                "kriging_peaks_full200", "kall_ellipsoids_tc05a", "steenbrf", "arki0019", 
                "elec200", "rocket100", "glider50", "ann_cumene_exp", "chain400"]

included = original_included[(!).(in.(original_included[:,1], Ref(exclude_list))),:]

for i in included[:,1]
    try
        include("./all_examples/$i.jl")
    catch
        println("$i failed")
    end
end

