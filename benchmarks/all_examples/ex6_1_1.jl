# ex6_1_1

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8]

# OBJECTIVE
obj = x1 * (log(x1) - log(x1 + x3)) + x3 * (log(x3) - log(x1 + x3)) + x2 * (log(x2) - log(x2 + x4)) + x4 * (log(x4) - log(x2 + x4)) + 0.925356626778358 * x1 * x7 + 0.746014540096753 * x3 * x5 + 0.925356626778358 * x2 * x8 + 0.746014540096753 * x4 * x6

# EQ (==0) CONSTRAINTS
eq = [
    x5 * (x1 + 0.159040857374844 * x3) - x1 - ( 0 ),
    x6 * (x2 + 0.159040857374844 * x4) - x2 - ( 0 ),
    x7 * (0.307941026821595 * x1 + x3) - x3 - ( 0 ),
    x8 * (0.307941026821595 * x2 + x4) - x4 - ( 0 ),
    x1 + x2 - ( 0.5 ),
    x3 + x4 - ( 0.5 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [1e-07, 1e-07, 1e-07, 1e-07, 0, 0, 0, 0]
uvbs = [0.5, 0.5, 0.5, 0.5, Inf, Inf, Inf, Inf]

# Create the problem
ex6_1_1 = Problem("ex6_1_1", vars, obj, lvbs, uvbs, eq=eq, leq=leq)