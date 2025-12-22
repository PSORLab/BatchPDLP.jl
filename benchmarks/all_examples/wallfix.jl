# wallfix

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = x1

# EQ (==0) CONSTRAINTS
eq = [
    x1 * x2 - ( 1 ),
    x3 / (x4 * x1) - ( 4.8 ),
    x5 / (x6 * x2) - ( 0.98 ),
    x6 * x4 - ( 1 ),
    x1 - x2 + 1e-07 * x3 - 1e-05 * x5 - ( 0 ),
    2 * x1 - 2 * x2 + 1e-07 * x3 - 0.01 * x4 - 1e-05 * x5 + 0.01 * x6 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, 0, 0, 0, 0, 0]
uvbs = [Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
wallfix = Problem("wallfix", vars, obj, lvbs, uvbs, eq=eq, leq=leq)