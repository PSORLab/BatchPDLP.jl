# chance

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = 24.55 * x1 + 26.75 * x2 + 39 * x3 + 40.5 * x4

# EQ (==0) CONSTRAINTS
eq = [
    x1 + x2 + x3 + x4 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    1.645 * sqrt(0.28 * x1^2 + 0.19 * x2^2 + 20.5 * x3^2 + 0.62 * x4^2) + 12 * x1 + 11.9 * x2 + 41.8 * x3 + 52.1 * x4 - ( -21 ),
    -2.3 * x1 + 5.6 * x2 + 11.1 * x3 + 1.3 * x4 - ( -5 ),
]

# BOUNDS
lvbs = [0, 0, 0, 0]
uvbs = [Inf, Inf, Inf, Inf]

# Create the problem
chance = Problem("chance", vars, obj, lvbs, uvbs, eq=eq, leq=leq)