# ex7_3_2

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = x4

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x1^4 * x2^4 - x1^4 - x3 * x2^4 - ( 0 ),
    -x1 - 0.25 * x4 - ( -1.4 ),
    x1 - 0.25 * x4 - ( 1.4 ),
    -x2 - 0.2 * x4 - ( -1.5 ),
    x2 - 0.2 * x4 - ( 1.5 ),
    -x3 - 0.2 * x4 - ( -0.8 ),
    x3 - 0.2 * x4 - ( 0.8 ),
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf]

# Create the problem
ex7_3_2 = Problem("ex7_3_2", vars, obj, lvbs, uvbs, eq=eq, leq=leq)