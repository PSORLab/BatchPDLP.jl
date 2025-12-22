# ex7_3_1

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = x4

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    10 * x2^2 * x3^3 + 10 * x2^3 * x3^2 + 200 * x2^2 * x3^2 + 100 * x2^3 * x3 + 100 * x3^3 * x2 + x1 * x2 * x3^2 + x1 * x2^2 * x3 + 1000 * x3^2 * x2 + 8 * x3^2 * x1 + 1000 * x2^2 * x3 + 8 * x2^2 * x1 + 6 * x1 * x2 * x3 - x1^2 + 60 * x1 * x3 + 60 * x1 * x2 - 200 * x1 - ( 0 ),
    -x1 - 800 * x4 - ( -800 ),
    x1 - 800 * x4 - ( 800 ),
    -x2 - 2 * x4 - ( -4 ),
    x2 - 2 * x4 - ( 4 ),
    -x3 - 3 * x4 - ( -6 ),
    x3 - 3 * x4 - ( 6 ),
]

# BOUNDS
lvbs = [0, 0, 0, 0]
uvbs = [Inf, Inf, Inf, Inf]

# Create the problem
ex7_3_1 = Problem("ex7_3_1", vars, obj, lvbs, uvbs, eq=eq, leq=leq)