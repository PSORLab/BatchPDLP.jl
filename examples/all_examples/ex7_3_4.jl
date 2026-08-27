# ex7_3_4

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12]

# OBJECTIVE
obj = x12

# EQ (==0) CONSTRAINTS
eq = [
    x10 * x11^4 - x8 * x11^2 + x6 - ( 0 ),
    x9 * x11^2 - x7 - ( 0 ),
    -54.387 * x3 * x2 + x6 - ( 0 ),
    -0.2 * (1364.67 * x3 * x2 - 147.15 * x4 * x3 * x2) + 5.544 * x5 + x7 - ( 0 ),
    -3 * (-9.81 * x2^2 * x3 - 9.81 * x3 * x1 * x2 - 4.312 * x3^2 * x2 + 264.896 * x3 * x2 + x4 * x5 - 9.274 * x5) + x8 - ( 0 ),
    -7 * x3^2 * x4 * x2 + 64.918 * x3^2 * x2 - 380.067 * x3 * x2 - 3 * x5 * x2 - 3 * x5 * x1 + x9 - ( 0 ),
    -x2 * x3^2 * (7 * x1 + 4 * x2) + x10 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    -x1 - x12 - ( -10 ),
    x1 - x12 - ( 10 ),
    x2 - 0.1 * x12 - ( 1 ),
    -x2 - 0.1 * x12 - ( -1 ),
    -x3 - 0.1 * x12 - ( -1 ),
    x3 - 0.1 * x12 - ( 1 ),
    -x4 - 0.01 * x12 - ( -0.2 ),
    x4 - 0.01 * x12 - ( 0.2 ),
    -x5 - 0.005 * x12 - ( -0.05 ),
    x5 - 0.005 * x12 - ( 0.05 ),
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 10, Inf]

# Create the problem
ex7_3_4 = Problem("ex7_3_4", vars, obj, lvbs, uvbs, eq=eq, leq=leq)