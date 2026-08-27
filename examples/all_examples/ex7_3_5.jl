# ex7_3_5

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13]

# OBJECTIVE
obj = x4

# EQ (==0) CONSTRAINTS
eq = [
    x13 * x3^8 - x11 * x3^6 + x9 * x3^4 - x7 * x3^2 + x5 - ( 0 ),
    x12 * x3^6 - x10 * x3^4 + x8 * x3^2 - x6 - ( 0 ),
    -4.53 * x1^2 + x5 - ( 0 ),
    -5.28 * x1^2 - 0.364 * x1 + x6 - ( 0 ),
    -5.72 * x1^2 * x2 - 1.13 * x1^2 - 0.425 * x1 + x7 - ( 0 ),
    -6.93 * x1^2 * x2 - 0.0911 * x1 + x8 - ( 0.00422 ),
    -1.45 * x1^2 * x2 - 0.168 * x1 * x2 + x9 - ( 0.000338 ),
    -1.56 * x1^2 * x2^2 - 0.00084 * x1^2 * x2 - 0.0135 * x1 * x2 + x10 - ( 1.35e-05 ),
    -0.125 * x1^2 * x2^2 - 1.68e-05 * x1^2 * x2 - 0.000539 * x1 * x2 + x11 - ( 2.7e-07 ),
    -0.005 * x1^2 * x2^2 - 1.08e-05 * x1 * x2 + x12 - ( 0 ),
    -0.0001 * x1^2 * x2^2 + x13 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    -x1 - 0.145 * x4 - ( -0.175 ),
    x1 - 0.145 * x4 - ( 0.175 ),
    -x2 - 0.15 * x4 - ( -0.2 ),
    x2 - 0.15 * x4 - ( 0.2 ),
]

# BOUNDS
lvbs = [-Inf, -Inf, 0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, 10, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
ex7_3_5 = Problem("ex7_3_5", vars, obj, lvbs, uvbs, eq=eq, leq=leq)