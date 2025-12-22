# hhfair

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29]

# OBJECTIVE
obj = -x24 * x25^0.944 * x26^0.891136

# EQ (==0) CONSTRAINTS
eq = [
    -0.01 * (0.5 * x5^0.5 + 0.5 * (1004.72366 - x8 - x15)^0.5)^2 + x24 - ( 0 ),
    -0.01 * (0.5 * x6^0.5 + 0.5 * (1004.72366 - x9 - x16)^0.5)^2 + x25 - ( 0 ),
    -0.01 * (0.5 * x7^0.5 + 0.5 * (1004.72366 - x10 - x17)^0.5)^2 + x26 - ( 0 ),
    -0.07 * x2 - x8 + x27 - ( 0 ),
    -0.07 * x3 - x9 + x28 - ( 0 ),
    -0.07 * x4 - x10 + x29 - ( 0 ),
    x21 - 0.2 * x27 - ( 0 ),
    x22 - 0.2 * x28 - ( 0 ),
    x23 - 0.2 * x29 - ( 0 ),
    x5 + x18 + x21 - x27 - ( 0 ),
    x6 + x19 + x22 - x28 - ( 0 ),
    x7 + x20 + x23 - x29 - ( 0 ),
    x1 - x2 + x11 - x12 + x18 - ( 0 ),
    x2 - x3 + x12 - x13 + x19 - ( 0 ),
    x3 - x4 + x13 - x14 + x20 - ( 0 ),
    x15 * (-0.255905 * x5 + x12) - ( 1 ),
    x16 * (-0.255905 * x6 + x13) - ( 1 ),
    x17 * (-0.255905 * x7 + x14) - ( 1 ),
    x4 + x14 - ( 1100 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    0.25846405 * x5 - x12 - ( 0 ),
    0.25846405 * x6 - x13 - ( 0 ),
    0.25846405 * x7 - x14 - ( 0 ),
    x8 + x15 - ( 904.251294 ),
    x9 + x16 - ( 904.251294 ),
    x10 + x17 - ( 904.251294 ),
]

# BOUNDS
lvbs = [1000, -Inf, -Inf, -Inf, 100, 100, 100, 100, 100, 100, 100, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0.01, 0.01, 0.01, -Inf, -Inf, -Inf]
uvbs = [1000, Inf, Inf, Inf, Inf, Inf, Inf, 400, 400, 400, 100, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
hhfair = Problem("hhfair", vars, obj, lvbs, uvbs, eq=eq, leq=leq)