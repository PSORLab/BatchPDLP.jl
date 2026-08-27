# ex7_2_3

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8]

# OBJECTIVE
obj = x1 + x2 + x3

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    833.33252 * x4 / x1 / x6 + 100 / x6 - 83333.333 / (x1 * x6) - ( 1 ),
    1250 * x5 / x2 / x7 + x4 / x7 - 1250 * x4 / x2 / x7 - ( 1 ),
    1.25e+06 / (x3 * x8) + x5 / x8 - 2500 * x5 / x3 / x8 - ( 1 ),
    0.0025 * x4 + 0.0025 * x6 - ( 1 ),
    -0.0025 * x4 + 0.0025 * x5 + 0.0025 * x7 - ( 1 ),
    -0.01 * x5 + 0.01 * x8 - ( 1 ),
]

# BOUNDS
lvbs = [100, 1000, 1000, 10, 10, 10, 10, 10]
uvbs = [1e+04, 1e+04, 1e+04, 1000, 1000, 1000, 1000, 1000]

# Create the problem
ex7_2_3 = Problem("ex7_2_3", vars, obj, lvbs, uvbs, eq=eq, leq=leq)