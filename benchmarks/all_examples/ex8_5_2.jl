# ex8_5_2

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = x1 * log(x1) + x2 * log(x2) + x3 * log(x3) + x6 / (x4 - x6) - log(x4 - x6) - 2 * x5 / x4 + 0.585616681390832 * x1 + 3.53797016206289 * x2 + 2.18345516206289 * x3

# EQ (==0) CONSTRAINTS
eq = [
    x4^3 - x4^2 * (1 + x6) + x5 * x4 - x5 * x6 - ( 0 ),
    -0.37943 * x1 * x1 - 0.75885 * x1 * x2 - 0.48991 * x1 * x3 - 0.75885 * x2 * x1 - 0.8836 * x2 * x2 - 0.23612 * x2 * x3 - 0.48991 * x3 * x1 - 0.23612 * x3 * x2 - 0.63263 * x3 * x3 + x5 - ( 0 ),
    -0.14998 * x1 - 0.14998 * x2 - 0.14998 * x3 + x6 - ( 0 ),
    x1 + x2 + x3 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
ex8_5_2 = Problem("ex8_5_2", vars, obj, lvbs, uvbs, eq=eq, leq=leq)