# ex8_5_3

# VARIABLES
vars = [x1, x2, x3, x4, x5]

# OBJECTIVE
obj = x1 * log(x1) + x2 * log(x2) - log(x3 - x5) + x3 - x4 * log(x5 / x3 + 1) / x5 + 5.0464317551216 * x1 + 0.366877055769689 * x2 - 1

# EQ (==0) CONSTRAINTS
eq = [
    x3^3 - x3^2 + x3 * (-x5^2 - x5 + x4) - x4 * x5 - ( 0 ),
    -1.04633 * x1 * x1 - 0.579822 * x1 * x2 - 0.579822 * x2 * x1 - 0.379615 * x2 * x2 + x4 - ( 0 ),
    -0.0771517 * x1 - 0.0765784 * x2 + x5 - ( 0 ),
    x1 + x2 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf]

# Create the problem
ex8_5_3 = Problem("ex8_5_3", vars, obj, lvbs, uvbs, eq=eq, leq=leq)