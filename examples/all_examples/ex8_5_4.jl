# ex8_5_4

# VARIABLES
vars = [x1, x2, x3, x4, x5]

# OBJECTIVE
obj = x1 * log(x1) + x2 * log(x2) - log(x3 - x5) + x3 - x4 * log(x5 / x3 + 1) / x5 + 0.362259780811985 * x1 + 3.27527428318836 * x2 - 1

# EQ (==0) CONSTRAINTS
eq = [
    x3^3 - x3^2 + x3 * (-x5^2 - x5 + x4) - x4 * x5 - ( 0 ),
    -0.352565 * x1 * x1 - 0.844083 * x1 * x2 - 0.844083 * x2 * x1 - 2.14335 * x2 * x2 + x4 - ( 0 ),
    -0.12932 * x1 - 0.271567 * x2 + x5 - ( 0 ),
    x1 + x2 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf]

# Create the problem
ex8_5_4 = Problem("ex8_5_4", vars, obj, lvbs, uvbs, eq=eq, leq=leq)