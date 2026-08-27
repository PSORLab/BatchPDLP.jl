# ex8_5_5

# VARIABLES
vars = [x1, x2, x3, x4, x5]

# OBJECTIVE
obj = x1 * log(x1) + x2 * log(x2) - log(x3 - x5) + x3 - 0.353553390593274 * log((x3 + 2.41421356237309 * x5) / (x3 - 0.414213562373095 * x5)) * x4 / x5 + 2.5746329124341 * x1 + 0.54639755131421 * x2 - 1

# EQ (==0) CONSTRAINTS
eq = [
    x3^3 - x3^2 * (1 - x5) + x3 * (-3 * x5^2 - 2 * x5 + x4) - x4 * x5 + x5^3 + x5^2 - ( 0 ),
    -0.884831 * x1 * x1 - 0.555442 * x1 * x2 - 0.555442 * x2 * x1 - 0.427888 * x2 * x2 + x4 - ( 0 ),
    -0.0885973 * x1 - 0.0890893 * x2 + x5 - ( 0 ),
    x1 + x2 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf]

# Create the problem
ex8_5_5 = Problem("ex8_5_5", vars, obj, lvbs, uvbs, eq=eq, leq=leq)