# ex8_5_6

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = x1 * log(x1) + x2 * log(x2) + x3 * log(x3) - log(x4 - x6) + x4 - 0.353553390593274 * log((x4 + 2.41421356237309 * x6) / (x4 - 0.414213562373095 * x6)) * x5 / x6 + 1.42876598488588 * x1 + 1.27098480432594 * x2 + 1.62663700075562 * x3 - 1

# EQ (==0) CONSTRAINTS
eq = [
    x4^3 - x4^2 * (1 - x6) + x4 * (-3 * x6^2 - 2 * x6 + x5) - x5 * x6 + x6^3 + x6^2 - ( 0 ),
    -0.142724 * x1 * x1 - 0.206577 * x1 * x2 - 0.342119 * x1 * x3 - 0.206577 * x2 * x1 - 0.323084 * x2 * x2 - 0.547748 * x2 * x3 - 0.342119 * x3 * x1 - 0.547748 * x3 * x2 - 0.968906 * x3 * x3 + x5 - ( 0 ),
    -0.0815247 * x1 - 0.0907391 * x2 - 0.13705 * x3 + x6 - ( 0 ),
    x1 + x2 + x3 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
ex8_5_6 = Problem("ex8_5_6", vars, obj, lvbs, uvbs, eq=eq, leq=leq)