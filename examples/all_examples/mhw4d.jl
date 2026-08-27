# mhw4d

# VARIABLES
vars = [x1, x2, x3, x4, x5]

# OBJECTIVE
obj = (-1 + x1)^2 + (x1 - x2)^2 + (x2 - x3)^3 + (x3 - x4)^4 + (x4 - x5)^4

# EQ (==0) CONSTRAINTS
eq = [
    x2^2 + x3^3 + x1 - ( 6.24264068711929 ),
    -x3^2 + x2 + x4 - ( 0.82842712474619 ),
    x1 * x5 - ( 2 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf]

# Create the problem
mhw4d = Problem("mhw4d", vars, obj, lvbs, uvbs, eq=eq, leq=leq)