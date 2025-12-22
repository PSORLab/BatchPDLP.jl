# ex14_1_5

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = x6

# EQ (==0) CONSTRAINTS
eq = [
    2 * x1 + x2 + x3 + x4 + x5 - ( 6 ),
    x1 + 2 * x2 + x3 + x4 + x5 - ( 6 ),
    x1 + x2 + 2 * x3 + x4 + x5 - ( 6 ),
    x1 + x2 + x3 + 2 * x4 + x5 - ( 6 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x1 * x2 * x3 * x4 * x5 - x6 - ( 1 ),
    -x1 * x2 * x3 * x4 * x5 - x6 - ( -1 ),
]

# BOUNDS
lvbs = [-2, -2, -2, -2, -2, -Inf]
uvbs = [2, 2, 2, 2, 2, Inf]

# Create the problem
ex14_1_5 = Problem("ex14_1_5", vars, obj, lvbs, uvbs, eq=eq, leq=leq)