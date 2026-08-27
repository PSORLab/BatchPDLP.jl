# st_e21

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = x1^0.6 + x2^0.6 + x3^0.4 - 4 * x3 + 2 * x4 + 5 * x5 - x6

# EQ (==0) CONSTRAINTS
eq = [
    -3 * x1 + x2 - 3 * x4 - ( 0 ),
    -2 * x2 + x3 - 2 * x5 - ( 0 ),
    4 * x4 - x6 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x1 + 2 * x4 - ( 4 ),
    x2 + x5 - ( 4 ),
    x3 + x6 - ( 6 ),
]

# BOUNDS
lvbs = [0, 0, 0, 0, 0, 0]
uvbs = [3, 4, 4, 2, 2, 6]

# Create the problem
st_e21 = Problem("st_e21", vars, obj, lvbs, uvbs, eq=eq, leq=leq)