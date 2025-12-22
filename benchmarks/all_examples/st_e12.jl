# st_e12

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = x1^0.6 + x2^0.6 - 6 * x1 - 4 * x3 + 3 * x4

# EQ (==0) CONSTRAINTS
eq = [
    -3 * x1 + x2 - 3 * x3 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x1 + 2 * x3 - ( 4 ),
    x2 + 2 * x4 - ( 4 ),
]

# BOUNDS
lvbs = [0, 0, 0, 0]
uvbs = [3, 4, 2, 1]

# Create the problem
st_e12 = Problem("st_e12", vars, obj, lvbs, uvbs, eq=eq, leq=leq)