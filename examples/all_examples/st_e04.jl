# st_e04

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = 400 * x1^0.9 + 22 * (-14.7 + x2)^1.2 + x3 + 1000

# EQ (==0) CONSTRAINTS
eq = [
    -exp(-3950 / (460 + x4) + 11.86) + x2 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    -x3 * x1 - 144 * x4 - ( 11520 ),
]

# BOUNDS
lvbs = [0, 14.7, 0, -459.67]
uvbs = [15.1, 94.2, 5371, 80]

# Create the problem
st_e04 = Problem("st_e04", vars, obj, lvbs, uvbs, eq=eq, leq=leq)