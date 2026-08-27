# ex7_2_2

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = -x4

# EQ (==0) CONSTRAINTS
eq = [
    0.09755988 * x1 * x5 + x1 - ( 1 ),
    0.0965842812 * x2 * x6 + x2 - x1 - ( 0 ),
    0.0391908 * x3 * x5 + x3 + x1 - ( 1 ),
    0.03527172 * x4 * x6 + x4 - x1 + x2 - x3 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x5^0.5 + x6^0.5 - ( 4 ),
]

# BOUNDS
lvbs = [0, 0, 0, 0, 1e-05, 1e-05]
uvbs = [1, 1, 1, 1, 16, 16]

# Create the problem
ex7_2_2 = Problem("ex7_2_2", vars, obj, lvbs, uvbs, eq=eq, leq=leq)