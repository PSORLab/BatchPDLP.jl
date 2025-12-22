# st_e41

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = 200 * x1^0.6 + 200 * x2^0.6 + 200 * x3^0.6 + 300 * x4^0.6

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x3 * (1 - x1)^2 * (1 - x4)^2 + (-x2 * (-(1 - x1) * (1 - x4) + 1) + 1)^2 * (1 - x3) - ( 0.1 ),
    -x3 * (1 - x1)^2 * (1 - x4)^2 - (-x2 * (-(1 - x1) * (1 - x4) + 1) + 1)^2 * (1 - x3) - ( 0 ),
]

# BOUNDS
lvbs = [0.5, 0.5, 0.5, 0.5]
uvbs = [1, 1, 1, 1]

# Create the problem
st_e41 = Problem("st_e41", vars, obj, lvbs, uvbs, eq=eq, leq=leq)