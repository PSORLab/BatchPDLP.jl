# st_e03

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]

# OBJECTIVE
obj = -0.063 * x4 * x7 + 5.04 * x1 + 0.035 * x2 + 10 * x3 + 3.36 * x5

# EQ (==0) CONSTRAINTS
eq = [
    x1 - 1.22 * x4 + x5 - ( 0 ),
    x9 + 0.222 * x10 - ( 35.82 ),
    3 * x7 - x10 - ( 133 ),
    0.038 * x8^2 - 1.098 * x8 - 0.325 * x6 + x7 - ( 57.425 ),
    x4 * x9 * x6 + 1000 * x3 * x6 - 98000 * x3 - ( 0 ),
    -x1 * x8 + x2 + x5 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    -0.13167 * x8 * x1 - 1.12 * x1 + 0.00667 * x8^2 * x1 + x4 - ( 0 ),
]

# BOUNDS
lvbs = [1, 1, 0, 1, 0, 85, 90, 3, 1.2, 145]
uvbs = [2000, 16000, 120, 5000, 2000, 93, 95, 12, 4, 162]

# Create the problem
st_e03 = Problem("st_e03", vars, obj, lvbs, uvbs, eq=eq, leq=leq)