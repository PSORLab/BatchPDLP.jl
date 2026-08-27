# ex7_2_4

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8]

# OBJECTIVE
obj = 0.4 * x1^0.67 / x7^0.67 + 0.4 * x2^0.67 / x8^0.67 - x1 - x2 + 10

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    0.0588 * x5 * x7 + 0.1 * x1 - ( 1 ),
    0.0588 * x6 * x8 + 0.1 * x1 + 0.1 * x2 - ( 1 ),
    4 * x3 / x5 + 2 / (x5 * x3^0.71) + (0.0588 * x7) / x3^1.3 - ( 1 ),
    4 * x4 / x6 + 2 / (x6 * x4^0.71) + 0.0588 * x4^1.3 * x8 - ( 1 ),
]

# BOUNDS
lvbs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
uvbs = [10, 10, 10, 10, 10, 10, 10, 10]

# Create the problem
ex7_2_4 = Problem("ex7_2_4", vars, obj, lvbs, uvbs, eq=eq, leq=leq)