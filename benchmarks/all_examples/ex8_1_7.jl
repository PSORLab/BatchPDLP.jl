# ex8_1_7

# VARIABLES
vars = [x1, x2, x3, x4, x5]

# OBJECTIVE
obj = (-1 + x1)^2 + (x1 - x2)^2 + (x2 - x3)^3 + (x3 - x4)^4 + (x4 - x5)^4

# EQ (==0) CONSTRAINTS
eq = [
    0.5 * x1 * x5 + 0.5 * x1 * x5 - ( 2 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x2^2 + x3^3 + x1 - ( 6.24264068711929 ),
    -x3^3 - x2^2 - x1 - ( -6.24264068711929 ),
    -x3^2 + x2 + x4 - ( 0.82842712474619 ),
    x3^2 - x2 - x4 - ( -0.82842712474619 ),
]

# BOUNDS
lvbs = [-5, -5, -5, -5, -5]
uvbs = [5, 5, 5, 5, 5]

# Create the problem
ex8_1_7 = Problem("ex8_1_7", vars, obj, lvbs, uvbs, eq=eq, leq=leq)