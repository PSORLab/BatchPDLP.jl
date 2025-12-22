# ex8_4_6

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14]

# OBJECTIVE
obj = ((-0.1622 + x1) / x1)^2 + ((-0.6791 + x2) / x2)^2 + ((-0.679 + x3) / x3)^2 + ((-0.3875 + x4) / x4)^2 + ((-0.1822 + x5) / x5)^2 + ((-0.1249 + x6) / x6)^2 + ((-0.0857 + x7) / x7)^2 + ((-0.0616 + x8) / x8)^2

# EQ (==0) CONSTRAINTS
eq = [
    x9 * exp(-4 * x12) + x10 * exp(-4 * x13) + x11 * exp(-4 * x14) - x1 - ( 0 ),
    x9 * exp(-8 * x12) + x10 * exp(-8 * x13) + x11 * exp(-8 * x14) - x2 - ( 0 ),
    x9 * exp(-12 * x12) + x10 * exp(-12 * x13) + x11 * exp(-12 * x14) - x3 - ( 0 ),
    x9 * exp(-24 * x12) + x10 * exp(-24 * x13) + x11 * exp(-24 * x14) - x4 - ( 0 ),
    x9 * exp(-48 * x12) + x10 * exp(-48 * x13) + x11 * exp(-48 * x14) - x5 - ( 0 ),
    x9 * exp(-72 * x12) + x10 * exp(-72 * x13) + x11 * exp(-72 * x14) - x6 - ( 0 ),
    x9 * exp(-94 * x12) + x10 * exp(-94 * x13) + x11 * exp(-94 * x14) - x7 - ( 0 ),
    x9 * exp(-118 * x12) + x10 * exp(-118 * x13) + x11 * exp(-118 * x14) - x8 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [0, 0, 0, 0, 0, 0, 0, 0, -10, -10, -10, 0, 0, 0]
uvbs = [1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 10, 0.5, 0.5, 0.5]

# Create the problem
ex8_4_6 = Problem("ex8_4_6", vars, obj, lvbs, uvbs, eq=eq, leq=leq)