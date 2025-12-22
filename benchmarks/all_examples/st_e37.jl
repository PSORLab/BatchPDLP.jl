# st_e37

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = (-1.9837 + x3 + x4)^2 + (x3 * exp(-x1) + x4 * exp(-x2) - 0.8393)^2 + (x3 * exp(-2 * x1) + x4 * exp(-2 * x2) - 0.4305)^2 + (x3 * exp(-3 * x1) + x4 * exp(-3 * x2) - 0.2441)^2 + (x3 * exp(-4 * x1) + x4 * exp(-4 * x2) - 0.1248)^2 + (x3 * exp(-5 * x1) + x4 * exp(-5 * x2) - 0.0981)^2 + (x3 * exp(-6 * x1) + x4 * exp(-6 * x2) - 0.0549)^2 + (x3 * exp(-7 * x1) + x4 * exp(-7 * x2) - 0.0174)^2 + (x3 * exp(-8 * x1) + x4 * exp(-8 * x2) - 0.0249)^2 + (x3 * exp(-9 * x1) + x4 * exp(-9 * x2) - 0.0154)^2 + (x3 * exp(-10 * x1) + x4 * exp(-10 * x2) - 0.0127)^2

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x1 - x2 - ( 0 ),
]

# BOUNDS
lvbs = [0, 0, 1, 1]
uvbs = [100, 100, 1, 1]

# Create the problem
st_e37 = Problem("st_e37", vars, obj, lvbs, uvbs, eq=eq, leq=leq)