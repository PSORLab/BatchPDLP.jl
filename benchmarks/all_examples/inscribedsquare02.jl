# inscribedsquare02

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8]

# OBJECTIVE
obj = (sqr(x7) + sqr(x8))

# EQ (==0) CONSTRAINTS
eq = [
    sin(x1)*cos(x1 - x1*x1) - x5;
    sin(x1)*x1 - x6;
    sin(x2)*cos(x2 - x2*x2) - x5 - x7;
    sin(x2)*x2 - x6 - x8;
    sin(x3)*cos(x3 - x3*x3) - x5 + x8;
    sin(x3)*x3 - x6 - x7;
    sin(x4)*cos(x4 - x4*x4) - x5 - x7 + x8;
    sin(x4)*x4 - x6 - x7 - x8;
]

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-3.14159265358979, -3.14159265358979, -3.14159265358979, -3.14159265358979, -Inf, -Inf, 0.0, 0.0]
uvbs = [3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, Inf, Inf, Inf, Inf]

# Create the problem
inscribedsquare02 = Problem("inscribedsquare02", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)