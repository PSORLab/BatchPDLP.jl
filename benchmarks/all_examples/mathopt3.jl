# mathopt3

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6]

# OBJECTIVE
obj = (sqr(x1 + x2) + sqr(x3 - x5) + sqr(x6 - x4) + 2*sqr(x1 + x3 - x4) + sqr(
     x2 - x1 + x3 - x4) + 10*sqr(sin(x1 + x5 - x6)))

# EQ (==0) CONSTRAINTS
eq = [
    sqr(x1) - sin(x2) - x4 + x5 + x6;
    x1*x3 - x2*x4*x1 - sin((-x1) - x3 + x6) - x5;
    x2*x6*cos(x5) - sin(x3*x4) + x2 - x5;
    x1*x2 - sqr(x3) - x4*x5 - sqr(x6);
]

# LEQ (<=0) CONSTRAINTS
leq = [
    2*x1 + 5*x2 + x3 + x4 - 1;
    3*x1 - 2*x2 + x3 - 4*x4;
    x1 + x2 + x3 + x4 + x5 + x6 - 2;
]

# GEQ (>=0) CONSTRAINTS
geq = []


# BOUNDS
lvbs = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
mathopt3 = Problem("mathopt3", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)