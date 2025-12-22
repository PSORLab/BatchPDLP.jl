# mathopt4

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
obj = (sqr(2*sqr(x1) - sqr(x2)) + sqr(x2 - 6*sqr(x1)))

# EQ (==0) CONSTRAINTS
eq = [
    x1 - (100*sin(2*x1 + 3*x2) + 10*x2)
]

# LEQ (<=0) CONSTRAINTS
leq = [
    x1 + x2 - 2
]

# GEQ (>=0) CONSTRAINTS
geq = []


# BOUNDS
lvbs = [-10.0, -10.0]
uvbs = [10.0, 10.0]

# Create the problem
mathopt4 = Problem("mathopt4", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)