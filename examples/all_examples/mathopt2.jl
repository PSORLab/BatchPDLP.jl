# mathopt2

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = sqr(2*sqr(x1) - x2) + sqr(x2 - 6*sqr(x1))

# EQ (==0) CONSTRAINTS
e2 = x1 - (x1*x2 + 10*x2)
e3 = x1 - 3*x2

# LEQ (<=0) CONSTRAINTS
e4 = x1 + x2 - 1
e5 = -x1 + x2 - 2

# BOUNDS
lvbs = [-Inf, -Inf]
uvbs = [Inf, Inf];

# Create the problem
mathopt2 = Problem("mathopt2", vars, e1, lvbs, uvbs, eq=[e2,e3], leq=[e4,e5])