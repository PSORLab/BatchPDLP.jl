# mathopt1

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = 10*sqr(sqr(x1) - x2) + sqr((-1) + x1)

# EQ (==0) CONSTRAINTS
e2 = x1 - x1*x2

# LEQ (<=0) CONSTRAINTS
e3 = 3*x1 + 4*x2 - 25

# BOUNDS
lvbs = [-10.0, -15.0]
uvbs = [20.0, 20.0];

# Create the problem
mathopt1 = Problem("mathopt1", vars, e1, lvbs, uvbs, eq=[e2], leq=[e3])