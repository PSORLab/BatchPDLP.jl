# ex4_1_8

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = (sqr(x2) - 7*x2) - 12*x1

# EQ (==0) CONSTRAINTS
e2 = -2*(x1^4) - x2 + 2

# BOUNDS
lvbs = [0.0, 0.0]
uvbs = [2.0, 3.0];

# Create the problem
ex4_1_8 = Problem("ex4_1_8", vars, e1, lvbs, uvbs, eq=[e2])