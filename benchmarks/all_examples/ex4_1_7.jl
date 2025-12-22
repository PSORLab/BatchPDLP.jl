# ex4_1_7

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = (x1^4) - 3*(x1^3) - 1.5*(x1^2) + 10*x1

# BOUNDS
lvbs = [-5.0]
uvbs = [5.0];

# Create the problem
ex4_1_7 = Problem("ex4_1_7", vars, e1, lvbs, uvbs)