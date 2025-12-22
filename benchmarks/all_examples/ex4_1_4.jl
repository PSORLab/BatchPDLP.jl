# ex4_1_4

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = 4*(x1^2) - 4*(x1^3) + (x1^4)

# BOUNDS
lvbs = [-5.0]
uvbs = [5.0];

# Create the problem
ex4_1_4 = Problem("ex4_1_4", vars, e1, lvbs, uvbs)