# ex4_1_6

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = (x1^6) - 15*(x1^4) + 27*(x1^2) - 250.0

# BOUNDS
lvbs = [-5.0]
uvbs = [5.0];

# Create the problem
ex4_1_6 = Problem("ex4_1_6", vars, e1, lvbs, uvbs)