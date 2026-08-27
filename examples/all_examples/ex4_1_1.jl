# ex4_1_1

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = x1^6 - 2.08*x1^5 + 0.4875*x1^4 + 7.1*x1^3 - 3.95*x1^2 - x1 - 0.1

# BOUNDS
lvbs = [-2.0]
uvbs = [11.0];

# Create the problem
ex4_1_1 = Problem("ex4_1_1", vars, e1, lvbs, uvbs)