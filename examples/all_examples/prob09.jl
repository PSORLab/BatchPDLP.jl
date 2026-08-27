# prob09

# VARIABLES
vars = [x2, x3]

# OBJECTIVE
e1 = 100*sqr(x3 - sqr(x2)) + sqr(1 - x2)

# BOUNDS
lvbs = [-2, -2]
uvbs = [2, 2];

# Create the problem
prob09 = Problem("prob09", vars, e1, lvbs, uvbs)