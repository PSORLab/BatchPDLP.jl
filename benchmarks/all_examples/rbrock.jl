# rbrock

# VARIABLES
vars = [x2, x3]

# OBJECTIVE
e1 = 100*sqr(x3 - sqr(x2)) + sqr(1 - x2)

# BOUNDS
lvbs = [-10.0, -10.0]
uvbs = [5.0, 10.0];

# Create the problem
rbrock = Problem("rbrock", vars, e1, lvbs, uvbs)