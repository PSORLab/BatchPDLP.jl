# mathopt5_8

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = 2*(x1^2) - x1 - 1.05*(x1^4) + 0.1666667*(x1^6)

# BOUNDS
lvbs = [-2.0]
uvbs = [2.5];

# Create the problem
mathopt5_8 = Problem("mathopt5_8", vars, e1, lvbs, uvbs)