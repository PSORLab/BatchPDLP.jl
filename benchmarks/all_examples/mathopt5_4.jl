# mathopt5_4

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = (3 + 18*(x1^2) - 10*x1 - 13*(x1^3) + 2*(x1^4))^2

# BOUNDS
lvbs = [-1.0]
uvbs = [4.0];

# Create the problem
mathopt5_4 = Problem("mathopt5_4", vars, e1, lvbs, uvbs)