# mathopt5_7

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = 0.01*(-0.0218343*(x1^2) - 8.9248e-5*x1 + 0.998266*(x1^3) - 1.6995*
     (x1^4) + 0.2*(x1^5))

# BOUNDS
lvbs = [0.0]
uvbs = [8.0];

# Create the problem
mathopt5_7 = Problem("mathopt5_7", vars, e1, lvbs, uvbs)