# ex4_1_3

# VARIABLES
vars = [x1]

# OBJECTIVE
e1 = 8.9248e-5*x1 - 0.0218343*(x1^2) + 0.998266*(x1^3) - 1.6995*(
     x1^4) + 0.2*(x1^5)

# BOUNDS
lvbs = [0.0]
uvbs = [10.0];

# Create the problem
ex4_1_3 = Problem("ex4_1_3", vars, e1, lvbs, uvbs)