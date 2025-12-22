# ex8_1_4

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = 12*sqr(x1) - 6.3*(x1^4) + (x1^6) - 6*x1*x2 + 6*sqr(x2)

# BOUNDS
lvbs = [-Inf, -Inf]
uvbs = [Inf, Inf];

# Create the problem
ex8_1_4 = Problem("ex8_1_4", vars, e1, lvbs, uvbs)