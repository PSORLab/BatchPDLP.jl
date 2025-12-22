# ex8_1_5

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = 4*sqr(x1) - 2.1*(x1^4) + 0.333333333333333*(x1^6) + x1*x2 - 4*
     sqr(x2) + 4*(x2^4)

# BOUNDS
lvbs = [-Inf, -Inf]
uvbs = [Inf, Inf];

# Create the problem
ex8_1_5 = Problem("ex8_1_5", vars, e1, lvbs, uvbs)