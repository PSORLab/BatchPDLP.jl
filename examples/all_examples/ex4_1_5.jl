# ex4_1_5

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = 2*sqr(x1) - 1.05*(x1^4) + 0.166666666666667*(x1^6) - x1*x2 + 
     sqr(x2)

# BOUNDS
lvbs = [-5.0, -Inf]
uvbs = [5.0, Inf];

# Create the problem
ex4_1_5 = Problem("ex4_1_5", vars, e1, lvbs, uvbs)