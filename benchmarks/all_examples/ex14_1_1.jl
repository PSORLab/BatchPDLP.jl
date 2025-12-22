# ex14_1_1

# VARIABLES
vars = [x1, x2, x3]

# OBJECTIVE
e1 = x3

# LEQ (<=0) CONSTRAINTS
e2 = 2*sqr(x2) + 4*x1*x2 - 42*x1 + 4*(x1^3) - x3 - 14
e3 = (-2*sqr(x2)) - 4*x1*x2 + 42*x1 - 4*(x1^3) - x3 + 14
e4 = 2*sqr(x1) + 4*x1*x2 - 26*x2 + 4*(x2^3) - x3 - 22
e5 = (-2*sqr(x1)) - 4*x1*x2 + 26*x2 - 4*(x2^3) - x3 + 22

# BOUNDS
lvbs = [-5.0, -5.0, -Inf]
uvbs = [5.0, 5.0, Inf];

# Create the problem
ex14_1_1 = Problem("ex14_1_1", vars, e1, lvbs, uvbs, leq=[e2,e3,e4,e5])