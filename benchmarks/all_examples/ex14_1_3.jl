# ex14_1_3

# VARIABLES
vars = [x1, x2, x3]

# OBJECTIVE
e1 = x3

# LEQ (<=0) CONSTRAINTS
e2 = 10000*x1*x2 - x3 - 1
e3 = -10000*x1*x2 - x3 + 1
e4 = exp(-x1) + exp(-x2) - x3 - 1.001
e5 = (-exp(-x1)) - exp(-x2) - x3 + 1.001

# BOUNDS
lvbs = [5.49E-6, 0.0021961, -Inf]
uvbs = [4.553, 18.21, Inf];

# Create the problem
ex14_1_3 = Problem("ex14_1_3", vars, e1, lvbs, uvbs, leq=[e2,e3,e4,e5])