# ex14_1_8

# VARIABLES
vars = [x1, x2, x3]

# OBJECTIVE
e1 = x3

# LEQ (<=0) CONSTRAINTS
e2 = exp(10*x1/(1 + 0.01*x1))*(0.0476666666666666 - 0.0649999999999999*x1) - x1 - x3
e3 = x1 - exp(10*x1/(1 + 0.01*x1))*(0.0476666666666666 - 0.0649999999999999*x1) - x3
e4 = exp(10*x2/(1 + 0.01*x2))*(0.143 - 0.13*x1 - 0.195*x2) + x1 - 3*x2 - x3
e5 = (-exp(10*x2/(1 + 0.01*x2))*(0.143 - 0.13*x1 - 0.195*x2)) - x1 + 3*x2 - x3

# BOUNDS
lvbs = [0.0, 0.0, -Inf]
uvbs = [1.0, 1.0, Inf];

# Create the problem
ex14_1_8 = Problem("ex14_1_8", vars, e1, lvbs, uvbs, leq=[e2,e3,e4,e5])