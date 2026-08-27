# ex4_1_9

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = -(x1 + x2)

# LEQ (<=0) CONSTRAINTS
e2 = 8*(x1^3) - 2*(x1^4) - 8*sqr(x1) + x2 - 2
e3 = 32*(x1^3) - 4*(x1^4) - 88*sqr(x1) + 96*x1 + x2 - 36

# BOUNDS
lvbs = [0.0, 0.0]
uvbs = [3.0, 4.0];

# Create the problem
ex4_1_9 = Problem("ex4_1_9", vars, e1, lvbs, uvbs, leq=[e2,e3])