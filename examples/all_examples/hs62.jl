# hs62

# VARIABLES
vars = [x2, x3, x4]

# OBJECTIVE
e1 = -32.174*(255*log((0.03 + x2 + x3 + x4)/(0.03 + 0.09*x2 + x3 + x4)) + 280*
     log((0.03 + x3 + x4)/(0.03 + 0.07*x3 + x4)) + 290*log((0.03 + x4)/(0.03 + 
     0.13*x4)))

# EQ (==0) CONSTRAINTS
e2 = 20*sqr((-1) + x2 + x3 + x4)

# BOUNDS
lvbs = [0.0, 0.0, 0.0]
uvbs = [Inf, Inf, Inf];

# Create the problem
hs62 = Problem("hs62", vars, e1, lvbs, uvbs, eq=[e2])