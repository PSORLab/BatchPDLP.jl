# st_e19

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e3 = (x1^4) - 14*sqr(x1) + 24*x1 - sqr(x2)

# LEQ (<=0) CONSTRAINTS
e1 = -x1 + x2 - 8
e2 = (-sqr(x1)) - 2*x1 + x2 + 2

# BOUNDS
lvbs = [-8.0, 0.0]
uvbs = [10.0, 10.0];

# Create the problem
st_e19 = Problem("st_e19", vars, e3, lvbs, uvbs, leq=[e1,e2])