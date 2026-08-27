# st_e06

# VARIABLES
vars = [x1, x2, x3]

# OBJECTIVE
e4 = Num(0)

# EQ (==0) CONSTRAINTS
e1 = x3*x3 - 0.000169*(x2^3)*x1
e2 = x1 + x2 + x3 - 50
e3 = -3*x1 + x2

# BOUNDS
lvbs = [0.0, 0.0, 0.0]
uvbs = [12.5, 37.5, 50.0];

# Create the problem
st_e06 = Problem("st_e06", vars, e4, lvbs, uvbs, eq=[e1,e2,e3])