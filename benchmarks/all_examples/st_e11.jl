# st_e11

# VARIABLES
vars = [x1, x2, x3]

# OBJECTIVE
e3 = 35*x1^0.6 + 35*x2^0.6

# EQ (==0) CONSTRAINTS
e1 = 600*x1 - x1*x3 - 50*x3 + 5000
e2 = 600*x2 + 50*x3 - 15000

# BOUNDS
lvbs = [0.0, 0.0, 0.0]
uvbs = [34, 17, 300];

# Create the problem
st_e11 = Problem("st_e11", vars, e3, lvbs, uvbs, eq=[e1,e2])