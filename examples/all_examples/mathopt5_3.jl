# mathopt5_3

# VARIABLES
vars = [x1]

# OBJECTIVE
obj = sin(x1)*sqr(cos(x1) - sin(x1))

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [0.0]
uvbs = [10.0]

# Create the problem
mathopt5_3 = Problem("mathopt5_3", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)