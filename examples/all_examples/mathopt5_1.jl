# mathopt5_1

# VARIABLES
vars = [x1]

# OBJECTIVE
obj = (exp(-x1) - ((sin(x1))^3)) 

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-5.0]
uvbs = [10.0]

# Create the problem
mathopt5_1 = Problem("mathopt5_1", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)