# mathopt5_2

# VARIABLES
vars = [x1]

# OBJECTIVE
obj = (sqr(x1) - cos(18*x1)) 

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-5.0]
uvbs = [5.0]

# Create the problem
mathopt5_2 = Problem("mathopt5_2", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)