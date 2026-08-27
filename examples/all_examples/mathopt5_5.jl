# mathopt5_5

# VARIABLES
vars = [x1]

# OBJECTIVE
obj = (sin(1 + 2*x1) + 2*sin(2 + 3*x1) + 3*sin(3 + 4*x1) + 4*sin(4 + 5*x1) + 5*
     sin(5 + 6*x1))

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-10.0]
uvbs = [10.0]

# Create the problem
mathopt5_5 = Problem("mathopt5_5", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)