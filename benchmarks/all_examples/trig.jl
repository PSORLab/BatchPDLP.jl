# trig

# VARIABLES
vars = [x1]

# OBJECTIVE
obj = (sin(11*x1) + cos(13*x1) - sin(17*x1) - cos(19*x1))

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = [5*sin(x1) - x1]

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-2.0]
uvbs = [5.0]

# Create the problem
trig = Problem("trig", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)