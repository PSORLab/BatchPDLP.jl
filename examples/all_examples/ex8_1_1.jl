# ex8_1_1

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
obj = (cos(x1)*sin(x2) - x1/(1 + sqr(x2)))

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-1.0, -1.0]
uvbs = [2.0, 1.0]

# Create the problem
ex8_1_1 = Problem("ex8_1_1", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)