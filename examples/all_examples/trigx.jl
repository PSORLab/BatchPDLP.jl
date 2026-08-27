# trigx

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
obj = (x1*x1 + x2*x2)

# EQ (==0) CONSTRAINTS
eq = [
    x1 - sin(2*x1 + 3*x2) - cos(3*x1 - 5*x2);
    x2 - sin(x1 - 2*x2) + cos(x1 + 3*x2)
]

# LEQ (<=0) CONSTRAINTS
leq = []

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [-Inf, -Inf]
uvbs = [Inf, Inf]

# Create the problem
trigx = Problem("trigx", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)