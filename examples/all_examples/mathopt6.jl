# mathopt6

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
obj = (exp(sin(50*x1)) + sin(60*exp(x2)) + sin(70*sin(x1)) + sin(sin(80*x2)) - 
     sin(10*x1 + 10*x2) + 0.25*(sqr(x1) + sqr(x2)))

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# GEQ (>=0) CONSTRAINTS
geq = []


# BOUNDS
lvbs = [-3.0, -3.0]
uvbs = [3.0, 3.0]

# Create the problem
mathopt6 = Problem("mathopt6", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)