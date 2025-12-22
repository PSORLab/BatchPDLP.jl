# ex14_1_4

# VARIABLES
vars = [x1, x2, x3]

# OBJECTIVE
obj = x3

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
    0.5*sin(x1*x2) - 0.5*x1 - 0.0795774703703634*x2 - x3;
    0.920422529629637*exp(2*x1) - 5.4365636*x1 + 0.865255957591193*x2 - x3 - 2.5019678106022;
    0.5*x1 - 0.5*sin(x1*x2) + 0.0795774703703634*x2 - x3;
    5.4365636*x1 - 0.920422529629637*exp(2*x1) - 0.865255957591193*x2 - x3 + 2.5019678106022;
]

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [0.25, 1.5, -Inf]
uvbs = [1.0, 6.28, Inf]

# Create the problem
ex14_1_4 = Problem("ex14_1_4", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)