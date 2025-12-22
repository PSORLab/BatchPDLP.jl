# ex3_1_1

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8]

# OBJECTIVE
obj = x1 + x2 + x3

# LEQ (<=0) CONSTRAINTS
leq = [0.0025*x4 + 0.0025*x6 - 1,
    -0.0025*x4 + 0.0025*x5 + 0.0025*x7 - 1,
    -0.01*x5 + 0.01*x8 - 1,
    -x1*x6 + 100*x1 + 833.33252*x4 - 83333.333,
    x2*x4 - x2*x7 - 1250*x4 + 1250*x5,
    x3*x5 - x3*x8 - 2500*x5 + 1250000]

# BOUNDS
lvbs = [100.0, 1000.0, 1000.0, 10.0, 10.0, 10.0, 10.0, 10.0]
uvbs = [10000.0, 10000.0, 10000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0];

# Create the problem
ex3_1_1 = Problem("ex3_1_1", vars, obj, lvbs, uvbs, leq=leq)


