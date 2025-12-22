# ex6_2_14

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = x1 * (log(x1 / (x1 + x3)) + log(x1 / (x1 + 0.095173 * x3))) + x3 * (log(x3 / (x1 + x3)) + log(x3 / (0.30384 * x1 + x3))) + log(x1 + 2.6738 * x3) * (x1 + 2.6738 * x3) + log(0.374 * x1 + x3) * (0.374 * x1 + x3) + 2.6738 * log(x3 / (x1 + 2.6738 * x3)) * x3 + 0.374 * log(x1 / (0.374 * x1 + x3)) * x1 + x2 * (log(x2 / (x2 + x4)) + log(x2 / (x2 + 0.095173 * x4))) + x4 * (log(x4 / (x2 + x4)) + log(x4 / (0.30384 * x2 + x4))) + log(x2 + 2.6738 * x4) * (x2 + 2.6738 * x4) + log(0.374 * x2 + x4) * (0.374 * x2 + x4) + 2.6738 * log(x4 / (x2 + 2.6738 * x4)) * x4 + 0.374 * log(x2 / (0.374 * x2 + x4)) * x2 - 3.6838 * log(x1) * x1 - 1.59549 * log(x3) * x3 - 3.6838 * log(x2) * x2 - 1.59549 * log(x4) * x4

# EQ (==0) CONSTRAINTS
eq = [
    x1 + x2 - ( 0.5 ),
    x3 + x4 - ( 0.5 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [1e-07, 1e-07, 1e-07, 1e-07]
uvbs = [0.5, 0.5, 0.5, 0.5]

# Create the problem
ex6_2_14 = Problem("ex6_2_14", vars, obj, lvbs, uvbs, eq=eq, leq=leq)