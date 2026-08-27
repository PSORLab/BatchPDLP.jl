# chem

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11]

# OBJECTIVE
obj = x1 * (log(x1 / x11) - 6.05576803624071) + x2 * (log(x2 / x11) - 17.1307680362407) + x3 * (log(x3 / x11) - 34.0207680362407) + x4 * (log(x4 / x11) - 5.88076803624071) + x5 * (log(x5 / x11) - 24.6877680362407) + x6 * (log(x6 / x11) - 14.9527680362407) + x7 * (log(x7 / x11) - 24.0667680362407) + x8 * (log(x8 / x11) - 10.6747680362407) + x9 * (log(x9 / x11) - 26.6287680362407) + x10 * (log(x10 / x11) - 22.1447680362407)

# EQ (==0) CONSTRAINTS
eq = [
    x1 + 2 * x2 + 2 * x3 + x6 + x10 - ( 2 ),
    x4 + 2 * x5 + x6 + x7 - ( 1 ),
    x3 + x7 + x8 + 2 * x9 + x10 - ( 1 ),
    -x1 - x2 - x3 - x4 - x5 - x6 - x7 - x8 - x9 - x10 + x11 - ( 0 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.01]
uvbs = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
chem = Problem("chem", vars, obj, lvbs, uvbs, eq=eq, leq=leq)