# ex6_1_2

# VARIABLES
vars = [x1, x2, x3, x4]

# OBJECTIVE
obj = x1 * (log(x1) + 0.06391) + x2 * (log(x2) - 0.02875) + 0.925356626778358 * x1 * x4 + 0.746014540096753 * x2 * x3

# EQ (==0) CONSTRAINTS
eq = [
    x3 * (x1 + 0.159040857374844 * x2) - x1 - ( 0 ),
    x4 * (0.307941026821595 * x1 + x2) - x2 - ( 0 ),
    x1 + x2 - ( 1 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [1e-06, 1e-06, 0, 0]
uvbs = [1, 1, Inf, Inf]

# Create the problem
ex6_1_2 = Problem("ex6_1_2", vars, obj, lvbs, uvbs, eq=eq, leq=leq)