# ann_fermentation_exp

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12]

# OBJECTIVE
obj = -91.91176470588235 * x9 / x1

# EQ (==0) CONSTRAINTS
eq = [
    2 / (exp(2 * x11) + 1) + x7 - ( 1 ),
    2 / (exp(2 * x12) + 1) + x8 - ( 1 ),
    -0.0949474332688833 * x7 + 0.968637250639063 * x8 + x10 - ( 0.002499597315649 ),
    x9 - 86.324 * x10 - ( 92.74 ),
    -0.025 * x1 + x4 - ( -3.5 ),
    -x2 + x5 - ( -2 ),
    -0.04 * x3 + x6 - ( -1.4 ),
    -24.4380718077469 * x4 - 22.0304402344789 * x5 + 147.921509281049 * x6 + x11 - ( 101.235018055261 ),
    1.48567642304727 * x4 - 0.0532843142008436 * x5 + 0.910590580134437 * x6 + x12 - ( 0.18771256886977 ),
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [100, 1, 10, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [180, 3, 60, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
ann_fermentation_exp = Problem("ann_fermentation_exp", vars, obj, lvbs, uvbs, eq=eq, leq=leq)