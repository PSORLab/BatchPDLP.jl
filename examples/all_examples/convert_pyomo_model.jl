# convert_pyomo_model

# VARIABLES
vars = []

# EQ (==0) CONSTRAINTS
eq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = []
uvbs = []

# Create the problem
convert_pyomo_model = Problem("convert_pyomo_model", vars, obj, lvbs, uvbs, eq=eq, leq=leq)