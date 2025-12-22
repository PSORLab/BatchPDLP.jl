# ex8_1_6

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 =  -1/(0.1 + ((-4) + x1)^2 + ((-4) + x2)^2) - 1/(0.2 + ((-1) + x1)^2 + 
     ((-1) + x2)^2) - 1/(0.2 + ((-8) + x1)^2 + ((-8) + x2)^2)
    
# BOUNDS
lvbs = [-Inf, -Inf]
uvbs = [Inf, Inf]

# Create the problem
ex8_1_6 = Problem("ex8_1_6", vars, e1, lvbs, uvbs)