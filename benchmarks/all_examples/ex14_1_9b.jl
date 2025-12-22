# ex14_1_9

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = x2

# LEQ (<=0) CONSTRAINTS
e2 = 4510067.11409396*exp(-7548.11926028431/x1)*x1 + 0.00335570469798658*x1 - 
     2020510067.11409*exp(-7548.11926028431/x1) - x2 - 1
e3 = (-4510067.11409396*exp(-7548.11926028431/x1)*x1) - 0.00335570469798658*x1
      + 2020510067.11409*exp(-7548.11926028431/x1) - x2 + 1

# BOUNDS
lvbs = [100, -Inf]
uvbs = [1000, Inf];

# Create the problem
ex14_1_9 = Problem("ex14_1_9", vars, e1, lvbs, uvbs, leq=[e2,e3])