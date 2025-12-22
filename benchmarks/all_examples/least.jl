# least

# VARIABLES
vars = [x2, x3, x4]

# OBJECTIVE
e1 = sqr(127 - exp(-5*x4)*x3 - x2) + sqr(151 - exp(-3*x4)*x3 - x2) + sqr(379
      - exp(-x4)*x3 - x2) + sqr(421 - exp(5*x4)*x3 - x2) + sqr(460 - exp(3*x4)*
     x3 - x2) + sqr(426 - exp(x4)*x3 - x2)

# BOUNDS
lvbs = [-Inf, -Inf, -5]
uvbs = [Inf, Inf, 5];

# Create the problem
least = Problem("least", vars, e1, lvbs, uvbs)