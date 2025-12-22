# ex8_1_3

# VARIABLES
vars = [x1, x2]

# OBJECTIVE
e1 = ((1 + x1 + x2)^2 * (3 * x1^2 -
    14 * x1 + 6 * x1 * x2 - 14 * x2 + 3 * x2^2 + 19) + 1) * ((2 *
    x1 - 3 * x2)^2 * (12 * x1^2 - 32 * x1 - 36 * x1 * x2 + 48 *
    x2 + 27 * x2^2 + 18) + 30)


# BOUNDS
lvbs = [-Inf, -Inf]
uvbs = [Inf, Inf];

# Create the problem
ex8_1_3 = Problem("ex8_1_3", vars, e1, lvbs, uvbs)