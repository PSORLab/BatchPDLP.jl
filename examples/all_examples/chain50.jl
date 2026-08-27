# chain50

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50, x51, x52, x53, x54, x55, x56, x57, x58, x59, x60, x61, x62, x63, x64, x65, x66, x67, x68, x69, x70, x71, x72, x73, x74, x75, x76, x77, x78, x79, x80, x81, x82, x83, x84, x85, x86, x87, x88, x89, x90, x91, x92, x93, x94, x95, x96, x97, x98, x99, x100, x101, x102]

# OBJECTIVE
obj = 0.01 * (x1 * sqrt(x52^2 + 1) + x2 * sqrt(x53^2 + 1) + x2 * sqrt(x53^2 + 1) + x3 * sqrt(x54^2 + 1) + x3 * sqrt(x54^2 + 1) + x4 * sqrt(x55^2 + 1) + x4 * sqrt(x55^2 + 1) + x5 * sqrt(x56^2 + 1) + x5 * sqrt(x56^2 + 1) + x6 * sqrt(x57^2 + 1) + x6 * sqrt(x57^2 + 1) + x7 * sqrt(x58^2 + 1) + x7 * sqrt(x58^2 + 1) + x8 * sqrt(x59^2 + 1) + x8 * sqrt(x59^2 + 1) + x9 * sqrt(x60^2 + 1) + x9 * sqrt(x60^2 + 1) + x10 * sqrt(x61^2 + 1) + x10 * sqrt(x61^2 + 1) + x11 * sqrt(x62^2 + 1) + x11 * sqrt(x62^2 + 1) + x12 * sqrt(x63^2 + 1) + x12 * sqrt(x63^2 + 1) + x13 * sqrt(x64^2 + 1) + x13 * sqrt(x64^2 + 1) + x14 * sqrt(x65^2 + 1) + x14 * sqrt(x65^2 + 1) + x15 * sqrt(x66^2 + 1) + x15 * sqrt(x66^2 + 1) + x16 * sqrt(x67^2 + 1) + x16 * sqrt(x67^2 + 1) + x17 * sqrt(x68^2 + 1) + x17 * sqrt(x68^2 + 1) + x18 * sqrt(x69^2 + 1) + x18 * sqrt(x69^2 + 1) + x19 * sqrt(x70^2 + 1) + x19 * sqrt(x70^2 + 1) + x20 * sqrt(x71^2 + 1) + x20 * sqrt(x71^2 + 1) + x21 * sqrt(x72^2 + 1) + x21 * sqrt(x72^2 + 1) + x22 * sqrt(x73^2 + 1) + x22 * sqrt(x73^2 + 1) + x23 * sqrt(x74^2 + 1) + x23 * sqrt(x74^2 + 1) + x24 * sqrt(x75^2 + 1) + x24 * sqrt(x75^2 + 1) + x25 * sqrt(x76^2 + 1) + x25 * sqrt(x76^2 + 1) + x26 * sqrt(x77^2 + 1) + x26 * sqrt(x77^2 + 1) + x27 * sqrt(x78^2 + 1) + x27 * sqrt(x78^2 + 1) + x28 * sqrt(x79^2 + 1) + x28 * sqrt(x79^2 + 1) + x29 * sqrt(x80^2 + 1) + x29 * sqrt(x80^2 + 1) + x30 * sqrt(x81^2 + 1) + x30 * sqrt(x81^2 + 1) + x31 * sqrt(x82^2 + 1) + x31 * sqrt(x82^2 + 1) + x32 * sqrt(x83^2 + 1) + x32 * sqrt(x83^2 + 1) + x33 * sqrt(x84^2 + 1) + x33 * sqrt(x84^2 + 1) + x34 * sqrt(x85^2 + 1) + x34 * sqrt(x85^2 + 1) + x35 * sqrt(x86^2 + 1) + x35 * sqrt(x86^2 + 1) + x36 * sqrt(x87^2 + 1) + x36 * sqrt(x87^2 + 1) + x37 * sqrt(x88^2 + 1) + x37 * sqrt(x88^2 + 1) + x38 * sqrt(x89^2 + 1) + x38 * sqrt(x89^2 + 1) + x39 * sqrt(x90^2 + 1) + x39 * sqrt(x90^2 + 1) + x40 * sqrt(x91^2 + 1) + x40 * sqrt(x91^2 + 1) + x41 * sqrt(x92^2 + 1) + x41 * sqrt(x92^2 + 1) + x42 * sqrt(x93^2 + 1) + x42 * sqrt(x93^2 + 1) + x43 * sqrt(x94^2 + 1) + x43 * sqrt(x94^2 + 1) + x44 * sqrt(x95^2 + 1) + x44 * sqrt(x95^2 + 1) + x45 * sqrt(x96^2 + 1) + x45 * sqrt(x96^2 + 1) + x46 * sqrt(x97^2 + 1) + x46 * sqrt(x97^2 + 1) + x47 * sqrt(x98^2 + 1) + x47 * sqrt(x98^2 + 1) + x48 * sqrt(x99^2 + 1) + x48 * sqrt(x99^2 + 1) + x49 * sqrt(x100^2 + 1) + x49 * sqrt(x100^2 + 1) + x50 * sqrt(x101^2 + 1) + x50 * sqrt(x101^2 + 1) + x51 * sqrt(x102^2 + 1))

# EQ (==0) CONSTRAINTS
eq = [
    -x1 + x2 - 0.01 * x52 - 0.01 * x53 - ( 0 ),
    -x2 + x3 - 0.01 * x53 - 0.01 * x54 - ( 0 ),
    -x3 + x4 - 0.01 * x54 - 0.01 * x55 - ( 0 ),
    -x4 + x5 - 0.01 * x55 - 0.01 * x56 - ( 0 ),
    -x5 + x6 - 0.01 * x56 - 0.01 * x57 - ( 0 ),
    -x6 + x7 - 0.01 * x57 - 0.01 * x58 - ( 0 ),
    -x7 + x8 - 0.01 * x58 - 0.01 * x59 - ( 0 ),
    -x8 + x9 - 0.01 * x59 - 0.01 * x60 - ( 0 ),
    -x9 + x10 - 0.01 * x60 - 0.01 * x61 - ( 0 ),
    -x10 + x11 - 0.01 * x61 - 0.01 * x62 - ( 0 ),
    -x11 + x12 - 0.01 * x62 - 0.01 * x63 - ( 0 ),
    -x12 + x13 - 0.01 * x63 - 0.01 * x64 - ( 0 ),
    -x13 + x14 - 0.01 * x64 - 0.01 * x65 - ( 0 ),
    -x14 + x15 - 0.01 * x65 - 0.01 * x66 - ( 0 ),
    -x15 + x16 - 0.01 * x66 - 0.01 * x67 - ( 0 ),
    -x16 + x17 - 0.01 * x67 - 0.01 * x68 - ( 0 ),
    -x17 + x18 - 0.01 * x68 - 0.01 * x69 - ( 0 ),
    -x18 + x19 - 0.01 * x69 - 0.01 * x70 - ( 0 ),
    -x19 + x20 - 0.01 * x70 - 0.01 * x71 - ( 0 ),
    -x20 + x21 - 0.01 * x71 - 0.01 * x72 - ( 0 ),
    -x21 + x22 - 0.01 * x72 - 0.01 * x73 - ( 0 ),
    -x22 + x23 - 0.01 * x73 - 0.01 * x74 - ( 0 ),
    -x23 + x24 - 0.01 * x74 - 0.01 * x75 - ( 0 ),
    -x24 + x25 - 0.01 * x75 - 0.01 * x76 - ( 0 ),
    -x25 + x26 - 0.01 * x76 - 0.01 * x77 - ( 0 ),
    -x26 + x27 - 0.01 * x77 - 0.01 * x78 - ( 0 ),
    -x27 + x28 - 0.01 * x78 - 0.01 * x79 - ( 0 ),
    -x28 + x29 - 0.01 * x79 - 0.01 * x80 - ( 0 ),
    -x29 + x30 - 0.01 * x80 - 0.01 * x81 - ( 0 ),
    -x30 + x31 - 0.01 * x81 - 0.01 * x82 - ( 0 ),
    -x31 + x32 - 0.01 * x82 - 0.01 * x83 - ( 0 ),
    -x32 + x33 - 0.01 * x83 - 0.01 * x84 - ( 0 ),
    -x33 + x34 - 0.01 * x84 - 0.01 * x85 - ( 0 ),
    -x34 + x35 - 0.01 * x85 - 0.01 * x86 - ( 0 ),
    -x35 + x36 - 0.01 * x86 - 0.01 * x87 - ( 0 ),
    -x36 + x37 - 0.01 * x87 - 0.01 * x88 - ( 0 ),
    -x37 + x38 - 0.01 * x88 - 0.01 * x89 - ( 0 ),
    -x38 + x39 - 0.01 * x89 - 0.01 * x90 - ( 0 ),
    -x39 + x40 - 0.01 * x90 - 0.01 * x91 - ( 0 ),
    -x40 + x41 - 0.01 * x91 - 0.01 * x92 - ( 0 ),
    -x41 + x42 - 0.01 * x92 - 0.01 * x93 - ( 0 ),
    -x42 + x43 - 0.01 * x93 - 0.01 * x94 - ( 0 ),
    -x43 + x44 - 0.01 * x94 - 0.01 * x95 - ( 0 ),
    -x44 + x45 - 0.01 * x95 - 0.01 * x96 - ( 0 ),
    -x45 + x46 - 0.01 * x96 - 0.01 * x97 - ( 0 ),
    -x46 + x47 - 0.01 * x97 - 0.01 * x98 - ( 0 ),
    -x47 + x48 - 0.01 * x98 - 0.01 * x99 - ( 0 ),
    -x48 + x49 - 0.01 * x99 - 0.01 * x100 - ( 0 ),
    -x49 + x50 - 0.01 * x100 - 0.01 * x101 - ( 0 ),
    -x50 + x51 - 0.01 * x101 - 0.01 * x102 - ( 0 ),
    0.01 * (sqrt(x52^2 + 1) + sqrt(x53^2 + 1) + sqrt(x53^2 + 1) + sqrt(x54^2 + 1) + sqrt(x54^2 + 1) + sqrt(x55^2 + 1) + sqrt(x55^2 + 1) + sqrt(x56^2 + 1) + sqrt(x56^2 + 1) + sqrt(x57^2 + 1) + sqrt(x57^2 + 1) + sqrt(x58^2 + 1) + sqrt(x58^2 + 1) + sqrt(x59^2 + 1) + sqrt(x59^2 + 1) + sqrt(x60^2 + 1) + sqrt(x60^2 + 1) + sqrt(x61^2 + 1) + sqrt(x61^2 + 1) + sqrt(x62^2 + 1) + sqrt(x62^2 + 1) + sqrt(x63^2 + 1) + sqrt(x63^2 + 1) + sqrt(x64^2 + 1) + sqrt(x64^2 + 1) + sqrt(x65^2 + 1) + sqrt(x65^2 + 1) + sqrt(x66^2 + 1) + sqrt(x66^2 + 1) + sqrt(x67^2 + 1) + sqrt(x67^2 + 1) + sqrt(x68^2 + 1) + sqrt(x68^2 + 1) + sqrt(x69^2 + 1) + sqrt(x69^2 + 1) + sqrt(x70^2 + 1) + sqrt(x70^2 + 1) + sqrt(x71^2 + 1) + sqrt(x71^2 + 1) + sqrt(x72^2 + 1) + sqrt(x72^2 + 1) + sqrt(x73^2 + 1) + sqrt(x73^2 + 1) + sqrt(x74^2 + 1) + sqrt(x74^2 + 1) + sqrt(x75^2 + 1) + sqrt(x75^2 + 1) + sqrt(x76^2 + 1) + sqrt(x76^2 + 1) + sqrt(x77^2 + 1) + sqrt(x77^2 + 1) + sqrt(x78^2 + 1) + sqrt(x78^2 + 1) + sqrt(x79^2 + 1) + sqrt(x79^2 + 1) + sqrt(x80^2 + 1) + sqrt(x80^2 + 1) + sqrt(x81^2 + 1) + sqrt(x81^2 + 1) + sqrt(x82^2 + 1) + sqrt(x82^2 + 1) + sqrt(x83^2 + 1) + sqrt(x83^2 + 1) + sqrt(x84^2 + 1) + sqrt(x84^2 + 1) + sqrt(x85^2 + 1) + sqrt(x85^2 + 1) + sqrt(x86^2 + 1) + sqrt(x86^2 + 1) + sqrt(x87^2 + 1) + sqrt(x87^2 + 1) + sqrt(x88^2 + 1) + sqrt(x88^2 + 1) + sqrt(x89^2 + 1) + sqrt(x89^2 + 1) + sqrt(x90^2 + 1) + sqrt(x90^2 + 1) + sqrt(x91^2 + 1) + sqrt(x91^2 + 1) + sqrt(x92^2 + 1) + sqrt(x92^2 + 1) + sqrt(x93^2 + 1) + sqrt(x93^2 + 1) + sqrt(x94^2 + 1) + sqrt(x94^2 + 1) + sqrt(x95^2 + 1) + sqrt(x95^2 + 1) + sqrt(x96^2 + 1) + sqrt(x96^2 + 1) + sqrt(x97^2 + 1) + sqrt(x97^2 + 1) + sqrt(x98^2 + 1) + sqrt(x98^2 + 1) + sqrt(x99^2 + 1) + sqrt(x99^2 + 1) + sqrt(x100^2 + 1) + sqrt(x100^2 + 1) + sqrt(x101^2 + 1) + sqrt(x101^2 + 1) + sqrt(x102^2 + 1)) - ( 4 ),
]

# GEQ (>=0) CONSTRAINTS
geq = [
]

# LEQ (<=0) CONSTRAINTS
leq = [
]

# BOUNDS
lvbs = [1, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 3, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf]
uvbs = [1, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, 3, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]

# Create the problem
chain50 = Problem("chain50", vars, obj, lvbs, uvbs, eq=eq, geq=geq, leq=leq)