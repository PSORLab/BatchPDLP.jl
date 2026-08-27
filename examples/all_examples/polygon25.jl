# polygon25

# VARIABLES
vars = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50]

# OBJECTIVE
obj = -0.5*(x2*x1*sin(x27 - x26) + x3*x2*sin(x28 - x27) + x4*x3*sin(x29 - x28)
      + x5*x4*sin(x30 - x29) + x6*x5*sin(x31 - x30) + x7*x6*sin(x32 - x31) + x8
     *x7*sin(x33 - x32) + x9*x8*sin(x34 - x33) + x10*x9*sin(x35 - x34) + x11*
     x10*sin(x36 - x35) + x12*x11*sin(x37 - x36) + x13*x12*sin(x38 - x37) + x14
     *x13*sin(x39 - x38) + x15*x14*sin(x40 - x39) + x16*x15*sin(x41 - x40) + 
     x17*x16*sin(x42 - x41) + x18*x17*sin(x43 - x42) + x19*x18*sin(x44 - x43)
      + x20*x19*sin(x45 - x44) + x21*x20*sin(x46 - x45) + x22*x21*sin(x47 - x46
     ) + x23*x22*sin(x48 - x47) + x24*x23*sin(x49 - x48) + x25*x24*sin(x50 - 
     x49))

# EQ (==0) CONSTRAINTS
eq = []

# LEQ (<=0) CONSTRAINTS
leq = [
    sqr(x1) + sqr(x2) - 2*x1*x2*cos(x27 - x26) - 1;
    sqr(x1) + sqr(x3) - 2*x1*x3*cos(x28 - x26) - 1;
    sqr(x1) + sqr(x4) - 2*x1*x4*cos(x29 - x26) - 1;
    sqr(x1) + sqr(x5) - 2*x1*x5*cos(x30 - x26) - 1;
    sqr(x1) + sqr(x6) - 2*x1*x6*cos(x31 - x26) - 1;
    sqr(x1) + sqr(x7) - 2*x1*x7*cos(x32 - x26) - 1;
    sqr(x1) + sqr(x8) - 2*x1*x8*cos(x33 - x26) - 1;
    sqr(x1) + sqr(x9) - 2*x1*x9*cos(x34 - x26) - 1;
    sqr(x1) + sqr(x10) - 2*x1*x10*cos(x35 - x26) - 1;
    sqr(x1) + sqr(x11) - 2*x1*x11*cos(x36 - x26) - 1;
    sqr(x1) + sqr(x12) - 2*x1*x12*cos(x37 - x26) - 1;
    sqr(x1) + sqr(x13) - 2*x1*x13*cos(x38 - x26) - 1;
    sqr(x1) + sqr(x14) - 2*x1*x14*cos(x39 - x26) - 1;
    sqr(x1) + sqr(x15) - 2*x1*x15*cos(x40 - x26) - 1;
    sqr(x1) + sqr(x16) - 2*x1*x16*cos(x41 - x26) - 1;
    sqr(x1) + sqr(x17) - 2*x1*x17*cos(x42 - x26) - 1;
    sqr(x1) + sqr(x18) - 2*x1*x18*cos(x43 - x26) - 1;
    sqr(x1) + sqr(x19) - 2*x1*x19*cos(x44 - x26) - 1;
    sqr(x1) + sqr(x20) - 2*x1*x20*cos(x45 - x26) - 1;
    sqr(x1) + sqr(x21) - 2*x1*x21*cos(x46 - x26) - 1;
    sqr(x1) + sqr(x22) - 2*x1*x22*cos(x47 - x26) - 1;
    sqr(x1) + sqr(x23) - 2*x1*x23*cos(x48 - x26) - 1;
    sqr(x1) + sqr(x24) - 2*x1*x24*cos(x49 - x26) - 1;
    sqr(x1) + sqr(x25) - 2*x1*x25*cos(x50 - x26) - 1;
    sqr(x2) + sqr(x3) - 2*x2*x3*cos(x28 - x27) - 1;
    sqr(x2) + sqr(x4) - 2*x2*x4*cos(x29 - x27) - 1;
    sqr(x2) + sqr(x5) - 2*x2*x5*cos(x30 - x27) - 1;
    sqr(x2) + sqr(x6) - 2*x2*x6*cos(x31 - x27) - 1;
    sqr(x2) + sqr(x7) - 2*x2*x7*cos(x32 - x27) - 1;
    sqr(x2) + sqr(x8) - 2*x2*x8*cos(x33 - x27) - 1;
    sqr(x2) + sqr(x9) - 2*x2*x9*cos(x34 - x27) - 1;
    sqr(x2) + sqr(x10) - 2*x2*x10*cos(x35 - x27) - 1;
    sqr(x2) + sqr(x11) - 2*x2*x11*cos(x36 - x27) - 1;
    sqr(x2) + sqr(x12) - 2*x2*x12*cos(x37 - x27) - 1;
    sqr(x2) + sqr(x13) - 2*x2*x13*cos(x38 - x27) - 1;
    sqr(x2) + sqr(x14) - 2*x2*x14*cos(x39 - x27) - 1;
    sqr(x2) + sqr(x15) - 2*x2*x15*cos(x40 - x27) - 1;
    sqr(x2) + sqr(x16) - 2*x2*x16*cos(x41 - x27) - 1;
    sqr(x2) + sqr(x17) - 2*x2*x17*cos(x42 - x27) - 1;
    sqr(x2) + sqr(x18) - 2*x2*x18*cos(x43 - x27) - 1;
    sqr(x2) + sqr(x19) - 2*x2*x19*cos(x44 - x27) - 1;
    sqr(x2) + sqr(x20) - 2*x2*x20*cos(x45 - x27) - 1;
    sqr(x2) + sqr(x21) - 2*x2*x21*cos(x46 - x27) - 1;
    sqr(x2) + sqr(x22) - 2*x2*x22*cos(x47 - x27) - 1;
    sqr(x2) + sqr(x23) - 2*x2*x23*cos(x48 - x27) - 1;
    sqr(x2) + sqr(x24) - 2*x2*x24*cos(x49 - x27) - 1;
    sqr(x2) + sqr(x25) - 2*x2*x25*cos(x50 - x27) - 1;
    sqr(x3) + sqr(x4) - 2*x3*x4*cos(x29 - x28) - 1;
    sqr(x3) + sqr(x5) - 2*x3*x5*cos(x30 - x28) - 1;
    sqr(x3) + sqr(x6) - 2*x3*x6*cos(x31 - x28) - 1;
    sqr(x3) + sqr(x7) - 2*x3*x7*cos(x32 - x28) - 1;
    sqr(x3) + sqr(x8) - 2*x3*x8*cos(x33 - x28) - 1;
    sqr(x3) + sqr(x9) - 2*x3*x9*cos(x34 - x28) - 1;
    sqr(x3) + sqr(x10) - 2*x3*x10*cos(x35 - x28) - 1;
    sqr(x3) + sqr(x11) - 2*x3*x11*cos(x36 - x28) - 1;
    sqr(x3) + sqr(x12) - 2*x3*x12*cos(x37 - x28) - 1;
    sqr(x3) + sqr(x13) - 2*x3*x13*cos(x38 - x28) - 1;
    sqr(x3) + sqr(x14) - 2*x3*x14*cos(x39 - x28) - 1;
    sqr(x3) + sqr(x15) - 2*x3*x15*cos(x40 - x28) - 1;
    sqr(x3) + sqr(x16) - 2*x3*x16*cos(x41 - x28) - 1;
    sqr(x3) + sqr(x17) - 2*x3*x17*cos(x42 - x28) - 1;
    sqr(x3) + sqr(x18) - 2*x3*x18*cos(x43 - x28) - 1;
    sqr(x3) + sqr(x19) - 2*x3*x19*cos(x44 - x28) - 1;
    sqr(x3) + sqr(x20) - 2*x3*x20*cos(x45 - x28) - 1;
    sqr(x3) + sqr(x21) - 2*x3*x21*cos(x46 - x28) - 1;
    sqr(x3) + sqr(x22) - 2*x3*x22*cos(x47 - x28) - 1;
    sqr(x3) + sqr(x23) - 2*x3*x23*cos(x48 - x28) - 1;
    sqr(x3) + sqr(x24) - 2*x3*x24*cos(x49 - x28) - 1;
    sqr(x3) + sqr(x25) - 2*x3*x25*cos(x50 - x28) - 1;
    sqr(x4) + sqr(x5) - 2*x4*x5*cos(x30 - x29) - 1;
    sqr(x4) + sqr(x6) - 2*x4*x6*cos(x31 - x29) - 1;
    sqr(x4) + sqr(x7) - 2*x4*x7*cos(x32 - x29) - 1;
    sqr(x4) + sqr(x8) - 2*x4*x8*cos(x33 - x29) - 1;
    sqr(x4) + sqr(x9) - 2*x4*x9*cos(x34 - x29) - 1;
    sqr(x4) + sqr(x10) - 2*x4*x10*cos(x35 - x29) - 1;
    sqr(x4) + sqr(x11) - 2*x4*x11*cos(x36 - x29) - 1;
    sqr(x4) + sqr(x12) - 2*x4*x12*cos(x37 - x29) - 1;
    sqr(x4) + sqr(x13) - 2*x4*x13*cos(x38 - x29) - 1;
    sqr(x4) + sqr(x14) - 2*x4*x14*cos(x39 - x29) - 1;
    sqr(x4) + sqr(x15) - 2*x4*x15*cos(x40 - x29) - 1;
    sqr(x4) + sqr(x16) - 2*x4*x16*cos(x41 - x29) - 1;
    sqr(x4) + sqr(x17) - 2*x4*x17*cos(x42 - x29) - 1;
    sqr(x4) + sqr(x18) - 2*x4*x18*cos(x43 - x29) - 1;
    sqr(x4) + sqr(x19) - 2*x4*x19*cos(x44 - x29) - 1;
    sqr(x4) + sqr(x20) - 2*x4*x20*cos(x45 - x29) - 1;
    sqr(x4) + sqr(x21) - 2*x4*x21*cos(x46 - x29) - 1;
    sqr(x4) + sqr(x22) - 2*x4*x22*cos(x47 - x29) - 1;
    sqr(x4) + sqr(x23) - 2*x4*x23*cos(x48 - x29) - 1;
    sqr(x4) + sqr(x24) - 2*x4*x24*cos(x49 - x29) - 1;
    sqr(x4) + sqr(x25) - 2*x4*x25*cos(x50 - x29) - 1;
    sqr(x5) + sqr(x6) - 2*x5*x6*cos(x31 - x30) - 1;
    sqr(x5) + sqr(x7) - 2*x5*x7*cos(x32 - x30) - 1;
    sqr(x5) + sqr(x8) - 2*x5*x8*cos(x33 - x30) - 1;
    sqr(x5) + sqr(x9) - 2*x5*x9*cos(x34 - x30) - 1;
    sqr(x5) + sqr(x10) - 2*x5*x10*cos(x35 - x30) - 1;
    sqr(x5) + sqr(x11) - 2*x5*x11*cos(x36 - x30) - 1;
    sqr(x5) + sqr(x12) - 2*x5*x12*cos(x37 - x30) - 1;
    sqr(x5) + sqr(x13) - 2*x5*x13*cos(x38 - x30) - 1;
    sqr(x5) + sqr(x14) - 2*x5*x14*cos(x39 - x30) - 1;
    sqr(x5) + sqr(x15) - 2*x5*x15*cos(x40 - x30) - 1;
    sqr(x5) + sqr(x16) - 2*x5*x16*cos(x41 - x30) - 1;
    sqr(x5) + sqr(x17) - 2*x5*x17*cos(x42 - x30) - 1;
    sqr(x5) + sqr(x18) - 2*x5*x18*cos(x43 - x30) - 1;
    sqr(x5) + sqr(x19) - 2*x5*x19*cos(x44 - x30) - 1;
    sqr(x5) + sqr(x20) - 2*x5*x20*cos(x45 - x30) - 1;
    sqr(x5) + sqr(x21) - 2*x5*x21*cos(x46 - x30) - 1;
    sqr(x5) + sqr(x22) - 2*x5*x22*cos(x47 - x30) - 1;
    sqr(x5) + sqr(x23) - 2*x5*x23*cos(x48 - x30) - 1;
    sqr(x5) + sqr(x24) - 2*x5*x24*cos(x49 - x30) - 1;
    sqr(x5) + sqr(x25) - 2*x5*x25*cos(x50 - x30) - 1;
    sqr(x6) + sqr(x7) - 2*x6*x7*cos(x32 - x31) - 1;
    sqr(x6) + sqr(x8) - 2*x6*x8*cos(x33 - x31) - 1;
    sqr(x6) + sqr(x9) - 2*x6*x9*cos(x34 - x31) - 1;
    sqr(x6) + sqr(x10) - 2*x6*x10*cos(x35 - x31) - 1;
    sqr(x6) + sqr(x11) - 2*x6*x11*cos(x36 - x31) - 1;
    sqr(x6) + sqr(x12) - 2*x6*x12*cos(x37 - x31) - 1;
    sqr(x6) + sqr(x13) - 2*x6*x13*cos(x38 - x31) - 1;
    sqr(x6) + sqr(x14) - 2*x6*x14*cos(x39 - x31) - 1;
    sqr(x6) + sqr(x15) - 2*x6*x15*cos(x40 - x31) - 1;
    sqr(x6) + sqr(x16) - 2*x6*x16*cos(x41 - x31) - 1;
    sqr(x6) + sqr(x17) - 2*x6*x17*cos(x42 - x31) - 1;
    sqr(x6) + sqr(x18) - 2*x6*x18*cos(x43 - x31) - 1;
    sqr(x6) + sqr(x19) - 2*x6*x19*cos(x44 - x31) - 1;
    sqr(x6) + sqr(x20) - 2*x6*x20*cos(x45 - x31) - 1;
    sqr(x6) + sqr(x21) - 2*x6*x21*cos(x46 - x31) - 1;
    sqr(x6) + sqr(x22) - 2*x6*x22*cos(x47 - x31) - 1;
    sqr(x6) + sqr(x23) - 2*x6*x23*cos(x48 - x31) - 1;
    sqr(x6) + sqr(x24) - 2*x6*x24*cos(x49 - x31) - 1;
    sqr(x6) + sqr(x25) - 2*x6*x25*cos(x50 - x31) - 1;
    sqr(x7) + sqr(x8) - 2*x7*x8*cos(x33 - x32) - 1;
    sqr(x7) + sqr(x9) - 2*x7*x9*cos(x34 - x32) - 1;
    sqr(x7) + sqr(x10) - 2*x7*x10*cos(x35 - x32) - 1;
    sqr(x7) + sqr(x11) - 2*x7*x11*cos(x36 - x32) - 1;
    sqr(x7) + sqr(x12) - 2*x7*x12*cos(x37 - x32) - 1;
    sqr(x7) + sqr(x13) - 2*x7*x13*cos(x38 - x32) - 1;
    sqr(x7) + sqr(x14) - 2*x7*x14*cos(x39 - x32) - 1;
    sqr(x7) + sqr(x15) - 2*x7*x15*cos(x40 - x32) - 1;
    sqr(x7) + sqr(x16) - 2*x7*x16*cos(x41 - x32) - 1;
    sqr(x7) + sqr(x17) - 2*x7*x17*cos(x42 - x32) - 1;
    sqr(x7) + sqr(x18) - 2*x7*x18*cos(x43 - x32) - 1;
    sqr(x7) + sqr(x19) - 2*x7*x19*cos(x44 - x32) - 1;
    sqr(x7) + sqr(x20) - 2*x7*x20*cos(x45 - x32) - 1;
    sqr(x7) + sqr(x21) - 2*x7*x21*cos(x46 - x32) - 1;
    sqr(x7) + sqr(x22) - 2*x7*x22*cos(x47 - x32) - 1;
    sqr(x7) + sqr(x23) - 2*x7*x23*cos(x48 - x32) - 1;
    sqr(x7) + sqr(x24) - 2*x7*x24*cos(x49 - x32) - 1;
    sqr(x7) + sqr(x25) - 2*x7*x25*cos(x50 - x32) - 1;
    sqr(x8) + sqr(x9) - 2*x8*x9*cos(x34 - x33) - 1;
    sqr(x8) + sqr(x10) - 2*x8*x10*cos(x35 - x33) - 1;
    sqr(x8) + sqr(x11) - 2*x8*x11*cos(x36 - x33) - 1;
    sqr(x8) + sqr(x12) - 2*x8*x12*cos(x37 - x33) - 1;
    sqr(x8) + sqr(x13) - 2*x8*x13*cos(x38 - x33) - 1;
    sqr(x8) + sqr(x14) - 2*x8*x14*cos(x39 - x33) - 1;
    sqr(x8) + sqr(x15) - 2*x8*x15*cos(x40 - x33) - 1;
    sqr(x8) + sqr(x16) - 2*x8*x16*cos(x41 - x33) - 1;
    sqr(x8) + sqr(x17) - 2*x8*x17*cos(x42 - x33) - 1;
    sqr(x8) + sqr(x18) - 2*x8*x18*cos(x43 - x33) - 1;
    sqr(x8) + sqr(x19) - 2*x8*x19*cos(x44 - x33) - 1;
    sqr(x8) + sqr(x20) - 2*x8*x20*cos(x45 - x33) - 1;
    sqr(x8) + sqr(x21) - 2*x8*x21*cos(x46 - x33) - 1;
    sqr(x8) + sqr(x22) - 2*x8*x22*cos(x47 - x33) - 1;
    sqr(x8) + sqr(x23) - 2*x8*x23*cos(x48 - x33) - 1;
    sqr(x8) + sqr(x24) - 2*x8*x24*cos(x49 - x33) - 1;
    sqr(x8) + sqr(x25) - 2*x8*x25*cos(x50 - x33) - 1;
    sqr(x9) + sqr(x10) - 2*x9*x10*cos(x35 - x34) - 1;
    sqr(x9) + sqr(x11) - 2*x9*x11*cos(x36 - x34) - 1;
    sqr(x9) + sqr(x12) - 2*x9*x12*cos(x37 - x34) - 1;
    sqr(x9) + sqr(x13) - 2*x9*x13*cos(x38 - x34) - 1;
    sqr(x9) + sqr(x14) - 2*x9*x14*cos(x39 - x34) - 1;
    sqr(x9) + sqr(x15) - 2*x9*x15*cos(x40 - x34) - 1;
    sqr(x9) + sqr(x16) - 2*x9*x16*cos(x41 - x34) - 1;
    sqr(x9) + sqr(x17) - 2*x9*x17*cos(x42 - x34) - 1;
    sqr(x9) + sqr(x18) - 2*x9*x18*cos(x43 - x34) - 1;
    sqr(x9) + sqr(x19) - 2*x9*x19*cos(x44 - x34) - 1;
    sqr(x9) + sqr(x20) - 2*x9*x20*cos(x45 - x34) - 1;
    sqr(x9) + sqr(x21) - 2*x9*x21*cos(x46 - x34) - 1;
    sqr(x9) + sqr(x22) - 2*x9*x22*cos(x47 - x34) - 1;
    sqr(x9) + sqr(x23) - 2*x9*x23*cos(x48 - x34) - 1;
    sqr(x9) + sqr(x24) - 2*x9*x24*cos(x49 - x34) - 1;
    sqr(x9) + sqr(x25) - 2*x9*x25*cos(x50 - x34) - 1;
    sqr(x10) + sqr(x11) - 2*x10*x11*cos(x36 - x35) - 1;
    sqr(x10) + sqr(x12) - 2*x10*x12*cos(x37 - x35) - 1;
    sqr(x10) + sqr(x13) - 2*x10*x13*cos(x38 - x35) - 1;
    sqr(x10) + sqr(x14) - 2*x10*x14*cos(x39 - x35) - 1;
    sqr(x10) + sqr(x15) - 2*x10*x15*cos(x40 - x35) - 1;
    sqr(x10) + sqr(x16) - 2*x10*x16*cos(x41 - x35) - 1;
    sqr(x10) + sqr(x17) - 2*x10*x17*cos(x42 - x35) - 1;
    sqr(x10) + sqr(x18) - 2*x10*x18*cos(x43 - x35) - 1;
    sqr(x10) + sqr(x19) - 2*x10*x19*cos(x44 - x35) - 1;
    sqr(x10) + sqr(x20) - 2*x10*x20*cos(x45 - x35) - 1;
    sqr(x10) + sqr(x21) - 2*x10*x21*cos(x46 - x35) - 1;
    sqr(x10) + sqr(x22) - 2*x10*x22*cos(x47 - x35) - 1;
    sqr(x10) + sqr(x23) - 2*x10*x23*cos(x48 - x35) - 1;
    sqr(x10) + sqr(x24) - 2*x10*x24*cos(x49 - x35) - 1;
    sqr(x10) + sqr(x25) - 2*x10*x25*cos(x50 - x35) - 1;
    sqr(x11) + sqr(x12) - 2*x11*x12*cos(x37 - x36) - 1;
    sqr(x11) + sqr(x13) - 2*x11*x13*cos(x38 - x36) - 1;
    sqr(x11) + sqr(x14) - 2*x11*x14*cos(x39 - x36) - 1;
    sqr(x11) + sqr(x15) - 2*x11*x15*cos(x40 - x36) - 1;
    sqr(x11) + sqr(x16) - 2*x11*x16*cos(x41 - x36) - 1;
    sqr(x11) + sqr(x17) - 2*x11*x17*cos(x42 - x36) - 1;
    sqr(x11) + sqr(x18) - 2*x11*x18*cos(x43 - x36) - 1;
    sqr(x11) + sqr(x19) - 2*x11*x19*cos(x44 - x36) - 1;
    sqr(x11) + sqr(x20) - 2*x11*x20*cos(x45 - x36) - 1;
    sqr(x11) + sqr(x21) - 2*x11*x21*cos(x46 - x36) - 1;
    sqr(x11) + sqr(x22) - 2*x11*x22*cos(x47 - x36) - 1;
    sqr(x11) + sqr(x23) - 2*x11*x23*cos(x48 - x36) - 1;
    sqr(x11) + sqr(x24) - 2*x11*x24*cos(x49 - x36) - 1;
    sqr(x11) + sqr(x25) - 2*x11*x25*cos(x50 - x36) - 1;
    sqr(x12) + sqr(x13) - 2*x12*x13*cos(x38 - x37) - 1;
    sqr(x12) + sqr(x14) - 2*x12*x14*cos(x39 - x37) - 1;
    sqr(x12) + sqr(x15) - 2*x12*x15*cos(x40 - x37) - 1;
    sqr(x12) + sqr(x16) - 2*x12*x16*cos(x41 - x37) - 1;
    sqr(x12) + sqr(x17) - 2*x12*x17*cos(x42 - x37) - 1;
    sqr(x12) + sqr(x18) - 2*x12*x18*cos(x43 - x37) - 1;
    sqr(x12) + sqr(x19) - 2*x12*x19*cos(x44 - x37) - 1;
    sqr(x12) + sqr(x20) - 2*x12*x20*cos(x45 - x37) - 1;
    sqr(x12) + sqr(x21) - 2*x12*x21*cos(x46 - x37) - 1;
    sqr(x12) + sqr(x22) - 2*x12*x22*cos(x47 - x37) - 1;
    sqr(x12) + sqr(x23) - 2*x12*x23*cos(x48 - x37) - 1;
    sqr(x12) + sqr(x24) - 2*x12*x24*cos(x49 - x37) - 1;
    sqr(x12) + sqr(x25) - 2*x12*x25*cos(x50 - x37) - 1;
    sqr(x13) + sqr(x14) - 2*x13*x14*cos(x39 - x38) - 1;
    sqr(x13) + sqr(x15) - 2*x13*x15*cos(x40 - x38) - 1;
    sqr(x13) + sqr(x16) - 2*x13*x16*cos(x41 - x38) - 1;
    sqr(x13) + sqr(x17) - 2*x13*x17*cos(x42 - x38) - 1;
    sqr(x13) + sqr(x18) - 2*x13*x18*cos(x43 - x38) - 1;
    sqr(x13) + sqr(x19) - 2*x13*x19*cos(x44 - x38) - 1;
    sqr(x13) + sqr(x20) - 2*x13*x20*cos(x45 - x38) - 1;
    sqr(x13) + sqr(x21) - 2*x13*x21*cos(x46 - x38) - 1;
    sqr(x13) + sqr(x22) - 2*x13*x22*cos(x47 - x38) - 1;
    sqr(x13) + sqr(x23) - 2*x13*x23*cos(x48 - x38) - 1;
    sqr(x13) + sqr(x24) - 2*x13*x24*cos(x49 - x38) - 1;
    sqr(x13) + sqr(x25) - 2*x13*x25*cos(x50 - x38) - 1;
    sqr(x14) + sqr(x15) - 2*x14*x15*cos(x40 - x39) - 1;
    sqr(x14) + sqr(x16) - 2*x14*x16*cos(x41 - x39) - 1;
    sqr(x14) + sqr(x17) - 2*x14*x17*cos(x42 - x39) - 1;
    sqr(x14) + sqr(x18) - 2*x14*x18*cos(x43 - x39) - 1;
    sqr(x14) + sqr(x19) - 2*x14*x19*cos(x44 - x39) - 1;
    sqr(x14) + sqr(x20) - 2*x14*x20*cos(x45 - x39) - 1;
    sqr(x14) + sqr(x21) - 2*x14*x21*cos(x46 - x39) - 1;
    sqr(x14) + sqr(x22) - 2*x14*x22*cos(x47 - x39) - 1;
    sqr(x14) + sqr(x23) - 2*x14*x23*cos(x48 - x39) - 1;
    sqr(x14) + sqr(x24) - 2*x14*x24*cos(x49 - x39) - 1;
    sqr(x14) + sqr(x25) - 2*x14*x25*cos(x50 - x39) - 1;
    sqr(x15) + sqr(x16) - 2*x15*x16*cos(x41 - x40) - 1;
    sqr(x15) + sqr(x17) - 2*x15*x17*cos(x42 - x40) - 1;
    sqr(x15) + sqr(x18) - 2*x15*x18*cos(x43 - x40) - 1;
    sqr(x15) + sqr(x19) - 2*x15*x19*cos(x44 - x40) - 1;
    sqr(x15) + sqr(x20) - 2*x15*x20*cos(x45 - x40) - 1;
    sqr(x15) + sqr(x21) - 2*x15*x21*cos(x46 - x40) - 1;
    sqr(x15) + sqr(x22) - 2*x15*x22*cos(x47 - x40) - 1;
    sqr(x15) + sqr(x23) - 2*x15*x23*cos(x48 - x40) - 1;
    sqr(x15) + sqr(x24) - 2*x15*x24*cos(x49 - x40) - 1;
    sqr(x15) + sqr(x25) - 2*x15*x25*cos(x50 - x40) - 1;
    sqr(x16) + sqr(x17) - 2*x16*x17*cos(x42 - x41) - 1;
    sqr(x16) + sqr(x18) - 2*x16*x18*cos(x43 - x41) - 1;
    sqr(x16) + sqr(x19) - 2*x16*x19*cos(x44 - x41) - 1;
    sqr(x16) + sqr(x20) - 2*x16*x20*cos(x45 - x41) - 1;
    sqr(x16) + sqr(x21) - 2*x16*x21*cos(x46 - x41) - 1;
    sqr(x16) + sqr(x22) - 2*x16*x22*cos(x47 - x41) - 1;
    sqr(x16) + sqr(x23) - 2*x16*x23*cos(x48 - x41) - 1;
    sqr(x16) + sqr(x24) - 2*x16*x24*cos(x49 - x41) - 1;
    sqr(x16) + sqr(x25) - 2*x16*x25*cos(x50 - x41) - 1;
    sqr(x17) + sqr(x18) - 2*x17*x18*cos(x43 - x42) - 1;
    sqr(x17) + sqr(x19) - 2*x17*x19*cos(x44 - x42) - 1;
    sqr(x17) + sqr(x20) - 2*x17*x20*cos(x45 - x42) - 1;
    sqr(x17) + sqr(x21) - 2*x17*x21*cos(x46 - x42) - 1;
    sqr(x17) + sqr(x22) - 2*x17*x22*cos(x47 - x42) - 1;
    sqr(x17) + sqr(x23) - 2*x17*x23*cos(x48 - x42) - 1;
    sqr(x17) + sqr(x24) - 2*x17*x24*cos(x49 - x42) - 1;
    sqr(x17) + sqr(x25) - 2*x17*x25*cos(x50 - x42) - 1;
    sqr(x18) + sqr(x19) - 2*x18*x19*cos(x44 - x43) - 1;
    sqr(x18) + sqr(x20) - 2*x18*x20*cos(x45 - x43) - 1;
    sqr(x18) + sqr(x21) - 2*x18*x21*cos(x46 - x43) - 1;
    sqr(x18) + sqr(x22) - 2*x18*x22*cos(x47 - x43) - 1;
    sqr(x18) + sqr(x23) - 2*x18*x23*cos(x48 - x43) - 1;
    sqr(x18) + sqr(x24) - 2*x18*x24*cos(x49 - x43) - 1;
    sqr(x18) + sqr(x25) - 2*x18*x25*cos(x50 - x43) - 1;
    sqr(x19) + sqr(x20) - 2*x19*x20*cos(x45 - x44) - 1;
    sqr(x19) + sqr(x21) - 2*x19*x21*cos(x46 - x44) - 1;
    sqr(x19) + sqr(x22) - 2*x19*x22*cos(x47 - x44) - 1;
    sqr(x19) + sqr(x23) - 2*x19*x23*cos(x48 - x44) - 1;
    sqr(x19) + sqr(x24) - 2*x19*x24*cos(x49 - x44) - 1;
    sqr(x19) + sqr(x25) - 2*x19*x25*cos(x50 - x44) - 1;
    sqr(x20) + sqr(x21) - 2*x20*x21*cos(x46 - x45) - 1;
    sqr(x20) + sqr(x22) - 2*x20*x22*cos(x47 - x45) - 1;
    sqr(x20) + sqr(x23) - 2*x20*x23*cos(x48 - x45) - 1;
    sqr(x20) + sqr(x24) - 2*x20*x24*cos(x49 - x45) - 1;
    sqr(x20) + sqr(x25) - 2*x20*x25*cos(x50 - x45) - 1;
    sqr(x21) + sqr(x22) - 2*x21*x22*cos(x47 - x46) - 1;
    sqr(x21) + sqr(x23) - 2*x21*x23*cos(x48 - x46) - 1;
    sqr(x21) + sqr(x24) - 2*x21*x24*cos(x49 - x46) - 1;
    sqr(x21) + sqr(x25) - 2*x21*x25*cos(x50 - x46) - 1;
    sqr(x22) + sqr(x23) - 2*x22*x23*cos(x48 - x47) - 1;
    sqr(x22) + sqr(x24) - 2*x22*x24*cos(x49 - x47) - 1;
    sqr(x22) + sqr(x25) - 2*x22*x25*cos(x50 - x47) - 1;
    sqr(x23) + sqr(x24) - 2*x23*x24*cos(x49 - x48) - 1;
    sqr(x23) + sqr(x25) - 2*x23*x25*cos(x50 - x48) - 1;
    sqr(x24) + sqr(x25) - 2*x24*x25*cos(x50 - x49) - 1;
    x26 - x27;
    x27 - x28;
    x28 - x29;
    x29 - x30;
    x30 - x31;
    x31 - x32;
    x32 - x33;
    x33 - x34;
    x34 - x35;
    x35 - x36;
    x36 - x37;
    x37 - x38;
    x38 - x39;
    x39 - x40;
    x40 - x41;
    x41 - x42;
    x42 - x43;
    x43 - x44;
    x44 - x45;
    x45 - x46;
    x46 - x47;
    x47 - x48;
    x48 - x49;
    x49 - x50;
]

# GEQ (>=0) CONSTRAINTS
geq = []

# BOUNDS
lvbs = [
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.14159265358979,
]
uvbs = [
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 
    3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979,
    3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979,
    3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979,
    3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979,
    3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979,
]

# Create the problem
polygon25 = Problem("polygon25", vars, obj, lvbs, uvbs, eq=eq, leq=leq, geq=geq)