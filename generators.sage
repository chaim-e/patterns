print "loading generators ..."

s = sqrt

GEN = {
    
    (2,):[
        Matrix([[  s(1) ]]),
        Matrix([[  s(1) ]]),
    ],
    (1,1):[
        Matrix([[  -s(1) ]]),
        Matrix([[  -s(1) ]]),
    ],    

    (3,):[
        Matrix([[  s(1) ]]),
        Matrix([[  s(1) ]]),
    ],
    (2,1):[
        Matrix([[  -s(1/4),   s(3/4) ],
                [   s(3/4),   s(1/4) ]]),
        Matrix([[  -s(1/4),  -s(3/4) ],
                [   s(3/4),  -s(1/4) ]]),
    ],
    (1,1,1):[
        Matrix([[  -s(1) ]]),
        Matrix([[   s(1) ]]),
    ],

    (4,):[
        Matrix([[  s(1) ]]),
        Matrix([[  s(1) ]]),
    ],
    (3,1):[
        Matrix([[    s(1/25),     s(4/5),   -s(4/25) ],
                [     s(4/5),       s(0),     s(1/5) ],
                [   -s(4/25),     s(1/5),   s(16/25) ]]),
        Matrix([[  -s(16/25),    -s(1/5),   -s(4/25) ],
                [     s(1/5),       s(0),    -s(4/5) ],
                [   -s(4/25),     s(4/5),   -s(1/25) ]]),
    ],
    (2,2):[
        Matrix([[    -s(1),     s(0) ],
                [     s(0),     s(1) ]]),
        Matrix([[   s(1/4),   s(3/4) ],
                [   s(3/4),  -s(1/4) ]]),
    ],
    (2,1,1):[
        Matrix([[   s(0),   s(0),   s(1) ],
                [   s(0),  -s(1),   s(0) ],
                [   s(1),   s(0),   s(0) ]]),
        Matrix([[   s(0),   s(1),   s(0) ],
                [  -s(1),   s(0),   s(0) ],
                [   s(0),   s(0),   s(1) ]]),
    ],
    (1,1,1,1):[
        Matrix([[  -s(1) ]]),
        Matrix([[  -s(1) ]]),
    ],

    (5,):[
        Matrix([[  s(1) ]]),
        Matrix([[  s(1) ]]),
    ],
    (4,1):[
        Matrix([[     s(1/100),    -s(9/100),      s(9/28),    s(81/140) ],
                [    -s(9/100),    s(81/100),      s(1/28),     s(9/140) ],
                [      s(9/28),      s(1/28),    s(81/196),   -s(45/196) ],
                [    s(81/140),     s(9/140),   -s(45/196),    s(25/196) ]]),
        Matrix([[      -s(1/4),      -s(1/4),      s(9/28),     -s(5/28) ],
                [      -s(1/4),         s(0),      s(1/28),       s(5/7) ],
                [     -s(9/28),     -s(1/28),  -s(121/196),    -s(5/196) ],
                [      s(5/28),      -s(5/7),    -s(5/196),      s(4/49) ]]),
    ],
    (3,2):[
        Matrix([[   -s(36/49),        s(0),     -s(1/7),        s(0),    -s(6/49) ],
                [        s(0),        s(1),        s(0),        s(0),        s(0) ],
                [     -s(1/7),        s(0),        s(0),        s(0),      s(6/7) ],
                [        s(0),        s(0),        s(0),        s(1),        s(0) ],
                [    -s(6/49),        s(0),      s(6/7),        s(0),    -s(1/49) ]]),
        Matrix([[     s(4/49),     s(9/28),     -s(4/7),    -s(1/56),   -s(3/392) ],
                [    -s(9/28),     -s(1/4),     -s(1/4),     -s(1/8),    -s(3/56) ],
                [      s(4/7),     -s(1/4),        s(0),     -s(1/8),    -s(3/56) ],
                [    -s(1/56),      s(1/8),      s(1/8),    -s(1/16),  -s(75/112) ],
                [   -s(3/392),     s(3/56),     s(3/56),  -s(75/112),  s(169/784) ]]),
    ],
    (3,1,1):[
        Matrix([[       s(4/25),      -s(6/25),        s(3/5),          s(0),          s(0),          s(0) ],
                [      -s(6/25),       s(9/25),        s(2/5),          s(0),          s(0),          s(0) ],
                [        s(3/5),        s(2/5),          s(0),          s(0),          s(0),          s(0) ],
                [          s(0),          s(0),          s(0),          s(0),        s(1/7),       -s(6/7) ],
                [          s(0),          s(0),          s(0),        s(1/7),     -s(36/49),      -s(6/49) ],
                [          s(0),          s(0),          s(0),       -s(6/7),      -s(6/49),      -s(1/49) ]]),
        Matrix([[          s(1),          s(0),          s(0),          s(0),          s(0),          s(0) ],
                [          s(0),      -s(1/16),      -s(5/72),     -s(25/72),      s(25/56),     s(25/336) ],
                [          s(0),       s(5/72),       s(1/36),       s(5/36),      s(5/252),    s(125/168) ],
                [          s(0),      s(25/72),       s(5/36),       -s(4/9),      -s(4/63),      s(1/168) ],
                [          s(0),      s(25/56),     -s(5/252),       s(4/63),      s(16/49),  -s(169/1176) ],
                [          s(0),     s(25/336),   -s(125/168),     -s(1/168),  -s(169/1176),     s(25/784) ]]),
    ],
    (2,2,1):[
        Matrix([[    -s(1/25),     s(4/25),        s(0),     s(4/15),     s(8/15) ],
                [     s(4/25),   -s(16/25),        s(0),     s(1/15),     s(2/15) ],
                [        s(0),        s(0),        s(1),        s(0),        s(0) ],
                [     s(4/15),     s(1/15),        s(0),     -s(4/9),      s(2/9) ],
                [     s(8/15),     s(2/15),        s(0),      s(2/9),     -s(1/9) ]]),
        Matrix([[    s(9/100),   -s(1/100),    s(27/40),     s(3/20),     s(3/40) ],
                [   -s(1/100),   -s(16/25),    -s(3/40),     s(4/15),   -s(1/120) ],
                [   -s(27/40),     s(3/40),     s(1/16),      s(1/8),    -s(1/16) ],
                [    -s(3/20),    -s(4/15),      s(1/8),     -s(4/9),     s(1/72) ],
                [    -s(3/40),    s(1/120),    -s(1/16),     s(1/72),  s(121/144) ]]),
    ],
    (2,1,1,1):[
        Matrix([[   -s(1/36),    s(5/12),    s(5/36),    s(5/12) ],
                [    s(5/12),    -s(1/4),    s(1/12),     s(1/4) ],
                [    s(5/36),    s(1/12),  -s(25/36),    s(1/12) ],
                [    s(5/12),     s(1/4),    s(1/12),    -s(1/4) ]]),
        Matrix([[    -s(4/9),       s(0),    s(5/36),    s(5/12) ],
                [       s(0),       s(0),     s(3/4),    -s(1/4) ],
                [    s(5/36),    -s(3/4),    s(1/36),    s(1/12) ],
                [   -s(5/12),    -s(1/4),   -s(1/12),    -s(1/4) ]]),
    ],
    (1,1,1,1,1):[
        Matrix([[  -s(1) ]]),
        Matrix([[   s(1) ]]),
    ],

    (6,):[
        Matrix([[  s(1) ]]),
        Matrix([[  s(1) ]]),
    ],
    (5,1):[
        Matrix([[         s(16/49),         s(12/35),         s(6/245),        -s(12/49),          s(3/49) ],
                [         s(12/35),          s(1/25),        -s(8/175),         s(16/35),         -s(4/35) ],
                [         s(6/245),        -s(8/175),     s(1089/1225),         s(8/245),        -s(2/245) ],
                [        -s(12/49),         s(16/35),         s(8/245),          s(9/49),          s(4/49) ],
                [          s(3/49),         -s(4/35),        -s(2/245),          s(4/49),         s(36/49) ]]),
        Matrix([[       s(169/784),        s(15/112),        -s(30/49),       -s(27/784),         s(3/784) ],
                [       -s(15/112),     -s(121/1296),          -s(2/7),     s(2209/5040),   -s(2209/45360) ],
                [         s(30/49),          -s(2/7),          s(1/49),        s(18/245),        -s(2/245) ],
                [       -s(27/784),    -s(2209/5040),       -s(18/245),      -s(169/784),     s(1681/7056) ],
                [        -s(3/784),   -s(2209/45360),        -s(2/245),    -s(1681/7056),  -s(44521/63504) ]]),
    ],
    (4,2):[
        Matrix([[             s(0),             s(0),             s(0),             s(0),          s(1/21),        -s(32/63),             s(0),          -s(1/9),          -s(1/3) ],
                [             s(0),             s(1),             s(0),             s(0),             s(0),             s(0),             s(0),             s(0),             s(0) ],
                [             s(0),             s(0),             s(0),             s(0),        -s(25/42),         -s(1/63),             s(0),        -s(25/72),          s(1/24) ],
                [             s(0),             s(0),             s(0),             s(1),             s(0),             s(0),             s(0),             s(0),             s(0) ],
                [          s(1/21),             s(0),        -s(25/42),             s(0),          s(4/49),        s(25/294),             s(0),         -s(4/21),             s(0) ],
                [        -s(32/63),             s(0),         -s(1/63),             s(0),        s(25/294),         -s(4/49),             s(0),         -s(1/56),          s(7/24) ],
                [             s(0),             s(0),             s(0),             s(0),             s(0),             s(0),             s(1),             s(0),             s(0) ],
                [          -s(1/9),             s(0),        -s(25/72),             s(0),         -s(4/21),         -s(1/56),             s(0),           s(1/4),         -s(1/12) ],
                [          -s(1/3),             s(0),          s(1/24),             s(0),             s(0),          s(7/24),             s(0),         -s(1/12),          -s(1/4) ]]),
        Matrix([[       -s(25/144),          s(1/48),          s(1/18),      s(605/3024),        -s(1/336),       -s(1/1134),        -s(1/168),        -s(25/81),        s(25/108) ],
                [          s(1/48),         -s(1/16),           s(1/6),      s(125/1008),        -s(9/112),       s(121/378),         -s(9/56),          s(1/27),         -s(1/36) ],
                [          s(1/18),           s(1/6),          s(1/36),        -s(5/378),          s(2/21),      s(841/2268),          s(4/21),       -s(49/648),        -s(1/216) ],
                [     -s(605/3024),     -s(125/1008),         s(5/378),  -s(20449/63504),        -s(5/784),     s(605/23814),        -s(5/392),     -s(320/1701),       -s(35/324) ],
                [         s(1/336),         s(9/112),         -s(2/21),        -s(5/784),       s(169/784),         s(2/147),      -s(225/392),         -s(1/84),             s(0) ],
                [        s(1/1134),      -s(121/378),     -s(841/2268),     s(605/23814),         s(2/147),  s(20449/142884),         s(4/147),    -s(361/40824),      s(175/1944) ],
                [         s(1/168),          s(9/56),         -s(4/21),        -s(5/392),      -s(225/392),         s(4/147),        -s(1/196),         -s(1/42),             s(0) ],
                [         s(25/81),         -s(1/27),        s(49/648),     -s(320/1701),         -s(1/84),    -s(361/40824),         -s(1/42),    -s(361/11664),     s(1225/3888) ],
                [       -s(25/108),          s(1/36),         s(1/216),       -s(35/324),             s(0),      s(175/1944),             s(0),     s(1225/3888),      s(289/1296) ]]),
    ],
    (4,1,1):[
        Matrix([[           s(0),        s(1/21),           s(0),       s(32/63),           s(0),           s(0),           s(0),           s(0),         s(4/9),           s(0) ],
                [        s(1/21),        s(4/49),           s(0),     -s(25/294),       s(25/42),           s(0),           s(0),           s(0),        s(1/21),         s(1/7) ],
                [           s(0),           s(0),        s(9/25),           s(0),           s(0),        s(2/35),       s(12/35),       -s(6/25),           s(0),           s(0) ],
                [       s(32/63),     -s(25/294),           s(0),       -s(4/49),       -s(1/63),           s(0),           s(0),           s(0),        s(9/56),     -s(25/168) ],
                [           s(0),       s(25/42),           s(0),       -s(1/63),           s(0),           s(0),           s(0),           s(0),       -s(1/72),        -s(3/8) ],
                [           s(0),           s(0),        s(2/35),           s(0),           s(0),       s(36/49),       -s(6/49),        s(3/35),           s(0),           s(0) ],
                [           s(0),           s(0),       s(12/35),           s(0),           s(0),       -s(6/49),        s(1/49),       s(18/35),           s(0),           s(0) ],
                [           s(0),           s(0),       -s(6/25),           s(0),           s(0),        s(3/35),       s(18/35),        s(4/25),           s(0),           s(0) ],
                [         s(4/9),        s(1/21),           s(0),        s(9/56),       -s(1/72),           s(0),           s(0),           s(0),        -s(1/4),        s(1/12) ],
                [           s(0),         s(1/7),           s(0),     -s(25/168),        -s(3/8),           s(0),           s(0),           s(0),        s(1/12),         s(1/4) ]]),
        Matrix([[        s(1/64),    -s(25/1344),           s(0),       -s(1/56),        -s(1/8),     s(121/448),      -s(1/168),       s(5/192),      -s(25/64),      s(25/192) ],
                [     s(25/1344),    s(289/3136),           s(0),     -s(75/392),       s(1/168),    -s(27/3136),    -s(225/392),      s(45/448),       s(3/448),      -s(1/448) ],
                [           s(0),           s(0),           s(1),           s(0),           s(0),           s(0),           s(0),           s(0),           s(0),           s(0) ],
                [        s(1/56),     -s(75/392),           s(0),  s(4225/15876),       -s(1/28),   -s(121/3528),      -s(4/147),    s(605/1512),       s(2/567),       -s(1/42) ],
                [        -s(1/8),      -s(1/168),           s(0),        s(1/28),        -s(1/4),        s(1/56),       -s(4/21),       -s(5/24),           s(0),        -s(1/6) ],
                [     s(121/448),     s(27/3136),           s(0),    s(121/3528),        s(1/56),   -s(529/3136),      -s(3/392),   -s(125/1344),  -s(1369/4032),     -s(27/448) ],
                [      -s(1/168),     s(225/392),           s(0),       s(4/147),       -s(4/21),      -s(3/392),        s(4/49),        s(5/56),      -s(1/168),       -s(1/56) ],
                [      -s(5/192),      s(45/448),           s(0),    s(605/1512),        s(5/24),    s(125/1344),       -s(5/56),      -s(1/576),     -s(5/1728),        s(5/64) ],
                [       s(25/64),       s(3/448),           s(0),       s(2/567),           s(0),   s(1369/4032),       s(1/168),     -s(5/1728),    s(625/5184),     -s(25/192) ],
                [     -s(25/192),      -s(1/448),           s(0),       -s(1/42),         s(1/6),      s(27/448),        s(1/56),        s(5/64),     -s(25/192),      -s(25/64) ]]),
    ],
    (3,3):[
        Matrix([[     s(36/49),       s(1/7),         s(0),     -s(6/49),         s(0) ],
                [       s(1/7),         s(0),         s(0),       s(6/7),         s(0) ],
                [         s(0),         s(0),         s(1),         s(0),         s(0) ],
                [     -s(6/49),       s(6/7),         s(0),      s(1/49),         s(0) ],
                [         s(0),         s(0),         s(0),         s(0),        -s(1) ]]),
        Matrix([[    -s(1/196),     -s(9/28),       s(3/7),     -s(3/98),     -s(3/14) ],
                [      s(9/28),       s(1/4),       s(1/3),      s(2/21),         s(0) ],
                [       s(3/7),      -s(1/3),      -s(1/9),      s(8/63),         s(0) ],
                [     -s(3/98),     -s(2/21),      s(8/63),  s(169/7056),    s(81/112) ],
                [     -s(3/14),         s(0),         s(0),    s(81/112),     -s(1/16) ]]),
    ],
    (3,2,1):[
        Matrix([[      -s(625/36864),        -s(15/2048),      s(1681/61440),         -s(25/512),          s(1/3584),        s(289/2240),       -s(105/2048),               s(0),               s(0),          s(35/576),      s(4205/28672),        -s(361/896),               s(0),       -s(405/4096),           s(5/512),               s(0) ],
                [        -s(15/2048),        s(225/1024),        -s(81/2048),          s(15/256),        -s(15/1792),         -s(27/224),        s(841/7168),               s(0),               s(0),          s(75/224),       s(675/14336),         -s(15/448),               s(0),          s(3/2048),          -s(3/256),               s(0) ],
                [      s(1681/61440),        -s(81/2048),    -s(4761/102400),        -s(27/2560),      -s(507/17920),       s(867/11200),         s(63/2048),               s(0),               s(0),           s(7/192),     -s(5043/28672),         -s(15/896),               s(0),       -s(147/4096),        -s(243/512),               s(0) ],
                [         -s(25/512),          s(15/256),        -s(27/2560),           s(25/64),            s(7/64),            s(7/40),        -s(15/1792),               s(0),               s(0),           -s(5/56),         s(45/3584),          -s(1/112),               s(0),           s(5/512),           -s(5/64),               s(0) ],
                [          s(1/3584),        -s(15/1792),      -s(507/17920),            s(7/64),           -s(1/64),            s(1/40),      s(2535/12544),               s(0),               s(0),          -s(5/392),     -s(1805/25088),          -s(9/784),               s(0),       -s(845/3584),         s(125/448),               s(0) ],
                [        s(289/2240),         -s(27/224),       s(867/11200),            s(7/40),            s(1/40),            s(1/25),        -s(27/1568),               s(0),               s(0),           s(16/49),         -s(1/3136),            s(5/98),               s(0),           s(9/448),            s(1/56),               s(0) ],
                [       -s(105/2048),        s(841/7168),         s(63/2048),        -s(15/1792),      s(2535/12544),        -s(27/1568),     -s(4761/50176),               s(0),               s(0),         s(75/1568),   -s(29403/100352),         s(15/3136),               s(0),     -s(1875/14336),          s(3/1792),               s(0) ],
                [               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(1),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0) ],
                [               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(1),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0) ],
                [          s(35/576),          s(75/224),           s(7/192),           -s(5/56),          -s(5/392),           s(16/49),         s(75/1568),               s(0),               s(0),           s(1/441),         -s(1/3136),            s(5/98),               s(0),           s(9/448),            s(1/56),               s(0) ],
                [      s(4205/28672),       s(675/14336),     -s(5043/28672),         s(45/3584),     -s(1805/25088),         -s(1/3136),   -s(29403/100352),               s(0),               s(0),         -s(1/3136),   -s(13225/200704),       -s(605/6272),               s(0),        s(175/4096),        s(169/3584),               s(0) ],
                [        -s(361/896),         -s(15/448),         -s(15/896),          -s(1/112),          -s(9/784),            s(5/98),         s(15/3136),               s(0),               s(0),            s(5/98),       -s(605/6272),          -s(1/196),               s(0),          s(35/128),           s(5/112),               s(0) ],
                [               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),              -s(1),               s(0),               s(0),               s(0) ],
                [       -s(405/4096),          s(3/2048),       -s(147/4096),           s(5/512),       -s(845/3584),           s(9/448),     -s(1875/14336),               s(0),               s(0),           s(9/448),        s(175/4096),          s(35/128),               s(0),       -s(529/4096),          -s(1/512),               s(0) ],
                [           s(5/512),          -s(3/256),        -s(243/512),           -s(5/64),         s(125/448),            s(1/56),          s(3/1792),               s(0),               s(0),            s(1/56),        s(169/3584),           s(5/112),               s(0),          -s(1/512),            s(1/64),               s(0) ],
                [               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),              -s(1) ]]),
        Matrix([[     -s(2809/16384),        -s(15/8192),    -s(17787/81920),          s(1/2048),       s(169/14336),     -s(4761/35840),       -s(15/57344),        s(125/1024),         -s(5/1024),        s(125/7168),      s(845/114688),       -s(529/3584),           s(9/128),     -s(1125/16384),         -s(5/2048),           s(3/128) ],
                [         s(15/8192),        s(361/4096),          s(9/8192),        -s(15/1024),        s(135/7168),        -s(27/3584),    -s(13689/28672),          -s(3/512),         -s(27/512),         -s(3/3584),       -s(189/8192),               s(0),          s(15/256),         -s(3/8192),        -s(75/1024),         -s(45/256) ],
                [     s(17787/81920),          s(9/8192),  -s(124609/409600),        s(27/10240),        -s(3/71680),    -s(2187/179200),         s(9/57344),        -s(75/1024),          s(3/1024),        -s(75/7168),   -s(14283/114688),          s(3/4480),         -s(15/512),      -s(243/16384),       -s(243/2048),          s(45/512) ],
                [          s(1/2048),         s(15/1024),       -s(27/10240),         -s(25/256),       -s(169/1792),         -s(1/4480),       -s(135/7168),           s(5/128),          s(45/128),           s(5/896),        s(45/14336),            -s(1/7),           -s(9/64),         s(45/2048),          -s(5/256),           -s(3/64) ],
                [      -s(169/14336),        s(135/7168),        -s(3/71680),        s(169/1792),    s(49729/112896),      s(1089/31360),    -s(1805/150528),          -s(5/896),       s(2645/8064),         -s(5/6272),    -s(1125/100352),         s(1/28224),               s(0),        -s(35/2048),         s(45/1792),               s(0) ],
                [      s(4761/35840),        -s(27/3584),    -s(2187/179200),          s(1/4480),      s(1089/31360),     s(14641/78400),       -s(27/25088),         s(225/448),          -s(9/448),        s(225/3136),      -s(169/50176),         s(81/7840),               s(0),          s(7/1024),          -s(9/896),               s(0) ],
                [       -s(15/57344),     s(13689/28672),        -s(9/57344),       -s(135/7168),     s(1805/150528),        s(27/25088),    s(29929/200704),          s(3/3584),      -s(529/10752),         s(3/25088),   -s(45387/401408),     -s(2645/37632),               s(0),      s(3267/57344),        s(363/7168),               s(0) ],
                [        s(125/1024),           s(3/512),         s(75/1024),           s(5/128),           s(5/896),        -s(225/448),          s(3/3584),            s(9/64),            s(1/64),         -s(25/448),         -s(1/7168),           s(5/224),               s(0),          s(9/1024),           s(1/128),               s(0) ],
                [          s(5/1024),         -s(27/512),          s(3/1024),         -s(45/128),       s(2645/8064),          -s(9/448),       s(529/10752),           -s(1/64),          -s(1/576),          -s(1/448),        s(361/7168),         -s(5/2016),               s(0),         s(49/1024),          -s(9/128),               s(0) ],
                [        s(125/7168),          s(3/3584),         s(75/7168),           s(5/896),          s(5/6272),       -s(225/3136),         s(3/25088),         -s(25/448),           s(1/448),       s(2601/3136),        -s(1/50176),          s(5/1568),               s(0),          s(9/7168),           s(1/896),               s(0) ],
                [     -s(845/114688),       -s(189/8192),   -s(14283/114688),       -s(45/14336),    -s(1125/100352),      -s(169/50176),    s(45387/401408),          s(1/7168),        s(361/7168),         s(1/50176),   -s(46225/802816),        s(845/6272),        s(405/3584),     s(2025/114688),        s(25/14336),      -s(1215/3584) ],
                [        s(529/3584),               s(0),          s(3/4480),             s(1/7),         s(1/28224),         s(81/7840),      s(2645/37632),          -s(5/224),         -s(5/2016),         -s(5/1568),        s(845/6272),   -s(27889/112896),         s(81/1792),        -s(45/3584),        -s(45/1792),       -s(243/1792) ],
                [          -s(9/128),          s(15/256),         -s(15/512),            s(9/64),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),        s(405/3584),         s(81/1792),          -s(1/256),          s(45/128),         -s(45/256),           s(3/256) ],
                [     -s(1125/16384),          s(3/8192),       s(243/16384),         s(45/2048),         s(35/2048),         -s(7/1024),      s(3267/57344),          s(9/1024),        -s(49/1024),          s(9/7168),    -s(2025/114688),         s(45/3584),         -s(45/128),     -s(2401/16384),       -s(225/2048),         -s(15/128) ],
                [         -s(5/2048),         s(75/1024),        s(243/2048),          -s(5/256),        -s(45/1792),           s(9/896),        s(363/7168),           s(1/128),           s(9/128),           s(1/896),       -s(25/14336),         s(45/1792),          s(45/256),       -s(225/2048),            -s(1/4),          s(15/256) ],
                [           s(3/128),          s(45/256),         -s(45/512),           -s(3/64),               s(0),               s(0),               s(0),               s(0),               s(0),               s(0),       s(1215/3584),        s(243/1792),          -s(3/256),         -s(15/128),          s(15/256),          -s(1/256) ]]),
    ],
    (3,1,1,1):[
        Matrix([[        s(4/25),           s(0),         s(1/7),           s(0),       -s(4/35),       s(16/35),       -s(1/25),        s(1/20),      -s(3/280),        s(1/40) ],
                [           s(0),      -s(9/400),     -s(12/175),        s(3/80),     s(121/420),       s(3/140),        s(1/48),       s(27/80),       -s(9/70),        s(3/40) ],
                [         s(1/7),     -s(12/175),   -s(256/1225),        s(4/35),      -s(4/245),      s(16/245),         s(1/7),      -s(9/140),     -s(3/1960),       -s(7/40) ],
                [           s(0),        s(3/80),        s(4/35),       -s(1/16),       -s(1/28),       -s(1/28),        s(5/16),        s(1/16),       -s(3/14),        -s(1/8) ],
                [       -s(4/35),     s(121/420),      -s(4/245),       -s(1/28),      s(25/441),       s(16/49),    s(121/1260),       -s(1/28),        s(3/98),           s(0) ],
                [       s(16/35),       s(3/140),      s(16/245),       -s(1/28),       s(16/49),       -s(1/49),       s(1/140),       -s(1/28),        s(3/98),           s(0) ],
                [       -s(1/25),        s(1/48),         s(1/7),        s(5/16),    s(121/1260),       s(1/140),   -s(361/3600),       -s(1/80),       -s(3/70),       -s(9/40) ],
                [        s(1/20),       s(27/80),      -s(9/140),        s(1/16),       -s(1/28),       -s(1/28),       -s(1/80),       -s(1/16),       -s(3/14),         s(1/8) ],
                [      -s(3/280),       -s(9/70),     -s(3/1960),       -s(3/14),        s(3/98),        s(3/98),       -s(3/70),       -s(3/14),      -s(16/49),           s(0) ],
                [        s(1/40),        s(3/40),       -s(7/40),        -s(1/8),           s(0),           s(0),       -s(9/40),         s(1/8),           s(0),        -s(1/4) ]]),
        Matrix([[      -s(9/100),      -s(3/400),       s(1/700),        s(5/16),           s(0),       -s(1/35),       s(9/400),       -s(5/16),           s(0),        s(9/40) ],
                [       s(3/400),   -s(169/6400),     s(27/2800),     s(27/1280),      -s(15/28),     -s(27/140),   -s(147/6400),    s(147/1280),     -s(45/896),       s(3/160) ],
                [      -s(1/700),     s(27/2800),   s(4489/4900),       s(1/560),           s(0),       s(9/245),    -s(81/2800),      -s(1/560),           s(0),      -s(1/280) ],
                [       -s(5/16),     s(27/1280),       s(1/560),     -s(49/256),        s(1/28),       -s(1/28),      -s(5/256),      s(25/256),       s(3/896),        s(9/32) ],
                [           s(0),       s(15/28),           s(0),       -s(1/28),       -s(9/49),           s(0),        s(5/28),       -s(1/28),        s(3/98),           s(0) ],
                [        s(1/35),     -s(27/140),       s(9/245),       -s(1/28),           s(0),        s(1/49),      s(81/140),        s(1/28),           s(0),        s(1/14) ],
                [       s(9/400),    s(147/6400),     s(81/2800),       s(5/256),        s(5/28),     -s(81/140),    s(361/6400),       s(5/256),      s(15/896),      -s(9/160) ],
                [       -s(5/16),   -s(147/1280),       s(1/560),     -s(25/256),       -s(1/28),       -s(1/28),       s(5/256),     -s(25/256),      -s(3/896),       -s(9/32) ],
                [           s(0),      s(45/896),           s(0),      -s(3/896),        s(3/98),           s(0),      s(15/896),      -s(3/896),  -s(2809/3136),           s(0) ],
                [        s(9/40),      -s(3/160),       s(1/280),       -s(9/32),           s(0),       -s(1/14),      -s(9/160),       -s(9/32),           s(0),        s(1/16) ]]),
    ],
    (2,2,2):[
        Matrix([[  -s(4/225),     s(5/9),   -s(1/15),   -s(9/25),       s(0) ],
                [     s(5/9),     s(1/9),     s(1/3),       s(0),       s(0) ],
                [   -s(1/15),     s(1/3),       s(0),     s(3/5),       s(0) ],
                [   -s(9/25),       s(0),     s(3/5),   -s(1/25),       s(0) ],
                [       s(0),       s(0),       s(0),       s(0),      -s(1) ]]),
        Matrix([[    -s(1/4),       s(0),   -s(3/20),       s(0),     s(3/5) ],
                [       s(0),       s(1),       s(0),       s(0),       s(0) ],
                [    s(3/20),       s(0),    -s(1/4),     s(3/5),       s(0) ],
                [       s(0),       s(0),    -s(3/5),    -s(1/4),   -s(3/20) ],
                [     s(3/5),       s(0),       s(0),   -s(3/20),     s(1/4) ]]),
    ],
    (2,2,1,1):[
        Matrix([[       s(1/25),       s(4/15),       -s(2/5),          s(0),      -s(2/15),       s(4/25),          s(0),          s(0),          s(0) ],
                [       s(4/15),          s(0),       s(8/27),          s(0),          s(0),     s(49/135),          s(0),      -s(2/27),          s(0) ],
                [       -s(2/5),       s(8/27),       s(4/81),          s(0),       s(1/27),      s(8/405),          s(0),     -s(16/81),          s(0) ],
                [          s(0),          s(0),          s(0),         -s(1),          s(0),          s(0),          s(0),          s(0),          s(0) ],
                [      -s(2/15),          s(0),       s(1/27),          s(0),          s(0),     s(32/135),          s(0),      s(16/27),          s(0) ],
                [       s(4/25),     s(49/135),      s(8/405),          s(0),     s(32/135),  -s(196/2025),          s(0),      s(10/81),          s(0) ],
                [          s(0),          s(0),          s(0),          s(0),          s(0),          s(0),         -s(1),          s(0),          s(0) ],
                [          s(0),      -s(2/27),     -s(16/81),          s(0),      s(16/27),      s(10/81),          s(0),      -s(1/81),          s(0) ],
                [          s(0),          s(0),          s(0),          s(0),          s(0),          s(0),          s(0),          s(0),         -s(1) ]]),
        Matrix([[       s(1/25),       s(1/15),       s(8/45),          s(0),      -s(8/15),      s(1/225),          s(0),       s(8/45),          s(0) ],
                [      -s(1/15),      -s(1/36),       s(1/54),       -s(1/6),      -s(1/18),    -s(49/540),        s(1/2),      -s(2/27),          s(0) ],
                [       s(8/45),      -s(1/54),  -s(289/1296),    -s(25/144),    -s(49/432),     -s(1/810),    -s(25/432),      -s(4/81),       s(5/27) ],
                [          s(0),       -s(1/6),     s(25/144),       s(1/16),      -s(1/48),      -s(5/18),      -s(3/16),       -s(1/9),          s(0) ],
                [       s(8/15),      -s(1/18),     s(49/432),      -s(1/48),     s(25/144),     -s(1/270),       s(1/16),       s(1/27),          s(0) ],
                [      s(1/225),     s(49/540),     -s(1/810),       s(5/18),      s(1/270),  -s(529/8100),       s(5/54),     -s(1/810),      s(25/54) ],
                [          s(0),       -s(1/2),    -s(25/432),       s(3/16),      -s(1/16),       s(5/54),       s(1/16),       s(1/27),          s(0) ],
                [       s(8/45),       s(2/27),      -s(4/81),        s(1/9),      -s(1/27),     -s(1/810),       s(1/27),  -s(289/1296),   -s(125/432) ],
                [          s(0),          s(0),       s(5/27),          s(0),          s(0),      s(25/54),          s(0),   -s(125/432),       s(1/16) ]]),
    ],
    (2,1,1,1,1):[
        Matrix([[     -s(1),      s(0),      s(0),      s(0),      s(0) ],
                [      s(0),     -s(1),      s(0),      s(0),      s(0) ],
                [      s(0),      s(0),      s(0),    s(1/3),    s(2/3) ],
                [      s(0),      s(0),    s(1/3),   -s(4/9),    s(2/9) ],
                [      s(0),      s(0),    s(2/3),    s(2/9),   -s(1/9) ]]),
        Matrix([[  -s(1/16),   s(3/16),  -s(9/16),  -s(3/16),      s(0) ],
                [  -s(3/16),  -s(1/16),   s(3/16),  -s(9/16),      s(0) ],
                [   s(9/16),   s(3/16),   s(1/16),  -s(3/16),      s(0) ],
                [  -s(3/16),   s(9/16),   s(3/16),   s(1/16),      s(0) ],
                [      s(0),      s(0),      s(0),      s(0),      s(1) ]]),
    ],
    (1,1,1,1,1,1):[
        Matrix([[  -s(1) ]]),
        Matrix([[  -s(1) ]]),
    ],

}
