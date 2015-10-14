/* The Image Registration Toolkit (IRTK)
 *
 * Copyright 2008-2015 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. */

#include <irtkTransformation.h>
#include <irtkFreeFormTransformationRungeKutta.h>


#define BT_R(...) {__VA_ARGS__}
#define BT_A(...) {__VA_ARGS__}
#define BT_B(b)   b
#define BT_C(c)   c


// ---------------------------------------------------------------------------
// Forward Euler method
IRTK_DEFINE_FFDRK_EXPLICIT(RKE1, 1, 1, false,
  BT_C(BT_R(0.0)),
  /* ------------- */
  BT_A(BT_R(0.0)),
  /* ------------- */
  BT_B(BT_R(1.0))
);

// ---------------------------------------------------------------------------
// Modified Euler method (Heun's method, explicit midpoint rule)
IRTK_DEFINE_FFDRK_EXPLICIT(RKE2, 2, 2, false,
  BT_C(BT_R(0.0, 0.5)),
  /* ----------------- */
  BT_A(BT_R(0.0, 0.0),
       BT_R(0.5, 0.0)),
  /* ----------------- */
  BT_B(BT_R(0.0, 1.0))
);

// ---------------------------------------------------------------------------
// Improved Euler method (Heun's method, explicit trapezoidal rule)
IRTK_DEFINE_FFDRK_EXPLICIT(RKH2, 2, 2, false,
  BT_C(BT_R(0.0, 0.5)),
  /* ----------------- */
  BT_A(BT_R(0.0, 0.0),
       BT_R(1.0, 0.0)),
  /* ----------------- */
  BT_B(BT_R(0.5, 0.5))
);

// ---------------------------------------------------------------------------
// Classical Runge-Kutta method
IRTK_DEFINE_FFDRK_EXPLICIT(RK4, 4, 4, false,
  BT_C(BT_R(0.0,     1.0/2.0, 1.0/2.0, 1.0)),
  /* ----------------------------------------- */
  BT_A(BT_R(0.0,     0.0,     0.0,     0.0),
       BT_R(1.0/2.0, 0.0,     0.0,     0.0),
       BT_R(0.0,     1.0/2.0, 0.0,     0.0),
       BT_R(0.0,     0.0,     1.0,     0.0)),
  /* ----------------------------------------- */
  BT_B(BT_R(1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0))
);


// ---------------------------------------------------------------------------
// Euler-Heun method of order 2(1)
IRTK_DEFINE_FFDRK_EMBEDDED(RKEH12, 2, 2, false,
  BT_C(BT_R(0.0, 1.0)),
  /* ----------------- */
  BT_A(BT_R(0.0, 0.0),
       BT_R(1.0, 0.0)),
  /* ----------------- */
  BT_B(BT_R(0.5, 0.5)),
  BT_B(BT_R(1.0, 0.0))
);

// ---------------------------------------------------------------------------
// Bogacki-Shampine method of order 3(2)
IRTK_DEFINE_FFDRK_EMBEDDED(RKBS23, 4, 3, true,
  BT_C(BT_R(0.0,      1.0/2.0,  3.0/4.0, 1.0    )),
  /* ---------------------------------------------- */
  BT_A(BT_R(0.0,      0.0,      0.0,     0.0    ),
       BT_R(1.0/2.0,  0.0,      0.0,     0.0    ),
       BT_R(0.0,      3.0/4.0,  0.0,     0.0    ),
       BT_R(2.0/9.0,  1.0/3.0,  4.0/9.0, 0.0    )),
  /* ---------------------------------------------- */
  BT_B(BT_R(2.0/9.0,  2.0/3.0,  4.0/9.0, 0.0    )),
  BT_B(BT_R(7.0/24.0, 1.0/4.0,  1.0/3.0, 1.0/8.0))
);

// ---------------------------------------------------------------------------
// Fehlberg method of order 5(4)
IRTK_DEFINE_FFDRK_EMBEDDED(RKF45, 6, 5, false,
  BT_C(BT_R(   0.0,            1.0/4.0,        3.0/8.0,        12.0/13.0,     1.0,       1.0/2.0 )),
  /* ----------------------------------------------------------------------------------------------- */
  BT_A(BT_R(   0.0,            0.0,            0.0,             0.0,          0.0,       0.0     ),
       BT_R(   1.0/4.0,        0.0,            0.0,             0.0,          0.0,       0.0     ),
       BT_R(   3.0/32.0,       9.0/32.0,       0.0,             0.0,          0.0,       0.0     ),
       BT_R(1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0,      0.0,          0.0,       0.0     ),
       BT_R( 439.0/216.0,     -8.0,         3680.0/513.0,    -845.0/4104.0,   0.0,       0.0     ),
       BT_R(  -8.0/27.0,       2.0,        -3544.0/2565.0,   1859.0/4104.0, -11.0/40.0,  0.0     )),
  /* ----------------------------------------------------------------------------------------------- */
  BT_B(BT_R(  25.0/216.0,      0.0,         1408.0/2565.0,   2197.0/4104.0,  -1.0/5.0,   0.0     )),
  BT_B(BT_R(  16.0/135.0,      0.0,         6656.0/12825.0, 28561.0/56430.0, -9.0/50.0,  2.0/55.0))
);

// ---------------------------------------------------------------------------
// Dormand–Prince method of order 5(4)
IRTK_DEFINE_FFDRK_EMBEDDED(RKDP45, 7, 5, true,
  BT_C(BT_R(     0.0,             1.0/5.0,        3.0/10.0,        4.0/5.0,        8.0/9.0,        1.0,        1.0     )),
  /* --------------------------------------------------------------------------------------------------------------------- */
  BT_A(BT_R(     0.0,             0.0,            0.0,             0.0,            0.0,            0.0,        0.0     ),
       BT_R(     1.0/5.0,         0.0,            0.0,             0.0,            0.0,            0.0,        0.0     ),
       BT_R(     3.0/40.0,        9.0/40.0,       0.0,             0.0,            0.0,            0.0,        0.0     ),
       BT_R(    44.0/45.0,      -56.0/15.0,      32.0/9.0,         0.0,            0.0,            0.0,        0.0     ),
       BT_R( 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0,   -212.0/729.0,      0.0,            0.0,        0.0     ),
       BT_R(  9017.0/3168.0,   -355.0/33.0,   46732.0/5247.0,     49.0/176.0,  -5103.0/18656.0,    0.0,        0.0     ),
       BT_R(    35.0/384.0,       0.0,          500.0/1113.0,    125.0/192.0,  -2187.0/6784.0,    11.0/84.0,   0.0     )),
  /* --------------------------------------------------------------------------------------------------------------------- */
  BT_B(BT_R(    35.0/384.0,       0.0,          500.0/1113.0,    125.0/192.0,  -2187.0/6784.0,    11.0/84.0,   0.0     )),
  BT_B(BT_R(  5179.0/57600.0,     0.0,         7571.0/16695.0,   393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0))
);

// ---------------------------------------------------------------------------
// Cash-Karp method of order 5(4)
IRTK_DEFINE_FFDRK_EMBEDDED(RKCK45, 6, 5, false,
  BT_C(BT_R(   0.0,           1.0/5.0,     3.0/10.0,        3.0/5.0,        1.0,           7.0/8.0   )),
  /* --------------------------------------------------------------------------------------------------- */
  BT_A(BT_R(   0.0,           0.0,         0.0,             0.0,            0.0,           0.0       ),
       BT_R(   1.0/5.0,       0.0,         0.0,             0.0,            0.0,           0.0       ),
       BT_R(   3.0/40.0,      9.0/40.0,    0.0,             0.0,            0.0,           0.0       ),
       BT_R(   3.0/10.0,     -9.0/10.0,    6.0/5.0,         0.0,            0.0,           0.0       ),
       BT_R( -11.0/54.0,      5.0/2.0,   -70.0/27.0,       35.0/27.0,       0.0,           0.0       ),
       BT_R(1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0,    0.0       )),
  /* --------------------------------------------------------------------------------------------------- */
  BT_B(BT_R(  37.0/378.0,     0.0,       250.0/621.0,     125.0/594.0,      0.0,         512.0/1771.0)),
  BT_B(BT_R(2825.0/27648.0,   0.0,     18575.0/48384.0, 13525.0/55296.0,  277.0/14336.0,   1.0/4.0   ))
);
