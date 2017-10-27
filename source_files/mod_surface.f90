!******************************************************************************
!> Iso surface extraction from volume data by means of marching cubes
!>
!> \author Ralf Schneider
!> \date 22.01.2010
Module iso_surface

!  Use kinds !** A precision module which should at least provide
            !** an integer kind ik and a real kind rk

  implicit none
INTEGER,parameter ::rk=8,ik=8
  !----------------------------------------------------------------------------
  !> Edge lookup table
  integer, dimension(256), Parameter :: edge_table = (/&
          0,  265,  515,  778, 1030, 1295, 1541, 1804,&
       2060, 2309, 2575, 2822, 3082, 3331, 3593, 3840,&
        400,  153,  915,  666, 1430, 1183, 1941, 1692,&
       2460, 2197, 2975, 2710, 3482, 3219, 3993, 3728,&
        560,  825,   51,  314, 1590, 1855, 1077, 1340,&
       2620, 2869, 2111, 2358, 3642, 3891, 3129, 3376,&
        928,  681,  419,  170, 1958, 1711, 1445, 1196,&
       2988, 2725, 2479, 2214, 4010, 3747, 3497, 3232,&
       1120, 1385, 1635, 1898,  102,  367,  613,  876,&
       3180, 3429, 3695, 3942, 2154, 2403, 2665, 2912,&
       1520, 1273, 2035, 1786,  502,  255, 1013,  764,&
       3580, 3317, 4095, 3830, 2554, 2291, 3065, 2800,&
       1616, 1881, 1107, 1370,  598,  863,   85,  348,&
       3676, 3925, 3167, 3414, 2650, 2899, 2137, 2384,&
       1984, 1737, 1475, 1226,  966,  719,  453,  204,&
       4044, 3781, 3535, 3270, 3018, 2755, 2505, 2240,&
       2240, 2505, 2755, 3018, 3270, 3535, 3781, 4044,&
        204,  453,  719,  966, 1226, 1475, 1737, 1984,&
       2384, 2137, 2899, 2650, 3414, 3167, 3925, 3676,&
        348,   85,  863,  598, 1370, 1107, 1881, 1616,&
       2800, 3065, 2291, 2554, 3830, 4095, 3317, 3580,&
        764, 1013,  255,  502, 1786, 2035, 1273, 1520,&
       2912, 2665, 2403, 2154, 3942, 3695, 3429, 3180,&
        876,  613,  367,  102, 1898, 1635, 1385, 1120,&
       3232, 3497, 3747, 4010, 2214, 2479, 2725, 2988,&
       1196, 1445, 1711, 1958,  170,  419,  681,  928,&
       3376, 3129, 3891, 3642, 2358, 2111, 2869, 2620,&
       1340, 1077, 1855, 1590,  314,   51,  825,  560,&
       3728, 3993, 3219, 3482, 2710, 2975, 2197, 2460,&
       1692, 1941, 1183, 1430,  666,  915,  153,  400,&
       3840, 3593, 3331, 3082, 2822, 2575, 2309, 2060,&
       1804, 1541, 1295, 1030,  778,  515,  265,    0 &
       /)

  !----------------------------------------------------------------------------
  !> Triangle lookup table
  integer, dimension(16,256), Parameter :: tri_table = reshape((/ &
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  8,  3,  9,  8,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  3,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         9,  2, 10,  0,  2,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         2,  8,  3,  2, 10,  8, 10,  9,  8, -1, -1, -1, -1, -1, -1, -1,&
         3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0, 11,  2,  8, 11,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  9,  0,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1, 11,  2,  1,  9, 11,  9,  8, 11, -1, -1, -1, -1, -1, -1, -1,&
         3, 10,  1, 11, 10,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0, 10,  1,  0,  8, 10,  8, 11, 10, -1, -1, -1, -1, -1, -1, -1,&
         3,  9,  0,  3, 11,  9, 11, 10,  9, -1, -1, -1, -1, -1, -1, -1,&
         9,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  3,  0,  7,  3,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  1,  9,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  1,  9,  4,  7,  1,  7,  3,  1, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 10,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         3,  4,  7,  3,  0,  4,  1,  2, 10, -1, -1, -1, -1, -1, -1, -1,&
         9,  2, 10,  9,  0,  2,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1,&
         2, 10,  9,  2,  9,  7,  2,  7,  3,  7,  9,  4, -1, -1, -1, -1,&
         8,  4,  7,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        11,  4,  7, 11,  2,  4,  2,  0,  4, -1, -1, -1, -1, -1, -1, -1,&
         9,  0,  1,  8,  4,  7,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1,&
         4,  7, 11,  9,  4, 11,  9, 11,  2,  9,  2,  1, -1, -1, -1, -1,&
         3, 10,  1,  3, 11, 10,  7,  8,  4, -1, -1, -1, -1, -1, -1, -1,&
         1, 11, 10,  1,  4, 11,  1,  0,  4,  7, 11,  4, -1, -1, -1, -1,&
         4,  7,  8,  9,  0, 11,  9, 11, 10, 11,  0,  3, -1, -1, -1, -1,&
         4,  7, 11,  4, 11,  9,  9, 11, 10, -1, -1, -1, -1, -1, -1, -1,&
         9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         9,  5,  4,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  5,  4,  1,  5,  0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         8,  5,  4,  8,  3,  5,  3,  1,  5, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 10,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         3,  0,  8,  1,  2, 10,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1,&
         5,  2, 10,  5,  4,  2,  4,  0,  2, -1, -1, -1, -1, -1, -1, -1,&
         2, 10,  5,  3,  2,  5,  3,  5,  4,  3,  4,  8, -1, -1, -1, -1,&
         9,  5,  4,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0, 11,  2,  0,  8, 11,  4,  9,  5, -1, -1, -1, -1, -1, -1, -1,&
         0,  5,  4,  0,  1,  5,  2,  3, 11, -1, -1, -1, -1, -1, -1, -1,&
         2,  1,  5,  2,  5,  8,  2,  8, 11,  4,  8,  5, -1, -1, -1, -1,&
        10,  3, 11, 10,  1,  3,  9,  5,  4, -1, -1, -1, -1, -1, -1, -1,&
         4,  9,  5,  0,  8,  1,  8, 10,  1,  8, 11, 10, -1, -1, -1, -1,&
         5,  4,  0,  5,  0, 11,  5, 11, 10, 11,  0,  3, -1, -1, -1, -1,&
         5,  4,  8,  5,  8, 10, 10,  8, 11, -1, -1, -1, -1, -1, -1, -1,&
         9,  7,  8,  5,  7,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         9,  3,  0,  9,  5,  3,  5,  7,  3, -1, -1, -1, -1, -1, -1, -1,&
         0,  7,  8,  0,  1,  7,  1,  5,  7, -1, -1, -1, -1, -1, -1, -1,&
         1,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         9,  7,  8,  9,  5,  7, 10,  1,  2, -1, -1, -1, -1, -1, -1, -1,&
        10,  1,  2,  9,  5,  0,  5,  3,  0,  5,  7,  3, -1, -1, -1, -1,&
         8,  0,  2,  8,  2,  5,  8,  5,  7, 10,  5,  2, -1, -1, -1, -1,&
         2, 10,  5,  2,  5,  3,  3,  5,  7, -1, -1, -1, -1, -1, -1, -1,&
         7,  9,  5,  7,  8,  9,  3, 11,  2, -1, -1, -1, -1, -1, -1, -1,&
         9,  5,  7,  9,  7,  2,  9,  2,  0,  2,  7, 11, -1, -1, -1, -1,&
         2,  3, 11,  0,  1,  8,  1,  7,  8,  1,  5,  7, -1, -1, -1, -1,&
        11,  2,  1, 11,  1,  7,  7,  1,  5, -1, -1, -1, -1, -1, -1, -1,&
         9,  5,  8,  8,  5,  7, 10,  1,  3, 10,  3, 11, -1, -1, -1, -1,&
         5,  7,  0,  5,  0,  9,  7, 11,  0,  1,  0, 10, 11, 10,  0, -1,&
        11, 10,  0, 11,  0,  3, 10,  5,  0,  8,  0,  7,  5,  7,  0, -1,&
        11, 10,  5,  7, 11,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  3,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         9,  0,  1,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  8,  3,  1,  9,  8,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1,&
         1,  6,  5,  2,  6,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  6,  5,  1,  2,  6,  3,  0,  8, -1, -1, -1, -1, -1, -1, -1,&
         9,  6,  5,  9,  0,  6,  0,  2,  6, -1, -1, -1, -1, -1, -1, -1,&
         5,  9,  8,  5,  8,  2,  5,  2,  6,  3,  2,  8, -1, -1, -1, -1,&
         2,  3, 11, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        11,  0,  8, 11,  2,  0, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1,&
         0,  1,  9,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1, -1, -1, -1,&
         5, 10,  6,  1,  9,  2,  9, 11,  2,  9,  8, 11, -1, -1, -1, -1,&
         6,  3, 11,  6,  5,  3,  5,  1,  3, -1, -1, -1, -1, -1, -1, -1,&
         0,  8, 11,  0, 11,  5,  0,  5,  1,  5, 11,  6, -1, -1, -1, -1,&
         3, 11,  6,  0,  3,  6,  0,  6,  5,  0,  5,  9, -1, -1, -1, -1,&
         6,  5,  9,  6,  9, 11, 11,  9,  8, -1, -1, -1, -1, -1, -1, -1,&
         5, 10,  6,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  3,  0,  4,  7,  3,  6,  5, 10, -1, -1, -1, -1, -1, -1, -1,&
         1,  9,  0,  5, 10,  6,  8,  4,  7, -1, -1, -1, -1, -1, -1, -1,&
        10,  6,  5,  1,  9,  7,  1,  7,  3,  7,  9,  4, -1, -1, -1, -1,&
         6,  1,  2,  6,  5,  1,  4,  7,  8, -1, -1, -1, -1, -1, -1, -1,&
         1,  2,  5,  5,  2,  6,  3,  0,  4,  3,  4,  7, -1, -1, -1, -1,&
         8,  4,  7,  9,  0,  5,  0,  6,  5,  0,  2,  6, -1, -1, -1, -1,&
         7,  3,  9,  7,  9,  4,  3,  2,  9,  5,  9,  6,  2,  6,  9, -1,&
         3, 11,  2,  7,  8,  4, 10,  6,  5, -1, -1, -1, -1, -1, -1, -1,&
         5, 10,  6,  4,  7,  2,  4,  2,  0,  2,  7, 11, -1, -1, -1, -1,&
         0,  1,  9,  4,  7,  8,  2,  3, 11,  5, 10,  6, -1, -1, -1, -1,&
         9,  2,  1,  9, 11,  2,  9,  4, 11,  7, 11,  4,  5, 10,  6, -1,&
         8,  4,  7,  3, 11,  5,  3,  5,  1,  5, 11,  6, -1, -1, -1, -1,&
         5,  1, 11,  5, 11,  6,  1,  0, 11,  7, 11,  4,  0,  4, 11, -1,&
         0,  5,  9,  0,  6,  5,  0,  3,  6, 11,  6,  3,  8,  4,  7, -1,&
         6,  5,  9,  6,  9, 11,  4,  7,  9,  7, 11,  9, -1, -1, -1, -1,&
        10,  4,  9,  6,  4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4, 10,  6,  4,  9, 10,  0,  8,  3, -1, -1, -1, -1, -1, -1, -1,&
        10,  0,  1, 10,  6,  0,  6,  4,  0, -1, -1, -1, -1, -1, -1, -1,&
         8,  3,  1,  8,  1,  6,  8,  6,  4,  6,  1, 10, -1, -1, -1, -1,&
         1,  4,  9,  1,  2,  4,  2,  6,  4, -1, -1, -1, -1, -1, -1, -1,&
         3,  0,  8,  1,  2,  9,  2,  4,  9,  2,  6,  4, -1, -1, -1, -1,&
         0,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         8,  3,  2,  8,  2,  4,  4,  2,  6, -1, -1, -1, -1, -1, -1, -1,&
        10,  4,  9, 10,  6,  4, 11,  2,  3, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  2,  2,  8, 11,  4,  9, 10,  4, 10,  6, -1, -1, -1, -1,&
         3, 11,  2,  0,  1,  6,  0,  6,  4,  6,  1, 10, -1, -1, -1, -1,&
         6,  4,  1,  6,  1, 10,  4,  8,  1,  2,  1, 11,  8, 11,  1, -1,&
         9,  6,  4,  9,  3,  6,  9,  1,  3, 11,  6,  3, -1, -1, -1, -1,&
         8, 11,  1,  8,  1,  0, 11,  6,  1,  9,  1,  4,  6,  4,  1, -1,&
         3, 11,  6,  3,  6,  0,  0,  6,  4, -1, -1, -1, -1, -1, -1, -1,&
         6,  4,  8, 11,  6,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         7, 10,  6,  7,  8, 10,  8,  9, 10, -1, -1, -1, -1, -1, -1, -1,&
         0,  7,  3,  0, 10,  7,  0,  9, 10,  6,  7, 10, -1, -1, -1, -1,&
        10,  6,  7,  1, 10,  7,  1,  7,  8,  1,  8,  0, -1, -1, -1, -1,&
        10,  6,  7, 10,  7,  1,  1,  7,  3, -1, -1, -1, -1, -1, -1, -1,&
         1,  2,  6,  1,  6,  8,  1,  8,  9,  8,  6,  7, -1, -1, -1, -1,&
         2,  6,  9,  2,  9,  1,  6,  7,  9,  0,  9,  3,  7,  3,  9, -1,&
         7,  8,  0,  7,  0,  6,  6,  0,  2, -1, -1, -1, -1, -1, -1, -1,&
         7,  3,  2,  6,  7,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         2,  3, 11, 10,  6,  8, 10,  8,  9,  8,  6,  7, -1, -1, -1, -1,&
         2,  0,  7,  2,  7, 11,  0,  9,  7,  6,  7, 10,  9, 10,  7, -1,&
         1,  8,  0,  1,  7,  8,  1, 10,  7,  6,  7, 10,  2,  3, 11, -1,&
        11,  2,  1, 11,  1,  7, 10,  6,  1,  6,  7,  1, -1, -1, -1, -1,&
         8,  9,  6,  8,  6,  7,  9,  1,  6, 11,  6,  3,  1,  3,  6, -1,&
         0,  9,  1, 11,  6,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         7,  8,  0,  7,  0,  6,  3, 11,  0, 11,  6,  0, -1, -1, -1, -1,&
         7, 11,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         3,  0,  8, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  1,  9, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         8,  1,  9,  8,  3,  1, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1,&
        10,  1,  2,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 10,  3,  0,  8,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1,&
         2,  9,  0,  2, 10,  9,  6, 11,  7, -1, -1, -1, -1, -1, -1, -1,&
         6, 11,  7,  2, 10,  3, 10,  8,  3, 10,  9,  8, -1, -1, -1, -1,&
         7,  2,  3,  6,  2,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         7,  0,  8,  7,  6,  0,  6,  2,  0, -1, -1, -1, -1, -1, -1, -1,&
         2,  7,  6,  2,  3,  7,  0,  1,  9, -1, -1, -1, -1, -1, -1, -1,&
         1,  6,  2,  1,  8,  6,  1,  9,  8,  8,  7,  6, -1, -1, -1, -1,&
        10,  7,  6, 10,  1,  7,  1,  3,  7, -1, -1, -1, -1, -1, -1, -1,&
        10,  7,  6,  1,  7, 10,  1,  8,  7,  1,  0,  8, -1, -1, -1, -1,&
         0,  3,  7,  0,  7, 10,  0, 10,  9,  6, 10,  7, -1, -1, -1, -1,&
         7,  6, 10,  7, 10,  8,  8, 10,  9, -1, -1, -1, -1, -1, -1, -1,&
         6,  8,  4, 11,  8,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         3,  6, 11,  3,  0,  6,  0,  4,  6, -1, -1, -1, -1, -1, -1, -1,&
         8,  6, 11,  8,  4,  6,  9,  0,  1, -1, -1, -1, -1, -1, -1, -1,&
         9,  4,  6,  9,  6,  3,  9,  3,  1, 11,  3,  6, -1, -1, -1, -1,&
         6,  8,  4,  6, 11,  8,  2, 10,  1, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 10,  3,  0, 11,  0,  6, 11,  0,  4,  6, -1, -1, -1, -1,&
         4, 11,  8,  4,  6, 11,  0,  2,  9,  2, 10,  9, -1, -1, -1, -1,&
        10,  9,  3, 10,  3,  2,  9,  4,  3, 11,  3,  6,  4,  6,  3, -1,&
         8,  2,  3,  8,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1,&
         0,  4,  2,  4,  6,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  9,  0,  2,  3,  4,  2,  4,  6,  4,  3,  8, -1, -1, -1, -1,&
         1,  9,  4,  1,  4,  2,  2,  4,  6, -1, -1, -1, -1, -1, -1, -1,&
         8,  1,  3,  8,  6,  1,  8,  4,  6,  6, 10,  1, -1, -1, -1, -1,&
        10,  1,  0, 10,  0,  6,  6,  0,  4, -1, -1, -1, -1, -1, -1, -1,&
         4,  6,  3,  4,  3,  8,  6, 10,  3,  0,  3,  9, 10,  9,  3, -1,&
        10,  9,  4,  6, 10,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  9,  5,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  3,  4,  9,  5, 11,  7,  6, -1, -1, -1, -1, -1, -1, -1,&
         5,  0,  1,  5,  4,  0,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1,&
        11,  7,  6,  8,  3,  4,  3,  5,  4,  3,  1,  5, -1, -1, -1, -1,&
         9,  5,  4, 10,  1,  2,  7,  6, 11, -1, -1, -1, -1, -1, -1, -1,&
         6, 11,  7,  1,  2, 10,  0,  8,  3,  4,  9,  5, -1, -1, -1, -1,&
         7,  6, 11,  5,  4, 10,  4,  2, 10,  4,  0,  2, -1, -1, -1, -1,&
         3,  4,  8,  3,  5,  4,  3,  2,  5, 10,  5,  2, 11,  7,  6, -1,&
         7,  2,  3,  7,  6,  2,  5,  4,  9, -1, -1, -1, -1, -1, -1, -1,&
         9,  5,  4,  0,  8,  6,  0,  6,  2,  6,  8,  7, -1, -1, -1, -1,&
         3,  6,  2,  3,  7,  6,  1,  5,  0,  5,  4,  0, -1, -1, -1, -1,&
         6,  2,  8,  6,  8,  7,  2,  1,  8,  4,  8,  5,  1,  5,  8, -1,&
         9,  5,  4, 10,  1,  6,  1,  7,  6,  1,  3,  7, -1, -1, -1, -1,&
         1,  6, 10,  1,  7,  6,  1,  0,  7,  8,  7,  0,  9,  5,  4, -1,&
         4,  0, 10,  4, 10,  5,  0,  3, 10,  6, 10,  7,  3,  7, 10, -1,&
         7,  6, 10,  7, 10,  8,  5,  4, 10,  4,  8, 10, -1, -1, -1, -1,&
         6,  9,  5,  6, 11,  9, 11,  8,  9, -1, -1, -1, -1, -1, -1, -1,&
         3,  6, 11,  0,  6,  3,  0,  5,  6,  0,  9,  5, -1, -1, -1, -1,&
         0, 11,  8,  0,  5, 11,  0,  1,  5,  5,  6, 11, -1, -1, -1, -1,&
         6, 11,  3,  6,  3,  5,  5,  3,  1, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 10,  9,  5, 11,  9, 11,  8, 11,  5,  6, -1, -1, -1, -1,&
         0, 11,  3,  0,  6, 11,  0,  9,  6,  5,  6,  9,  1,  2, 10, -1,&
        11,  8,  5, 11,  5,  6,  8,  0,  5, 10,  5,  2,  0,  2,  5, -1,&
         6, 11,  3,  6,  3,  5,  2, 10,  3, 10,  5,  3, -1, -1, -1, -1,&
         5,  8,  9,  5,  2,  8,  5,  6,  2,  3,  8,  2, -1, -1, -1, -1,&
         9,  5,  6,  9,  6,  0,  0,  6,  2, -1, -1, -1, -1, -1, -1, -1,&
         1,  5,  8,  1,  8,  0,  5,  6,  8,  3,  8,  2,  6,  2,  8, -1,&
         1,  5,  6,  2,  1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  3,  6,  1,  6, 10,  3,  8,  6,  5,  6,  9,  8,  9,  6, -1,&
        10,  1,  0, 10,  0,  6,  9,  5,  0,  5,  6,  0, -1, -1, -1, -1,&
         0,  3,  8,  5,  6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        10,  5,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        11,  5, 10,  7,  5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        11,  5, 10, 11,  7,  5,  8,  3,  0, -1, -1, -1, -1, -1, -1, -1,&
         5, 11,  7,  5, 10, 11,  1,  9,  0, -1, -1, -1, -1, -1, -1, -1,&
        10,  7,  5, 10, 11,  7,  9,  8,  1,  8,  3,  1, -1, -1, -1, -1,&
        11,  1,  2, 11,  7,  1,  7,  5,  1, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  3,  1,  2,  7,  1,  7,  5,  7,  2, 11, -1, -1, -1, -1,&
         9,  7,  5,  9,  2,  7,  9,  0,  2,  2, 11,  7, -1, -1, -1, -1,&
         7,  5,  2,  7,  2, 11,  5,  9,  2,  3,  2,  8,  9,  8,  2, -1,&
         2,  5, 10,  2,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1,&
         8,  2,  0,  8,  5,  2,  8,  7,  5, 10,  2,  5, -1, -1, -1, -1,&
         9,  0,  1,  5, 10,  3,  5,  3,  7,  3, 10,  2, -1, -1, -1, -1,&
         9,  8,  2,  9,  2,  1,  8,  7,  2, 10,  2,  5,  7,  5,  2, -1,&
         1,  3,  5,  3,  7,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  7,  0,  7,  1,  1,  7,  5, -1, -1, -1, -1, -1, -1, -1,&
         9,  0,  3,  9,  3,  5,  5,  3,  7, -1, -1, -1, -1, -1, -1, -1,&
         9,  8,  7,  5,  9,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         5,  8,  4,  5, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1,&
         5,  0,  4,  5, 11,  0,  5, 10, 11, 11,  3,  0, -1, -1, -1, -1,&
         0,  1,  9,  8,  4, 10,  8, 10, 11, 10,  4,  5, -1, -1, -1, -1,&
        10, 11,  4, 10,  4,  5, 11,  3,  4,  9,  4,  1,  3,  1,  4, -1,&
         2,  5,  1,  2,  8,  5,  2, 11,  8,  4,  5,  8, -1, -1, -1, -1,&
         0,  4, 11,  0, 11,  3,  4,  5, 11,  2, 11,  1,  5,  1, 11, -1,&
         0,  2,  5,  0,  5,  9,  2, 11,  5,  4,  5,  8, 11,  8,  5, -1,&
         9,  4,  5,  2, 11,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         2,  5, 10,  3,  5,  2,  3,  4,  5,  3,  8,  4, -1, -1, -1, -1,&
         5, 10,  2,  5,  2,  4,  4,  2,  0, -1, -1, -1, -1, -1, -1, -1,&
         3, 10,  2,  3,  5, 10,  3,  8,  5,  4,  5,  8,  0,  1,  9, -1,&
         5, 10,  2,  5,  2,  4,  1,  9,  2,  9,  4,  2, -1, -1, -1, -1,&
         8,  4,  5,  8,  5,  3,  3,  5,  1, -1, -1, -1, -1, -1, -1, -1,&
         0,  4,  5,  1,  0,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         8,  4,  5,  8,  5,  3,  9,  0,  5,  0,  3,  5, -1, -1, -1, -1,&
         9,  4,  5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4, 11,  7,  4,  9, 11,  9, 10, 11, -1, -1, -1, -1, -1, -1, -1,&
         0,  8,  3,  4,  9,  7,  9, 11,  7,  9, 10, 11, -1, -1, -1, -1,&
         1, 10, 11,  1, 11,  4,  1,  4,  0,  7,  4, 11, -1, -1, -1, -1,&
         3,  1,  4,  3,  4,  8,  1, 10,  4,  7,  4, 11, 10, 11,  4, -1,&
         4, 11,  7,  9, 11,  4,  9,  2, 11,  9,  1,  2, -1, -1, -1, -1,&
         9,  7,  4,  9, 11,  7,  9,  1, 11,  2, 11,  1,  0,  8,  3, -1,&
        11,  7,  4, 11,  4,  2,  2,  4,  0, -1, -1, -1, -1, -1, -1, -1,&
        11,  7,  4, 11,  4,  2,  8,  3,  4,  3,  2,  4, -1, -1, -1, -1,&
         2,  9, 10,  2,  7,  9,  2,  3,  7,  7,  4,  9, -1, -1, -1, -1,&
         9, 10,  7,  9,  7,  4, 10,  2,  7,  8,  7,  0,  2,  0,  7, -1,&
         3,  7, 10,  3, 10,  2,  7,  4, 10,  1, 10,  0,  4,  0, 10, -1,&
         1, 10,  2,  8,  7,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  9,  1,  4,  1,  7,  7,  1,  3, -1, -1, -1, -1, -1, -1, -1,&
         4,  9,  1,  4,  1,  7,  0,  8,  1,  8,  7,  1, -1, -1, -1, -1,&
         4,  0,  3,  7,  4,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         4,  8,  7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         9, 10,  8, 10, 11,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         3,  0,  9,  3,  9, 11, 11,  9, 10, -1, -1, -1, -1, -1, -1, -1,&
         0,  1, 10,  0, 10,  8,  8, 10, 11, -1, -1, -1, -1, -1, -1, -1,&
         3,  1, 10, 11,  3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  2, 11,  1, 11,  9,  9, 11,  8, -1, -1, -1, -1, -1, -1, -1,&
         3,  0,  9,  3,  9, 11,  1,  2,  9,  2, 11,  9, -1, -1, -1, -1,&
         0,  2, 11,  8,  0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         3,  2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         2,  3,  8,  2,  8, 10, 10,  8,  9, -1, -1, -1, -1, -1, -1, -1,&
         9, 10,  2,  0,  9,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         2,  3,  8,  2,  8, 10,  0,  1,  8,  1, 10,  8, -1, -1, -1, -1,&
         1, 10,  2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         1,  3,  8,  9,  1,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  9,  1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
         0,  3,  8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,&
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/),&
        (/16,256/))

  !----------------------------------------------------------------------------
  !> Linked list types for dynamic data growth
  TYPE tField

     INTEGER(Kind=ik)          :: nx    ! No. of pixels in x-direction
     INTEGER(Kind=ik)          :: ny    ! No. of pixels in y-direction
     INTEGER(Kind=ik)          :: nz    ! No. of pixels in z-direction

     REAL(Kind=rk)             :: dx    ! Pixels spacing in x-direction
     REAL(Kind=rk)             :: dy    ! Pixels spacing in y-direction
     REAL(Kind=rk)             :: dz    ! Pixels spacing in z-direction

  End type tField

  type tEdge

     integer(kind=rk), dimension(:,:), allocatable :: edge_ls
     type(tedge)     , pointer                     :: edge_next => null()

  End type tEdge

  type tSort

     real(kind=rk), dimension(:,:), allocatable :: node_ls_sort
     type(tsort)  , pointer                     :: sort_next => null()

  End type tSort

  type tRef

     integer(kind=ik), dimension(:), allocatable :: node_ls_ref
     type(tref)      , pointer                   :: ref_next => null()

  End type tRef

  !----------------------------------------------------------------------------
  !> Linked list chunk size
  Integer, Parameter :: fsize = 300000

!******************************************************************************
Contains

  !----------------------------------------------------------------------------
  !> Subroutine which extracts an iso-surface from volume data by means of the
  !> marching cubes algorithm
  Subroutine MARCHING_CUBES(p_phi, phi, rho, node_ls, cell_ls)

    Type(tField)                                , Intent(In)  :: p_phi
    Real(kind=rk), Dimension(:,:,:)             , Intent(In)  :: phi
    Integer(kind=4)                            , Intent(In)  :: rho

    Real(kind=rk)   , Dimension(:,:)  , Allocatable, Intent(Out) :: node_ls
    Integer(kind=ik), Dimension(:,:)  , Allocatable, Intent(Out) :: cell_ls
    Logical         , Dimension(:)    , Allocatable              :: nalloc

    Integer(kind=ik), Dimension(:), Allocatable :: edge_ls_ref

    type(tedge), pointer :: ll_edge, edge_start
    type(tsort), pointer :: ll_sort, sort_start
    type(tref) , pointer :: ll_ref , ref_start

    Real(kind=rk)    , Dimension(3,12)  :: intpol_val
    Integer(kind=ik) , Dimension(15)    :: nodes

    Integer(kind=ik)     :: i, j, k, m, r, s, t, u, v
    Integer(kind=ik)     :: vertex_code, edge_code, edge_num, tri_vertex
    Integer(kind=ik)     :: edge_cnt, ll_edge_cnt, node_cnt, ll_sort_cnt
    Integer(kind=ik)     :: ll_cnt, ll_ref_cnt
    Integer(kind=ik)     :: mval, msize
    Real(kind=rk)        :: lin_fac, ii, jj, kk
    !==========================================================================

    i = 0
    j = 0
    k = 0
    m = 0
    r = 0
    s = 1
    t = 1
    u = 0
    v = 0
    mval = 0
    msize = 0
    ii = 0_rk
    jj = 0_rk
    kk = 0_rk
    node_cnt = 0
    edge_num = 1
    edge_cnt = 1
    ll_edge_cnt = 1
    ll_sort_cnt = 1
    ll_ref_cnt = 0
    ll_cnt = 1

    Allocate(node_ls(3,3*p_phi%nx*p_phi%ny*p_phi%nz))
    node_ls = huge(0.0_rk)
    Allocate(nalloc(3*p_phi%nx*p_phi%ny*p_phi%nz))
    nalloc = .FALSE.

    Allocate(ll_edge)
    Allocate(ll_edge%edge_ls(2,fsize))
    ll_edge%edge_ls = 0_ik

    Allocate(ll_sort)
    Allocate(ll_sort%node_ls_sort(3,fsize))
    ll_sort%node_ls_sort = huge(0.0_rk)

    Allocate(ll_ref)
    Allocate(ll_ref%node_ls_ref(fsize))
    ll_ref%node_ls_ref = 0

    call init_ll_edge(ll_edge, edge_start)
    call init_ll_sort(ll_sort, sort_start)
    call init_ll_ref (ll_ref , ref_start)

    Do k = 1,p_phi%nz - 1
       Do j = 1,p_phi%ny - 1
          Do i = 1,p_phi%nx - 1

             !*****************************************************************
             !** building vertex code *****************************************
             !** representing vertexes above and below choosen boundary value *
             !*****************************************************************

             vertex_code = 0

             If (phi(i,j,k) >= rho) Then
                vertex_code = Ibset(vertex_code,0)
             End If

             If (phi(i+1,j,k) >= rho) Then
                vertex_code = Ibset(vertex_code,1)
             End If

             If (phi(i+1,j+1,k) >= rho) Then
                vertex_code = Ibset(vertex_code,2)
             End If

             If (phi(i,j+1,k) >= rho) Then
                vertex_code = Ibset(vertex_code,3)
             End If

             If (phi(i,j,k+1) >= rho) Then
                vertex_code = Ibset(vertex_code,4)
             End If

             If (phi(i+1,j,k+1) >= rho) Then
                vertex_code = Ibset(vertex_code,5)
             End If

             If (phi(i+1,j+1,k+1) >= rho) Then
                vertex_code = Ibset(vertex_code,6)
             End If

             If (phi(i,j+1,k+1) >= rho) Then
                vertex_code = Ibset(vertex_code,7)
             End If

             !** looking up intersected edges by isosurface in edge table *****
             edge_code = edge_table(vertex_code + 1)

             !*****************************************************************
             !** interpolating between two vertexes in local coordinates ******
             !** transforming into global coordinates *************************
             !*****************************************************************

             intpol_val = 0_rk

             lin_fac    = 0_rk
             If(Btest(edge_code, 0)) Then
                lin_fac = -(rho - phi(i,j,k))/(phi(i,j,k) - phi(i+1,j,k))
             End If
             intpol_val(1,1) = Real(i,rk) + lin_fac
             intpol_val(2,1) = Real(j,rk)
             intpol_val(3,1) = Real(k,rk)

             If(Btest(edge_code, 1)) Then
                lin_fac = -(rho - phi(i+1,j,k))/(phi(i+1,j,k) - phi(i+1,j+1,k))
             End If
             intpol_val(1,2) = Real(i+1,rk)
             intpol_val(2,2) = Real(j  ,rk) + lin_fac
             intpol_val(3,2) = Real(k  ,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 2)) Then
                lin_fac = -(rho - phi(i+1,j+1,k))/(phi(i+1,j+1,k) - phi(i,j+1,k))
             End If
             intpol_val(1,3) = Real(i+1,rk) - lin_fac
             intpol_val(2,3) = Real(j+1,rk)
             intpol_val(3,3) = Real(k  ,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 3)) Then
                lin_fac = -(rho - phi(i,j+1,k))/(phi(i,j+1,k) - phi(i,j,k))
             End If
             intpol_val(1,4) = Real(i  ,rk)
             intpol_val(2,4) = Real(j+1,rk) - lin_fac
             intpol_val(3,4) = Real(k  ,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 4)) Then
                lin_fac = -(rho - phi(i,j,k+1))/(phi(i,j,k+1) - phi(i+1,j,k+1))
             End If
             intpol_val(1,5) = Real(i  ,rk) + lin_fac
             intpol_val(2,5) = Real(j  ,rk)
             intpol_val(3,5) = Real(k+1,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 5)) Then
                lin_fac = -(rho - phi(i+1,j,k+1))/(phi(i+1,j,k+1) - phi(i+1,j+1,k+1))
             End If
             intpol_val(1,6) = Real(i+1,rk)
             intpol_val(2,6) = Real(j  ,rk) + lin_fac
             intpol_val(3,6) = Real(k+1,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 6)) Then
                lin_fac = -(rho - phi(i+1,j+1,k+1))/(phi(i+1,j+1,k+1) - phi(i,j+1,k+1))
             End If
             intpol_val(1,7) = Real(i+1,rk) - lin_fac
             intpol_val(2,7) = Real(j+1,rk)
             intpol_val(3,7) = Real(k+1,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 7)) Then
                lin_fac = -(rho - phi(i,j+1,k+1))/(phi(i,j+1,k+1) - phi(i,j,k+1))
             End If
             intpol_val(1,8) =  Real(i,rk)
             intpol_val(2,8) = Real(j+1,rk) - lin_fac
             intpol_val(3,8) = Real(k+1,rk)

             lin_fac = 0_rk
             If(Btest(edge_code, 8)) Then
                lin_fac = -(rho - phi(i,j,k))/(phi(i,j,k) - phi(i,j,k+1))
             End If
             intpol_val(1,9) = Real(i  ,rk)
             intpol_val(2,9) = Real(j  ,rk)
             intpol_val(3,9) = Real(k,rk) + lin_fac

             lin_fac = 0_rk
             If(Btest(edge_code, 9)) Then
                lin_fac = -(rho - phi(i+1,j,k))/(phi(i+1,j,k) - phi(i+1,j,k+1))
             End If
             intpol_val(1,10) = Real(i+1,rk)
             intpol_val(2,10) = Real(j  ,rk)
             intpol_val(3,10) = Real(k,rk) + lin_fac

             lin_fac = 0_rk
             If(Btest(edge_code, 10)) Then
                lin_fac = -(rho - phi(i+1,j+1,k))/(phi(i+1,j+1,k) - phi(i+1,j+1,k+1))
             End If
             intpol_val(1,11) = Real(i+1,rk)
             intpol_val(2,11) = Real(j+1,rk)
             intpol_val(3,11) = Real(k,rk) + lin_fac

             lin_fac = 0_rk
             If(Btest(edge_code, 11)) Then
                lin_fac = -(rho - phi(i,j+1,k))/(phi(i,j+1,k) - phi(i,j+1,k+1))
             End If
             intpol_val(1,12) = Real(i  ,rk)
             intpol_val(2,12) = Real(j+1,rk)
             intpol_val(3,12) = Real(k,rk) + lin_fac

             !*****************************************************************
             !** looking up vertexes of built triangles in the tri table     **
             !** linking globol node coordinates with local node coordinates **
             !*****************************************************************

             Call get_tri_edge(intpol_val, vertex_code, node_cnt, &
                  p_phi%nx, p_phi%ny, p_phi%nz, &
                  node_ls, i, j, k, nodes, tri_vertex , nalloc)

             !** writing global edge list **

             do m = 1,tri_vertex-3,3

                ll_edge%edge_ls(1,edge_num)   = nodes(m)
                ll_edge%edge_ls(2,edge_num)   = nodes(m+1)
                ll_edge%edge_ls(1,edge_num+1) = nodes(m+1)
                ll_edge%edge_ls(2,edge_num+1) = nodes(m+2)
                ll_edge%edge_ls(1,edge_num+2) = nodes(m+2)
                ll_edge%edge_ls(2,edge_num+2) = nodes(m)
                edge_num = edge_num + 3
                edge_cnt = edge_cnt + 3

                if(edge_num > fsize) then
                   call get_ll_edge(ll_edge, fsize)
                   edge_num = 1
                   ll_edge_cnt = ll_edge_cnt + 1
                end if

             end do

          End Do
       End Do
    End Do

    !**************************************************************************
    !** Pack Node and element lists *******************************************
    !**************************************************************************
    do r = 1,size(node_ls,2)

!!$       if(node_ls(1,r) /= huge(0.0_rk) .AND. node_ls(2,r) /= huge(0.0_rk)&
!!$            .AND. node_ls(3,r) /= huge(0.0_rk)) then
       if(nalloc(r)) then

          ll_sort%node_ls_sort(1,ll_sort_cnt) = node_ls(1,r)
          ll_sort%node_ls_sort(2,ll_sort_cnt) = node_ls(2,r)
          ll_sort%node_ls_sort(3,ll_sort_cnt) = node_ls(3,r)

          ll_ref%node_ls_ref(ll_sort_cnt) = r

          ll_sort_cnt = ll_sort_cnt + 1

          mval = maxval(ll_ref%node_ls_ref)

          if(ll_sort_cnt > fsize) then
             ll_cnt = ll_cnt + 1
             msize = mval

             call get_ll_sort(ll_sort,fsize)
             call get_ll_ref(ll_ref,fsize)
             ll_sort_cnt = 1

          end if
       end if
       msize = max(msize, mval)
    end do

    Allocate(edge_ls_ref(msize))
    edge_ls_ref = 0

    ll_ref => ref_start
    do s = 1, ll_cnt-1
       do t = 1, fsize
          edge_ls_ref(ll_ref%node_ls_ref(t)) = ll_ref_cnt*fsize+t
       end do
       ll_ref_cnt = ll_ref_cnt + 1
       ll_ref => ll_ref%ref_next
    end do
    do t = 1, ll_sort_cnt-1
       edge_ls_ref(ll_ref%node_ls_ref(t)) = ll_ref_cnt*fsize+t
    end do

    Deallocate(node_ls)

    !** Generate node list ****************************************************
    Allocate(node_ls(3,(ll_cnt-1)*fsize + (ll_sort_cnt-1)))

    i = 1

    ll_sort => sort_start
    Do u = 1, ll_cnt-1
       Do v = 1, fsize
          node_ls(:,i) = ll_sort%node_ls_sort(:,v)
          i = i + 1
       End Do
       ll_sort => ll_sort%sort_next
    End Do
    Do v = 1, ll_sort_cnt-1
       node_ls(:,i) = ll_sort%node_ls_sort(:,v)
       i = i + 1
    End Do

    node_ls(1,:) = node_ls(1,:) * p_phi%dx
    node_ls(2,:) = node_ls(2,:) * p_phi%dy
    node_ls(3,:) = node_ls(3,:) * p_phi%dz

    !** Generate cell list ****************************************************
    Allocate(cell_ls(3,(edge_cnt-1)/3))

    i = 1

    ll_edge => edge_start
    Do u = 1, ll_edge_cnt-1
       Do v = 0, fsize/3-1
          cell_ls(1,i) = edge_ls_ref(ll_edge%edge_ls(1,3*v+1))
          cell_ls(2,i) = edge_ls_ref(ll_edge%edge_ls(1,3*v+2))
          cell_ls(3,i) = edge_ls_ref(ll_edge%edge_ls(1,3*v+3))
          i = i + 1
       End Do
       ll_edge => ll_edge%edge_next
    End Do
    Do v = 0, (edge_num-1)/3-1
          cell_ls(1,i) = edge_ls_ref(ll_edge%edge_ls(1,3*v+1))
          cell_ls(2,i) = edge_ls_ref(ll_edge%edge_ls(1,3*v+2))
          cell_ls(3,i) = edge_ls_ref(ll_edge%edge_ls(1,3*v+3))
          i = i + 1
    End Do

  End Subroutine marching_cubes

  !****************************************************************************
  !** Transform local to global coordinates ***********************************
  !****************************************************************************
  subroutine get_tri_edge(intpol_val, vertex_code, node_cnt, x_dim, y_dim, z_dim,&
       node_ls, i, j, k, nodes, tri_vertex, nalloc)

    integer(kind=ik), dimension(15) :: nodes
    integer(kind=ik), intent(in) :: vertex_code, x_dim, y_dim, z_dim
    real(kind = rk), dimension(3,12), intent(in) :: intpol_val
    real(kind = rk), dimension(3,3*x_dim*y_dim*z_dim), intent(inout):: node_ls
    Logical        , dimension(3*x_dim*y_dim*z_dim)  , intent(inout):: nalloc

    integer(kind=ik) :: tri_vertex, node_num, i, j, k, tri_edge
    integer(kind=ik), intent(inout) :: node_cnt

    tri_vertex = 1
    tri_edge = 0
    node_num = 0
    nodes = 0


    do while (tri_table(tri_vertex, vertex_code + 1) /= -1)
       tri_edge = tri_table(tri_vertex, vertex_code + 1)

       !** transforming local node number in global node number **

       select case(tri_edge)

       case(0)
          node_num = i + (j-1) * (x_dim-1) + (k-1) * (y_dim) * (x_dim-1)
       case(1)
          node_num = (i+1) + (j-1) * x_dim + (k-1) * (y_dim-1) * x_dim + (x_dim-1) * y_dim * z_dim
       case(2)
          node_num = i + j * (x_dim-1) + (k-1) * (y_dim) * (x_dim-1)
       case(3)
          node_num = i + (j-1) * x_dim + (k-1) * (y_dim-1) * x_dim + (x_dim-1) * y_dim * z_dim
       case(4)
          node_num = i + (j-1) * (x_dim-1) + k * (y_dim) * (x_dim-1)
       case(5)
          node_num = (i+1) + (j-1) * x_dim + k * (y_dim-1) * x_dim + (x_dim-1) * y_dim * z_dim
       case(6)
          node_num = i + j * (x_dim-1) + k * (y_dim) * (x_dim-1)
       case(7)
          node_num = i + (j-1) * x_dim + k * (y_dim-1) * x_dim + (x_dim-1) * y_dim * z_dim
       case(8)
          node_num = i + (j-1) * x_dim + (k-1) * y_dim * x_dim + (x_dim-1) * y_dim * z_dim + x_dim * (y_dim-1) * z_dim
       case(9)
          node_num = (i+1) + (j-1) * x_dim + (k-1) * y_dim * x_dim + (x_dim-1) * y_dim * z_dim + x_dim * (y_dim-1) * z_dim
       case(10)
          node_num = (i+1) + j * x_dim + (k-1) * y_dim * x_dim + (x_dim-1) * y_dim * z_dim + x_dim * (y_dim-1) * z_dim
       case(11)
          node_num = i + j * x_dim + (k-1) * y_dim * x_dim + (x_dim-1) * y_dim * z_dim + x_dim * (y_dim-1) * z_dim
       case default

          write(*,*)"Tri_edge case ",tri_edge," not supported"
          stop

       end select

       !** writing global node list **

       node_ls(1,node_num) = intpol_val(1,tri_edge + 1)
       node_ls(2,node_num) = intpol_val(2,tri_edge + 1)
       node_ls(3,node_num) = intpol_val(3,tri_edge + 1)

       nalloc(node_num) = .TRUE.

       nodes(tri_vertex) = node_num

       node_cnt = node_cnt + 1
       tri_vertex = tri_vertex + 1

    end do

  end subroutine get_tri_edge

  !****************************************************************************
  !** Linked list subroutines *************************************************
  !****************************************************************************
  subroutine init_ll_edge(ll_edge, edge_start)

    type(tedge), pointer, intent(inout) :: ll_edge, edge_start

    nullify(ll_edge%edge_next)
    edge_start => ll_edge

  end subroutine init_ll_edge

  subroutine init_ll_sort(ll_sort, sort_start)

    type(tsort), pointer, intent(inout) :: ll_sort, sort_start

    nullify(ll_sort%sort_next)
    sort_start => ll_sort

  end subroutine init_ll_sort

  subroutine init_ll_ref(ll_ref, ref_start)

    type(tref), pointer, intent(inout):: ll_ref, ref_start

    nullify(ll_ref%ref_next)
    ref_start => ll_ref

  end subroutine init_ll_ref

  subroutine get_ll_edge(ll_edge, fsize)

    type(tedge), pointer, intent(inout) :: ll_edge

    integer :: fsize

    allocate(ll_edge%edge_next)
    ll_edge => ll_edge%edge_next
    allocate(ll_edge%edge_ls(2,fsize))
    nullify(ll_edge%edge_next)

  end subroutine get_ll_edge

  subroutine get_ll_sort(ll_sort, fsize)

    type(tsort), pointer, intent(inout) :: ll_sort

    integer :: fsize

    allocate(ll_sort%sort_next)
    ll_sort => ll_sort%sort_next
    allocate(ll_sort%node_ls_sort(3,fsize))
    nullify(ll_sort%sort_next)

  end subroutine get_ll_sort

  subroutine get_ll_ref(ll_ref, fsize)

    type(tref), pointer, intent(inout) :: ll_ref

    integer :: fsize

    allocate(ll_ref%ref_next)
    ll_ref => ll_ref%ref_next
    allocate(ll_ref%node_ls_ref(fsize))
    ll_ref%node_ls_ref = -1
    nullify(ll_ref%ref_next)

  end subroutine get_ll_ref

End Module iso_surface

!        p_phi%nx=ct%pixel_x
!        p_phi%ny=ct%pixel_y
!        p_phi%nz=ct%slices
!
!        p_phi%dx=ct%dx
!        p_phi%dy=ct%dy
!        p_phi%dz=ct%dz
!        rho=2600
!        CALL MARCHING_CUBES(p_phi,ct_data,1000,node_ls,cell_ls)
!        no_nodes=Int(size(node_ls(1,:)),4)
!        no_elems=Int(size(cell_ls(1,:)),4)
!        write(file_name,'(A,I0,A)')'contour_1000.vtk'
!        CALL write_vtk_unstructured_grid(node_ls,cell_ls,no_nodes,no_elems,trim(file_name))
!        !========================================
!         type(tField)                                             ::p_phi
!         INTEGER(rk8)                                              ::rho
!         Real(kind=rk8),Dimension(:,:),Allocatable                 ::node_ls
!         Integer(rk8), Dimension(:,:), Allocatable                 ::cell_ls
!         Integer(rk4)                                              ::no_nodes,no_elems
!
!==============================================================
