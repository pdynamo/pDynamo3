# ifndef _MNDODEFINITIONS
# define _MNDODEFINITIONS

# include "Units.h"

/*==================================================================================================================================
! . Parameters.
!=================================================================================================================================*/
/* . Orbital identities - molecular frame. */
# define S     0
# define PZ    1
# define PX    2
# define PY    3
# define DZ2   4
# define DXZ   5
# define DYZ   6
# define DX2Y2 7
# define DXY   8

/* . Local frame. */
# define PSIGMA      1
# define PPIPLUS     2
# define PPIMINUS    3
# define DSIGMA      4
# define DPIPLUS     5
# define DPIMINUS    6
# define DDELTAPLUS  7
# define DDELTAMINUS 8

/* . Charge distributions. */
# define SS           0
# define PZS          1
# define PZPZ         2
# define PXS          3
# define PXPZ         4
# define PXPX         5
# define PYS          6
# define PYPZ         7
# define PYPX         8
# define PYPY         9
# define DZ2S        10
# define DZ2PZ       11
# define DZ2PX       12
# define DZ2PY       13
# define DZ2DZ2      14
# define DXZS        15
# define DXZPZ       16
# define DXZPX       17
# define DXZPY       18
# define DXZDZ2      19
# define DXZDXZ      20
# define DYZS        21
# define DYZPZ       22
# define DYZPX       23
# define DYZPY       24
# define DYZDZ2      25
# define DYZDXZ      26
# define DYZDYZ      27
# define DX2Y2S      28
# define DX2Y2PZ     29
# define DX2Y2PX     30
# define DX2Y2PY     31
# define DX2Y2DZ2    32
# define DX2Y2DXZ    33
# define DX2Y2DYZ    34
# define DX2Y2DX2Y2  35
# define DXYS        36
# define DXYPZ       37
# define DXYPX       38
# define DXYPY       39
# define DXYDZ2      40
# define DXYDXZ      41
# define DXYDYZ      42
# define DXYDX2Y2    43
# define DXYDXY      44

/*==================================================================================================================================
! . Integral definitions.
!=================================================================================================================================*/
/* . Coefficients of charge interaction terms in the local frame. */
/* . Maximum number of CH terms in an interaction. */
# define MAXCHTERMS 10

/* . Increments into the CHINDICES array - each charge distribution has 15 values: l (0, 1, 2) and m (-2, -1, 0, 1, 2). */
# define CHINCREMENT1 15
# define CHINCREMENT2  5
# define CHINCREMENT3  2

/* . Other parameters. */
/* . Constants. */
# define EXPONENT_TOLERANCE         25.00e+00 /* . This value may be important for certain molecules, e.g. H2. */
# define MNDO_AU_TO_EV              27.21e+00
# define MNDO_ZERO_CONTACT_ENERGY    1.0e+05
# define PDDG_EXPONENT            ( 10.0e+00 * Units_Length_Bohrs_To_Angstroms * Units_Length_Bohrs_To_Angstroms )
# define PM6_UNPOLARIZABLECORE    (  1.0e-08 / MNDO_AU_TO_EV ) /* . Assumed to be in eV A^12. */
# define SMALL_RIJ                   5.0e-03
# define SMALL_RIJ2                  2.5e-05

/* . The maximum number of unique one-center two-electron integrals - 1 (s), 16 (sp), 155 (spd). */
# define N1CTEIS 155

/* . The maximum number of two-center two-electron integrals - 1 (s), 100 (sp), 2025 (spd). */
# define N2CTEIS 2025

# endif
