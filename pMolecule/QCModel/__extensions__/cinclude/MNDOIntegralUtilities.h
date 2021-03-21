# ifndef _MNDOINTEGRALUTILITIES
# define _MNDOINTEGRALUTILITIES

# include "Boolean.h"
# include "Integer.h"
# include "MNDOParameters.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type definition for function pointer. */
typedef Real ( * ChargeInteractionFunction ) ( const Real r, const Integer l1, const Integer l2, const Integer m, const Real da, const Real db, const Real add ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Displacement - requires rJI. */
# define MNDOIntegralUtilities_GetDisplacement( rI, rJ, r, x, y, z ) \
    x = rJ[0] - rI[0] ; \
    y = rJ[1] - rI[1] ; \
    z = rJ[2] - rI[2] ; \
    r = sqrt ( x * x + y * y + z * z ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDOIntegralUtilities_GetTransformationMatrices ( const Integer                   ni               ,
                                                              const Integer                   nj               ,
                                                              const Real                      r                ,
                                                              const Real                      x                ,
                                                              const Real                      y                ,
                                                              const Real                      z                ,
                                                                    RealArray2D             **iTransformation  ,
                                                                    RealArray2D             **jTransformation  ,
                                                                    RealArray2D             **iTransformationX ,
                                                                    RealArray2D             **iTransformationY ,
                                                                    RealArray2D             **iTransformationZ ,
                                                                    RealArray2D             **jTransformationX ,
                                                                    RealArray2D             **jTransformationY ,
                                                                    RealArray2D             **jTransformationZ ) ;
extern void MNDOIntegralUtilities_LocalFrame2CTEIs          ( const MNDOParameters           *iData            ,
                                                              const MNDOParameters           *jData            ,
                                                              const Real                      r                ,
                                                                    RealArray2D              *lfteis           ,
                                                                    RealArray1D              *core1b           ,
                                                                    RealArray1D              *core2a           ,
                                                                    RealArray2D              *dlfteis          ,
                                                                    RealArray1D              *dcore1b          ,
                                                                    RealArray1D              *dcore2a          ) ;
extern void MNDOIntegralUtilities_LocalFrame2COEIsSP        ( const MNDOParameters           *iData            ,
                                                              const Real                      jPO8             ,
                                                              const Real                      r                ,
                                                              const Boolean                   swapped          ,
                                                                    RealArray1D              *core             ,
                                                                    RealArray1D              *dcore            ) ;
extern void MNDOIntegralUtilities_LocalFrame2CTEIsSP        ( const MNDOParameters           *iData            ,
                                                              const MNDOParameters           *jData            ,
                                                              const Real                      r                ,
                                                                    RealArray2D              *lfteis           ,
                                                                    RealArray2D              *dlfteis          ) ;
extern Real MNDOIntegralUtilities_LocalFrame2CTEI           ( const ChargeInteractionFunction Evaluate         ,
                                                              const MNDOParameters           *iData            ,
                                                              const MNDOParameters           *jData            ,
                                                              const Integer                   ij               ,
                                                              const Integer                   kl               ,
                                                              const Integer                   i                ,
                                                              const Integer                   j                ,
                                                              const Integer                   k                ,
                                                              const Integer                   l                ,
                                                              const Integer                   c                ,
                                                              const Real                      r                ) ;

/* . Charge interaction functions. */
extern Real MNDOIntegralUtilities_2CChargeInteraction       ( const Real r, const Integer  l1, const Integer  l2, const Integer  m, const Real da, const Real db, const Real add ) ;
extern Real MNDOIntegralUtilities_2CChargeInteractionD      ( const Real r, const Integer  l1, const Integer  l2, const Integer  m, const Real da, const Real db, const Real add ) ;

# endif
