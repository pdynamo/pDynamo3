# ifndef _MNDOINTEGRALSMM
# define _MNDOINTEGRALSMM

# include "CubicSpline.h"
# include "Integer.h"
# include "MNDOParameters.h"
# include "Real.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDOIntegralsMM_CoreCharge     ( const MNDOParameters *qData       ,
                                             const Real            qM          ,
                                             const Real            R           ,
                                                   Real           *fCore0      ,
                                                   Real           *fCore1      ,
                                                   Real           *gCore0      ,
                                                   Real           *gCore1      ) ;
extern void MNDOIntegralsMM_FromSpline     ( const MNDOParameters *qData       ,
                                             const CubicSpline    *qSpline     ,
                                             const Real            qM          ,
                                             const Real            R           ,
                                                   Real           *fCore       ,
                                                   Real           *gCore       ,
                                                   RealArray1D    *integrals   ,
                                                   RealArray1D    *gIntegrals  ) ;
extern void MNDOIntegralsMM_LocalFrame     ( const MNDOParameters *qData       ,
                                             const Real            R           ,
                                                   RealArray1D    *integrals   ,
                                                   RealArray1D    *gIntegrals  ) ;
extern void MNDOIntegralsMM_MolecularFrame ( const Integer         nQ          ,
                                             const Real            r           ,
                                             const Real            x           ,
                                             const Real            y           ,
                                             const Real            z           ,
                                             const RealArray1D    *iLocal      ,
                                             const RealArray1D    *gLocal      ,
                                                   RealArray1D    *iMolecular  ,
                                                   RealArray1D    *gMolecularX ,
                                                   RealArray1D    *gMolecularY ,
                                                   RealArray1D    *gMolecularZ ) ;
extern void MNDOIntegralsMM_Values         ( const MNDOParameters *qData       ,
                                             const Real            R           ,
                                             const Integer         index       ,
                                                   Real           *f           ,
                                                   Real           *g           ) ;
# endif
