# ifndef _MNDOMULTIPOLES
# define _MNDOMULTIPOLES

# include "IntegerArray1D.h"
# include "MNDOParametersContainer.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Enumerations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Multipole types. */
typedef enum {
    MultipoleRepresentation_Buckingham = 1 , /* . Cartesian traceless. */
    MultipoleRepresentation_Cartesian  = 2 , /* . Cartesian traced. */
    MultipoleRepresentation_Spherical  = 3
} MultipoleRepresentation ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDO_AtomicMultipoles     ( const MNDOParametersContainer *parameters              ,
                                        const IntegerArray1D          *basisIndices            ,
                                        const SymmetricMatrix         *density                 ,
                                        const MultipoleRepresentation  multipoleRepresentation ,
                                        const Integer                  multipoleOrder          ,
                                              RealArray1D             *multipoles              ) ;
extern void MNDO_AtomicMultipolesFock ( const MNDOParametersContainer *parameters              ,
                                        const IntegerArray1D          *basisIndices            ,
                                        const RealArray1D             *potentials              ,
                                        const Integer                  multipoleOrder          ,
                                              SymmetricMatrix         *fock                    ) ;
extern void MNDO_BondOrders           ( const IntegerArray1D          *basisIndices            ,
                                        const SymmetricMatrix         *density                 ,
                                              SymmetricMatrix         *bondOrders              ) ;
# endif
