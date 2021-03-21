# ifndef _GAUSSIANBASISCONTAINER_INTEGRALS_B0E0N2
# define _GAUSSIANBASISCONTAINER_INTEGRALS_B0E0N2

# include "Coordinates3.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real GaussianBasisContainer_NuclearNuclearEnergy     ( const RealArray1D  *chargesI      ,
                                                              const RealArray1D  *chargesJ      ,
                                                              const Coordinates3 *coordinates3I ,
                                                              const Coordinates3 *coordinates3J ,
                                                                    Selection    *selectionI    ,
                                                                    Selection    *selectionJ    ,
                                                              const RealArray1D  *widthsEI      ,
                                                              const RealArray1D  *widthsEJ      ,
                                                              const RealArray1D  *widthsNI      ,
                                                              const RealArray1D  *widthsNJ      ,
                                                                    Coordinates3 *gradients3I   ,
                                                                    Coordinates3 *gradients3J   ) ;
extern void GaussianBasisContainer_NuclearNuclearPotentials ( const RealArray1D  *chargesJ      ,
                                                              const Coordinates3 *coordinates3I ,
                                                              const Coordinates3 *coordinates3J ,
                                                                    Selection    *selectionI    ,
                                                                    Selection    *selectionJ    ,
                                                              const RealArray1D  *widthsEI      ,
                                                              const RealArray1D  *widthsEJ      ,
                                                              const RealArray1D  *widthsNI      ,
                                                              const RealArray1D  *widthsNJ      ,
                                                                    RealArray1D  *potentialsI   ) ;                                  
# endif
