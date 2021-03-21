/*==================================================================================================================================
! . This module handles DFT integration.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DFTGrid.h"
# include "DFTGridWeights.h"
# include "DFTIntegrator.h"
# include "DFTIntegratorDataBlock.h"
# include "GaussianBasisContainerIntegrals_b1e0n1.h"
# include "GridFunctionDataBlock.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Memory.h"

# define _DFTGRIDWEIGHTDERIVATIVES

/* . Arrays are defined as number of basis functions (columns) * number of grid points (rows). */

/*
!
! . Formulae for rho-dependent terms:
!    Rp     = Sum_mn Bmp Bnp Pmn
!    Del Rp = Sum_mn ( Xmp Bnp + Bmp Xnp ) Pmn, etc.
!    Sp     = Del Rp . Del Rp
!    Lp     = Sum_mn ( XXmp Bnp + 2 Xmp Xnp + Bmp XXnp ) Pmn + ...
!    Tp     = 1/2 Sum_mn ( Xmp Xnp + Ymp Ynp + Zmp Znp ) Pmn
!
! . Formulae for derivatives (all v-terms weighted by Wp):
!    Rp      -     Sum_p Bmp Bnp vRp
!    Sp      - 2 * Sum_p ( Xmp Bnp + Bmp Xnp ) * dRhoXp         * vSp + ...
!    Sp (ab) -     Sum_p ( Xmp Bnp + Bmp Xnp ) * dRhoXp (other) * vSp + ...
!    Lp      -     Sum_p ( XXmp Bnp + 2 Xmp Xnp + Bmp XXnp ) vLp + ...
!    Tp      - 1/2 Sum_p ( Xmp Xnp + Ymp Ynp + Zmp Znp ) vTp
!
! . Key:
!    B                      - basis function values
!    P                      - density matrix
!    X, Y, Z                - derivative values of basis functions
!    XX, XY, XZ, YY, YZ, ZZ - second derivative values of basis functions
!
!    m, n refer to basis functions, p to points.
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_Fock                        ( const Boolean                hasSigma       ,
                                                        const Boolean                hasLaplacian   ,
                                                        const Boolean                hasTau         ,
                                                        const GridFunctionDataBlock *basisData      ,
                                                        const RealArray1D           *dRhoX          ,
                                                        const RealArray1D           *dRhoY          ,
                                                        const RealArray1D           *dRhoZ          ,
                                                        const RealArray1D           *vRho           ,
                                                        const RealArray1D           *vSigma         ,
                                                        const RealArray1D           *vLaplacianRho  ,
                                                        const RealArray1D           *vTau           ,
                                                              SymmetricMatrix       *fock           ,
                                                              RealArray2D           *work2D         ) ;
static void DFTIntegrator_FockSigma                   ( const Boolean                hasSigma       ,
                                                        const Boolean                scale          ,
                                                        const GridFunctionDataBlock *basisData      ,
                                                        const RealArray1D           *dRhoX          ,
                                                        const RealArray1D           *dRhoY          ,
                                                        const RealArray1D           *dRhoZ          ,
                                                        const RealArray1D           *vSigma         ,
                                                              SymmetricMatrix       *fock           ,
                                                              RealArray2D           *work2D         ) ;
static void DFTIntegrator_FormReducedDensity          ( const IntegerArray1D        *indices        ,
                                                        const SymmetricMatrix       *density        ,
                                                              RealArray2D          **reducedDensity ,
                                                              Status                *status         ) ;
static void DFTIntegrator_Gradients                   ( const Boolean                hasSigma       ,  
                                                        const Boolean                hasLaplacian   ,  
                                                        const Boolean                hasTau         ,  
                                                        const IntegerArray1D        *atomIndices    ,  
                                                        const GridFunctionDataBlock *basisData      ,  
                                                        const RealArray1D           *dRhoX          ,  
                                                        const RealArray1D           *dRhoY          ,  
                                                        const RealArray1D           *dRhoZ          ,  
                                                        const RealArray1D           *vRho           ,  
                                                        const RealArray1D           *vSigma         ,  
                                                        const RealArray1D           *vLaplacianRho  ,  
                                                        const RealArray1D           *vTau           ,  
                                                        const RealArray2D           *density        ,  
                                                        const Integer                gridAtom       ,  
                                                              Coordinates3          *gradients3     ,  
                                                              RealArray2D           *temp2D         ,  
                                                              RealArray2D           *work2D         ) ;
static void DFTIntegrator_GradientsSigma              ( const Boolean                hasSigma       ,  
                                                        const Boolean                isSelf         ,  
                                                        const IntegerArray1D        *atomIndices    ,  
                                                        const GridFunctionDataBlock *basisData      ,  
                                                        const RealArray1D           *dRhoX          ,  
                                                        const RealArray1D           *dRhoY          ,  
                                                        const RealArray1D           *dRhoZ          ,  
                                                        const RealArray1D           *vSigma         ,  
                                                        const RealArray2D           *density        ,  
                                                        const Integer                gridAtom       ,  
                                                              Coordinates3          *gradients3     ,  
                                                              RealArray2D           *temp2D         ,  
                                                              RealArray2D           *work2D         ) ;
static void DFTIntegrator_GridPointRho                ( const Boolean                hasSigma       ,
                                                        const Boolean                hasLaplacian   ,
                                                        const Boolean                hasTau         ,
                                                        const GridFunctionDataBlock *basisData      ,
                                                        const RealArray2D           *density        ,
                                                              RealArray1D           *rho            ,
                                                              RealArray1D           *dRhoX          ,
                                                              RealArray1D           *dRhoY          ,
                                                              RealArray1D           *dRhoZ          ,
                                                              RealArray1D           *sigma          ,
                                                              RealArray1D           *laplacianRho   ,
                                                              RealArray1D           *tau            ,
                                                              RealArray1D           *work1D         ,
                                                              RealArray2D           *work2D         ) ;
static void DFTIntegrator_GridPointSigma              ( const Boolean                hasSigma       ,  
                                                        const RealArray1D           *dRhoXa         ,  
                                                        const RealArray1D           *dRhoYa         ,  
                                                        const RealArray1D           *dRhoZa         ,  
                                                        const RealArray1D           *dRhoXb         ,  
                                                        const RealArray1D           *dRhoYb         ,  
                                                        const RealArray1D           *dRhoZb         ,  
                                                              RealArray1D           *sigma          ) ;
static void UGradientContributions                    ( const IntegerArray1D        *atomIndices    ,  
                                                        const IntegerArray1D        *indices        ,  
                                                        const RealArray2D           *a              ,  
                                                        const RealArray2D           *x              ,  
                                                        const RealArray2D           *y              ,  
                                                        const RealArray2D           *z              ,  
                                                        const Integer                gridAtom       ,  
                                                              Coordinates3          *gradients3     ) ;
static void URealArray2D_ColumnAddScaledArray         ( const RealArray2D           *a              ,
                                                        const RealArray1D           *weights        ,
                                                              RealArray2D           *b              ) ;
static void URealArray2D_ColumnDotProducts            ( const Boolean                initialize     ,
                                                        const RealArray2D           *a              ,
                                                        const RealArray2D           *b              ,
                                                              RealArray1D           *c              ) ;
static void URealArray2D_ColumnScale                  ( const RealArray1D           *a              ,
                                                              RealArray2D           *b              ) ;
static void USymmetricMatrix_DotProductIncrement      ( const IntegerArray1D        *indices        ,
                                                        const RealArray2D           *a              ,
                                                        const RealArray2D           *b              ,
                                                              SymmetricMatrix       *c              ) ;
static void USymmetricMatrix_IndexedCopyToRealArray2D ( const SymmetricMatrix       *self           ,
                                                        const IntegerArray1D        *indices        ,
                                                              RealArray2D           *target         ,
                                                              Status                *status         ) ;

/*==================================================================================================================================
! . Integrate over a grid.
!=================================================================================================================================*/
extern void DFTIntegrator_Integrate ( const DFTFunctionalModel     *functionalModel    ,
                                            DFTGrid                *grid               ,
                                      const GaussianBasisContainer *gaussianBases      ,
                                            Coordinates3           *qcCoordinates3     ,
                                      const SymmetricMatrix        *densityP           ,
                                      const SymmetricMatrix        *densityQ           ,
                                      const Boolean                 inCore             ,
                                      const Boolean                 isSpinUnrestricted ,
                                            Real                   *eQuad              ,
                                            Real                   *rhoQuad            ,
                                            SymmetricMatrix        *fockA              ,
                                            SymmetricMatrix        *fockB              ,
                                            Coordinates3           *gradients3         ,
                                            Status                 *status             )
{
    if ( eQuad   != NULL ) (*eQuad  ) = 0.0e+00 ;
    if ( rhoQuad != NULL ) (*rhoQuad) = 0.0e+00 ;
    if ( ( functionalModel != NULL ) &&
         ( grid            != NULL ) &&
         ( gaussianBases   != NULL ) &&
         ( qcCoordinates3  != NULL ) &&
         ( densityP        != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean determineFunctionData, doFock, doGradients, storeFunctionData ;        
        auto Integer gridAtom, numberOfBasisFunctions = 0, order, r ;                                                   
        auto Real    eXCTotal = 0.0e+00, rhoTotal = 0.0e+00 ;                               
        auto Status  localStatus = Status_OK ;                                              
        auto Coordinates3                  *coordinates3    = NULL ;
        auto DFTGridPointBlock             *block           = NULL ;
        auto DFTGridWeightsDerivativesWork *weightsWork     = NULL ;
        auto DFTIntegratorDataBlock        *rhoData         = NULL ;
        auto DFTIntegratorDataBlockView    *rhoDataP        = NULL , *rhoDataQ        = NULL ;
        auto GridFunctionDataBlock         *basisData       = NULL ;
        auto IntegerArray1D                *atomIndices     = NULL , *basisIndices    = NULL ;
        auto RealArray1D                   *weights         = NULL , *work1D          = NULL ;
        auto RealArray2D                   *reducedDensityP = NULL , *reducedDensityQ = NULL, *temp2D = NULL, *work2D = NULL ;
        /* . Initialization. */
        DFTGrid_MakeRecords ( grid, &localStatus ) ;
        doFock                 = ( fockA      != NULL ) ;
        doGradients            = ( gradients3 != NULL ) ;
        numberOfBasisFunctions = GaussianBasisContainer_NumberOfBasisFunctions ( gaussianBases, True ) ;
        order                  = functionalModel->order ;
        if ( doGradients )
        {
            order      += 1 ;
            atomIndices = IntegerArray1D_AllocateWithExtent ( numberOfBasisFunctions, &localStatus ) ;
            GaussianBasisContainer_MakeBasisAtomIndices ( gaussianBases, True, atomIndices, &localStatus ) ;
            DFTGrid_DeallocateFunctionData ( grid, &localStatus ) ;
        }
        if ( localStatus != Status_OK ) goto FinishUp ;
        determineFunctionData = ( ! inCore      ) || ( inCore && ( grid->records[0]->functionData == NULL ) ) ;
        storeFunctionData     = ( ! doGradients ) && ( inCore && ( grid->records[0]->functionData == NULL ) ) ;
        if ( determineFunctionData )
        {
            basisIndices = IntegerArray1D_AllocateWithExtent ( gaussianBases->capacity + 1, &localStatus ) ;
            GaussianBasisContainer_MakeBasisIndices ( gaussianBases, True, basisIndices, &localStatus ) ;
        }
        /* . Loop over the grid point blocks. */
        for ( r = 0 ; r < grid->numberOfRecords ; r++ )
        {
            /* . Get the block and associated data. */
            block        = grid->records[r]    ;
            gridAtom     = block->atom         ;
            coordinates3 = block->coordinates3 ;
            weights      = block->weights      ;
            /* . Determine basis function values and their derivatives at the grid points. */
            if ( determineFunctionData )
            {
                if ( ( basisData == NULL ) || ( basisData->numberOfPoints != block->numberOfPoints ) )
                {
                    GridFunctionDataBlock_Deallocate ( &basisData ) ;
                    basisData = GridFunctionDataBlock_Allocate ( numberOfBasisFunctions, block->numberOfPoints, order, &localStatus ) ;
                }
                else GridFunctionDataBlock_Resize ( basisData, numberOfBasisFunctions, &localStatus ) ;
                GaussianBasisContainerIntegrals_GridFunctionDataBlock ( gaussianBases        ,
                                                                        basisIndices         ,
                                                                        qcCoordinates3       ,
                                                                        coordinates3         ,
                                                                        True                 ,
                                                                        &(grid->bfTolerance) ,
                                                                        basisData            ,
                                                                        &localStatus         ) ;
            }
            /* . Retrieve function data. */
            else basisData = block->functionData ;
            if ( ( basisData == NULL ) || ( localStatus != Status_OK ) || ( basisData->numberOfFunctions <= 0 ) ) goto EndOfLoop ;
            /* . Ensure that there is an integration data block of the correct size. */
            if ( ( rhoData == NULL ) || ( rhoData->numberOfPoints != block->numberOfPoints ) )
            {
                DFTIntegratorDataBlock_Deallocate ( &rhoData ) ;
                rhoData = DFTIntegratorDataBlock_Allocate ( functionalModel->numberOfFunctionals ,
                                                            block->numberOfPoints                ,
                                                            functionalModel->hasSigma            ,
                                                            functionalModel->hasLaplacian        ,
                                                            functionalModel->hasTau              ,
                                                            functionalModel->isSpinRestricted    ,
                                                            &localStatus                         ) ;
                if ( rhoData == NULL ) goto EndOfLoop ;
                rhoDataP = &(rhoData->viewP) ;
                rhoDataQ = &(rhoData->viewQ) ;
            }
            /* . Allocate scratch space. */
            if ( (                  work2D   == NULL                         ) ||
                 ( View2D_Rows    ( work2D ) != basisData->numberOfFunctions ) ||
                 ( View2D_Columns ( work2D ) != rhoData->numberOfPoints      ) )
            {
                RealArray2D_Deallocate ( &work2D ) ;
                work2D = RealArray2D_AllocateWithExtents ( basisData->numberOfFunctions, rhoData->numberOfPoints, &localStatus ) ;
                if ( work2D == NULL ) goto EndOfLoop ;
            }
            if ( ( functionalModel->hasLaplacian || functionalModel->hasTau ) &&
                 ( ( work1D == NULL ) || ( View1D_Extent ( work1D ) != rhoData->numberOfPoints ) ) )
            {
                RealArray1D_Deallocate ( &work1D ) ;
                work1D = RealArray1D_AllocateWithExtent ( rhoData->numberOfPoints, &localStatus ) ;
                if ( work1D == NULL ) goto EndOfLoop ;
            }
            if ( doGradients && functionalModel->hasSigma &&
                 ( (                  temp2D   == NULL                         ) ||
                   ( View2D_Rows    ( temp2D ) != basisData->numberOfFunctions ) ||
                   ( View2D_Columns ( temp2D ) != rhoData->numberOfPoints      ) ) )
            {
                RealArray2D_Deallocate ( &temp2D ) ;
                temp2D = RealArray2D_AllocateWithExtents ( basisData->numberOfFunctions, rhoData->numberOfPoints, &localStatus ) ;
                if ( temp2D == NULL ) goto EndOfLoop ;
            }
            /* . Evaluate the densities and associated quantities at the grid points. */
            DFTIntegrator_FormReducedDensity ( basisData->indices ,
                                               densityP           ,
                                               &reducedDensityP   ,
                                               &localStatus       ) ;
            if ( reducedDensityP == NULL ) goto EndOfLoop ;
            DFTIntegrator_GridPointRho   ( functionalModel->hasSigma     ,
                                           functionalModel->hasLaplacian ,
                                           functionalModel->hasTau       ,
                                           basisData                     ,
                                           reducedDensityP               ,  
                                           &(rhoDataP->rho         )     ,  
                                           &(rhoDataP->dRhoX       )     ,
                                           &(rhoDataP->dRhoY       )     ,
                                           &(rhoDataP->dRhoZ       )     ,
                                           &(rhoDataP->sigma       )     ,
                                           &(rhoDataP->laplacianRho)     ,
                                           &(rhoDataP->tau         )     ,
                                           work1D                        ,  
                                           work2D                        ) ;
            if ( isSpinUnrestricted )
            {
                DFTIntegrator_FormReducedDensity ( basisData->indices ,
                                                   densityQ           ,
                                                   &reducedDensityQ   ,
                                                   &localStatus       ) ;
                if ( reducedDensityQ == NULL ) goto EndOfLoop ;
                DFTIntegrator_GridPointRho   ( functionalModel->hasSigma     ,
                                               functionalModel->hasLaplacian ,
                                               functionalModel->hasTau       ,
                                               basisData                     ,
                                               reducedDensityQ               ,  
                                               &(rhoDataQ->rho         )     ,  
                                               &(rhoDataQ->dRhoX       )     ,
                                               &(rhoDataQ->dRhoY       )     ,
                                               &(rhoDataQ->dRhoZ       )     ,
                                               &(rhoDataQ->sigma       )     ,
                                               &(rhoDataQ->laplacianRho)     ,
                                               &(rhoDataQ->tau         )     ,
                                               work1D                        ,  
                                               work2D                        ) ;
                DFTIntegrator_GridPointSigma ( functionalModel->hasSigma     ,
                                               &(rhoDataP->dRhoX )           ,
                                               &(rhoDataP->dRhoY )           ,
                                               &(rhoDataP->dRhoZ )           ,
                                               &(rhoDataQ->dRhoX )           ,
                                               &(rhoDataQ->dRhoY )           ,
                                               &(rhoDataQ->dRhoZ )           ,
                                               &(rhoData->sigmaPQ)           ) ;
            }
            /* . Skip the block if all densities are insignificant. */
            if ( RealArray2D_AbsoluteMaximum ( rhoData->rho ) <= grid->rhoTolerance )  goto EndOfLoop ;
            /* . Evaluate the functional terms. */
            DFTFunctionalModel_Evaluate ( functionalModel, rhoData ) ;
            /* . Accumulation and weighting of the integration data. */
            RealArray1D_Multiply ( &(rhoDataP->vRho), weights, NULL ) ;
            if ( functionalModel->hasLaplacian ) RealArray1D_Multiply ( &(rhoDataP->vLaplacianRho), weights, NULL ) ;
            if ( functionalModel->hasSigma     ) RealArray1D_Multiply ( &(rhoDataP->vSigma       ), weights, NULL ) ;
            if ( functionalModel->hasTau       ) RealArray1D_Multiply ( &(rhoDataP->vTau         ), weights, NULL ) ;
            if ( isSpinUnrestricted )
            {
                RealArray1D_Add      ( &(rhoDataP->rho ), 1.0e+00, &(rhoDataQ->rho), NULL ) ;
                RealArray1D_Multiply ( &(rhoDataQ->vRho),          weights         , NULL ) ;
                if ( functionalModel->hasLaplacian ) RealArray1D_Multiply ( &(rhoDataQ->vLaplacianRho), weights, NULL ) ;
                if ( functionalModel->hasSigma     ) RealArray1D_Multiply ( &(rhoDataQ->vSigma       ), weights, NULL ) ;
                if ( functionalModel->hasSigma     ) RealArray1D_Multiply ( &(rhoData->vSigmaPQ      ), weights, NULL ) ;
                if ( functionalModel->hasTau       ) RealArray1D_Multiply ( &(rhoDataQ->vTau         ), weights, NULL ) ;
            }
            /* . Total energy and density - eXC is multiplied by rhoP which is now the total density. */
            RealArray1D_Multiply ( rhoData->eXC, &(rhoDataP->rho), NULL ) ;
            eXCTotal += RealArray1D_Dot ( rhoData->eXC    , weights, NULL ) ;
            rhoTotal += RealArray1D_Dot ( &(rhoDataP->rho), weights, NULL ) ;
            /* . Fock terms. */
            if ( doFock )
            {
                DFTIntegrator_Fock ( functionalModel->hasSigma     ,
                                     functionalModel->hasLaplacian ,
                                     functionalModel->hasTau       ,
                                     basisData                     ,
                                     &(rhoDataP->dRhoX        )    ,
                                     &(rhoDataP->dRhoY        )    ,
                                     &(rhoDataP->dRhoZ        )    ,
                                     &(rhoDataP->vRho         )    ,
                                     &(rhoDataP->vSigma       )    ,
                                     &(rhoDataP->vLaplacianRho)    ,
                                     &(rhoDataP->vTau         )    ,
                                     fockA                         ,
                                     work2D                        ) ;
                if ( isSpinUnrestricted )
                {
                    DFTIntegrator_Fock      ( functionalModel->hasSigma     ,
                                              functionalModel->hasLaplacian ,
                                              functionalModel->hasTau       ,
                                              basisData                     ,
                                              &(rhoDataQ->dRhoX        )    ,
                                              &(rhoDataQ->dRhoY        )    ,
                                              &(rhoDataQ->dRhoZ        )    ,
                                              &(rhoDataQ->vRho         )    ,
                                              &(rhoDataQ->vSigma       )    ,
                                              &(rhoDataQ->vLaplacianRho)    ,
                                              &(rhoDataQ->vTau         )    ,
                                              fockB                         ,
                                              work2D                        ) ;
                    DFTIntegrator_FockSigma ( functionalModel->hasSigma     ,
                                              False                         ,
                                              basisData                     ,
                                              &(rhoDataQ->dRhoX  )          ,
                                              &(rhoDataQ->dRhoY  )          ,
                                              &(rhoDataQ->dRhoZ  )          ,
                                              &(rhoData->vSigmaPQ)          ,
                                              fockA                         ,
                                              work2D                        ) ;
                    DFTIntegrator_FockSigma ( functionalModel->hasSigma     ,
                                              False                         ,
                                              basisData                     ,
                                              &(rhoDataP->dRhoX  )          ,
                                              &(rhoDataP->dRhoY  )          ,
                                              &(rhoDataP->dRhoZ  )          ,
                                              &(rhoData->vSigmaPQ)          ,
                                              fockB                         ,
                                              work2D                        ) ;
                }
            }
            /* . Gradient terms. */
            if ( doGradients )
            {
                /* . Direct terms. */
                DFTIntegrator_Gradients ( functionalModel->hasSigma     ,
                                          functionalModel->hasLaplacian ,
                                          functionalModel->hasTau       ,
                                          atomIndices                   ,
                                          basisData                     ,
                                          &(rhoDataP->dRhoX        )    ,
                                          &(rhoDataP->dRhoY        )    ,
                                          &(rhoDataP->dRhoZ        )    ,
                                          &(rhoDataP->vRho         )    ,
                                          &(rhoDataP->vSigma       )    ,
                                          &(rhoDataP->vLaplacianRho)    ,
                                          &(rhoDataP->vTau         )    ,
                                          reducedDensityP               ,
                                          gridAtom                      ,
                                          gradients3                    ,
                                          temp2D                        ,
                                          work2D                        ) ;
                if ( isSpinUnrestricted )
                {   
                    DFTIntegrator_Gradients      ( functionalModel->hasSigma     ,
                                                   functionalModel->hasLaplacian ,
                                                   functionalModel->hasTau       ,
                                                   atomIndices                   ,
                                                   basisData                     ,
                                                   &(rhoDataQ->dRhoX        )    ,
                                                   &(rhoDataQ->dRhoY        )    ,
                                                   &(rhoDataQ->dRhoZ        )    ,
                                                   &(rhoDataQ->vRho         )    ,
                                                   &(rhoDataQ->vSigma       )    ,
                                                   &(rhoDataQ->vLaplacianRho)    ,
                                                   &(rhoDataQ->vTau         )    ,
                                                   reducedDensityQ               ,
                                                   gridAtom                      ,
                                                   gradients3                    ,
                                                   temp2D                        ,
                                                   work2D                        ) ;
                    DFTIntegrator_GradientsSigma ( functionalModel->hasSigma     ,
                                                   False                         ,
                                                   atomIndices                   ,
                                                   basisData                     ,
                                                   &(rhoDataQ->dRhoX )           ,
                                                   &(rhoDataQ->dRhoY )           ,
                                                   &(rhoDataQ->dRhoZ )           ,
                                                   &(rhoData->vSigmaPQ)          ,
                                                   reducedDensityP               ,
                                                   gridAtom                      ,
                                                   gradients3                    ,
                                                   temp2D                        ,
                                                   work2D                        ) ;
                    DFTIntegrator_GradientsSigma ( functionalModel->hasSigma     ,
                                                   False                         ,
                                                   atomIndices                   ,
                                                   basisData                     ,
                                                   &(rhoDataP->dRhoX )           ,
                                                   &(rhoDataP->dRhoY )           ,
                                                   &(rhoDataP->dRhoZ )           ,
                                                   &(rhoData->vSigmaPQ)          ,
                                                   reducedDensityQ               ,
                                                   gridAtom                      ,
                                                   gradients3                    ,
                                                   temp2D                        ,
                                                   work2D                        ) ;
                }
# ifdef _DFTGRIDWEIGHTDERIVATIVES
                /* . Weight terms. */
                if ( weightsWork == NULL )
                {
                    weightsWork = DFTGridWeightsDerivativesWork_Allocate ( grid->weights, &localStatus ) ;
                    if ( weightsWork == NULL ) goto EndOfLoop ;
                }
                DFTGridWeights_Derivatives ( grid->weights         ,
                                             gridAtom              ,
                                             block->numberOfPoints ,
                                             coordinates3          ,
                                             weights               ,
                                             rhoData->eXC          ,
                                             gradients3            ,
                                             weightsWork           ) ;
# endif
            }
            /* . End of loop. */
        EndOfLoop:
            /* . Store function data. */
            if ( storeFunctionData ) { block->functionData = basisData ; basisData = NULL ; }
            if ( localStatus != Status_OK ) break ;
        }
        /* . Deallocate space. */
    FinishUp:
        DFTGridWeightsDerivativesWork_Deallocate ( &weightsWork     ) ;
        DFTIntegratorDataBlock_Deallocate        ( &rhoData         ) ;
        IntegerArray1D_Deallocate                ( &atomIndices     ) ;
        IntegerArray1D_Deallocate                ( &basisIndices    ) ;
        RealArray1D_Deallocate                   ( &work1D          ) ;
        RealArray2D_Deallocate                   ( &reducedDensityP ) ;
        RealArray2D_Deallocate                   ( &reducedDensityQ ) ;
        RealArray2D_Deallocate                   ( &temp2D          ) ;
        RealArray2D_Deallocate                   ( &work2D          ) ;
        if ( doGradients || ( ! inCore ) ) GridFunctionDataBlock_Deallocate ( &basisData ) ;
        /* . Finish up. */
        if ( eQuad   != NULL ) (*eQuad  ) = eXCTotal ;
        if ( rhoQuad != NULL ) (*rhoQuad) = rhoTotal ;
        if ( localStatus != Status_OK ) Status_Set ( status, localStatus ) ;
    }
}

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Contributions to a Fock matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_Fock ( const Boolean                hasSigma     ,
                                 const Boolean                hasLaplacian ,
                                 const Boolean                hasTau       ,
                                 const GridFunctionDataBlock *basisData    ,
                                 const RealArray1D           *dRhoX        ,
                                 const RealArray1D           *dRhoY        ,
                                 const RealArray1D           *dRhoZ        ,
                                 const RealArray1D           *vRho         ,
                                 const RealArray1D           *vSigma       ,
                                 const RealArray1D           *vLaplacianRho,
                                 const RealArray1D           *vTau         ,
                                       SymmetricMatrix       *fock         ,
                                       RealArray2D           *work2D       )
{
    /* . Aliases. */
    IntegerArray1D *indices = basisData->indices ;
    RealArray2D *b   = basisData->f   ,
                *bX  = basisData->fX  ,
                *bY  = basisData->fY  ,
                *bZ  = basisData->fZ  ,
                *bXX = basisData->fXX ,
                *bYY = basisData->fYY ,
                *bZZ = basisData->fZZ ;
   /* . Rho. */
    RealArray2D_CopyTo       ( b, work2D, NULL ) ;
    URealArray2D_ColumnScale ( vRho, work2D ) ;
    USymmetricMatrix_DotProductIncrement ( indices, b, work2D, fock ) ;
    /* . Sigma. */
    DFTIntegrator_FockSigma ( hasSigma, True, basisData, dRhoX, dRhoY, dRhoZ, vSigma, fock, work2D ) ;
    /* . Laplacian. */
    if ( hasLaplacian )
    {
        /* . First derivative contribution. */
        RealArray2D_CopyTo       ( bX, work2D, NULL ) ;
        URealArray2D_ColumnScale ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale        ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bX, work2D, fock ) ;
        RealArray2D_CopyTo       ( bY, work2D, NULL ) ;
        URealArray2D_ColumnScale ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale        ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bY, work2D, fock ) ;
        RealArray2D_CopyTo       ( bZ, work2D, NULL ) ;
        URealArray2D_ColumnScale ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale        ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bZ, work2D, fock ) ;
        /* . Second derivative contribution. */
        RealArray2D_CopyTo         ( bXX, work2D, NULL ) ;
        RealArray2D_Add            ( work2D, 1.0e+00, bYY, NULL ) ;
        RealArray2D_Add            ( work2D, 1.0e+00, bZZ, NULL ) ;
        URealArray2D_ColumnScale   ( vLaplacianRho, work2D ) ;
        USymmetricMatrix_DotProductIncrement ( indices, b, work2D, fock ) ;
        USymmetricMatrix_DotProductIncrement ( indices, work2D, b, fock ) ;
    }
    /* . Tau. */
    if ( hasTau )
    {
        RealArray2D_CopyTo       ( bX, work2D, NULL ) ;
        URealArray2D_ColumnScale ( vTau, work2D ) ;
        RealArray2D_Scale        ( work2D, 0.5e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bX, work2D, fock ) ;
        RealArray2D_CopyTo       ( bY, work2D, NULL ) ;
        URealArray2D_ColumnScale ( vTau, work2D ) ;
        RealArray2D_Scale        ( work2D, 0.5e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bY, work2D, fock ) ;
        RealArray2D_CopyTo       ( bZ, work2D, NULL ) ;
        URealArray2D_ColumnScale ( vTau, work2D ) ;
        RealArray2D_Scale        ( work2D, 0.5e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, bZ, work2D, fock ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sigma contributions to a Fock matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_FockSigma ( const Boolean                hasSigma  ,
                                      const Boolean                scale     ,
                                      const GridFunctionDataBlock *basisData ,
                                      const RealArray1D           *dRhoX     ,
                                      const RealArray1D           *dRhoY     ,
                                      const RealArray1D           *dRhoZ     ,
                                      const RealArray1D           *vSigma    ,
                                            SymmetricMatrix       *fock      ,
                                            RealArray2D           *work2D    )
{
    if ( hasSigma )
    {
        /* . Aliases. */
        auto IntegerArray1D *indices = basisData->indices ;
        auto RealArray2D    *b  = basisData->f  ,
                            *bX = basisData->fX ,
                            *bY = basisData->fY ,
                            *bZ = basisData->fZ ;
        RealArray2D_CopyTo                ( bX, work2D, NULL ) ;
        URealArray2D_ColumnScale          ( dRhoX , work2D ) ;
        URealArray2D_ColumnAddScaledArray ( bY, dRhoY, work2D ) ;
        URealArray2D_ColumnAddScaledArray ( bZ, dRhoZ, work2D ) ;
        URealArray2D_ColumnScale          ( vSigma, work2D ) ;
        if ( scale ) RealArray2D_Scale    ( work2D, 2.0e+00 ) ;
        USymmetricMatrix_DotProductIncrement ( indices, b, work2D, fock ) ;
        USymmetricMatrix_DotProductIncrement ( indices, work2D, b, fock ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form a square reduced density matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_FormReducedDensity ( const IntegerArray1D   *indices        ,
                                               const SymmetricMatrix  *density        ,
                                                     RealArray2D     **reducedDensity ,
                                                     Status           *status         )
{
    Integer      n = View1D_Extent ( indices ) ;
    RealArray2D *new = NULL, *old = (*reducedDensity) ;
    if ( ( old != NULL ) && ( View2D_Rows ( old ) == n ) && ( View2D_Columns ( old ) == n ) ) new = old ;
    else
    {
        RealArray2D_Deallocate ( reducedDensity ) ;
        new = RealArray2D_AllocateWithExtents ( n, n, status ) ;
    }
    if ( new != NULL ) USymmetricMatrix_IndexedCopyToRealArray2D ( density, indices, new, NULL ) ;
    (*reducedDensity) = new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Contributions to the gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_Gradients ( const Boolean                hasSigma      ,
                                      const Boolean                hasLaplacian  ,
                                      const Boolean                hasTau        ,
                                      const IntegerArray1D        *atomIndices   ,
                                      const GridFunctionDataBlock *basisData     ,
                                      const RealArray1D           *dRhoX         ,
                                      const RealArray1D           *dRhoY         ,
                                      const RealArray1D           *dRhoZ         ,
                                      const RealArray1D           *vRho          ,
                                      const RealArray1D           *vSigma        ,
                                      const RealArray1D           *vLaplacianRho ,
                                      const RealArray1D           *vTau          ,
                                      const RealArray2D           *density       ,
                                      const Integer                gridAtom      ,
                                            Coordinates3          *gradients3    ,
                                            RealArray2D           *temp2D        ,
                                            RealArray2D           *work2D        )
{
    /* . Aliases. */
    auto IntegerArray1D *indices = basisData->indices ;
    auto RealArray2D *b   = basisData->f   ,
                     *bX  = basisData->fX  ,
                     *bY  = basisData->fY  ,
                     *bZ  = basisData->fZ  ,
                     *bXX = basisData->fXX ,
                     *bXY = basisData->fXY ,
                     *bXZ = basisData->fXZ ,
                     *bYY = basisData->fYY ,
                     *bYZ = basisData->fYZ ,
                     *bZZ = basisData->fZZ ;
    /* . Rho. */
    RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, b, 0.0e+00, work2D, NULL ) ;
    URealArray2D_ColumnScale   ( vRho, work2D ) ;
    RealArray2D_Scale          ( work2D, 2.0e+00 ) ;
    UGradientContributions     ( atomIndices, indices, work2D, bX, bY, bZ, gridAtom, gradients3 ) ;
    /* . Sigma. */
    if ( hasSigma )
    {
        DFTIntegrator_GradientsSigma ( hasSigma, True, atomIndices, basisData, dRhoX, dRhoY, dRhoZ, vSigma, density, gridAtom, gradients3, temp2D, work2D ) ;
    }
    /* . Laplacian. */
    if ( hasLaplacian )
    {
        /* . More aliases. */
        auto RealArray2D *bXXX = basisData->fXXX ,
                         *bXXY = basisData->fXXY ,
                         *bXXZ = basisData->fXXZ ,
                         *bXYY = basisData->fXYY ,
                         *bXYZ = basisData->fXYZ ,
                         *bXZZ = basisData->fXZZ ,
                         *bYYY = basisData->fYYY ,
                         *bYYZ = basisData->fYYZ ,
                         *bYZZ = basisData->fYZZ ,
                         *bZZZ = basisData->fZZZ ,
                         *sum, *sumX, *sumY, *sumZ ;
        /* . Aliases - some of the integrals are destroyed. */
        sum = bXYZ ; sumX = bXXX ; sumY = bYYY ; sumZ = bZZZ ;
        /* . First derivative contribution. */
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, bX, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale          ( work2D, 4.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXX, bXY, bXZ, gridAtom, gradients3 ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, bY, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale          ( work2D, 4.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXY, bYY, bYZ, gridAtom, gradients3 ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, bZ, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale          ( work2D, 4.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXZ, bYZ, bZZ, gridAtom, gradients3 ) ;
        /* . Second derivative contribution - 1. */
        RealArray2D_CopyTo         ( bXX, sum, NULL ) ;
        RealArray2D_Add            ( sum, 1.0e+00, bYY, NULL ) ;
        RealArray2D_Add            ( sum, 1.0e+00, bZZ, NULL ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, sum, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale          ( work2D, 2.0e+00 ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bX, bY, bZ, gridAtom, gradients3 ) ;
        /* . Second derivative contribution - 2. */
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, b, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vLaplacianRho, work2D ) ;
        RealArray2D_Scale          ( work2D, 2.0e+00 ) ;
        RealArray2D_Add            ( sumX, 1.0e+00, bXYY, NULL ) ;
        RealArray2D_Add            ( sumX, 1.0e+00, bXZZ, NULL ) ;
        RealArray2D_Add            ( sumY, 1.0e+00, bXXY, NULL ) ;
        RealArray2D_Add            ( sumY, 1.0e+00, bYZZ, NULL ) ;
        RealArray2D_Add            ( sumZ, 1.0e+00, bXXZ, NULL ) ;
        RealArray2D_Add            ( sumZ, 1.0e+00, bYYZ, NULL ) ;
        UGradientContributions     ( atomIndices, indices, work2D, sumX, sumY, sumZ, gridAtom, gradients3 ) ;
    }
    /* . Tau. */
    if ( hasTau )
    {
        /* . No scaling as factors of 2 and 1/2. */
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, bX, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vTau, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXX, bXY, bXZ, gridAtom, gradients3 ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, bY, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vTau, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXY, bYY, bYZ, gridAtom, gradients3 ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, bZ, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vTau, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXZ, bYZ, bZZ, gridAtom, gradients3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sigma contributions to the gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_GradientsSigma ( const Boolean                hasSigma    ,
                                           const Boolean                isSelf      ,
                                           const IntegerArray1D        *atomIndices ,
                                           const GridFunctionDataBlock *basisData   ,
                                           const RealArray1D           *dRhoX       ,
                                           const RealArray1D           *dRhoY       ,
                                           const RealArray1D           *dRhoZ       ,
                                           const RealArray1D           *vSigma      ,
                                           const RealArray2D           *density     ,
                                           const Integer                gridAtom    ,
                                                 Coordinates3          *gradients3  ,
                                                 RealArray2D           *temp2D      ,
                                                 RealArray2D           *work2D      )
{
    if ( hasSigma )
    {
        auto Real factor ;
        /* . Aliases. */
        auto IntegerArray1D *indices = basisData->indices ;
        auto RealArray2D *b   = basisData->f   ,
                         *bX  = basisData->fX  ,
                         *bY  = basisData->fY  ,
                         *bZ  = basisData->fZ  ,
                         *bXX = basisData->fXX ,
                         *bXY = basisData->fXY ,
                         *bXZ = basisData->fXZ ,
                         *bYY = basisData->fYY ,
                         *bYZ = basisData->fYZ ,
                         *bZZ = basisData->fZZ ;
        if ( isSelf ) factor = 4.0e+00 ;
        else          factor = 2.0e+00 ;
        /* . Contribution 1. */
        RealArray2D_CopyTo         ( bX, temp2D, NULL ) ;
        URealArray2D_ColumnScale   ( dRhoX, temp2D ) ;
        URealArray2D_ColumnAddScaledArray ( bY, dRhoY, temp2D ) ;
        URealArray2D_ColumnAddScaledArray ( bZ, dRhoZ, temp2D ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, temp2D, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( vSigma, work2D ) ;
        RealArray2D_Scale          ( work2D, factor ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bX, bY, bZ, gridAtom, gradients3 ) ;
        /* . Set up for contributions 2, 3 and 4. */
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, density, b, 0.0e+00, temp2D, NULL ) ;
        URealArray2D_ColumnScale   ( vSigma, temp2D ) ;
        RealArray2D_Scale          ( temp2D, factor ) ;
        /* . Contribution 2. */
        RealArray2D_CopyTo         ( temp2D, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( dRhoX, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXX, bXY, bXZ, gridAtom, gradients3 ) ;
        /* . Contribution 3. */
        RealArray2D_CopyTo         ( temp2D, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( dRhoY, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXY, bYY, bYZ, gridAtom, gradients3 ) ;
        /* . Contribution 4. */
        RealArray2D_CopyTo         ( temp2D, work2D, NULL ) ;
        URealArray2D_ColumnScale   ( dRhoZ, work2D ) ;
        UGradientContributions     ( atomIndices, indices, work2D, bXZ, bYZ, bZZ, gridAtom, gradients3 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the values of a single density and associated quantities at the grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_GridPointRho ( const Boolean                hasSigma     ,
                                         const Boolean                hasLaplacian ,
                                         const Boolean                hasTau       ,
                                         const GridFunctionDataBlock *basisData    ,
                                         const RealArray2D           *density      ,
                                               RealArray1D           *rho          ,
                                               RealArray1D           *dRhoX        ,
                                               RealArray1D           *dRhoY        ,
                                               RealArray1D           *dRhoZ        ,
                                               RealArray1D           *sigma        ,
                                               RealArray1D           *laplacianRho ,
                                               RealArray1D           *tau          ,
                                               RealArray1D           *work1D       ,
                                               RealArray2D           *work2D       )
{
    /* . Aliases. */
    RealArray2D *b   = basisData->f   ,
                *bX  = basisData->fX  ,
                *bY  = basisData->fY  ,
                *bZ  = basisData->fZ  ,
                *bXX = basisData->fXX ,
                *bYY = basisData->fYY ,
                *bZZ = basisData->fZZ ;
    /* . Rho. */
    RealArray2D_MatrixMultiply     ( False, False, 1.0e+00, density, b, 0.0e+00, work2D, NULL ) ;
    URealArray2D_ColumnDotProducts ( True, b, work2D, rho ) ;
    /* . Sigma. */
    if ( hasSigma )
    {
        URealArray2D_ColumnDotProducts ( True, bX, work2D, dRhoX ) ;
        URealArray2D_ColumnDotProducts ( True, bY, work2D, dRhoY ) ;
        URealArray2D_ColumnDotProducts ( True, bZ, work2D, dRhoZ ) ;
        DFTIntegrator_GridPointSigma ( True, dRhoX, dRhoY, dRhoZ, dRhoX, dRhoY, dRhoZ, sigma ) ;
        RealArray1D_Scale ( dRhoX, 2.0e+00 ) ;
        RealArray1D_Scale ( dRhoY, 2.0e+00 ) ;
        RealArray1D_Scale ( dRhoZ, 2.0e+00 ) ;
        RealArray1D_Scale ( sigma, 4.0e+00 ) ;
    }
    /* . Laplacian or tau. */
    if ( hasLaplacian || hasTau )
    {
        RealArray2D_MatrixMultiply     ( False, False, 1.0e+00, density, bX, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnDotProducts ( True , bX, work2D, work1D ) ;
        RealArray2D_MatrixMultiply     ( False, False, 1.0e+00, density, bY, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnDotProducts ( False, bY, work2D, work1D ) ;
        RealArray2D_MatrixMultiply     ( False, False, 1.0e+00, density, bZ, 0.0e+00, work2D, NULL ) ;
        URealArray2D_ColumnDotProducts ( False, bZ, work2D, work1D ) ;
    }
    /* . Laplacian. */
    if ( hasLaplacian )
    {
        RealArray1D_CopyTo ( work1D, laplacianRho, NULL ) ;
        URealArray2D_ColumnDotProducts ( False, bXX, work2D, laplacianRho ) ;
        URealArray2D_ColumnDotProducts ( False, bYY, work2D, laplacianRho ) ;
        URealArray2D_ColumnDotProducts ( False, bZZ, work2D, laplacianRho ) ;
        RealArray1D_Scale  ( laplacianRho, 2.0e+00 ) ;
    }
    /* . Tau. */
    if ( hasTau )
    {
        RealArray1D_CopyTo ( work1D, tau, NULL ) ;
        RealArray1D_Scale ( tau, 0.5e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the cross-sigma values for two densities at the grid points.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTIntegrator_GridPointSigma ( const Boolean      hasSigma ,
                                           const RealArray1D *dRhoXa   ,
                                           const RealArray1D *dRhoYa   ,
                                           const RealArray1D *dRhoZa   ,
                                           const RealArray1D *dRhoXb   ,
                                           const RealArray1D *dRhoYb   ,
                                           const RealArray1D *dRhoZb   ,
                                                 RealArray1D *sigma    )
{
    if ( hasSigma )
    {
        auto Integer p ;
        for ( p = 0 ; p < View1D_Extent ( sigma ) ; p++ )
        {
            Array1D_Item ( sigma, p ) = Array1D_Item ( dRhoXa, p ) * Array1D_Item ( dRhoXb, p ) +
                                        Array1D_Item ( dRhoYa, p ) * Array1D_Item ( dRhoYb, p ) +
                                        Array1D_Item ( dRhoZa, p ) * Array1D_Item ( dRhoZb, p ) ;
        }
    }
}

/*==================================================================================================================================
! . Local utilities - may be generalized and moved.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Utility procedure for determining gradient contributions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void UGradientContributions ( const IntegerArray1D *atomIndices ,
                                     const IntegerArray1D *indices     ,
                                     const RealArray2D    *a           ,
                                     const RealArray2D    *x           ,
                                     const RealArray2D    *y           ,
                                     const RealArray2D    *z           ,
                                     const Integer         gridAtom    ,
                                           Coordinates3   *gradients3  )
{
    Integer     c, i, m ;
    Real        gX, gY, gZ ;
    RealArray1D viewA, viewB ;
    for ( i = 0 ; i < View1D_Extent ( indices ) ; i++ )
    {
        m = Array1D_Item ( indices    , i ) ;
        c = Array1D_Item ( atomIndices, m ) ;
        /* . Calculation the contributions. */
        RealArray2D_RowView  ( a     , i, False, &viewA, NULL ) ;
        RealArray2D_RowView  ( x     , i, False, &viewB, NULL ) ;
        gX = RealArray1D_Dot ( &viewA, &viewB, NULL ) ;
        RealArray2D_RowView  ( y     , i, False, &viewB, NULL ) ;
        gY = RealArray1D_Dot ( &viewA, &viewB, NULL ) ;
        RealArray2D_RowView  ( z     , i, False, &viewB, NULL ) ;
        gZ = RealArray1D_Dot ( &viewA, &viewB, NULL ) ;
        /* . Add in the contributions. */
        Coordinates3_DecrementRow ( gradients3, c, gX, gY, gZ ) ;
# ifdef _DFTGRIDWEIGHTDERIVATIVES
        /* . This is the dE/drg term (i.e. the derivative with respect to the grid point which belongs to gridAtom). */
        Coordinates3_IncrementRow ( gradients3, gridAtom, gX, gY, gZ ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add scaled array by column for a 2-D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void URealArray2D_ColumnAddScaledArray ( const RealArray2D *a       ,
                                                const RealArray1D *weights ,
                                                      RealArray2D *b       )
{
    Integer     i ;
    Real        w ;
    RealArray1D viewA, viewB ;
    for ( i = 0 ; i < View1D_Extent ( weights ) ; i++ )
    {
        w = Array1D_Item ( weights, i ) ;
        RealArray2D_ColumnView ( a, i, False, &viewA, NULL ) ;
        RealArray2D_ColumnView ( b, i, False, &viewB, NULL ) ;
        RealArray1D_Add        ( &viewB, w, &viewA, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form a 1-D array from the dot products of the columns of 2 2-D arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void URealArray2D_ColumnDotProducts ( const Boolean      initialize ,
                                             const RealArray2D *a          ,
                                             const RealArray2D *b          ,
                                                   RealArray1D *c          )
{
    Integer p ;
    RealArray1D viewA, viewB ;
    if ( initialize ) RealArray1D_Set ( c, 0.0e+00 ) ;
    for ( p = 0 ; p < View1D_Extent ( c ) ; p++ )
    {
        RealArray2D_ColumnView ( a, p, False, &viewA, NULL ) ;
        RealArray2D_ColumnView ( b, p, False, &viewB, NULL ) ;
        Array1D_Item ( c, p ) += RealArray1D_Dot ( &viewA, &viewB, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the columns of a 2-D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void URealArray2D_ColumnScale ( const RealArray1D *a ,
                                             RealArray2D *b )
{
    Integer     i ;
    Real        w ;
    RealArray1D view ;
    for ( i = 0 ; i < View1D_Extent ( a ) ; i++ )
    {
        w = Array1D_Item       ( a, i ) ;
        RealArray2D_ColumnView ( b, i, False, &view, NULL ) ;
        RealArray1D_Scale      ( &view, w ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment selected entries of a symmetric matrix using the dot products of the rows of two 2-D arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void USymmetricMatrix_DotProductIncrement ( const IntegerArray1D  *indices ,
                                                   const RealArray2D     *a       ,
                                                   const RealArray2D     *b       ,
                                                         SymmetricMatrix *c       )
{
    Integer     i, j, m, n ;
    RealArray1D viewA, viewB ;
    for ( i = 0 ; i < View1D_Extent ( indices ) ; i++ )
    {
        m = Array1D_Item ( indices, i ) ;
        RealArray2D_RowView ( a, i, False, &viewA, NULL ) ;
        for ( j = 0 ; j <= i ; j++ )
        {
            n = Array1D_Item ( indices, j ) ;
            RealArray2D_RowView ( b, j, False, &viewB, NULL ) ;
            SymmetricMatrix_Item ( c, m, n ) += RealArray1D_Dot ( &viewA, &viewB, NULL ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy selected elements of the matrix to a square form.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void USymmetricMatrix_IndexedCopyToRealArray2D ( const SymmetricMatrix *self    ,
                                                        const IntegerArray1D  *indices ,
                                                              RealArray2D     *target  ,
                                                              Status          *status  )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( target != NULL ) )
    {
        auto Integer i, j, m, n ;
        auto Real    v ;
        for ( i = 0 ; i < View1D_Extent ( indices ) ; i++ )
        {
            m = Array1D_Item ( indices, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                n = Array1D_Item         ( indices, j ) ;
                v = SymmetricMatrix_Item ( self, m, n ) ;
                Array2D_Item ( target, i, j ) = v ;
                Array2D_Item ( target, j, i ) = v ;
            }
        }
    }
}
