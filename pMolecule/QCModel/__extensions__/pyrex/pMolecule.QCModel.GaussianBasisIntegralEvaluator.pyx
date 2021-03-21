"""Gaussian basis integral evaluator."""

from  pCore                     import logFile              , \
                                       LogFileActive        , \
                                       RawObjectConstructor
from  pScientific               import Units
from  pScientific.Arrays        import Array                , \
                                       StorageType
from  pScientific.LinearAlgebra import MatrixPowerInverse
from .QCDefinitions             import BasisRepresentation
from .QCModelError              import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisIntegralEvaluator:

    def CoreCoreEnergy ( self, target ):
        """The core-core energy."""
        cdef Coordinates3    coordinates3
        cdef Coordinates3    gradients3
        cdef RealArray1D     charges
        cdef CRealArray2D   *cGradients3 = NULL
        cdef CReal           energy      = 0.0
        charges      = target.qcState.nuclearCharges
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        gradients3   = scratch.Get ( "qcGradients3AU", None )
        if gradients3 is not None: cGradients3 = gradients3.cObject
        energy = GaussianBasisContainer_NuclearNuclearEnergy ( charges.cObject      ,
                                                               charges.cObject      ,
                                                               coordinates3.cObject ,
                                                               coordinates3.cObject ,
                                                               NULL                 ,
                                                               NULL                 ,
                                                               NULL                 ,
                                                               NULL                 ,
                                                               NULL                 ,
                                                               NULL                 ,
                                                               cGradients3          ,
                                                               cGradients3          )
        scratch.qcEnergyReport["QC Core-Core"] =   energy
        scratch.energyTerms   ["QC Core-Core"] = ( energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def DipoleIntegrals ( self, target, Vector3 center = None ):
        """The dipole integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef SymmetricMatrix        dX
        cdef SymmetricMatrix        dY
        cdef SymmetricMatrix        dZ
        cdef CInteger               n
        cdef CRealArray1D          *cCenter = NULL
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        if center is not None: cCenter = center.cObject
        n  = len ( bases )
        dX = SymmetricMatrix.WithExtent ( n )
        dY = SymmetricMatrix.WithExtent ( n )
        dZ = SymmetricMatrix.WithExtent ( n )
        GaussianBasisContainerIntegrals_Dipole ( bases.cObject        ,
                                                 basisIndices.cObject ,
                                                 coordinates3.cObject ,
                                                 cCenter              ,
                                                 dX.cObject           ,
                                                 dY.cObject           ,
                                                 dZ.cObject           ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating dipole integrals." )
        return ( dX, dY, dZ )

    def ElectronFitGradients ( self, target ):
        """The electron-fit gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer fBases
        cdef GaussianBasisContainer oBases
        cdef IntegerArray1D         fIndices
        cdef IntegerArray1D         oIndices
        cdef RealArray1D            fPotential
        cdef SymmetricMatrix        dTotal
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            fBases       = target.qcState.fitBases
            oBases       = target.qcState.orbitalBases
            fIndices     = fBases._centerFunctionPointers[self.basisRepresentation]
            oIndices     = oBases._centerFunctionPointers[self.basisRepresentation]
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            fPotential   = scratch.fitPotential
            gradients3   = scratch.qcGradients3AU
            GaussianBasisContainerIntegrals_ElectronFitD ( oBases.cObject       ,
                                                           oIndices.cObject     ,
                                                           fBases.cObject       ,
                                                           fIndices.cObject     ,
                                                           coordinates3.cObject ,
                                                           dTotal.cObject       ,
                                                           fPotential.cObject   ,
                                                           gradients3.cObject   ,
                                                           &cStatus             )
            if cStatus != CStatus_OK: raise QCModelError ( "Error calculating electron-fit gradients." )

    def ElectronFitIntegrals ( self, target ):
        """The electron-fit integrals."""
        cdef BlockStorage           fitIntegrals
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer fBases
        cdef GaussianBasisContainer oBases
        cdef IntegerArray1D         fIndices
        cdef IntegerArray1D         oIndices
        cdef CStatus                cStatus = CStatus_OK
        fBases       = target.qcState.fitBases
        oBases       = target.qcState.orbitalBases
        fIndices     = fBases._centerFunctionPointers[self.basisRepresentation]
        oIndices     = oBases._centerFunctionPointers[self.basisRepresentation]
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        fitIntegrals = scratch.Get ( "electronFitIntegrals", None )
        if fitIntegrals is None:
            fitIntegrals = BlockStorage ( )
            scratch.electronFitIntegrals = fitIntegrals
        else:
            fitIntegrals.Empty ( )
        GaussianBasisContainerIntegrals_ElectronFit ( oBases.cObject       ,
                                                      oIndices.cObject     ,
                                                      fBases.cObject       ,
                                                      fIndices.cObject     ,
                                                      coordinates3.cObject ,
                                                      fitIntegrals.cObject ,
                                                      &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating electron-fit integrals." )
        f        = float ( len ( fBases ) )
        o        = float ( len ( oBases ) )
        n        = fitIntegrals.count
        ( s, m ) = fitIntegrals.byteSize
        report = scratch.qcEnergyReport
        report["Electron-Fit Integrals"            ] = n
        report["Electron-Fit Integral Sparsity (%)"] = ( ( 1.0 - float ( n ) / ( f * ( o * ( o + 1.0 ) ) / 2.0 ) ) * 100.0, "{:.1f}" )
        report["Electron-Fit Integral Storage ({:s}B)".format ( m.symbol )] =  ( s, "{:.3f}" )

    def ElectronNuclearGradients ( self, target ):
        """The electron-nuclear gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef RealArray1D            charges
        cdef SymmetricMatrix        dTotal
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
            charges      = target.qcState.nuclearCharges
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            GaussianBasisContainerIntegrals_ElectronNuclearD ( bases.cObject        ,
                                                               basisIndices.cObject ,
                                                               charges.cObject      ,
                                                               NULL                 ,
                                                               NULL                 ,
                                                               coordinates3.cObject ,
                                                               coordinates3.cObject ,
                                                               NULL                 ,
                                                               dTotal.cObject       ,
                                                               gradients3.cObject   ,
                                                               gradients3.cObject   ,
                                                               &cStatus             )
            if cStatus != CStatus_OK: raise QCModelError ( "Error calculating electron-nuclear gradients." )

    def ElectronNuclearIntegrals ( self, target ):
        """The electron-nuclear integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef RealArray1D            charges
        cdef SymmetricMatrix        oneElectronMatrix
        cdef CStatus                cStatus = CStatus_OK
        bases             = target.qcState.orbitalBases
        basisIndices      = bases._centerFunctionPointers[self.basisRepresentation]
        charges           = target.qcState.nuclearCharges
        scratch           = target.scratch
        coordinates3      = scratch.qcCoordinates3AU
        oneElectronMatrix = scratch.oneElectronMatrix
        GaussianBasisContainerIntegrals_ElectronNuclear ( bases.cObject             ,
                                                          basisIndices.cObject      ,
                                                          charges.cObject           ,
                                                          NULL                      ,
                                                          NULL                      ,
                                                          coordinates3.cObject      ,
                                                          coordinates3.cObject      ,
                                                          NULL                      ,
                                                          oneElectronMatrix.cObject ,
                                                          &cStatus                  )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating electron-nuclear integrals." )

    def ElectronPotentials ( self, target, SymmetricMatrix density    not None ,
                                           Coordinates3    gridPoints not None ,
                                           RealArray1D     potentials not None ):
        """The potentials due to an electron density at a set of grid points."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef RealArray1D            charges
        cdef SymmetricMatrix        oneElectronMatrix
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
        coordinates3 = target.scratch.qcCoordinates3AU
        GaussianBasisContainerIntegrals_ElectronNuclearPotentials ( bases.cObject        ,
                                                                    basisIndices.cObject ,
                                                                    NULL                 ,
                                                                    NULL                 ,
                                                                    coordinates3.cObject ,
                                                                    gridPoints.cObject   ,
                                                                    NULL                 ,
                                                                    density.cObject      ,
                                                                    potentials.cObject   ,
                                                                    &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating electron potentials." )

    def FitFitGradients ( self, target ):
        """The fit-fit gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer fBases
        cdef IntegerArray1D         fIndices
        cdef RealArray1D            fPotential
        cdef SymmetricMatrix        dTotal
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            fBases       = target.qcState.fitBases
            fIndices     = fBases._centerFunctionPointers[self.basisRepresentation]
            coordinates3 = scratch.qcCoordinates3AU
            fPotential   = scratch.fitPotential
            gradients3   = scratch.qcGradients3AU
            GaussianBasisContainerIntegrals_2CoulombD ( fBases.cObject       ,
                                                        fIndices.cObject     ,
                                                        coordinates3.cObject ,
                                                        fPotential.cObject   ,
                                                        NULL                 ,
                                                        gradients3.cObject   ,
                                                        &cStatus             )
            if cStatus != CStatus_OK: raise QCModelError ( "Error calculating fit-fit gradients." )

    def FitFitIntegrals ( self, target, eigenValueTolerance = 1.0e-5 ):
        """Calculate the inverse fit matrix from the fit-fit and fit self-overlap integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer fBases
        cdef IntegerArray1D         fIndices
        cdef RealArray1D            fso
        cdef SymmetricMatrix        eri
        cdef SymmetricMatrix        inverse
        cdef CStatus                cStatus = CStatus_OK
        fBases       = target.qcState.fitBases
        fIndices     = fBases._centerFunctionPointers[self.basisRepresentation]
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        # . Create the fit matrix (n+1)x(n+1).
        n   = len ( fBases )
        eri = Array.WithExtent ( n+1, storageType = StorageType.Symmetric ) ; eri.Set ( 0.0 )
#NoConstraint        eri = Array.WithExtent ( n, storageType = StorageType.Symmetric ) ; eri.Set ( 0.0 )
        fso = Array.WithExtent ( n                                        ) ; fso.Set ( 0.0 )
        GaussianBasisContainerIntegrals_2Coulomb    ( fBases.cObject       ,
                                                      fIndices.cObject     ,
                                                      coordinates3.cObject ,
                                                      eri.cObject          ,
                                                      &cStatus             )
# . Remove self-overlap lines if no constraint.
        GaussianBasisContainerIntegrals_SelfOverlap ( fBases.cObject       ,
                                                      fIndices.cObject     ,
                                                      fso.cObject          ,
                                                      &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating fit-fit integrals." )
        for i in range ( n ): eri[i,n] = fso[i]
        scratch.fitMatrix = eri
        # . Use inverse fit matrix instead (less accurate).
        #scratch.fitMatrix = Array.WithExtent ( eri.extent, storageType = StorageType.Symmetric )
        #MatrixPowerInverse ( eri, 1.0, eigenValueTolerance, scratch.fitMatrix )

    def GridValues ( self, target, Coordinates3 gridPoints not None ,
                                   RealArray2D  values     not None ):
        """Calculate the values of the basis functions at a set of grid points."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef CStatus                cStatus = CStatus_OK
        bases                     = target.qcState.orbitalBases
        basisIndices              = bases._centerFunctionPointers[self.basisRepresentation]
        coordinates3              = target.scratch.qcCoordinates3AU
        GaussianBasisContainerIntegrals_Grid ( bases.cObject        ,
                                               basisIndices.cObject ,
                                               coordinates3.cObject ,
                                               gridPoints.cObject   ,
                                               values.cObject       ,
                                               &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating basis function grid point values." )

    def KineticOverlapGradients ( self, target ):
        """The kinetic energy and overlap gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef SymmetricMatrix        dTotal
        cdef SymmetricMatrix        wTotal
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
            charges      = target.qcState.nuclearCharges
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            wTotal       = scratch.weightedDensity
            GaussianBasisContainerIntegrals_Kinetic2OverlapD ( bases.cObject        ,
                                                               basisIndices.cObject ,
                                                               coordinates3.cObject ,
                                                               dTotal.cObject       ,
                                                               wTotal.cObject       ,
                                                               gradients3.cObject   ,
                                                               &cStatus             )
            if cStatus != CStatus_OK: raise QCModelError ( "Error calculating kinetic energy and overlap gradients." )

    def KineticOverlapIntegrals ( self, target ):
        """The kinetic energy and overlap integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef SymmetricMatrix        oneElectronMatrix
        cdef SymmetricMatrix        overlap
        cdef CStatus                cStatus = CStatus_OK
        bases             = target.qcState.orbitalBases
        basisIndices      = bases._centerFunctionPointers[self.basisRepresentation]
        scratch           = target.scratch
        coordinates3      = scratch.qcCoordinates3AU
        oneElectronMatrix = scratch.oneElectronMatrix
        overlap           = scratch.overlapMatrix
        GaussianBasisContainerIntegrals_Kinetic2Overlap ( bases.cObject             ,
                                                          basisIndices.cObject      ,
                                                          coordinates3.cObject      ,
                                                          oneElectronMatrix.cObject ,
                                                          overlap.cObject           ,
                                                          &cStatus                  )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating kinetic energy and overlap integrals." )

    def NuclearPotentials ( self, target, Coordinates3 gridPoints not None ,
                                          RealArray1D  potentials not None ):
        """The potentials due to the nuclei at a set of grid points."""
        cdef Coordinates3 coordinates3
        cdef RealArray1D  charges
        charges      = target.qcState.nuclearCharges
        coordinates3 = target.scratch.qcCoordinates3AU
        GaussianBasisContainer_NuclearNuclearPotentials ( charges.cObject      ,
                                                          gridPoints.cObject   ,
                                                          coordinates3.cObject ,
                                                          NULL                 ,
                                                          NULL                 ,
                                                          NULL                 ,
                                                          NULL                 ,
                                                          NULL                 ,
                                                          NULL                 ,
                                                          potentials.cObject   )          

    def OverlapIntegrals ( self, target ):
        """The overlap integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef SymmetricMatrix        overlap
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        overlap      = scratch.overlapMatrix
        GaussianBasisContainerIntegrals_2Overlap ( bases.cObject        ,
                                                   basisIndices.cObject ,
                                                   coordinates3.cObject ,
                                                   overlap.cObject      ,
                                                   &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating overlap integrals." )

    def TwoElectronGradients ( self, target, doCoulomb = True, CReal exchangeScaling = 1.0 ):
        """The two-electron gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef SymmetricMatrix        dSpin
        cdef SymmetricMatrix        dTotal
        cdef CBoolean               cDoCoulomb
        cdef CStatus                cStatus = CStatus_OK
        cdef CSymmetricMatrix      *cDSpin  = NULL
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
            scratch      = target.scratch
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            if hasattr ( scratch, "onePDMQ" ):
                dSpin  = scratch.onePDMQ.density
                cDSpin = dSpin.cObject
            if doCoulomb: cDoCoulomb = CTrue
            else:         cDoCoulomb = CFalse
            GaussianBasisContainerIntegrals_TEIsD ( bases.cObject        ,
                                                    basisIndices.cObject ,
                                                    coordinates3.cObject ,
                                                    dTotal.cObject       ,
                                                    cDSpin               ,
                                                    cDoCoulomb           ,
                                                    exchangeScaling      ,
                                                    gradients3.cObject   ,
                                                    &cStatus             )
            if cStatus != CStatus_OK: raise QCModelError ( "Error calculating two-electron gradients." )

    def TwoElectronIntegrals ( self, target ):
        """The two-electron integrals."""
        cdef BlockStorage           teis
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef IntegerArray1D         basisIndices
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        teis         = scratch.Get ( "twoElectronIntegrals", None )
        if teis is None:
            teis = BlockStorage ( )
            scratch.twoElectronIntegrals = teis
        else:
            teis.Empty ( )
        GaussianBasisContainerIntegrals_TEIs ( bases.cObject        ,
                                               basisIndices.cObject ,
                                               coordinates3.cObject ,
                                               teis.cObject         ,
                                               &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating two-electron integrals." )
        n        = float ( len ( bases ) )
        p        = ( n * ( n + 1.0 ) ) / 2.0
        n        = teis.count
        ( s, m ) = teis.byteSize
        report   = scratch.qcEnergyReport
        report["TEI Number"      ] = n
        report["TEI Sparsity (%)"] = ( ( 1.0 - float ( n ) / ( ( p * ( p + 1.0 ) ) / 2.0 ) ) * 100.0, "{:.1f}" )
        report["TEI Storage ({:s}B)".format ( m.symbol )] = ( s, "{:.3f}" )

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Work

