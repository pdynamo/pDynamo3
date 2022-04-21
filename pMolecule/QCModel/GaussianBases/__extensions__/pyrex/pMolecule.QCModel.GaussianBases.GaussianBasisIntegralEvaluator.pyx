"""Gaussian basis integral evaluator."""

from  pCore                     import logFile               , \
                                       LogFileActive         , \
                                       RawObjectConstructor
from  pScientific               import Units
from  pScientific.Arrays        import Array                 , \
                                       StorageType
from  pScientific.LinearAlgebra import MatrixPowerInverse
from .GaussianBasis             import GaussianBasisOperator
from .GaussianBasisError        import GaussianBasisError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisIntegralEvaluator:

    def f1Cp1V ( self, target, RealArray1D  fitCoefficients not None ,
                               Coordinates3 gridPoints      not None ,
                               RealArray1D  potentials      not None ,
                               fitBases = None ):
        """The potentials due to a fitted electron density at a set of grid points."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef CStatus                cStatus = CStatus_OK
        if fitBases is None: bases = target.qcState.fitBases
        else:                bases = fitBases
        coordinates3 = target.scratch.qcCoordinates3AU
        GaussianBasisContainerIntegrals_f1Cp1V ( bases.cObject           ,
                                                 NULL                    ,
                                                 NULL                    ,
                                                 coordinates3.cObject    ,
                                                 gridPoints.cObject      ,
                                                 NULL                    ,
                                                 fitCoefficients.cObject ,
                                                 potentials.cObject      ,
                                                 &cStatus                )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating one-basis electron potentials." )

    def f1Df1i ( self, target, Vector3 center = None ):
        """The dipole integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef SymmetricMatrix        dX
        cdef SymmetricMatrix        dY
        cdef SymmetricMatrix        dZ
        cdef CInteger               n
        cdef CRealArray1D          *cCenter = NULL
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        if center is not None: cCenter = center.cObject
        n  = len ( bases )
        dX = SymmetricMatrix.WithExtent ( n )
        dY = SymmetricMatrix.WithExtent ( n )
        dZ = SymmetricMatrix.WithExtent ( n )
        GaussianBasisContainerIntegrals_f1Df1i ( bases.cObject        ,
                                                 coordinates3.cObject ,
                                                 cCenter              ,
                                                 dX.cObject           ,
                                                 dY.cObject           ,
                                                 dZ.cObject           ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating dipole integrals." )
        return ( dX, dY, dZ )

    def f1Di ( self, target, fitBases = None, Vector3 center = None ):
        """The dipole integrals for a fit basis."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef RealArray1D            dX
        cdef RealArray1D            dY
        cdef RealArray1D            dZ
        cdef CInteger               n
        cdef CRealArray1D          *cCenter = NULL
        cdef CStatus                cStatus = CStatus_OK
        if fitBases is None: bases = target.qcState.fitBases
        else:                bases = fitBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        if center is not None: cCenter = center.cObject
        n  = len ( bases )
        dX = RealArray1D.WithExtent ( n )
        dY = RealArray1D.WithExtent ( n )
        dZ = RealArray1D.WithExtent ( n )
        GaussianBasisContainerIntegrals_f1Di ( bases.cObject        ,
                                               coordinates3.cObject ,
                                               cCenter              ,
                                               dX.cObject           ,
                                               dY.cObject           ,
                                               dZ.cObject           ,
                                               &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating fit dipole integrals." )
        return ( dX, dY, dZ )

    def f1KOf1i ( self, target ):
        """The kinetic energy and overlap integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef SymmetricMatrix        oneElectronMatrix
        cdef SymmetricMatrix        overlap
        cdef CStatus                cStatus = CStatus_OK
        bases             = target.qcState.orbitalBases
        scratch           = target.scratch
        coordinates3      = scratch.qcCoordinates3AU
        oneElectronMatrix = scratch.oneElectronMatrix
        overlap           = scratch.overlapMatrix
        GaussianBasisContainerIntegrals_f1KOf1i ( bases.cObject             ,
                                                  coordinates3.cObject      ,
                                                  oneElectronMatrix.cObject ,
                                                  overlap.cObject           ,
                                                  &cStatus                  )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating kinetic energy and overlap integrals." )

    def f1KOf1R1 ( self, target ):
        """The kinetic energy and overlap gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer bases
        cdef SymmetricMatrix        dTotal
        cdef SymmetricMatrix        wTotal
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            charges      = target.qcState.nuclearCharges
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            wTotal       = scratch.weightedDensity
            GaussianBasisContainerIntegrals_f1KOf1R1 ( bases.cObject        ,
                                                       coordinates3.cObject ,
                                                       dTotal.cObject       ,
                                                       wTotal.cObject       ,
                                                       gradients3.cObject   ,
                                                       &cStatus             )
            if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating kinetic energy and overlap gradients." )

    def f1Of1i ( self, target ):
        """The overlap integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef SymmetricMatrix        overlap
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        overlap      = scratch.overlapMatrix
        GaussianBasisContainerIntegrals_f1Of1i ( bases.cObject        ,
                                                 coordinates3.cObject ,
                                                 overlap.cObject      ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating overlap integrals." )

    def f1Oi ( self, target, fitBases = None ):
        """The overlap integrals for a fit basis."""
        cdef GaussianBasisContainer bases
        cdef RealArray1D            overlap
        cdef CStatus                cStatus = CStatus_OK
        if fitBases is None: bases = target.qcState.fitBases
        else:                bases = fitBases
        overlap = RealArray1D.WithExtent ( len ( bases ) )
        GaussianBasisContainerIntegrals_f1Oi ( bases.cObject   ,
                                               overlap.cObject ,
                                               &cStatus        )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating fit overlap integrals." )
        return overlap

    def f1Op1i ( self, target, Coordinates3 gridPoints not None ,
                                   RealArray2D  values     not None ):
        """The values of the basis functions at a set of grid points."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        coordinates3 = target.scratch.qcCoordinates3AU
        GaussianBasisContainerIntegrals_f1Op1i ( bases.cObject        ,
                                                 coordinates3.cObject ,
                                                 gridPoints.cObject   ,
                                                 values.cObject       ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating basis function grid point values." )

    def f1Qf1i ( self, target, Vector3 center = None ):
        """The quadrupole integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef SymmetricMatrix        qXX
        cdef SymmetricMatrix        qXY
        cdef SymmetricMatrix        qXZ
        cdef SymmetricMatrix        qYY
        cdef SymmetricMatrix        qYZ
        cdef SymmetricMatrix        qZZ
        cdef CInteger               n
        cdef CRealArray1D          *cCenter = NULL
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        if center is not None: cCenter = center.cObject
        n   = len ( bases )
        qXX = SymmetricMatrix.WithExtent ( n )
        qXY = SymmetricMatrix.WithExtent ( n )
        qXZ = SymmetricMatrix.WithExtent ( n )
        qYY = SymmetricMatrix.WithExtent ( n )
        qYZ = SymmetricMatrix.WithExtent ( n )
        qZZ = SymmetricMatrix.WithExtent ( n )
        GaussianBasisContainerIntegrals_f1Qf1i ( bases.cObject        ,
                                                 coordinates3.cObject ,
                                                 cCenter              ,
                                                 qXX.cObject          ,
                                                 qYY.cObject          ,
                                                 qZZ.cObject          ,
                                                 qXY.cObject          ,
                                                 qXZ.cObject          ,
                                                 qYZ.cObject          ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating quadrupole integrals." )
        return ( qXX, qYY, qZZ, qXY, qXZ, qYZ )

    def f1Qi ( self, target, fitBases = None, Vector3 center = None ):
        """The quadrupole integrals for a fit basis."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef RealArray1D            qXX
        cdef RealArray1D            qXY
        cdef RealArray1D            qXZ
        cdef RealArray1D            qYY
        cdef RealArray1D            qYZ
        cdef RealArray1D            qZZ
        cdef CInteger               n
        cdef CRealArray1D          *cCenter = NULL
        cdef CStatus                cStatus = CStatus_OK
        if fitBases is None: bases = target.qcState.fitBases
        else:                bases = fitBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        if center is not None: cCenter = center.cObject
        n   = len ( bases )
        qXX = RealArray1D.WithExtent ( n )
        qXY = RealArray1D.WithExtent ( n )
        qXZ = RealArray1D.WithExtent ( n )
        qYY = RealArray1D.WithExtent ( n )
        qYZ = RealArray1D.WithExtent ( n )
        qZZ = RealArray1D.WithExtent ( n )
        GaussianBasisContainerIntegrals_f1Qi ( bases.cObject        ,
                                               coordinates3.cObject ,
                                               cCenter              ,
                                               qXX.cObject          ,
                                               qYY.cObject          ,
                                               qZZ.cObject          ,
                                               qXY.cObject          ,
                                               qXZ.cObject          ,
                                               qYZ.cObject          ,
                                               &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating fit quadrupole integrals." )
        return ( qXX, qYY, qZZ, qXY, qXZ, qYZ )

    def f1Xf1i_f1Oi  ( self, target, attribute = "fitMatrix", fitBases = None, operator = GaussianBasisOperator.Coulomb, withConstraints = True ):
        """The fit-fit and fit self-overlap integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer fBases
        cdef RealArray1D            fso
        cdef SymmetricMatrix        eri
        cdef SymmetricMatrix        inverse
        cdef CStatus                cStatus = CStatus_OK
        if fitBases is None: fBases = target.qcState.fitBases
        else:                fBases = fitBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        # . Create the fit matrix (n+1)x(n+1).
        n   = len ( fBases )
        eri = Array.WithExtent ( n+1, storageType = StorageType.Symmetric )
        eri.Set ( 0.0 )
        # . Fit-fit integrals.
        if operator is GaussianBasisOperator.AntiCoulomb:
            GaussianBasisContainerIntegrals_f1Af1i ( fBases.cObject       ,
                                                     coordinates3.cObject ,
                                                     eri.cObject          ,
                                                     &cStatus             )
        elif operator is GaussianBasisOperator.Coulomb:
            GaussianBasisContainerIntegrals_f1Cf1i ( fBases.cObject       ,
                                                     coordinates3.cObject ,
                                                     eri.cObject          ,
                                                     &cStatus             )
        else:
            GaussianBasisContainerIntegrals_f1Of1i ( fBases.cObject       ,
                                                     coordinates3.cObject ,
                                                     eri.cObject          ,
                                                     &cStatus             )
        # . With constraints.
        if withConstraints:
            fso = Array.WithExtent ( n )
            fso.Set ( 0.0 )
            GaussianBasisContainerIntegrals_f1Oi ( fBases.cObject ,
                                                   fso.cObject    ,
                                                   &cStatus       )
            for i in range ( n ): eri[i,n] = fso[i]
        # . Finish up.
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating fit-fit integrals." )
        scratch.Set ( attribute, eri )

    def f1Xf1R1 ( self, target, attributeA = "fitCoefficients", attributeX = "fitGradientVectorM", fitBases = None, operator = GaussianBasisOperator.Coulomb ):
        """The fit-fit gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer fBases
        cdef RealArray1D            aVector
        cdef RealArray1D            xVector
        cdef SymmetricMatrix        dTotal
        cdef CGaussianBasisOperator cOperator
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            if fitBases is None: fBases = target.qcState.fitBases
            else:                fBases = fitBases
            cOperator    = operator.value
            coordinates3 = scratch.qcCoordinates3AU
            aVector      = scratch.Get ( attributeA, None )
            gradients3   = scratch.qcGradients3AU
            xVector      = scratch.Get ( attributeX, None )
            GaussianBasisContainerIntegrals_f1Xf1R1 ( fBases.cObject       ,
                                                      coordinates3.cObject ,
                                                      aVector.cObject      ,
                                                      xVector.cObject      ,
                                                      cOperator            ,
                                                      gradients3.cObject   ,
                                                      &cStatus             )
            if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating fit-fit gradients." )

    def f1Xg2i ( self, target, attribute = "fitIntegrals", fitBases = None, operator = GaussianBasisOperator.Coulomb, reportTag = "Fit" ):
        """The electron-fit integrals."""
        cdef BlockStorage           fitIntegrals
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer fBases
        cdef GaussianBasisContainer oBases
        cdef CGaussianBasisOperator cOperator
        cdef CStatus                cStatus = CStatus_OK
        if fitBases is None: fBases = target.qcState.fitBases
        else:                fBases = fitBases
        cOperator    = operator.value
        oBases       = target.qcState.orbitalBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        fitIntegrals = scratch.Get ( attribute, None )
        if fitIntegrals is None:
            fitIntegrals = BlockStorage ( )
            scratch.Set ( attribute, fitIntegrals )
        else:
            fitIntegrals.Empty ( )
        GaussianBasisContainerIntegrals_f1Xg2i ( oBases.cObject       ,
                                                 fBases.cObject       ,
                                                 coordinates3.cObject ,
                                                 cOperator            ,
                                                 fitIntegrals.cObject ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating fit integrals." )
        f        = float ( len ( fBases ) )
        o        = float ( len ( oBases ) )
        n        = fitIntegrals.count
        ( s, m ) = fitIntegrals.byteSize
        report = scratch.qcEnergyReport
        tag    = reportTag + " Integral"
        report[tag + "s"            ] = n
        report[tag + " Sparsity (%)"] = ( ( 1.0 - float ( n ) / ( f * ( o * ( o + 1.0 ) ) / 2.0 ) ) * 100.0, "{:.1f}" )
        report[tag + " Storage ({:s}B)".format ( m.symbol )] =  ( s, "{:.3f}" )

    def f1Xg2R1 ( self, target, attribute = "fitGradientVectorB", fitBases = None, operator = GaussianBasisOperator.Coulomb ):
        """The electron-fit gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer fBases
        cdef GaussianBasisContainer oBases
        cdef RealArray1D            xVector
        cdef SymmetricMatrix        dTotal
        cdef CGaussianBasisOperator cOperator
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            if fitBases is None: fBases = target.qcState.fitBases
            else:                fBases = fitBases
            cOperator    = operator.value
            oBases       = target.qcState.orbitalBases
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            xVector      = scratch.Get ( attribute, None ) # . fitGradientVectorB or qcmmFitGradientVectorB.
            GaussianBasisContainerIntegrals_f1Xg2R1 ( oBases.cObject       ,
                                                      fBases.cObject       ,
                                                      coordinates3.cObject ,
                                                      dTotal.cObject       ,
                                                      xVector.cObject      ,
                                                      cOperator            ,
                                                      gradients3.cObject   ,
                                                      &cStatus             )
            if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating electron-fit gradients." )

    def f2Cf2i ( self, target ):
        """The two-electron integrals."""
        cdef BlockStorage           teis
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        teis         = scratch.Get ( "twoElectronIntegrals", None )
        if teis is None:
            teis = BlockStorage ( )
            scratch.twoElectronIntegrals = teis
        else:
            teis.Empty ( )
        GaussianBasisContainerIntegrals_f2Cf2i ( bases.cObject        ,
                                                 coordinates3.cObject ,
                                                 teis.cObject         ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating two-electron integrals." )
        n        = float ( len ( bases ) )
        p        = ( n * ( n + 1.0 ) ) / 2.0
        n        = teis.count
        ( s, m ) = teis.byteSize
        report   = scratch.qcEnergyReport
        report["TEI Number"      ] = n
        report["TEI Sparsity (%)"] = ( ( 1.0 - float ( n ) / ( ( p * ( p + 1.0 ) ) / 2.0 ) ) * 100.0, "{:.1f}" )
        report["TEI Storage ({:s}B)".format ( m.symbol )] = ( s, "{:.3f}" )

    def f2Cf2R1 ( self, target, doCoulomb = True, CReal exchangeScaling = 1.0 ):
        """The two-electron gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer bases
        cdef SymmetricMatrix        dSpin
        cdef SymmetricMatrix        dTotal
        cdef CBoolean               cDoCoulomb
        cdef CStatus                cStatus = CStatus_OK
        cdef CSymmetricMatrix      *cDSpin  = NULL
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            scratch      = target.scratch
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            if hasattr ( scratch, "onePDMQ" ):
                dSpin  = scratch.onePDMQ.density
                cDSpin = dSpin.cObject
            if doCoulomb: cDoCoulomb = CTrue
            else:         cDoCoulomb = CFalse
            GaussianBasisContainerIntegrals_f2Cf2R1 ( bases.cObject        ,
                                                      coordinates3.cObject ,
                                                      dTotal.cObject       ,
                                                      cDSpin               ,
                                                      cDoCoulomb           ,
                                                      exchangeScaling      ,
                                                      gradients3.cObject   ,
                                                      &cStatus             )
            if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating two-electron gradients." )

    def f2Cm1R1 ( self, target ):
        """The electron-nuclear gradients."""
        cdef Coordinates3           coordinates3
        cdef Coordinates3           gradients3
        cdef GaussianBasisContainer bases
        cdef RealArray1D            charges
        cdef SymmetricMatrix        dTotal
        cdef CStatus                cStatus = CStatus_OK
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            charges      = target.qcState.nuclearCharges
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            GaussianBasisContainerIntegrals_f2Cm1R1 ( bases.cObject        ,
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
            if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating electron-nuclear gradients." )

    def f2Cm1V ( self, target ):
        """The electron-nuclear integrals."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef RealArray1D            charges
        cdef SymmetricMatrix        oneElectronMatrix
        cdef CStatus                cStatus = CStatus_OK
        bases             = target.qcState.orbitalBases
        charges           = target.qcState.nuclearCharges
        scratch           = target.scratch
        coordinates3      = scratch.qcCoordinates3AU
        oneElectronMatrix = scratch.oneElectronMatrix
        GaussianBasisContainerIntegrals_f2Cm1V ( bases.cObject             ,
                                                 charges.cObject           ,
                                                 NULL                      ,
                                                 NULL                      ,
                                                 coordinates3.cObject      ,
                                                 coordinates3.cObject      ,
                                                 NULL                      ,
                                                 oneElectronMatrix.cObject ,
                                                 &cStatus                  )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating electron-nuclear integrals." )

    def f2Cp1V ( self, target, SymmetricMatrix density    not None ,
                               Coordinates3    gridPoints not None ,
                               RealArray1D     potentials not None ):
        """The potentials due to an electron density at a set of grid points."""
        cdef Coordinates3           coordinates3
        cdef GaussianBasisContainer bases
        cdef CStatus                cStatus = CStatus_OK
        bases        = target.qcState.orbitalBases
        coordinates3 = target.scratch.qcCoordinates3AU
        GaussianBasisContainerIntegrals_f2Cp1V ( bases.cObject        ,
                                                 NULL                 ,
                                                 NULL                 ,
                                                 coordinates3.cObject ,
                                                 gridPoints.cObject   ,
                                                 NULL                 ,
                                                 density.cObject      ,
                                                 potentials.cObject   ,
                                                 &cStatus             )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error calculating two-basis electron potentials." )

    def m1Cm1ER1 ( self, target ):
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
        energy = GaussianBasisContainer_m1Cn1ER1 ( charges.cObject      ,
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

    def m1Cp1V ( self, target, Coordinates3 gridPoints not None ,
                                          RealArray1D  potentials not None ):
        """The potentials due to the nuclei at a set of grid points."""
        cdef Coordinates3 coordinates3
        cdef RealArray1D  charges
        charges      = target.qcState.nuclearCharges
        coordinates3 = target.scratch.qcCoordinates3AU
        GaussianBasisContainer_m1Cp1V ( charges.cObject      ,
                                        gridPoints.cObject   ,
                                        coordinates3.cObject ,
                                        NULL                 ,
                                        NULL                 ,
                                        NULL                 ,
                                        NULL                 ,
                                        NULL                 ,
                                        NULL                 ,
                                        potentials.cObject   )          

