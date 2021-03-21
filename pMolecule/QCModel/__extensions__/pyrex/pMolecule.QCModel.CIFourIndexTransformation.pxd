from pCore.CPrimitiveTypes                    cimport CBoolean               , \
                                                      CInteger               , \
                                                      CFalse                 , \
                                                      CTrue                  , \
                                                      CInteger               , \
                                                      CReal              
from pCore.Status                             cimport CStatus, CStatus_OK
from pMolecule.QCModel.BlockStorage           cimport CBlockStorage          , \
                                                      BlockStorage        
from pScientific.Arrays.DoubleSymmetricMatrix cimport CDoubleSymmetricMatrix , \
                                                      DoubleSymmetricMatrix        
from pScientific.Arrays.RealArray2D           cimport CRealArray2D           , \
                                                      RealArray2D     
from pScientific.Arrays.RealArrayND           cimport CRealArrayND           , \
                                                      RealArrayND     

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CIFourIndexTransformation.h":

    cdef void CIFourIndexTransformation ( CRealArray2D           *activeMOs            ,
                                          CBlockStorage          *twoElectronIntegrals ,
                                          CRealArray2D           *moTEI34              ,
                                          CRealArrayND           *moTEI234             ,
                                          CDoubleSymmetricMatrix *moTEIs               ) ;
