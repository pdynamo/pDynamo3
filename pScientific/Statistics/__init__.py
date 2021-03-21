"""A package for various statistical analyses."""

from .Correlation import Correlation_AutoCorrelation            , \
                         Correlation_CrossCorrelation           , \
                         Correlation_DotProductAutoCorrelation  , \
                         Correlation_DotProductCrossCorrelation
from .Histogram   import RegularHistogram                       , \
                         RegularHistogramBinMidPointIterator    , \
                         RegularHistogramDimension
from .Statistics  import Statistics                             , \
                         StatisticsAccumulator
