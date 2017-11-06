source /softraid/bockjoo/combine_LCG85_swan3/HiggsAnalysis/CombinedLimit/env_standalone.sh
python plotBand.py 2>&1 | tee plotBand.log

plotBand.py <- plots.cxx
plots.py <- plots.cxx
makeBands.py <- makeBands.cxx
bandUtils.py <- bandUtils.cxx
pdefs: the global variables that need to be imported by plotBand.py, plots.py, bandUtils.py, and makeBands.py
    