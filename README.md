#Script to run on files produced using the BristolAnalysisTools

#Produce a cutflow from one sample:
root cutFlowOneSamp.C+

#Produce a cutflow for more than one sample:
root cutFlow.C+

#Scripts for binning:
root do2DMET.C+   #For MET binning is already defined from last year
root do2DPlots.C+  #Can swap to find binning or use predefined

#Make various control plots:
root doPlots(Obj).C+

