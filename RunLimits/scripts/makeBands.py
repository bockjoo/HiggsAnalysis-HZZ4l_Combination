# ported from https://alachua.ihepa.ufl.edu:9999/tree/raid8/bockjoo/optimize/remote/Run2/Submit/limits/HZZ4l_Combination/CreateDatacards/CMSSW_6_1_1/src/HCG/makeBands_hzz4l.cc

from ROOT import *
import ROOT
import os
#%jsroot on
#ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
#ROOT.gSystem.AddIncludePath("-Iinclude/")
ROOT.gSystem.Load("/cvmfs/sft.cern.ch/lcg/releases/LCG_85swan3/vdt/0.3.6/x86_64-slc6-gcc49-opt/lib/libvdt.so")
#ROOT.gSystem.Load("../../lib/libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.Load("/softraid/bockjoo/combine_LCG85_swan3/HiggsAnalysis/CombinedLimit/lib/libHiggsAnalysisCombinedLimit.so")
#ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
#ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")
ROOT.gSystem.Load("/softraid/bockjoo/combine_LCG85_swan3/HiggsAnalysis/CombinedLimit/HiggsCSandWidth/include/HiggsCSandWidth_cc.so")
ROOT.gSystem.Load("/softraid/bockjoo/combine_LCG85_swan3/HiggsAnalysis/CombinedLimit/HiggsCSandWidth/include/HiggsCSandWidthSM4_cc.so")
# Setup python 
# import CMS tdrStyle
from tdrStyle import *
setTDRStyle()

# import some more python modules
import sys,glob
from array import array
import string
from scipy.special import erf
import math
#from bandUtils import *

import pdefs
from plots import *
from bandUtils import *

prefix = ""
nfcbands = 4
fcbandnam = ["68", "95", "99", "9973" ]
fcbandwid = [ 0.68, 0.95, 0.99, 0.9973 ]


def cleanCLs(graph):
    foundBad = False
    while (foundBad) :
        foundBad = False
        for i in xrange(n) : # for (int i = 0, n = graph.GetN(); i < n; ++i):
            if (graph.GetX()[i] == 0): 
                graph.RemovePoint(i)
                foundBad = True
                break
def doIt(bands, name, rootfile):
    #import bandUtils
    print " name " , name , " access rootfile " , rootfile , " result = " , ROOT.gSystem.AccessPathName(rootfile)
    #//if (/*error=*/gSystem.AccessPathName(rootfile)) print " if " << std::endl
    if (ROOT.gSystem.AccessPathName(rootfile)) : return
    pdefs.zero_is_valid = ("pval" in name or "sig" in name or "smcls" in name or "ml_" in name or "mlz_" in name)
    pdefs.do_bands_asimov = "pval" in name
    pdefs.use_precomputed_quantiles = "cls" in name or "ml_" in name or "mlz_" in name
    pdefs.do_bands_95 = not ("ml_" in name or "mlz_" in name)
    #//if ("pval" in name and name.Contains("comb")) do_bands_cputime = True
    makeBands(bands, name, rootfile, 0)
    print "Trying to do " , name
    if name.startswith("ml") : # if (name.index("ml") == 0):
        if (bands.Get(name+"_median") == 0):  print "Missing median??" ; return;  
        obs =  bands.Get(name+"_median").Clone()
        if (obs == 0) : return
        bands.Delete(name+"_obs;*")
        bands.Delete(name+"_median;*")
        obs.SetName(name+"_obs")
        bands.WriteTObject(obs, name+"_obs")
        printLineAErr(bands, name+"_obs", "results/"+prefix+name+".txt")
    else:
        #print "DEBUG results/"+prefix+name+".txt"
        printBand(bands, name, "results/"+prefix+name+".txt", False)
        if ("PLPE" in rootfile):
            printLine(bands, name+"_asimov", "results/"+prefix+name+"_expected.txt", False)
        
    
    pdefs.zero_is_valid = False
    pdefs.use_precomputed_quantiles = False
    pdefs.do_bands_95 = True
    pdefs.do_bands_asimov = False
    pdefs.do_bands_cputime = False; pdefs.do_bands_realtime = False

def doFc(bands, name, rootfile):
    if (ROOT.gSystem.AccessPathName(rootfile)) : return
    rfile = ROOT.TFile.Open(rootfile)
    for i in xrange(nfcbands) : # for (int i = 0; i < nfcbands; ++i):
        print "Make " , name
        fcb = theFcBelt(rfile, 1, 0, Observed, fcbandwid[i])
        if (fcb == 0) : continue
        fcb.SetName(name+"_"+fcbandnam[i])
        bands.WriteTObject(fcb)
    
    printFcBand(bands, name, name+".txt", nfcbands, fcbandnam)

def importCloned(bands, oldfile, oldname, nwhat, what,  nwho, who):
    oldBands = ROOT.TFile.Open(oldfile)
    for i in xrange(-2,nwho) : # for (int i = -2; i < nwho; ++i):        
        if i == -1: whoI = "comb"; who0 = "comb"
        elif i == -2: whoI = "hzz2l2q"; who0 = "hzz2l2q"
        else : whoI = who[i]; who0 = who[i]
        
        if (oldname == "lp11" and whoI == "hww") : whoI = "hwwc"
        for j in xrange(nwhat) : # for (int j = 0; j < nwhat; ++j):
            name  ="%s_%s" % (what[j], who0)
            nameI = "%s_%s" % (what[j], whoI)
            oldg_obs    =  oldBands.Get(nameI+"_obs")
            oldg_median =  oldBands.Get(nameI+"_median")
            oldg_medi95 =  oldBands.Get(nameI+"_median_95")
            if (oldg_obs):
                cleanCLs(oldg_obs)
                bands.WriteTObject(oldg_obs.Clone(name+"_"+oldname+"_obs"), name+"_"+oldname+"_obs")
                print "Imported " , oldname , " result " , name , "_obs" 
            
            if (oldg_median):
                cleanCLs(oldg_median)
                bands.WriteTObject(oldg_median.Clone(name+"_"+oldname+"_median"), name+"_"+oldname+"_median")
                print "Imported " , oldname , " result " , name << "_median"
            
            if (oldg_medi95):
                cleanCLs(oldg_medi95)
                bands.WriteTObject(oldg_medi95.Clone(name+"_"+oldname+"_median_95"), name+"_"+oldname+"_median_95")
                print "Imported " , oldname , " result " , name , "_median_95"
            
        
    
    oldBands.Close()

def importLandS(bands):
    importLandS(bands, "mlz_comb_lands",    "lands/comb_maxllfit.root", 1, 0)
    importLandS(bands, "mlz_combe_lands",   "lands/comb_maxllfit_ebereso.root", 1, 0)
    importLandS(bands, "pvala_comb_lands",  "lands/comb_pvalue.root",   1, 0)
    importLandS(bands, "pvala_comb_lands",  "lands/comb_pvalueSB.root", 0, 1)
    importLandS(bands, "pvala_combe_lands", "lands/comb_pvalue_ebereso.root", 1, 0)

def makebands( which=0):
    global prefix
    #gROOT.LoadMacro("$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/bandUtils.cxx+")

    if which == 1: prefix = "SM4_"
    elif which == 2: prefix = "FF_"
    

    bands = ROOT.TFile(("results/%sbands.root" % prefix),"RECREATE")

    print "Will save bands to " , bands.GetName()
    #from bandUtils import ObsAvgMode
    #ObsAvgMode = Enum(['MeanObs', 'LogMeanObs', 'MedianObs'])
    obs_avg_mode = ObsAvgMode.MedianObs
    pdefs.do_bands_nosyst = False
    pdefs.do_bands_mean   = False
    pdefs.do_bands_ntoys  = False
    pdefs.do_bands_asimov = False
    pdefs.halfint_masses = True

    nwhat = 13
    what = [ "pla", "pvala", "pvala", "pval", "acls" , "ml", "mlz", "cls",  "bayes", "smcls", "smacls", "acls90", "acls99"]
    WHAT = [ "PLC", "PLP",   "PLPE",  "PVAL", "ASCLS", "ML", "MLZ", "FREQ", "BAYES", "SMCLS", "SMASCLS", "ASCLS90", "ASCLS99"]

    
    nwho = 20
    who = [ "comb", "htt", "vhtt", "httm", "htt0", "hzz", "hzz2l2t", "hzz2l2nu", "hzz2l2q_all", "hzz4l", "hww", "vhww3l", "hww2l", "hgg_novbf", "hgg_vbf", "hgg", "vhbb", "combs", "combp_low", "combl_low" ]
    WHO = [ "COMB", "HTT", "VHTT", "HTTM", "HTT0", "HZZ", "HZZ2L2T", "HZZ2L2NU", "HZZ2L2Q",     "HZZ4L", "HWW", "VHWW3L", "HWW2L", "HGG_NOVBF", "HGG_VBF", "HGG", "VHBB", "COMBS", "COMBP", "COMBL" ]

    for i in xrange(nwho) : # for (int i = 0; i < nwho; ++i):
        #//if (!TString(who[i]).Contains("comb")) continue
        if (not "hzz4l" in who[i]) : continue
        #//doFc(bands, TString::Format("fc_%s",who[i]), TString::Format("higgsCombine%s_FC.root", WHO[i]))
        #//continue
        for j in xrange(nwhat) : # for (int j = 0; j < nwhat; ++j):
            if ( what[j] != "acls" and what[j] != "pvala" ) : continue
            print "doIt what = " , what[j] , " who = " , who[i] 
            doIt(bands, ("%s_%s" % (what[j], who[i])),
                        ("results/higgsCombine%s%s_%s.root" % (prefix, WHO[i], WHAT[j])))
        
        #//break
        #if 0
        #if (not ROOT.gSystem.AccessPathName(("timecls_%s.txt" %who[i]))):
        #    time = ROOT.TGraphErrors(("timecls_%s.txt" % who[i]), "Mass: %lg Time: %lg +/- %lg s")
        #    time.SetName(("timecls_%s_obs" % who[i]))
        #    bands.WriteTObject(time)
        #
        #endif
    

    #if 1 == 0 : ##if 0
    #const int npvflav = 9
    #pvflav[npvflav] = { "M2M", "S", "SZ3", "BZ1", "BZ", "TQ0", "TQ03", "TQ05N", "TQ05" 
    #for (int i = 0; i < npvflav; ++i):
    #    doIt(bands, TString::Format("pvala_%s_comb", pvflav[i]),
    #                TString::Format("results/higgsCombine%sCOMB_PLP.%s.root", prefix.Data(), pvflav[i]))
    
    ##endif
    ##if 0
    #if (which == 0):
    #    importLine(bands, "tevatron_obs",    "tevatron.obs.txt");    
    #    importLine(bands, "tevatron_median", "tevatron.exp.txt");    
    #    std::cout  << "Imported Tevatron results" << std::endl
    #    importCloned(bands, "results/bands.paper.root", "paper", nwhat, what, nwho, who)
    
    ##endif

    ##if 0
    #for (int j = 0; j < nwhat; ++j):
    #    w = what[j]
    #    cutBands(bands, w+"_hzz2l2q_all", w+"_hzz2l2q_low", 130, 164)
    #    cutBands(bands, w+"_hzz2l2q_all", w+"_hzz2l2q",     200, 600)
    #    cutBands(bands, w+"_hzz", w+"_combp_high", 150.5, 600)
    #    cutBands(bands, w+"_hww", w+"_combl_high", 145.5, 600)
    #    pasteBands(bands, w+"_combl_low", w+"_combl_high", w+"_combl")
    #    pasteBands(bands, w+"_combp_low", w+"_combp_high", w+"_combp")
    
    ##endif
    if (which == 0):
        cutBands(bands, "smcls_comb",       "smcls_comb_new", 110, 200)
        cutBands(bands, "smcls_comb_paper", "smcls_comb_old", 200, 600)
        pasteBands(bands, "smcls_comb_new", "smcls_comb_old", "smcls_comb_patch")
        cutBands(bands, "cls_comb",       "cls_comb_new", 110, 200)
        cutBands(bands, "cls_comb_paper", "cls_comb_old", 200, 600)
        pasteBands(bands, "cls_comb_new", "cls_comb_old", "cls_comb_patch")

#root -b -l -q makeBands.cxx
makebands()
