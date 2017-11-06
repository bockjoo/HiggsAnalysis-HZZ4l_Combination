#!/usr/env/bin python
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
import pdefs

# Testing 
#print "plot definitions ", pdefs.pdefs
#pdef.pdef = "it is indeed"
#print "plot definitinos after update ", pdefs.pdefs


#import plots
from plots import *
from bandUtils import *
which = 0
#def plots(which):
if 1 == 1:
    #print 1
    #gROOT.LoadMacro("$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/bandUtils.cxx+")
    #ROOT.gSystem.Load("/cvmfs/sft.cern.ch/lcg/releases/LCG_85swan3/vdt/0.3.6/x86_64-slc6-gcc49-opt/lib/libvdt.so")
    #os.environ['LD_LIBRARY_PATH']='/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/lib64:/cvmfs/sft.cern.ch/lcg/views/LCG_85swan3/x86_64-slc6-gcc49-opt/lib:/cvmfs/sft.cern.ch/lcg/views/LCG_85swan3/x86_64-slc6-gcc49-opt/lib64:$LD_LIBRARY_PATH'
    #os.environ['LD_LIBRARY_PATH']='/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib64:/cvmfs/sft.cern.ch/lcg/views/LCG_85swan3/x86_64-slc6-gcc49-opt/lib:/cvmfs/sft.cern.ch/lcg/views/LCG_85swan3/x86_64-slc6-gcc49-opt/lib64:$LD_LIBRARY_PATH'
    #ROOT.gSystem.Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/gcc/4.9.3/lib64/libstdc++.so")
    #ROOT.gROOT.LoadMacro("../HiggsAnalysis/CombinedLimit/test/plotting/bandUtils.cxx+")
    #theBand()
    #print 2
    inputPrefix = ""
    #switch (which)  {
    if which == 1:
        inputPrefix = "SM4_" ; pdefs.SM="SM4"; 
        pdefs.LEESPAM_2L = "#splitline{Global p-value ~0.5}{for 110-600 GeV range}"
        #        break
    elif which == 2:
        pdefs.inputPrefix = "FF_"  ; pdefs.SM="FP";  
        pdefs.LEESPAM_2L = "#splitline{Global significance 1.1#sigma}{for 110-300 GeV range}"
        #        break
    #}
    #print 4

    loadMasses("masses.txt")
    #print 5
    globalPrefix0 = inputPrefix+"plots/"
    pdefs.globalPrefix =  globalPrefix0+"/"
    ROOT.gStyle.SetOptTitle(0)
    print "cp results/"+inputPrefix+"bands.root "+inputPrefix+"bands.work.root"
    ROOT.gSystem.Exec("/bin/cp results/"+inputPrefix+"bands.root "+inputPrefix+"bands.work.root")
    #print 6
  
    #//TFile::Open(inputPrefix+"bands.work.root", "UPDATE")
    workband = ROOT.TFile.Open(inputPrefix+"bands.work.root", "UPDATE")
    #print 7

    #print inputPrefix
    #// WRITE OUT THE GRIDS TO RUN CLS AND BAYESIAN
    #////writeGrids(inputPrefix) return; 
    writeGrid("hzz4l","")
    #workband.Close()

    #workband.ls()

    #if not workband.Get("pval_comb"):
    #    print "NOTOK"
    #else:
    #    print "OK"
    #if workband.Get("pval_comb") is None:
    #    print "nill"
    #else:
    #    print "OK"
    ##workband.Get("pvala_coms")

    #from bandUtils import *
    #////return 
    #switch (which):
    #print 7.1
    if which == 0:
        print "gFile = " , ROOT.gFile.GetName() #, " gFile.GetN() ", ROOT.gFile.GetN()
        #//*test = (*) gFile.Get("acls_hzz4l_obs")
        #//bandIn="pval_comb"
        #//bandOut="pval_comb125"
        #//selectedPointsBands("pval_comb")
        selectedPointsBands(ROOT.gFile, "pval_comb", "pval_comb125", 125.0)
        #selectedPointsBands(workband, "pval_comb", "pval_comb125", 125.0)
        #//selectedPointsBands(gFile, "pval_comb", "pval_comb125", 125.0)
        #//selectedPointsBands(gFile, "pval_comb", "pval_comb12x", 125.0, 124.5, 124.0, 125.5, 126.0)
        selectedPointsBands(ROOT.gFile, "pval_comb", "pval_comb12x", 125.0, 124.5, 124.0, 125.5, 126.0)
        #selectedPointsBands(workband, "pval_comb", "pval_comb12x", 125.0, 124.5, 124.0, 125.5, 126.0)
        #//selectedPointsBands(gFile, "pval_comb", "pval_comb12x", 125.0)
        #break
    elif which == 1:
        selectedPointsBands(ROOT.gFile, "pval_comb", "pval_comb125", 119.0)
        selectedPointsBands(ROOT.gFile, "pval_comb", "pval_comb12x", 119.0, 120.0, 118.0, 121.0, 117.0)
        #break
    elif which == 2:
        selectedPointsBands(ROOT.gFile, "pval_comb", "pval_comb125", 126.0)
        selectedPointsBands(ROOT.gFile, "pval_comb", "pval_comb12x", 126.0, 124.0, 125.0, 127.0, 128.0)
        #break
    
    print 8
    #pdefs.c1 = ROOT.TCanvas("c1","c1")
    rectangleCanvas()
    #print 9

    #//#if 0
    #// PAPER PLOTS, keep it simple
    #squareCanvas(/*grid x,y=*/0,0); x_zoom = 145;
    squareCanvas(0,0) ; pdefs.x_zoom = 145;
    pdefs.globalPrefix = globalPrefix0+"/pre-paper/sqr_"
    pdefs.noLineStyles = True #// always use solid lines

    #// fig 1
    #//drawSMCLs(TString::Format("smcls_%s", "comb"))

    #// fig 2
    #//drawOneCLs(TString::Format("cls_%s",  "comb"))

    #// fig 3
    #//bockjoo drawCombObs("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann, nchann)
    pdefs.lineAt1Style = 2
    print "chann ", pdefs.chann, " nchann ",pdefs.nchann
    #drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+SM+"}", chann, nchann, /*obs=*/False, /*combined=*/False)
    drawCombBoth("acls", "95% CL limit on #sigma/#sigma_{"+pdefs.SM+"}", pdefs.chann, pdefs.nchann, False, False)
    pdefs.lineAt1Style = 1
    #// fig 4
    pdefs.forceYmin = 0.08 ; pdefs.forceYmax = 12;
    pdefs.channelSpamOnRightHandSide = True
    pdefs.CLs_debug_apriori_grid = True
    drawOneCLs("acls_%s" %  "hzz4l")
    #//drawOneCLs(TString::Format("acls_%s",  "combl_full"))
    pdefs.channelSpamOnRightHandSide = False
    pdefs.forceYmin = 0 ; pdefs.forceYmax = 0;
