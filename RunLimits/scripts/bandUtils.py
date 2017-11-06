
# coding: utf-8

# In[1]:


# ported from https://alachua.ihepa.ufl.edu:9999/tree/raid8/bockjoo/optimize/remote/Run2/Submit/limits/HZZ4l_Combination/CreateDatacards/CMSSW_6_1_1/src/HCG/makeBands_hzz4l.cc


# In[2]:


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

# In[3]:


#pdefs.band_safety_crop = 0
#pdefs.use_precomputed_quantiles = False
#pdefs.precomputed_median_only = False
#pdefs.zero_is_valid = False
#pdefs.seed_is_channel = False
#pdefs.halfint_masses  = False #// find the halfling!

class Enum(tuple): __getattr__ = tuple.index
ObsAvgMode = Enum(['MeanObs', 'LogMeanObs', 'MedianObs'])    
obs_avg_mode = ObsAvgMode.MeanObs
BandType = Enum(['Mean', 'Median', 'Quantile', 'Observed', 'Asimov', 'CountToys', 'MeanCPUTime', 'MeanRealTime', 'AdHoc', 'ObsQuantile']) 


# In[4]:


#// Maritz-Jarrett, JASA vol. 73 (1978)
#// http:#//www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/quantse.htm
def quantErr(n, vals, q):
    m = floor(q*n+0.5)
    a = m-1
    b = n-m
    ninv = 1.0/n
    c1 = 0
    c2 = 0
    last_cdf = 0
    for i in xrange(n): #(size_t i = 0; i < n; ++i):
        this_cdf = ROOT.Math.beta_cdf((i+1) * ninv, a, b)
        w = this_cdf - last_cdf
        c1 += w * vals[i]
        c2 += w * vals[i] * vals[i]
        last_cdf = this_cdf
    
    return sqrt(c2 - c1*c1)


# In[5]:


def polyFit(x0, y0,  npar,  n, xi, yi):
    #//std::cout << "smoothWithPolyFit(x = " << x <<", npar = " << npar << ", n = " << n << ", xi = {" << xi[0] << ", " << xi[1] << ", ...}, yi = {" << yi[0] << ", " << yi[1] << ", ...})" << std::endl
    mat=ROOT.TMatrixDSym (npar)
    vec=ROOT.TVectorD    (npar)
    for j in xrange(npar): #(size_t i = 0; i < n; ++i): #for (int j = 0; j < npar; ++j):
        for j2 in xrange(j,npar): #       for (int j2 = j; j2 < npar; ++j2):
            #cor[i][j]
            mat[j][j2] = 0 #mat(j,j2) = 0
        
        vec[j] = 0 # vec(j) = 0
    
    for i in xrange(n) : # for (int i = 0; i < n; ++i):
        for j in xrange(npar): #for (int j = 0; j < npar; ++j):
            for j2 in xrange(j,npar): # for (int j2 = j; j2 < npar; ++j2):
                mat[j][j2] += pow(xi[i]-x0, j+j2)
            
            vec[j] += (yi[i]-y0)*pow(xi[i]-x0, j)
        
    
    bk = ROOT.TDecompBK (mat)
    bk.Solve(vec)
    return vec


# In[6]:

def theBand(file, doSyst, whichChannel, type, width=0.68):
    if not file : return 0
    t = file.Get("limit")
    if not t : t = file.Get("test") #// backwards compatibility
    if not t: 
        print "TFile " , file.GetName() , " does not contain the tree" 
        return 0 

    dataset = {} # std::map<int,std::vector<double> >  dataset
    errors = {} # std::map<int,std::vector<double> >  errors    
    obsValues = {} # std::map<int,double>                obsValues
    
    n = t.GetEntries()
    #print "DEBUG printing tree "
    #t.Print()
    for i in xrange(n): # for (size_t i = 0, n = t.GetEntries(); i < n; ++i):
        t.GetEntry(i)

        mass = float(t.mh) ; limit = float(t.limit)
        #print "mass = ", mass , " vs t.mh = ",t.mh, " limit = ",limit, " t.limit = ",t.limit
        if t.GetBranch("limitErr") : limitErr = float(t.limitErr)
        if (t.GetBranch("t_cpu") != 0):
            t_cpu = float(t.t_cpu)
            t_real = float(t.t_real)
        if (pdefs.use_precomputed_quantiles):
            if (t.GetBranch("quantileExpected") == 0): 
                print "TFile " , file.GetName() , " does not have precomputed quantiles" 
                #return 0
            quant = float(t.quantileExpected)
        syst = t.syst
        if pdefs.seed_is_channel : iChannel = t.iSeed
        else : iChannel = t.iChannel
        iToy = t.iToy
        
        iMass = int(mass*100)
        #//printf("%6d mh=%.1f  limit=%8.3f +/- %8.3f toy=%5d quant=% .3f\n", i, mass, limit, limitErr, iToy, quant)
        if (syst != doSyst) :           continue
        if (iChannel != whichChannel) : continue
        if      (type == BandType.Asimov) :
            if (iToy != -1) : continue 
        elif (type == BandType.Observed): 
            if (iToy !=  0) : continue
        elif (type == BandType.ObsQuantile and iToy == 0): 
            #obsValues[iMass] = limit
            if iMass in obsValues:
                obsValues[iMass] = limit
            else:
                obsValues[iMass] = []
                obsValues[iMass] = limit
            continue 
        elif (iToy <= 0 and not pdefs.use_precomputed_quantiles) : continue
        if (limit == 0 and not pdefs.zero_is_valid) : continue
        if (type == BandType.MeanCPUTime): 
            if (limit < 0) : continue
            limit = t_cpu # limit = t_cpu
        
        if (type == BandType.MeanRealTime): 
            if (limit < 0) : continue
            limit = t_real
        
        if (pdefs.use_precomputed_quantiles):
            if (type == BandType.CountToys)   :return 0
            if (type == BandType.Mean)        :return 0
            #//std::cout << "Quantiles. What should I do " << (type == Observed ? " obs" : " exp") << std::endl
            if (type == BandType.Observed and quant > 0) :continue
            if (type == BandType.Median):
                if (abs(quant - 0.5) > 0.005 and abs(quant - (1-width)/2) > 0.005 and abs(quant - (1+width)/2) > 0.005):
                    #//std::cout << " don't care about " << quant << std::endl
                    continue
                 #else
                 #   #//std::cout << " will use " << quant << std::endl
                
            
        
        if iMass in dataset:
            dataset[iMass].append(limit)
            errors[iMass].append(limitErr)
        else:
            dataset[iMass] = []
            dataset[iMass].append(limit)
            errors[iMass] = []
            errors[iMass].append(limitErr)
        #print " appending limit to dataset iMass = ",iMass, " limit = ",limit," dataset[iMass] = ", dataset[iMass]
    #print "DEBUG dataset ",len(dataset), dataset
    tge = ROOT.TGraphAsymmErrors () # new () # TGraphAsymmErrors
    ip = 0
    #for value in values:
    #   # Print each row's length and its elements.
    #   print(len(value), value)
    # for (std::map<int,std::vector<double> >::iterator it = dataset.begin(), ed = dataset.end(); it != ed; ++it):
    #for idx, data in enumerate(dataset) :
    #for idx, data in dataset.iteritems():
    #    #data is it.second # std::vector<double> &data = it.second
    #    nd=len(data) # int nd = data.size()
    #    print "DEBUG theBand( nd = ",nd
    #    #data.sort() # std::sort(data.begin(), data.end())
    
    #for idx, data in dataset.iteritems():
    for idx in sorted(dataset):
        data = dataset[idx] #data is it.second # std::vector<double> &data = it.second
        nd=len(data) # int nd = data.size()
        print "DEBUG theBand( idx= ",idx, " nd = ",nd
        for j in xrange(nd): print "DEBUG theBand nd[",j,"]=",data[j]
        data.sort() # std::sort(data.begin(), data.end())
        if nd % 2 == 0 :
            #median = 0.5*(data[nd/2]+data[nd/2+1])
            median = 0.5*(data[nd/2-1]+data[nd/2])
        else :
            median = data[nd/2]
        #print "DEBUG theBand( idx = ",idx," 1 median = ",median
            
        if (pdefs.band_safety_crop > 0):
            data2 = [] # std::vector<double> data2
            for j in xrange(nd): # for (int j = 0; j < nd; ++j):
                if (data[j] > median*pdefs.band_safety_crop and data[j] < median/pdefs.band_safety_crop):
                    data2.append(data[j]) # data2.push_back(data[j])
                
            
            data2.swap(data)
            datatmp=data2, data
            data=datatmp[0] 
            data2=datatmp[1]
            nd = len(data) #data.size()
            #if nd % 2 == 0 : median = 0.5*(data[nd/2]+data[nd/2+1])
            if nd % 2 == 0 : median = 0.5*(data[nd/2-1]+data[nd/2])
            else : median =  data[nd/2]
        print "DEBUG theBand( nd = ",nd," idx = ",idx," 1 median = ",median
        mean = 0
        for j in xrange(nd):
            mean += data[j]
        mean /= nd
        summer68 = data[int(math.floor(nd * 0.5*(1-width)+0.5))]
        winter68 =  data[min(int(math.floor(nd * 0.5*(1+width)+0.5)), nd-1)]
        if (pdefs.use_precomputed_quantiles and type == BandType.Median):
            if (pdefs.precomputed_median_only and len(data) == 1):
                mean = median = summer68 = winter68 = data[0]
            elif (len(data) != 3): 
                print "Error for expected quantile for mass " , idx , ": size of data is " , len(data)
                continue
            else :
                mean = median = data[1]; summer68 = data[0]; winter68 = data[2]
            
        
        x = mean
        #switch (type):
        if type == BandType.MeanCPUTime: pass
        elif type == BandType.MeanRealTime: pass
        elif type == BandType.Mean: x = mean
        elif type == BandType.Median: x = median
        elif type == BandType.CountToys: x = summer68 = winter68 = nd
        elif type == BandType.Asimov: #// mean (in case we did it more than once), with no band
            if obs_avg_mode == mean : winter68 = mean
            else : winter68 = median
            x = summer68 = winter68 # = (obs_avg_mode == mean ? mean : median)
        elif type == BandType.Observed:
            x = mean
            if (nd == 1):
                if (len(errors[idx]) == 1):
                    summer68 = mean - errors[idx][0]
                    winter68 = mean + errors[idx][0]
                else :
                    #// could happen if limitErr is not available
                    summer68 = winter68 = mean

            else : #// if we have multiple, average and report rms (useful e.g. for MCMC)
                #switch (obs_avg_mode):
                if obs_avg_mode == ObsAvgMode.MeanObs:   x = mean
                elif obs_avg_mode ==  ObsAvgMode.MedianObs: x = median
                elif obs_avg_mode ==  ObsAvgMode.LogMeanObs:
                    x = 0.0
                    for j in xrange(nd) : x += log(data[j]) #for (int j = 0; j < nd; ++j): x += log(data[j]); 
                    x = exp(x/nd)

                rms = 0.0
                for j in xrange(nd) : rms += (x-data[j])*(x-data[j]) # for (int j = 0; j < nd; ++j): rms += (x-data[j])*(x-data[j]); 
                rms = math.sqrt(rms/(nd*(nd-1)))
                summer68 = mean - rms
                winter68 = mean + rms

        elif type == BandType.AdHoc:
            x = summer68 = winter68 = mean
        elif type == BandType.Quantile: #// get the quantile equal to width, and it's uncertainty
            x = data[int(math.floor(nd*width+0.5))]
            summer68 = x - quantErr(nd, data, width)
            winter68 = x + (x-summer68)
        elif type == BandType.ObsQuantile:   
            #if ( idx == ( len(obsValues) - 1)) : continue # ????? #if (obsValues.find(it.first) == obsValues.end()) continue
            # we need to check if idx exist in obsValues
            if not idx in obsValues : continue
            passed = 0
            failed = 0
            for i in xrange(nd): # for (int i = 0; i < nd and data[i] <= obsValues[it.first]; ++i):
                if data[i] > obsValues[idx] : continue
                failed += 1

            passed = nd - failed; x = double(passed)/nd
            alpha = (1.0 - .68540158589942957)/2.0
            if passed == 0 : summer68 = 0.0
            else : summer68 = ROOT.Math.beta_quantile(   alpha, passed,   failed+1 )
            if failed == 0 : winter68 = 1.0
            else : winter68 = ROOT.Math.beta_quantile( 1-alpha, passed+1, failed   )
            

        #// end switch
        tge.Set(ip+1)
        tge.SetPoint(ip, float(idx)*0.01, x)
        tge.SetPointError(ip, 0, 0, x-summer68, winter68-x)
        ip += 1
    
    return tge

def theBandOld(file, doSyst, whichChannel, type, width=0.68):
    if not file : return 0
    t = file.Get("limit")
    if not t : t = file.Get("test") #// backwards compatibility
    if not t: 
        print "TFile " , file.GetName() , " does not contain the tree" 
        return 0 
    #Double_t mass, limit, 
    #limitErr = 0
    #; Float_t t_cpu, t_real; Int_t syst, iChannel, iToy, iMass; Float_t 
    #                quant = -1
    #                mass = array('d',[])
    #                limit = array('d',[])
    #                limitErr = array('d',[0])
    #                t_cpu = array('d',[])
    #                t_real = array('d',[])
    #                quant = array('i',[])
    #                syst = array('i',[])
    #                iChannel = array('i',[])
    #                iToy = array('i',[])
    #                t.SetBranchAddress("mh", mass)
    #                t.SetBranchAddress("limit", limit)
    #                if (t.GetBranch("limitErr")) : t.SetBranchAddress("limitErr", limitErr)
    #                if (t.GetBranch("t_cpu") != 0):
    #                    t.SetBranchAddress("t_cpu", t_cpu)
    #                    t.SetBranchAddress("t_real", t_real)#

    #                if (pdefs.use_precomputed_quantiles):
    #                    if (t.GetBranch("quantileExpected") == 0): 
    #                        print "TFile " , file.GetName() , " does not have precomputed quantiles" 
    #                        return 0; 
    #                    t.SetBranchAddress("quantileExpected", quant)

    #                t.SetBranchAddress("syst", syst)
    #                if pdefs.seed_is_channel : t.SetBranchAddress("iSeed", iChannel)
    #                else : t.SetBranchAddress("iChannel", iChannel)
    #                t.SetBranchAddress("iToy", iToy)

    dataset = [[]] # std::map<int,std::vector<double> >  dataset
    errors = [[]] # std::map<int,std::vector<double> >  errors
    #obsValues # std::map<int,double>                obsValues
    obsValues = [[]]
    n = t.GetEntries()
    #print "DEBUG printing tree "
    #t.Print()
    for i in xrange(n): # for (size_t i = 0, n = t.GetEntries(); i < n; ++i):
        t.GetEntry(i)

        mass = t.mh ; limit = t.limit
        if t.GetBranch("limitErr") : limitErr = t.limitErr
        if (t.GetBranch("t_cpu") != 0):
            t_cpu = t.t_cpu
            t_real = t.t_real
        if (pdefs.use_precomputed_quantiles):
            if (t.GetBranch("quantileExpected") == 0): 
                print "TFile " , file.GetName() , " does not have precomputed quantiles" 
                #return 0
            quant = t.quantileExpected
        syst = t.syst
        if pdefs.seed_is_channel : iChannel = t.iSeed
        else : iChannel = t.iChannel
        iToy = t.iToy
        
        iMass = int(mass*100)
        #//printf("%6d mh=%.1f  limit=%8.3f +/- %8.3f toy=%5d quant=% .3f\n", i, mass, limit, limitErr, iToy, quant)
        if (syst != doSyst) :           continue
        if (iChannel != whichChannel) : continue
        if      (type == BandType.Asimov) :
            if (iToy != -1) : continue 
        elif (type == BandType.Observed): 
            if (iToy !=  0) : continue
        elif (type == BandType.ObsQuantile and iToy == 0): 
            obsValues[iMass] = limit
            continue 
        elif (iToy <= 0 and not pdefs.use_precomputed_quantiles) : continue
        if (limit == 0 and not pdefs.zero_is_valid) : continue
        if (type == BandType.MeanCPUTime): 
            if (limit < 0) : continue
            limit = t_cpu # limit = t_cpu
        
        if (type == BandType.MeanRealTime): 
            if (limit < 0) : continue
            limit = t_real
        
        if (pdefs.use_precomputed_quantiles):
            if (type == BandType.CountToys)   :return 0
            if (type == BandType.Mean)        :return 0
            #//std::cout << "Quantiles. What should I do " << (type == Observed ? " obs" : " exp") << std::endl
            if (type == BandType.Observed and quant > 0) :continue
            if (type == BandType.Median):
                if (fabs(quant - 0.5) > 0.005 and fabs(quant - (1-width)/2) > 0.005 and fabs(quant - (1+width)/2) > 0.005):
                    #//std::cout << " don't care about " << quant << std::endl
                    continue
                 #else
                 #   #//std::cout << " will use " << quant << std::endl
                
            
        
        dataset[iMass].append(limit) # .push_back(limit)
        errors[iMass].append(limitErr) # .push_back(limitErr)
    
    tge = ROOT.TGraphAsymmErrors () # new () # TGraphAsymmErrors
    ip = 0
    #for value in values:
    #   # Print each row's length and its elements.
    #   print(len(value), value)
    # for (std::map<int,std::vector<double> >::iterator it = dataset.begin(), ed = dataset.end(); it != ed; ++it):
    for idx, data in enumerate(dataset) :
        #data is it.second # std::vector<double> &data = it.second
        nd=len(data) # int nd = data.size()
        data.sort() # std::sort(data.begin(), data.end())
        if nd % 2 == 0 :
            median = 0.5*(data[nd/2]+data[nd/2+1])
        else :
            median = data[nd/2]
            
        if (pdefs.band_safety_crop > 0):
            data2 = [] # std::vector<double> data2
            for j in xrange(nd): # for (int j = 0; j < nd; ++j):
                if (data[j] > median*pdefs.band_safety_crop and data[j] < median/pdefs.band_safety_crop):
                    data2.append(data[j]) # data2.push_back(data[j])
                
            
            data2.swap(data)
            datatmp=data2, data
            data=datatmp[0] 
            data2=datatmp[1]
            nd = len(data) #data.size()
            if nd % 2 == 0 : median = 0.5*(data[nd/2]+data[nd/2+1])
            else : median =  data[nd/2]
        mean = 0
        for j in xrange(nd):
            mean += data[j]
        mean /= nd
        summer68 = data[floor(nd * 0.5*(1-width)+0.5)]
        winter68 =  data[min(int(floor(nd * 0.5*(1+width)+0.5)), nd-1)]
        if (pdefs.use_precomputed_quantiles and type == BandType.Median):
            if (pdefs.precomputed_median_only and len(data) == 1):
                mean = median = summer68 = winter68 = data[0]
            elif (len(data) != 3): 
                print "Error for expected quantile for mass " , idx , ": size of data is " , len(data)
                continue
            else :
                mean = median = data[1]; summer68 = data[0]; winter68 = data[2]
            
        
        x = mean
        #switch (type):
        if type == BandType.MeanCPUTime: pass
        elif type == BandType.MeanRealTime: pass
        elif type == BandType.Mean: x = mean; break
        elif type == BandType.Median: x = median; break
        elif type == BandType.CountToys: x = summer68 = winter68 = nd; break
        elif type == BandType.Asimov: #// mean (in case we did it more than once), with no band
            if obs_avg_mode == mean : winter68 = mean
            else : winter68 = median
            x = summer68 = winter68 # = (obs_avg_mode == mean ? mean : median)
            break
        elif type == BandType.Observed:
            x = mean
            if (nd == 1):
                if (len(errors[idx]) == 1):
                    summer68 = mean - errors[idx][0]
                    winter68 = mean + errors[idx][0]
                else :
                    #// could happen if limitErr is not available
                    summer68 = winter68 = mean

            else : #// if we have multiple, average and report rms (useful e.g. for MCMC)
                #switch (obs_avg_mode):
                if obs_avg_mode == ObsAvgMode.MeanObs:   x = mean; break
                elif obs_avg_mode ==  ObsAvgMode.MedianObs: x = median; break
                elif obs_avg_mode ==  ObsAvgMode.LogMeanObs:
                    x = 0
                    for j in xrange(nd) : x += log(data[j]) #for (int j = 0; j < nd; ++j): x += log(data[j]); 
                    x = exp(x/nd)
                    break

                rms = 0
                for j in xrange(nd) : rms += (x-data[j])*(x-data[j]) # for (int j = 0; j < nd; ++j): rms += (x-data[j])*(x-data[j]); 
                rms = sqrt(rms/(nd*(nd-1)))
                summer68 = mean - rms
                winter68 = mean + rms

            break
        elif type == BandType.AdHoc:
            x = summer68 = winter68 = mean
            break
        elif type == BandType.Quantile: #// get the quantile equal to width, and it's uncertainty
            x = data[floor(nd*width+0.5)]
            summer68 = x - quantErr(nd, data, width)
            winter68 = x + (x-summer68)
            break
        elif type == BandType.ObsQuantile:   
            if ( idx == ( len(obsValues) - 1)) : continue # ????? #if (obsValues.find(it.first) == obsValues.end()) continue
            passed = 0
            failed = 0
            for i in xrange(nd): # for (int i = 0; i < nd and data[i] <= obsValues[it.first]; ++i):
                if data[i] > obsValues[idx] : continue
                failed += 1

            passed = nd - failed; x = double(passed)/nd
            alpha = (1.0 - .68540158589942957)/2
            if passed == 0 : summer68 = 0.0
            else : summer68 = ROOT.Math.beta_quantile(   alpha, passed,   failed+1 )
            if failed == 0 : winter68 = 1.0
            else : winter68 = ROOT.Math.beta_quantile( 1-alpha, passed+1, failed   )
            break

        #// end switch
        tge.Set(ip+1)
        tge.SetPoint(ip, it.first*0.01, x)
        tge.SetPointError(ip, 0, 0, x-summer68, winter68-x)
        ip += 1
    
    return tge


# In[7]:


def theFcBelt(file,  doSyst,  whichChannel,  type, width=0.68):
    if (file is None) : return 0
    t = file.Get("limit")
    if (t is None) : t = file.Get("test") #// backwards compatibility
    if (t is None): 
        print "TFile " , file.GetName() , " does not contain the tree" 
        return 0 

    #quant = -1
    #mass = array('d',[])
    #limit = array('d',[])
    #limitErr = array('d',[0])
    #t_cpu = array('d',[])
    #t_real = array('d',[])
    #quant = array('i',[])
    #syst = array('i',[])
    #iChannel = array('i',[])
    #iToy = array('i',[])
    #t.SetBranchAddress("mh", mass)
    #t.SetBranchAddress("limit", limit)
    #t.SetBranchAddress("limitErr", limitErr)
    #t.SetBranchAddress("quantileExpected", quant)
    #t.SetBranchAddress("syst", syst)
    #if pdefs.seed_is_channel : t.SetBranchAddress("iSeed", iChannel)
    #else : t.SetBranchAddress("iChannel", iChannel)
    #t.SetBranchAddress("iToy", iToy)

    fitExp = ROOT.TF1 ("fitExp","[0]*exp([1]*(x-[2]))", 0, 1)
    fitErf = ROOT.TF1 ("fitErf","[0]*TMath::Erfc([1]*abs(x-[2]))", 0, 1)
    dataset = [] # std::map<int,TGraphErrors*>  dataset
    for i in xrange(n): # for (size_t i = 0, n = t.GetEntries(); i < n; ++i):
        t.GetEntry(i)
        mass = t.mh ; limit = t.limit ; limitErr=t.limitErr ; quant = t.quantileExpected ; syst=t.syst
        if pdefs.seed_is_channel : iChannel = t.iSeed
        else : iChannel = t.iChannel
        iToy = t.iToy
        iMass = int(mass*10)
        #//printf("%6d mh=%.1f  limit=%8.3f +/- %8.3f toy=%5d quant=% .3f\n", i, mass, limit, limitErr, iToy, quant)
        if (syst != doSyst) :           continue
        if (iChannel != whichChannel) : continue
        if      (type == BandType.Asimov) :
            if (iToy != -1) : continue 
        elif (type == BandType.Observed): 
            if (iToy !=  0) : continue

        if (quant < 0) : continue
        graph = dataset[iMass] # ????? this is strange where does dataset come from
        if (graph is None) : graph = ROOT.TGraphErrors()
        elif (graph == 0) : graph = ROOT.TGraphErrors()
        ipoint = graph.GetN()
        graph.Set(ipoint+1)
        graph.SetPoint(ipoint, limit, quant)
        graph.SetPointError(ipoint, 0, limitErr)
    
    #//std::cout << "Loaded " << dataset.size() << " masses " << std::endl
    #*tge = new (); int ip = 0
    tge = ROOT.TGraphAsymmErrors () # new () # TGraphAsymmErrors
    ip = 0
    for idx, graph in enumerate(datasets) : # for (std::map<int,TGraphErrors*>::iterator it = dataset.begin(), ed = dataset.end(); it != ed; ++it):
        #TGraphErrors *graph = it.second; 
        graph.sort()
        n = graph.GetN()
        if (n < 3) : continue
        #//std::cout << "For mass " << it.first/10 << " I have " << n << " points" << std::endl

        #blow, bhigh, bmid

        imax = 0
        ymax = graph.GetY()[0]
        for i in xrange(n) : # for (int i = 0; i < n; ++i):
            #//printf(" i = %2d mH = %.1f, r = %6.3f, pval = %8.6f +/- %8.6f\n", i, it.first/10., graph.GetX()[i], graph.GetY()[i], graph.GetEY()[i])
            if (graph.GetY()[i] > ymax):
                imax = i; ymax = graph.GetY()[i]
            
        
        if (imax == 0):
            bmid = graph.GetX()[0]
        elif (imax == n-1):
            bmid = graph.GetX()[n-1]
        else :
            #// ad hoc method
            sumxmax = 0 ; sumwmax = 0
            ib=max(0,imax-5)
            ie=max(n-1,imax+5)
            for j in xrange(ib,ie) : #for (int i = std::max<int>(0, imax-5); i < std::max<int>(n-1,imax+5); ++i):
                y4 = pow(graph.GetY()[j],4)
                sumxmax += graph.GetX()[j] * y4; sumwmax += y4
            
            bmid = sumxmax/sumwmax
        

        #//std::cout << "band center for " << it.first/10 << " is at " << bmid << " (imax = " << imax << ")\n" << std::endl

        if (graph.GetY()[0] > 1-width or imax == 0):
            blow = graph.GetX()[0]
        else:
            ilo = 0 ;  ihi = 0
            for ilo in xrange (1,imax) : # for (ilo = 1;  ilo < imax; ++ilo):
                if (graph.GetEY()[ilo] == 0) : continue
                if (graph.GetY()[ilo]  >= 0.05*(1-width)) : break
            
            ilo -= 1
            for ihi in xrange(imax, ilo+1,-1) : # for (ihi = imax; ihi > ilo+1; --ihi):
                if (graph.GetEY()[ihi] == 0)  : continue
                if (graph.GetY()[ihi] <= 3*(1-width)) : break
            
            xmin = graph.GetX()[ilo]
            xmax = graph.GetX()[ihi]
            if (ilo <= 1) : xmin = 0.0001
            fitErf.SetRange(xmin,xmax); fitErf.SetNpx(1000)
            fitErf.SetParameters(0.6,bmid,2.0/bmid)
            graph.Fit(fitErf,"WNR EX0","",xmin,xmax)
            fitErf.SetNpx(4)
            blow = fitErf.GetX(1-width,xmin,xmax)
            if (blow <= 2*0.0001) : blow = 0
            #//std::cout << width << " band low end " << it.first/10 << " is at " << blow << " (xmin = " << xmin << ", xmax = " << xmax << ")\n" << std::endl
        
        
        if (graph.GetY()[n-1] > 1-width or imax == n-1):
            bhigh = graph.GetX()[n-1]
        elif (imax == 0 and graph.GetY()[1] < 1-width):
            xmin = graph.GetX()[1], xmax = graph.GetX()[2]
            for i in xrange(3,max(5,n)+1) : # for (int i = 3; i <= std::max<int>(5,n); ++i):
                 if (graph.GetY()[i] < 0.5*(1-width)) : break
                 xmax = graph.GetX()[i]
            
            fitExp.SetRange(xmin,xmax)
            fitExp.SetNpx(1000)
            fitExp.SetParameters(1-width, -2.0/(xmax-xmin), 0.5*(xmax-xmin))
            fitExp.FixParameter(0,1-width)
            graph.Fit(fitExp,"WNR EX0","",xmin,xmax)
            bhigh = fitExp.GetParameter(2)
            if (bhigh < graph.GetX()[0]):
                bhigh = graph.GetX()[0] + ((1-width)-graph.GetY()[0])*(graph.GetX()[1]-graph.GetX()[0])/(graph.GetY()[1]-graph.GetY()[0])
                #//std::cout << width << " band high end forces stupid linear interpolation" << std::endl
            
            #//std::cout << width << " band high end " << it.first/10 << " is at " << bhigh << " (xmin = " << xmin << ", xmax = " << xmax << ")\n" << std::endl
        else :
            #int ilo = 0, ihi = 0
            for ilo in xrange(imax+1, n-2) : # for (ilo = imax+1;  ilo < n-2; ++ilo):
                if (graph.GetEY()[ilo] == 0) : continue
                if (graph.GetY()[ilo]  <= 3*(1-width)) : break
            
            if (ilo > 0 and graph.GetEY()[ilo-1] != 0) : ilo -= 1
            for ihi in xrange(ilo+1,n) : # for (ihi = ilo+1; ihi < n; ++ihi):
                if (graph.GetEY()[ihi] == 0): ihi -= 1 ; break 
                if (graph.GetY()[ihi] >= 0.05*(1-width)) : break
            
            if (ihi - ilo <= 1): 
                xmin = graph.GetX()[ilo] ;  xmax = graph.GetX()[ihi]
                bhigh = 0.5*(xmin+xmax)
                #//std::cout << width << " band high end " << it.first/10 << " is " << bhigh << " (xmin = " << xmin << ", xmax = " << xmax << ", no fit)\n" << std::endl
            else :
                xmin = graph.GetX()[ilo] ; xmax = graph.GetX()[ihi]
                fitExp.SetRange(xmin,xmax); fitExp.SetNpx(1000)
                fitExp.SetParameters(1-width, -2.0/(xmax-xmin), 0.5*(xmax-xmin))
                fitExp.FixParameter(0,1-width)
                graph.Fit(fitExp,"WNR EX0","",xmin,xmax)
                bhigh = fitExp.GetParameter(2)
                #//std::cout << width << " band high end " << it.first/10 << " is at " << bhigh << " (xmin = " << xmin << ", xmax = " << xmax << ")\n" << std::endl
            
        

        tge.Set(ip+1)
        tge.SetPoint(ip, idx*0.1, bmid)
        tge.SetPointError(ip, 0, 0, bmid-blow, bhigh-bmid)
        ip += 1

        continue
        #delete graph
    
    return tge


# In[8]:


#def theBand():
#    pass
do_bands_95 = True


# In[9]:


def makeBand(bands, name,  file=ROOT.TFile(),  doSyst=0,  whichChannel=0,  type=BandType.Median):
    suffix = ""
    #switch (type):
    if type == BandType.Asimov:    suffix = "_asimov" #; break
    elif type == BandType.Observed:  suffix = "_obs" #; break
    elif type == BandType.Mean:      suffix = "_mean" #; break
    elif type == BandType.Median:    suffix = "_median" #; break
    elif type == BandType.CountToys: suffix = "_ntoys" #; break
    elif type == BandType.MeanCPUTime: suffix = "_cputime" # ; break
    elif type == BandType.MeanRealTime: suffix = "_realtime" #; break
    elif type == BandType.AdHoc:       suffix = "" # ; break
    elif type == BandType.Quantile:    suffix = "" #; break
    elif type == BandType.ObsQuantile:    suffix = "_qobs" #; break
   
    print "DEBUG makeBand BandType ",type
    
    if (not doSyst and (type != BandType.AdHoc)): suffix = "_nosyst"+suffix
    if (type == BandType.Median or type == BandType.Mean):
        band68 = theBand(file, doSyst, whichChannel, type, 0.68)
        if pdefs.do_bands_95 : band95 = theBand(file, doSyst, whichChannel, type, 0.95)
        else : band95 = 0
        if (band68 != 0 and band68.GetN() > 0):
            band68.SetName(name+suffix)
            bands.WriteTObject(band68, name+suffix)
            if (pdefs.do_bands_95):
                band95.SetName(name+suffix+"_95")
                bands.WriteTObject(band95, name+suffix+"_95")
            
        else :
            print "Band " , name+suffix , " missing 1"
        
    else :
        band = theBand(file, doSyst, whichChannel, type)
        if band and band.GetN() > 0:
            band.SetName(name+suffix)
            bands.WriteTObject(band, name+suffix)
            print "Band " , name+suffix , " found (" , band.GetN() , " points)"
        else :
            print "Band " , name+suffix , " missing 2"


# In[10]:


def makeBandfromFilename(bands, name, filename,  doSyst,  whichChannel,  type):
    input = ROOT.TFile.Open(filename)
    if not input: print "input not defined for ",filename
    if input.isZombie() : print "Filename " , filename , " missing" #== 0): std::cerr << "Filename " << filename << " missing" << std::endl; return; 
    makeBand(bands, name, input, doSyst, whichChannel, type)


# In[11]:


def makeLine( bands, name, filename,   doSyst,  whichChannel):
    input = ROOT.TFile.Open(filename)
    if input.isZombie() : 
        print "Filename '" , filename , "' missing"
        return 
    makeBand(bands, name, input, doSyst, whichChannel, BandType.AdHoc)
    input.Close()


# In[12]:


def makeBands(bands, name, filename,  channel=0, quantiles=False):
    #pdefs.do_bands_nosyst = False
    #pdefs.do_bands_mean = False
    #pdefs.do_bands_median = True
    #pdefs.do_bands_ntoys = False
    #pdefs.do_bands_asimov = True
    #pdefs.do_bands_cputime = False
    #pdefs.do_bands_realtime = False
    #pdefs.do_bands_95 = True
    input = ROOT.TFile.Open(filename)
    if not input: print "input not defined"; return
    print "DEBUG checking ",filename
    print "DEBUG makeBands 0 do_bands_nosyst " , pdefs.do_bands_nosyst
    print "DEBUG makeBands 0 do_bands_mean " , pdefs.do_bands_mean
    print "                1 do_bands_median " , pdefs.do_bands_median
    print "                2 do_bands_quant " , ""
    print "                3 do_bands_obs ", ""
    print "                4 do_bands_asimov ", pdefs.do_bands_asimov
    print "                5 do_bands_ntoys ", pdefs.do_bands_ntoys
    print "                  do_bands_cputime ",""
    print "                  do_bands_95 ",pdefs.do_bands_95
    #input.ls()
    #print "DEBUG end of checking ",filename
    #if (input.isZombie()): print "Filename " , filename , " missing"; return; 
    if pdefs.do_bands_nosyst : isb = 0
    else : isb = 1
    for s in xrange(isb, 2): # for (int s = do_bands_nosyst ? 0 : 1; s <= 1; ++s):
        print "DEBUG makeBands s " , s ," do_bands_nosyst " , pdefs.do_bands_nosyst , " asimov ",pdefs.do_bands_asimov
        if (pdefs.do_bands_mean) : makeBand(bands, name, input, s, channel, BandType.Mean)
        makeBand(bands, name, input, s, channel, BandType.Median)
        makeBand(bands, name, input, s, channel, BandType.Observed)
        if (pdefs.do_bands_ntoys)  : makeBand(bands, name, input, s, channel, BandType.CountToys)
	#if (pdefs.do_bands_asimov) : print "DEBUG makeBands doing  makeBand(bands, name, in, s, channel, Asimov)"
        if (pdefs.do_bands_asimov) : makeBand(bands, name, input, s, channel, BandType.Asimov)
        if (pdefs.do_bands_cputime) : makeBand(bands, name, input, s, channel, BandType.MeanCPUTime)
        if (pdefs.do_bands_realtime) : makeBand(bands, name, input, s, channel, BandType.MeanRealTime)
    
    if (quantiles):
        quants = [ 0.025, 0.16, 0.5, 0.84, 0.975 ] 
        for i in xrange(5) : # for (int i = 0; i < 5; ++i):
            for s in xrange(0,2) : # for (int s = 0; s <= 1; ++s):
                band = theBand(input, s, channel, Quantile, quants[i])
                if s == 0 : qname = ROOT.TString.Format("%s%s_quant%03d", str(name), "_nosyst", int(1000*quants[i]))
                else : qname = ROOT.TString.Format("%s%s_quant%03d", str(name), "", int(1000*quants[i]))
                if (band != 0 and band.GetN() != 0):
                    band.SetName(qname)
                    bands.WriteTObject(band, qname)
                else :
                    print "Band " , qname , " missing"
                
            
        
    
    input.Close()


# In[13]:


def findBin(g, x, tolerance):
    if not g : return -1
    for i in xrange(g.GetN()) : # for (int i = 0; i < g.GetN(); ++i):
        xi = g.GetX()[i]
        if (abs(xi -  x) < tolerance):
            return i
        
    
    return -1


# In[14]:


def findBin(g, x):
    #if not g: print "DEBUG findBin  return -1 for ",x
    if not g: return -1
    for i in xrange(g.GetN()) :
        xi = g.GetX()[i]
        #print "DEBUG findBin " , xi - g.GetErrorXlow(i) , " <= " , x , " <= " , xi + g.GetErrorXhigh(i), " xi ", xi, " low ",g.GetErrorXlow(i), " high ",g.GetErrorXhigh(i)
        if (  ( xi - g.GetErrorXlow(i) <= x ) and ( x <= xi + g.GetErrorXhigh(i) ) ):
            #print "DEBUG findBin  return ", i, " for ",x
            return i
        
    #print "DEBUG findBin  return -1 for not finding bin for ",x
    return -1

# In[15]:


def significanceToPVal( bands, inName, outName):
    b1 = bands.Get(inName)
    if (b1 is None) : return
    if (b1 == 0 or b1.GetN() == 0) : return
    n = b1.GetN()
    b2 = ROOT.TGraphAsymmErrors (n)
    for i in xrange(n) : # for (int i = 0; i < n; ++i):
        x = b1.GetX()[i]
        s = b1.GetY()[i]
        slo = s - b1.GetErrorYlow(i) ; shi = s + b1.GetErrorYhigh(i)
        pval = ROOT.Math.normal_cdf_c(s)
        phi  = ROOT.Math.normal_cdf_c(slo)
        plo  = ROOT.Math.normal_cdf_c(shi)
        b2.SetPoint(i, x, pval)
        b2.SetPointError(i, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), pval - plo, phi - pval)
    
    b2.SetName(outName)
    bands.WriteTObject(b2, outName)


# In[16]:


def pvalToSignificance(bands, inName, outName):
    b1 = bands.Get(inName)
    if (b1 is None) : return
    if (b1 == 0 or b1.GetN() == 0) : return
    n = b1.GetN()
    b2 = ROOT.TGraphAsymmErrors (n)
    for i in xrange(n) : # for (int i = 0; i < n; ++i):
        x = b1.GetX()[i]; s = b1.GetY()[i]
        slo = s - b1.GetErrorYlow(i); shi = s + b1.GetErrorYhigh(i)
        pval = ROOT.Math.normal_quantile_c(s,   1.0)
        phi  = ROOT.Math.normal_quantile_c(slo, 1.0)
        plo  = ROOT.Math.normal_quantile_c(shi, 1.0)
        b2.SetPoint(i, x, pval)
        b2.SetPointError(i, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), pval - plo, phi - pval)
    
    b2.SetName(outName)
    bands.WriteTObject(b2, outName)



# In[17]:


def testStatToPVal(bands, inName, outName):
    b1 = bands.Get(inName)
    if (b1 is None) : return
    if (b1 == 0 or b1.GetN() == 0): return
    n = b1.GetN()
    b2 = ROOT.TGraphAsymmErrors (n)
    for i in xrange(n) : # for (int i = 0; i < n; ++i):
        x = b1.GetX()[i]; s = b1.GetY()[i]
        slo = s - b1.GetErrorYlow(i); shi = s + b1.GetErrorYhigh(i)
        if s > 0 : s = sqrt(s)
        else : s   = 0 # (s   > 0 ? sqrt(s  ) : 0)
        if slo > 0 : slo = sqrt(slo)
        else : slo = 0 # slo = (slo > 0 ? sqrt(slo) : 0)
        if shi > 0 : shi = sqrt(shi)
        else : shi = 0 # shi = (shi > 0 ? sqrt(shi) : 0)
        pval = ROOT.Math.normal_cdf_c(s)
        phi  = ROOT.Math.normal_cdf_c(slo)
        plo  = ROOT.Math.normal_cdf_c(shi)
        b2.SetPoint(i, x, pval)
        b2.SetPointError(i, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), pval - plo, phi - pval)

    
    b2.SetName(outName)
    bands.WriteTObject(b2, outName)


# In[18]:


def significanceToPVals(bands, inName, outName):
    significanceToPVal(bands, inName+"_obs",   outName+"_obs")
    significanceToPVal(bands, inName+"_mean",   outName+"_mean")
    significanceToPVal(bands, inName+"_median", outName+"_median")
    significanceToPVal(bands, inName+"_mean_95",   outName+"_mean_95")
    significanceToPVal(bands, inName+"_median_95", outName+"_median_95")
    significanceToPVal(bands, inName+"_asimov",    outName+"_asimov")

    significanceToPVal(bands, inName+"_nosyst_obs",   outName+"_nosyst_obs")
    significanceToPVal(bands, inName+"_nosyst_mean",   outName+"_nosyst_mean")
    significanceToPVal(bands, inName+"_nosyst_median", outName+"_nosyst_median")
    significanceToPVal(bands, inName+"_nosyst_mean_95",   outName+"_nosyst_mean_95")
    significanceToPVal(bands, inName+"_nosyst_asimov",    outName+"_nosyst_asimov")
    significanceToPVal(bands, inName+"_nosyst_ntoys",     outName+"_nosyst_ntoys")


# In[19]:


def testStatToPVals(bands, inName, outName):
    testStatToPVal(bands, inName+"_obs",   outName+"_obs")
    testStatToPVal(bands, inName+"_mean",   outName+"_mean")
    testStatToPVal(bands, inName+"_median", outName+"_median")
    testStatToPVal(bands, inName+"_mean_95",   outName+"_mean_95")
    testStatToPVal(bands, inName+"_median_95", outName+"_median_95")
    testStatToPVal(bands, inName+"_asimov",    outName+"_asimov")

    testStatToPVal(bands, inName+"_nosyst_obs",   outName+"_nosyst_obs")
    testStatToPVal(bands, inName+"_nosyst_mean",   outName+"_nosyst_mean")
    testStatToPVal(bands, inName+"_nosyst_median", outName+"_nosyst_median")
    testStatToPVal(bands, inName+"_nosyst_mean_95",   outName+"_nosyst_mean_95")
    testStatToPVal(bands, inName+"_nosyst_asimov",    outName+"_nosyst_asimov")
    testStatToPVal(bands, inName+"_nosyst_ntoys",     outName+"_nosyst_ntoys")


# In[20]:


def cutBand(bands, inName, outName, mMin, mMax):
    b1 =  bands.Get(inName)
    if not b1 : return
    if b1 == 0 : return
    if b1.GetN() == 0 : return
    b2 = ROOT.TGraphAsymmErrors()
    n = b1.GetN(); m = 0
    for i in xrange(n) : # for (int i = 0; i < n; ++i):
        if (mMin <= b1.GetX()[i] and b1.GetX()[i] <= mMax):
            b2.Set(m+1)
            b2.SetPoint(m, b1.GetX()[i], b1.GetY()[i])
            b2.SetPointError(m, b1.GetErrorXlow(i), b1.GetErrorXhigh(i),
                                 b1.GetErrorYlow(i), b1.GetErrorYhigh(i))
            m += 1
        
    
    b2.SetName(outName)
    bands.WriteTObject(b2, outName)
    #//bands.Add(b2)


# In[21]:


def cutBands(bands, inName, outName,  mMin,  mMax):
    cutBand(bands, inName+"_obs",   outName+"_obs",    mMin, mMax)
    cutBand(bands, inName+"_mean",   outName+"_mean",    mMin, mMax)
    cutBand(bands, inName+"_median", outName+"_median",  mMin, mMax)
    cutBand(bands, inName+"_mean_95",   outName+"_mean_95",    mMin, mMax)
    cutBand(bands, inName+"_median_95", outName+"_median_95",  mMin, mMax)
    cutBand(bands, inName+"_asimov",    outName+"_asimov",     mMin, mMax)
    cutBand(bands, inName+"_ntoys",     outName+"_ntoys",      mMin, mMax)

    cutBand(bands, inName+"_nosyst_obs",   outName+"_nosyst_obs",    mMin, mMax)
    cutBand(bands, inName+"_nosyst_mean",   outName+"_nosyst_mean",    mMin, mMax)
    cutBand(bands, inName+"_nosyst_median", outName+"_nosyst_median",  mMin, mMax)
    cutBand(bands, inName+"_nosyst_mean_95",   outName+"_nosyst_mean_95",    mMin, mMax)
    cutBand(bands, inName+"_nosyst_asimov",    outName+"_nosyst_asimov",     mMin, mMax)
    cutBand(bands, inName+"_nosyst_ntoys",     outName+"_nosyst_ntoys",      mMin, mMax)



# In[22]:


def cutFcBands(bands, inName, outName,  mMin,  mMax,  npostfix, postfixes):
    for i in xrange(postfix) : # for (int i = 0; i < npostfix; ++i):
        cutBand(bands, inName+"_"+postfixes[i],   outName+"_"+postfixes[i],  mMin, mMax)


# In[23]:


def combineBand(input, band1, band2, comb):
    b1 =  input.Get(band1)
    b2 =  input.Get(band2)
    if (b1 is None) : return
    if (b2 is None) : return
    if (b1 == 0 or b1.GetN() == 0): return
    if (b2 == 0 or b2.GetN() == 0): return
    n = b1.GetN(); m = b2.GetN()
    first = n; last = 0
    for i in xrange(n) : #for (int i = 0; i < n; ++i):
        for j in xrange(m) : # for (int j = 0; j < m; ++j):
            if (int(b1.GetX()[i]) == int(b2.GetX()[j])):
                if (i < first) : first = i
                last = i
            
        
    
    bc = ROOT.TGraphAsymmErrors((first-1) + m + (n-last-1))
    bc.SetName(comb)
    k = 0
    for i in xrange(first) : # for (int i = 0; i < first; ++i, ++k):
        k += 1
        bc.SetPoint(k, b1.GetX()[i], b1.GetY()[i])
        bc.SetPointError(k, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), 
                             b1.GetErrorYlow(i), b1.GetErrorYhigh(i))
    
    for i in xrange(m) : # for (int i = 0; i < m; ++i, ++k):
        k += 1
        bc.SetPoint(k, b2.GetX()[i], b2.GetY()[i])
        bc.SetPointError(k, b2.GetErrorXlow(i), b2.GetErrorXhigh(i), 
                             b2.GetErrorYlow(i), b2.GetErrorYhigh(i))
    
    for i in xrange(last+1,n) : # for (int i = last+1; i < n; ++i, ++k):
        k += 1
        bc.SetPoint(k, b1.GetX()[i], b1.GetY()[i])
        bc.SetPointError(k, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), 
                             b1.GetErrorYlow(i), b1.GetErrorYhigh(i))
    
    bc.SetName(comb)
    input.WriteTObject(bc, comb)


# In[24]:


def combineBands(input, band1, band2, comb):
    combineBand(input, band1+"_mean",   band2+"_mean",   comb+"_mean")
    combineBand(input, band1+"_median", band2+"_median", comb+"_median")
    combineBand(input, band1+"_mean_95",   band2+"_mean_95",   comb+"_mean_95")
    combineBand(input, band1+"_median_95", band2+"_median_95", comb+"_median_95")
    combineBand(input, band1+"_asimov",    band2+"_asimov",    comb+"_asimov")
    combineBand(input, band1+"_ntoys",    band2+"_ntoys",    comb+"_ntoys")

    combineBand(input, band1+"_nosyst_mean",   band2+"_nosyst_mean",   comb+"_nosyst_mean")
    combineBand(input, band1+"_nosyst_median", band2+"_nosyst_median", comb+"_nosyst_median")
    combineBand(input, band1+"_nosyst_mean_95",   band2+"_nosyst_mean_95",   comb+"_nosyst_mean_95")
    combineBand(input, band1+"_nosyst_median_95", band2+"_nosyst_median_95", comb+"_nosyst_median_95")
    combineBand(input, band1+"_nosyst_asimov",    band2+"_nosyst_asimov",    comb+"_nosyst_asimov")
    combineBand(input, band1+"_nosyst_ntoys",    band2+"_nosyst_ntoys",    comb+"_nosyst_ntoys")


# In[25]:


def mergeBand(input, band1, band2, comb):
    b1 =  input.Get(band1)
    b2 =  input.Get(band2)
    if (b1 is None) : return
    if (b2 is None) : return
    if (b1 == 0 or b1.GetN() == 0): return
    if (b2 == 0 or b2.GetN() == 0): return
    n = b1.GetN(); m = b2.GetN()
    bc = ROOT.TGraphAsymmErrors(n)
    bc.SetName(comb)
    k = 0
    for i in xrange(n) : # for (int i = 0; i < n; ++i, ++k):
        k += 1
        bc.SetPoint(k, b1.GetX()[i], b1.GetY()[i])
        bc.SetPointError(k, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), 
                             b1.GetErrorYlow(i), b1.GetErrorYhigh(i))
    
    for i in xrange(m) : # for (int i = 0; i < m; ++i):
        if (findBin(b1, b2.GetX()[i], 0.001) != -1) : continue
        bc.Set(k)
        bc.SetPoint(k, b2.GetX()[i], b2.GetY()[i])
        bc.SetPointError(k, b2.GetErrorXlow(i), b2.GetErrorXhigh(i), 
                             b2.GetErrorYlow(i), b2.GetErrorYhigh(i))
        k +=1
    
    bc.Sort()
    bc.SetName(comb)
    input.WriteTObject(bc, comb, "Overwrite")


# In[26]:


def mergeBands(input, band1, band2, comb):
    mergeBand(input, band1+"_obs",   band2+"_obs",   comb+"_obs")
    mergeBand(input, band1+"_mean",   band2+"_mean",   comb+"_mean")
    mergeBand(input, band1+"_median", band2+"_median", comb+"_median")
    mergeBand(input, band1+"_mean_95",   band2+"_mean_95",   comb+"_mean_95")
    mergeBand(input, band1+"_median_95", band2+"_median_95", comb+"_median_95")
    mergeBand(input, band1+"_asimov",    band2+"_asimov",    comb+"_asimov")
    mergeBand(input, band1+"_ntoys",    band2+"_ntoys",    comb+"_ntoys")

    mergeBand(input, band1+"_nosyst_mean",   band2+"_nosyst_mean",   comb+"_nosyst_mean")
    mergeBand(input, band1+"_nosyst_median", band2+"_nosyst_median", comb+"_nosyst_median")
    mergeBand(input, band1+"_nosyst_mean_95",   band2+"_nosyst_mean_95",   comb+"_nosyst_mean_95")
    mergeBand(input, band1+"_nosyst_median_95", band2+"_nosyst_median_95", comb+"_nosyst_median_95")
    mergeBand(input, band1+"_nosyst_asimov",    band2+"_nosyst_asimov",    comb+"_nosyst_asimov")
    mergeBand(input, band1+"_nosyst_ntoys",    band2+"_nosyst_ntoys",    comb+"_nosyst_ntoys")


# In[27]:


def pasteBand(input, band1, band2, comb):
    b1 =  input.Get(band1)
    b2 =  input.Get(band2)
    if not b1 : return
    if b1 == 0: return
    if b1.GetN() == 0 : return
    if not b2 : return
    if b2 == 0: return
    if b2.GetN() == 0 : return
    bc = ROOT.TGraphAsymmErrors(b1.GetN()+b2.GetN())
    bc.SetName(comb)
    k = 0 ; n = b1.GetN(); m = b2.GetN()
    for i in xrange(n): # for (int i = 0; i < n; ++i, ++k):
        k +=1
        bc.SetPoint(k, b1.GetX()[i], b1.GetY()[i])
        bc.SetPointError(k, b1.GetErrorXlow(i), b1.GetErrorXhigh(i), 
                             b1.GetErrorYlow(i), b1.GetErrorYhigh(i))
    
    for i in xrange(m) : # for (int i = 0; i < m; ++i, ++k):
        k +=1
        bc.SetPoint(k, b2.GetX()[i], b2.GetY()[i])
        bc.SetPointError(k, b2.GetErrorXlow(i), b2.GetErrorXhigh(i), 
                             b2.GetErrorYlow(i), b2.GetErrorYhigh(i))
    
    bc.Sort()
    input.WriteTObject(bc, comb)


# In[28]:


def pasteBands(input, band1, band2, comb):
    pasteBand(input, band1+"_obs",   band2+"_obs",   comb+"_obs")
    pasteBand(input, band1+"_mean",   band2+"_mean",   comb+"_mean")
    pasteBand(input, band1+"_median", band2+"_median", comb+"_median")
    pasteBand(input, band1+"_mean_95",   band2+"_mean_95",   comb+"_mean_95")
    pasteBand(input, band1+"_median_95", band2+"_median_95", comb+"_median_95")
    pasteBand(input, band1+"_asimov",    band2+"_asimov",    comb+"_asimov")
    pasteBand(input, band1+"_ntoys",    band2+"_ntoys",    comb+"_ntoys")

    pasteBand(input, band1+"_nosyst_obs",   band2+"_nosyst_obs",   comb+"_nosyst_obs")
    pasteBand(input, band1+"_nosyst_mean",   band2+"_nosyst_mean",   comb+"_nosyst_mean")
    pasteBand(input, band1+"_nosyst_median", band2+"_nosyst_median", comb+"_nosyst_median")
    pasteBand(input, band1+"_nosyst_mean_95",   band2+"_nosyst_mean_95",   comb+"_nosyst_mean_95")
    pasteBand(input, band1+"_nosyst_median_95", band2+"_nosyst_median_95", comb+"_nosyst_median_95")
    pasteBand(input, band1+"_nosyst_asimov",    band2+"_nosyst_asimov",    comb+"_nosyst_asimov")
    pasteBand(input, band1+"_nosyst_ntoys",    band2+"_nosyst_ntoys",    comb+"_nosyst_ntoys")


# In[29]:


def pasteFcBands(bands, band1, band2, comb,  npostfix,  postfixes):
    for i in xrange(npostfix) : # for (int i = 0; i < npostfix; ++i):
        pasteBand(bands, band1+"_"+postfixes[i],   band2+"_"+postfixes[i],   comb+"_"+postfixes[i])


# In[30]:


def stripPoint( band,  m):
    n = band.GetN()
    for i in xrange(n): # for (int i = 0, n = band.GetN(); i < n; ++i):
        if (float(band.GetX()[i]) == m):
            band.RemovePoint(i)
            return
        
    
    if ((band.GetN() > 0) and
        (band.GetX()[0] <= m) and
        (band.GetX()[band.GetN()-1] >= m)): pass


# In[31]:


def stripBand(input, band1,  m1,  m2=0,  m3=0,  m4=0,  m5=0):
    band =  input.Get(band1)
    if (band is None) : return
    if (band == 0 or band.GetN() == 0) : return
    if (m1): stripPoint(band,m1)
    if (m2): stripPoint(band,m2)
    if (m3): stripPoint(band,m3)
    if (m4): stripPoint(band,m4)
    if (m5): stripPoint(band,m5)
    input.WriteTObject(band, band.GetName(), "Overwrite")


# In[32]:


def stripBands(input, band,   m1,  m2=0,  m3=0,  m4=0,  m5=0):
    stripBand(input, band+"_obs",       m1,m2,m3,m4,m5)
    stripBand(input, band+"_mean",      m1,m2,m3,m4,m5)
    stripBand(input, band+"_median",    m1,m2,m3,m4,m5)
    stripBand(input, band+"_mean_95",   m1,m2,m3,m4,m5)
    stripBand(input, band+"_median_95", m1,m2,m3,m4,m5)
    stripBand(input, band+"_asimov",    m1,m2,m3,m4,m5)
    stripBand(input, band+"_nosyst_obs",       m1,m2,m3,m4,m5)
    stripBand(input, band+"_nosyst_mean",      m1,m2,m3,m4,m5)
    stripBand(input, band+"_nosyst_median",    m1,m2,m3,m4,m5)
    stripBand(input, band+"_nosyst_mean_95",   m1,m2,m3,m4,m5)
    stripBand(input, band+"_nosyst_median_95", m1,m2,m3,m4,m5)
    stripBand(input, band+"_nosyst_asimov",    m1,m2,m3,m4,m5)


# In[33]:


def copyPoint(fromwhat,  m, to,  idx=-1):
    j = findBin(fromwhat, m)
    if (j == -1) :return
    if (idx == -1): idx = to.GetN(); to.Set(idx+1); 
    to.SetPoint(idx, fromwhat.GetX()[j], fromwhat.GetY()[j])
    to.SetPointError(idx, fromwhat.GetErrorXlow(j), fromwhat.GetErrorXhigh(j), fromwhat.GetErrorYlow(j), fromwhat.GetErrorYhigh(j))


# In[34]:


def selectedPointsBand(input, band1, band2,  m1,  m2=0,  m3=0,  m4=0,  m5=0,  m6=0,  m7=0):
    band =  input.Get(band1)
    if not band: print "band is nil" ; return
 
    ret = ROOT.TGraphAsymmErrors()
    copyPoint(band, m1, ret)
    if (m2): copyPoint(band, m2, ret)
    if (m3): copyPoint(band, m3, ret)
    if (m4): copyPoint(band, m4, ret)
    if (m5): copyPoint(band, m5, ret)
    if (m6): copyPoint(band, m6, ret)
    if (m7): copyPoint(band, m7, ret)
    ret.SetName(band2)
    input.WriteTObject(ret, ret.GetName(), "Overwrite")


# In[35]:


def selectedPointsBands(bandIn): #// Directory *in): #//, bandIn, bandOut):
  #//std::cout << "Input = " << in.GetName() << std::endl #//" bandIn = " << bandIn << " bandOut = " << bandOut << std::endl #// << " m1 = " << m1 << std::endl
  print " bandIn = " , bandIn


# In[36]:


def selectedPointsBands(input): #//, bandIn, bandOut):
  print "Input = " , input.GetName() #//" bandIn = " << bandIn << " bandOut = " << bandOut << std::endl #// << " m1 = " << m1 << std::endl
  #//std::cout << " bandIn = " << bandIn << std::endl


#//void selectedPointsBands(in): #//, bandIn, bandOut):
#//  std::cout << "Input = " << in.GetName() << std::endl #//" bandIn = " << bandIn << " bandOut = " << bandOut << std::endl #// << " m1 = " << m1 << std::endl
#//


# In[37]:


def selectedPointsBands(input, bandIn, bandOut,   m1,  m2=0,  m3=0,  m4=0,  m5=0,  m6=0,  m7=0):
    print "selectedPointsBands "
    selectedPointsBand(input, bandIn+"_obs",       bandOut+"_obs",       m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_mean",      bandOut+"_mean",      m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_median",    bandOut+"_median",    m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_mean_95",   bandOut+"_mean_95",   m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_median_95", bandOut+"_median_95", m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_asimov",    bandOut+"_asimov",    m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_nosyst_obs",       bandOut+"_nosyst_obs",       m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_nosyst_mean",      bandOut+"_nosyst_mean",      m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_nosyst_median",    bandOut+"_nosyst_median",    m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_nosyst_mean_95",   bandOut+"_nosyst_mean_95",   m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_nosyst_median_95", bandOut+"_nosyst_median_95", m1,m2,m3,m4,m5,m6,m7)
    selectedPointsBand(input, bandIn+"_nosyst_asimov",    bandOut+"_nosyst_asimov",    m1,m2,m3,m4,m5,m6,m7)


# In[38]:


#def printLine(bands, who, fout, header="value"):
def printLineWithFileDescriptor(bands, who, fout, header="value"):
    mean = bands.Get(who)
    if not mean: 
        print "MISSING " , who
        return
    fout.write( "%4s \t %7s\n" % ( "mass",  header) )
    fout.write( "%5s\t %7s\n" % ( "-----", "-----") )
    fmt1="%4.0f "
    if pdefs.halfint_masses : fmt1="%5.1f"
    fmt2="%7.3f"
    if ("pval" in who  or "smcls" in who) : fmt2="%7.5f"
    fmtstring = fmt1 + "\t " + fmt2 + "\n"
    n=mean.GetN()
    for i in xrange(n): # for (int i = 0, n = mean.GetN(); i < n; ++i):
        fout.write( fmtstring % (  mean.GetX()[i], mean.GetY()[i]) )


# In[39]:


def printLine(bands, who, fileName, header="value"):
    mean =  bands.Get(who)
    if (mean == 0): 
        print "MISSING " , who
        return 
    fout = open(fileName, "w")
    printLineWithFileDescriptor(bands,who,fout,header)
    fout.close


# In[40]:


#def printLineErr(bands, who, fout, header="value"):
def printLineErrWithFileDescriptor(bands, who, fout, header="value"):
    mean = bands.Get(who)
    if not mean: 
        print "MISSING " , who
        return 
    fout.write( "%4s \t %7s +/- %6s\n" % ( "mass",  header," error"))
    fout.write( "%5s\t %7s-----%6s-\n" % ("-----", " ------","------"))
    fmt1="%4.0f "
    if pdefs.halfint_masses : fmt1="%5.1f"
    fmt2="%7.3f +/- %6.3f"
    if ("pval" in who  or "smcls" in who) : fmt2="%7.5f +/- %7.5f"
    #fmtstring = TString() +
    #                    (pdefs.halfint_masses ? "%5.1f" : "%4.0f ") + 
    #                    "\t " +
    #                    (who.Contains("pval")  or who.Contains("smcls")  ? "%7.5f +/- %7.5f" : "%7.3f +/- %6.3f") + 
    #                    "\n"
    fmtstring = fmt1 + "\t " + fmt2 + "\n"
    n=mean.GetN()
    for i in xrange(n): # for (int i = 0, n = mean.GetN(); i < n; ++i):
        fout.write(fmtstring % ( mean.GetX()[i],  mean.GetY()[i], ROOT.TMath.Max(mean.GetErrorYlow(i),mean.GetErrorYhigh(i))))


# In[41]:


def printLineErr(bands, who, fileName, header="value"):
    mean =  bands.Get(who)
    if (mean == 0): 
        print "MISSING " , who
        return 
    fout = open(fileName, "w")
    if not fout:
        print "CANNOT WRITE TO " , fileName
        return 
    if (fout is None) :
        print "CANNOT WRITE TO " , fileName
        return 
    printLineErrWithFileDescriptor(bands,who,fout,header)
    fout.close()


# In[42]:


def printLineAErr(bands, who, fout, header="value"):
    mean = bands.Get(who)
    if (mean == 0): 
        print "MISSING " , who
        return
    
    fprintf(fout, "%4s \t %7s  -%6s   +%6s\n",  "mass",  header,"error"," error")
    fmt1="%4.0f "
    if pdefs.halfint_masses : fmt1="%5.1f"
    fmt2="%7.3f  -%6.3f / +%6.3f"
    if ("pval" in who  or "smcls" in who) : fmt2="%7.5f  -%7.5f / +%7.5f" 
    #fmtstring = TString() +
    #                    (pdefs.halfint_masses ? "%5.1f" : "%4.0f ") + 
    #                    "\t " +
    #                    (who.Contains("pval") or who.Contains("smcls")   ? "%7.5f  -%7.5f / +%7.5f" : "%7.3f  -%6.3f / +%6.3f") + 
    #                    "\n"
    fmtstring = TString() + fmt1 + "\t " + fmt2 + "\n"
    fprintf(fout,  "%5s\t %7s------%6s----%6s-\n", "-----", " ------","------","------")
    n = mean.GetN()
    for i in xrange(n): # for (int i = 0, n = mean.GetN(); i < n; ++i):
        fprintf(fout, fmtstring.Data(),  
            mean.GetX()[i], 
            mean.GetY()[i], 
            mean.GetErrorYlow(i),mean.GetErrorYhigh(i))


# In[43]:


def printLineAErr(bands, who, fileName, header="value"):
    mean =  bands.Get(who)
    if (mean == 0): print "MISSING " , who; return; 
    fout = fopen(fileName.Data(), "w")
    printLineAErr(bands,who,fout,header)
    fclose(fout)


# In[49]:


#def printBand(bands, who, fout, mean=False):
def printBandWithFileDescriptor(bands, who, fout, mean=False):
    obs    =  bands.Get(who+"_obs")
    if mean : mean68 =  bands.Get(who+"_mean")
    else : mean68 =  bands.Get(who+"_median")
    if mean : mean95 =  bands.Get(who+"_mean_95")
    else : mean95 =  bands.Get(who+"_median_95")
    if not mean68 and not obs: print "MISSING " , who , "_mean and " , who , "_obs" ; return 
    if not mean68 : printLineErrWithFileDescriptor(bands, who+"_obs", fout); return 
    if mean: fout.write( "%4s \t %8s  %8s  %8s  %8s  %8s  %8s\n" % ( "mass", " obs ", "-95%", "-68%", "mean", "+68%", "+95%"))
    else : fout.write( "%4s \t %8s  %8s  %8s  %8s  %8s  %8s\n" % ( "mass", " obs ", "-95%", "-68%", "median", "+68%", "+95%"))
    fout.write(  "%5s\t %8s  %8s  %8s  %8s  %8s  %8s\n" % ( "-----","-----",  "-----", "-----", "-----", "-----", "-----"))
    fmt1="%4.0f "
    if pdefs.halfint_masses : fmt1="%5.1f"
    fmt2="%8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f"
    if ("pval" in who  or "smcls" in who) : fmt2="%8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f" 
    fmtstring = fmt1 + "\t " + fmt2 + "\n"
    n = mean68.GetN()
    for i in xrange(n) : # for (int i = 0, n = mean68.GetN(); i < n; ++i):
        if obs: j  = findBin(obs,    mean68.GetX()[i])
        else : j = -1
        if mean95: j2 = findBin(mean95, mean68.GetX()[i])
        else : j = -1
        #obsy = obs.GetY()[j]
        if j == -1 : obsy = float("nan")
        else : obsy = obs.GetY()[j]
        #mean95y = mean95.GetY()[j2]-mean95.GetErrorYlow(j2)
        if j2 == -1 : mean95y = float("nan")
        else : mean95y = mean95.GetY()[j2]-mean95.GetErrorYlow(j2)
        #mean95yh = mean95.GetY()[j2]+mean95.GetErrorYhigh(j2)
        if j2 == -1 : mean95yh = float("nan")
        else : mean95yh = mean95.GetY()[j2]+mean95.GetErrorYhigh(j2)
        fout.write( fmtstring % ( 
            mean68.GetX()[i],  
            obsy,
            mean95y, 
            mean68.GetY()[i]-mean68.GetErrorYlow(i), 
            mean68.GetY()[i],
            mean68.GetY()[i]+mean68.GetErrorYhigh(i),
            mean95yh))
        
    

def printFcBand(bands, who, fout,  npostfix, postfixes):
    bs=[]; names=[]; nbands = 0
    for i in xrange(npostfix) : # for (int i = 0; i < npostfix; ++i):
        bs[nbands] =  bands.Get(who+"_"+postfixes[i])
        if (bs[nbands] != 0): names[nbands] = postfixes[i]; nbands +=1 ; 
    
    if (nbands == 0) : return
    printf("Found %d bands\n", nbands)

    fprintf(fout, "%4s \t ", "mass")
    for i in xrange(nbands) : fprintf(fout, " -%-5s  ", names[nbands-i-1]) # for (int i = 0; i < nbands; ++i) fprintf(fout, " -%-5s  ", names[nbands-i-1])
    fprintf(fout, "  %8s  ", "  mid.  ")
    for i in xrange(nbands) : fprintf(fout, "   +%-5s", names[i]) # for (int i = 0; i < nbands; ++i) fprintf(fout, "   +%-5s", names[i])
    fprintf(fout, "\n")
    n = bs[0].GetN()
    for i in xrange(n) : # for (int i = 0, n = bs[0].GetN(); i < n; ++i):
        xi = bs[0].GetX()[i]
        if pdefs.halfint_masses : fprintf(fout, "%5.1f\t", xi)
        else : fprintf(fout, "%4.0f \t", xi)

        for j in xrange(nbands-1, -1,-1): # for (int j = nbands-1; j >= 0; --j):
            ij = findBin(bs[j], xi)
            y = bs[j].GetY()[ij] - bs[j].GetErrorYlow(ij)
            if ij == -1 : y = NAN
            fprintf(fout, "%7.5f  ", y)
        

        fprintf(fout, "  %7.5f  ", bs[0].GetY()[i])

        for j in xrange(nbands) : #for (int j = 0; j < nbands; ++j):
            ij = findBin(bs[j], xi)
            y = bs[j].GetY()[ij] + bs[j].GetErrorYhigh(ij)
            if ij == -1 : y = NAN
            fprintf(fout, "  %7.5f", y)
        

        fprintf(fout, "\n")
    

def printFcBand(bands, who, fileName,  npostfix, *postfixes):
    first =  bands.Get(who+"_"+postfixes[0])
    if (first == 0): print "MISSING " , who , "_" , postfixes[0]; return; 
    fout = fopen(fileName.Data(), "w")
    printFcBand(bands, who, fout, npostfix, postfixes)
    fclose(fout)


def printQuantiles(bands, who, fout):
    quants = [ 0.025, 0.16, 0.5, 0.84, 0.975 ]
    graphs=[]
    for i in xrange(5) : # for (int i = 0; i < 5; ++i):
        graphs[i] =  bands.Get(who+ROOT.TString.Format("_quant%03d", int(1000*quants[i])))
        if (graphs[i] == 0): print "Missing quantile band for p = " , quants[i]; return; 
    
    fprintf(fout, "%4s \t %6s %5s   %6s %5s   %6s %5s   %6s %5s   %6s %5s\n", "mass", "-95%","err", "-68%","err", "median","err", "+68%","err", "+95%","err")
    fprintf(fout, "%4s \t %6s %5s   %6s %5s   %6s %5s   %6s %5s   %6s %5s\n", "-----", "-----", "-----", "-----", "-----", "-----","-----", "-----", "-----", "-----", "-----")
    n = graphs[0].GetN()
    for i in xrange(n) : # for (int i = 0, n = graphs[0].GetN(); i < n; ++i):
        fprintf(fout, "%4d \t ", int(graphs[0].GetX()[i]))
        for i in xrange(5) : # for (int j = 0; j < 5; ++j):
            fprintf(fout, "%6.2f %5.2f   ", graphs[j].GetY()[i], graphs[j].GetErrorYlow(i))
        
        fprintf(fout, "\n")
    

def printQuantiles(bands, who, fileName):
    mean68 =  bands.Get(who+"_quant025")
    if (mean68 == 0): print "MISSING " , who , "_quant025"; return; 
    fout = fopen(fileName.Data(), "w")
    printQuantiles(bands,who,fout)
    fclose(fout)


def printBand(bands, who, fileName, mean=False):
    if mean : mean68 =  bands.Get(who+"_mean")
    else : mean68 =  bands.Get(who+"_median")
    obs  =  bands.Get(who+"_obs")
    if not mean68 and not obs : 
        print "MISSING " , who , "_mean and " , who , "_obs"
        return
    
    fout = open(fileName, "w")
    #printBand(bands,who,fout,mean)
    printBandWithFileDescriptor(bands,who,fout,mean)
    fout.close()




# In[50]:


def importLine(bands, name, fileName):
    input = fopen(fileName, "r")
    if (input == 0): print "Cannot open " , fileName ; return; 
    fclose(input)
    inObs = ROOT.TGraphAsymmErrors(); inObs.SetName(name)
    #float mH, yObs
    #for (int n = 0; fscanf(in,"%f %f", &mH, &yObs) == 2; ++n):
    #    inObs.SetPoint(n, mH, yObs)
    n=0
    with open(input) as f:
        for line in f:
            mH, yObs = line.split()
            mH = float(mH)
            yObs = float(yObs)
            inObs.SetPoint(n, mH, yObs)
            n +=1
    bands.WriteTObject(inObs)
    #fclose(input)
    input.close()


# In[52]:


def importBands(bands, name, fileName, hasObs = False, has95 = True): # ????? needs a test because fgetc, fgets, ungetc
    input = fopen(fileName, "r")
    if (input == 0): print "Cannot open " , fileName; return; 
    fclose(input)
    inObs = ROOT.TGraphAsymmErrors(); inObs.SetName(name+"_obs")
    in68  = ROOT.TGraphAsymmErrors();  in68.SetName(name+"_median")
    in95  = ROOT.TGraphAsymmErrors();  in95.SetName(name+"_median_95")
    #float mH, yObs, yLL, yLo, y, yHi, yHH
    #buff=[] # char buff[1025]
    #while (True): #do {
    #    c = ord(input.read(1)[0]) # fgetc(in)        
    #    if (c == 'm' or c == '-'):
    #        buff = input.read(1024) # fgets(buff,1024,in)
    #    else :
    #        ungetc(c,input)
    #        break
        
    # while(True)
    if (hasObs):
        n=0
        #fscanf(input,"%f %f %f %f %f %f %f", &mH, &yObs, &yLL, &yLo, &y, &yHi, &yHH) == 7 #for (int n = 0; fscanf(input,"%f %f %f %f %f %f %f", &mH, &yObs, &yLL, &yLo, &y, &yHi, &yHH) == 7; ++n):
        with open(input) as f:
            for line in f:
                if len(line.split()) == 7:
                    mH, yObs, yLL, yLo, y, yHi, yHH = line.split()
                    inObs.SetPoint(n, mH, yObs)
                    in68.SetPoint(n, mH, y); in68.SetPointError(n, 0, 0, y-yLo, yHi-y)
                    in95.SetPoint(n, mH, y); in95.SetPointError(n, 0, 0, y-yLL, yHH-y)
                    n += 1
        
    else :
        if (has95):
            n=0
            with open(input) as f:
                for line in f:
                    if len(line.split()) == 6:
                            #for (int n = 0; fscanf(input,"%f %f %f %f %f %f", &mH, &yLL, &yLo, &y, &yHi, &yHH) == 6; ++n):
                            mH, yLL, yLo, y, yHi, yHH = line.split()
                            in68.SetPoint(n, mH, y); in68.SetPointError(n, 0, 0, y-yLo, yHi-y)
                            in95.SetPoint(n, mH, y); in95.SetPointError(n, 0, 0, y-yLL, yHH-y)
                            n += 1
        else :
            #for (int n = 0; fscanf(input,"%f %f %f %f", &mH, &yLo, &y, &yHi) == 4; ++n):
            #    in68.SetPoint(n, mH, y); in68.SetPointError(n, 0, 0, y-yLo, yHi-y)
            
            n=0
            with open(input) as f:
                for line in f:
                    if len(line.split()) == 4:
                            #for (int n = 0; fscanf(input,"%f %f %f %f %f %f", &mH, &yLL, &yLo, &y, &yHi, &yHH) == 6; ++n):
                            mH, yLo, y, yHi = line.split()
                            in68.SetPoint(n, mH, y); in68.SetPointError(n, 0, 0, y-yLo, yHi-y)
                            #in68.SetPoint(n, mH, y); in68.SetPointError(n, 0, 0, y-yLo, yHi-y)
                            #in95.SetPoint(n, mH, y); in95.SetPointError(n, 0, 0, y-yLL, yHH-y)
                            n += 1
        
    
    bands.WriteTObject(in68)
    if (has95) :bands.WriteTObject(in95)
    if (hasObs) :bands.WriteTObject(inObs)
    #fclose(input)
    input.close()


# In[53]:


def importLandS(bands, name, thefile, doObserved=True, doExpected=True):
    if (thefile == 0) : return 
    t = thefile.Get("T")
    if (t == 0): print "TFile " , thefile.GetName() , " does not contain the tree"; return; 
    isML = (name.Index("ml") == 0); isPVal = (name.Index("pval") == 0)
    #Double_t mass, limit, limitErr, rmedian, rm1s, rp1s, rm2s, rp2s

    mass = array('d',[])
    rmedian = array('d',[])
    rm1s = array('d',[0])
    rm2s = array('d',[])
    rp2s = array('d',[])
    rp1s = array('i',[])
    t.SetBranchAddress("mH", mass)
    t.SetBranchAddress("rmedian", rmedian)
    t.SetBranchAddress("rm1s", rm1s)
    t.SetBranchAddress("rm2s", rm2s)
    t.SetBranchAddress("rp2s", rp2s)
    t.SetBranchAddress("rp1s", rp1s)
    what = "limit"
    if (isPVal) :what = "pvalue"
    if (isML)   :what = "rmean"
    print "For " , name , " will read " , what.Data()
    limit = array('d',[])
    limitErr = array('d',[])
    t.SetBranchAddress(what.Data(), limit)
    t.SetBranchAddress("limitErr", limitErr)
    obs       = ROOT.TGraphAsymmErrors(); obs.SetName(name+"_obs");              nobs = 0
    median    = ROOT.TGraphAsymmErrors(); median.SetName(name+"_median");        nmedian = 0
    median_95 = ROOT.TGraphAsymmErrors(); median_95.SetName(name+"_median_95");  nmedian_95 = 0
    n = t.GetEntries()
    for i in xrange(n): # for (size_t i = 0, n = t.GetEntries(); i < n; ++i):
        t.GetEntry(i)
        if (doObserved):
            if (isML):
                obs.Set(nobs+1)
                obs.SetPoint(nobs, mass, limit)
                obs.SetPointError(nobs, 0, 0, limit-rm1s, rp1s-limit)
                nobs +=1
            elif (limit != 0):
                obs.Set(nobs+1)
                obs.SetPoint(nobs, mass, limit)
                obs.SetPointError(nobs, 0, 0, limitErr, limitErr)
                nobs +=1
            
        
        if ( not isML and doExpected):
            if (isPVal): 
                if (limit != 0):
                    median.Set(nmedian+1)
                    median.SetPoint(nmedian, mass, limit)
                    median.SetPointError(nmedian, 0, 0, 0, 0)
                    nmedian += 1
                 
            else:
                if (limit != 0):
                    median.Set(nmedian+1)
                    median.SetPoint(nmedian, mass, rmedian)
                    if (rm1s != 0 and rp1s != 0):
                        median.SetPointError(nmedian, 0, 0, rmedian - rm1s, rp1s - rmedian)
                    else:
                        median.SetPointError(nmedian, 0, 0, 0, 0)
                    
                    nmedian +=1
                
                if (limit != 0 and rm2s != 0 and rp2s != 0):
                    median_95.Set(nmedian_95+1)
                    median_95.SetPoint(nmedian_95, mass, rmedian)
                    median_95.SetPointError(nmedian_95, 0, 0, rmedian - rm2s, rp2s - rmedian)
                    nmedian_95 +=1
                
            
        
    
    if (obs.GetN()): obs.Sort(); bands.WriteTObject(obs); print " imported " , obs.GetName() , " with " , obs.GetN() , " points."  
    if (median.GetN()): median.Sort(); bands.WriteTObject(median); print " imported " , median.GetName() , " with " , median.GetN() , " points."  
    if (median_95.GetN()): median_95.Sort(); bands.WriteTObject(median_95); print " imported " , median_95.GetName() , " with " , median_95.GetN() , " points."  

def importLandS(bands, name, fileName, doObserved=True, doExpected=True):
    input = ROOT.TFile.Open(fileName)
    if (input == 0): print "Cannot open " , fileName ; return;  
    importLandS(bands, name, input, doObserved, doExpected)
    input.Close()


# In[54]:


def smoothWithPolyFit(x,  npar,  n, xi, yi):
    fitRes = polyFit(x, yi[n/2], npar, n, xi, yi)
    return fitRes(0)+yi[n/2]


# In[59]:


def printValueFromScan1D(bands, name, out):
    graph =  bands.Get(name)
    if (graph is None) : return
    if (graph == 0) : return
    x = graph.GetX()
    y = graph.GetY()
    imin = 0; n = graph.GetN()
    for i in xrange(n) : # for (int i = 1; i < n; ++i):
        if (y[i] < y[imin]) : imin = i
    
    t1 = 1; t2 = 3.84
    hi68ok = False; hi95ok = False; lo68ok = False; lo95ok = False
    hi68 = x[n-1]; hi95 = x[n-1]; lo68 = x[0]; lo95 = x[0]
    for i in xrange(n-1) : # for (int i = 0; i < n-1; ++i):
        if (y[i] > t1 and y[i+1] < t1):
            d1 = fabs(y[i] - t1); d2 = fabs(y[i+1] - t1)
            lo68 = (x[i]*d2 + x[i+1]*d1)/(d1+d2); lo68ok = True
        elif (y[i] < t1 and y[i+1] > t1):
            d1 = fabs(y[i] - t1); d2 = fabs(y[i+1] - t1)
            hi68 = (x[i]*d2 + x[i+1]*d1)/(d1+d2); hi68ok = True
        
        if (y[i] > t2 and y[i+1] < t2):
            d1 = fabs(y[i] - t2); d2 = fabs(y[i+1] - t2)
            lo95 = (x[i]*d2 + x[i+1]*d1)/(d1+d2); lo95ok = True
        elif (y[i] < t2 and y[i+1] > t2):
            d1 = fabs(y[i] - t2); d2 = fabs(y[i+1] - t2)
            hi95 = (x[i]*d2 + x[i+1]*d1)/(d1+d2); hi95ok = True
        
    
    log = fopen(out.Data(), "w")
    fprintf(log, "Lowest point :  % 8.4f \n", x[imin])
    if (lo68ok) :fprintf(log, "Crossing at 1.00 from left:  % 8.4f \n", lo68)
    if (hi68ok) :fprintf(log, "Crossing at 1.00 from right: % 8.4f \n", hi68)
    if (lo95ok) :fprintf(log, "Crossing at 3.84 from left:  % 8.4f \n", lo95)
    if (hi95ok) :fprintf(log, "Crossing at 3.84 from right: % 8.4f \n", hi95)
    fclose(log)


# In[55]:


#def array_sort(*begin, *end): std::sort(begin, end); 
#def array_sort(float *begin, float *end): std::sort(begin, end); 
#def array_sort(int *begin, int *end): std::sort(begin, end); 
#def array_sort(&begin, &end): std::sort(&begin, &end); 
#def array_sort(float &begin, float &end): std::sort(&begin, &end); 
#def array_sort(int &begin, int &end): std::sort(&begin, &end); 
def array_sort (thelist) :
    thelist.sort()
    return thelist
def bandUtils(): pass

