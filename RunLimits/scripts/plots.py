
# coding: utf-8

# In[1]:


# ported from https://alachua.ihepa.ufl.edu:9999/tree/raid8/bockjoo/optimize/remote/Run2/Submit/limits/HZZ4l_Combination/CreateDatacards/CMSSW_6_1_1/src/HCG/makeBands_hzz4l.cc


# In[2]:


from ROOT import *
import ROOT
import os
#get_ipython().magic(u'jsroot on')
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
import time # for sleep

# In[3]:

import pdefs

from bandUtils import *


# In[4]:

# globals used to be here

# In[5]:


def drawExclusions(xmin0, xmax0, ymin, ymax):
    if (pdefs.cms_excluded):
        xmin = 127.5
        xmax = ROOT.TMath.Min(xmax0,600.)
        box = ROOT.TBox(xmin, ymin, xmax, ymax)
        box.SetLineStyle(0)
        box.SetFillStyle(3345)
        box.SetFillColor(95)
        box.DrawClone()
        line = ROOT.TLine(xmin, ymin, xmin, ymax)
        line.SetLineWidth(2)
        line.SetLineColor(95)
        line.DrawClone()
        if (xmax0 > 600) : line.DrawLine(xmax, ymin, xmax, ymax)
        if (pdefs.fakeCMS == 0):
            pdefs.fakeCMS = ROOT.TH1F("fake_CMS","fake_CMS",1,0.,1.)
            pdefs.fakeCMS.SetLineStyle(0)
            pdefs.fakeCMS.SetFillStyle(3345)
            pdefs.fakeCMS.SetFillColor(95)
            pdefs.fakeCMS.SetLineWidth(2)
            pdefs.fakeCMS.SetLineColor(95)
        
    
    if (pdefs.tev_excluded and xmax0 > 156):
        xmin = 156
        xmax = ROOT.TMath.Min(177., xmax0)
        box = ROOT.TBox(xmin, ymin, xmax, ymax)
        box.SetLineStyle(0)
        box.SetFillStyle(3354)
        box.SetFillColor(213)
        box.DrawClone()
        line = ROOT.TLine(xmin, ymin, xmin, ymax)
        line.SetLineWidth(2)
        line.SetLineColor(213)
        line.DrawClone()
        if (xmax0 > 177) : line.DrawLine(xmax, ymin, xmax, ymax)
        xmin = 100 ; xmax = 109
        box2 = ROOT.TBox(xmin, ymin, xmax, ymax)
        box2.SetLineStyle(0)
        box2.SetFillStyle(3354)
        box2.SetFillColor(213)
        if (xmin0 < 109 and pdefs.pdefs.tev_excluded_alsolow) : box2.DrawClone()
        line2 = ROOT.TLine(xmin, ymin, xmin, ymax)
        line2.SetLineWidth(2)
        line2.SetLineColor(213)
        if (xmin0 < 109 and pdefs.pdefs.tev_excluded_alsolow) : line2.DrawClone()
        if (xmin0 < 109 and pdefs.pdefs.tev_excluded_alsolow) : line2.DrawLine(xmax, ymin, xmax, ymax)
        if (pdefs.fakeTEV == 0):
            pdefs.fakeTEV = ROOT.TH1F("fake_tev","fake_tev",1,0.,1.)
            pdefs.fakeTEV.SetLineStyle(0)
            pdefs.fakeTEV.SetFillStyle(3354)
            pdefs.fakeTEV.SetFillColor(213)
            pdefs.fakeTEV.SetLineWidth(2)
            pdefs.fakeTEV.SetLineColor(213)
     
    if (pdefs.lep_excluded):
        box=ROOT.TBox(xmin0, ymin, 114.4, ymax)
        box.SetLineStyle(0)
        box.SetFillStyle(3345)
        box.SetFillColor(209)
        box.DrawClone()
        line=ROOT.TLine(114.4, ymin, 114.4, ymax)
        line.SetLineWidth(2)
        line.SetLineColor(209)
        line.DrawClone()
        if (pdefs.fakeLEP == 0):
            pdefs.fakeLEP = ROOT.TH1F("fake_lep","fake_lep",1,0.,1.)
            pdefs.fakeLEP.SetLineStyle(0)
            pdefs.fakeLEP.SetFillStyle(3345)
            pdefs.fakeLEP.SetFillColor(209)
            pdefs.fakeLEP.SetLineWidth(2)
            pdefs.fakeLEP.SetLineColor(209)


# In[6]:


def setCanvas(first, title, min=0.1, max=30, ytitle="Limit (#sigma_{95%}/#sigma_{SM})"): 
    if (first is None) : print "Error for ", title
    first.SetTitle(title)
    first.GetYaxis().SetTitle(ytitle)
    pdefs.c1.SetLogy(min > 0 or ("pval" in title))    
    first.GetYaxis().SetRangeUser(min,max)
    if (min <= 0) : first.GetYaxis().SetDecimals(True)
    first.GetYaxis().SetTitleOffset(1.00+0.2*pdefs.isSquareCanvas)
    first.GetXaxis().SetTitle("Higgs boson mass ("+pdefs.massUnits+")")
    if (pdefs.SM == "SM4"):
        first.GetXaxis().SetTitle("SM4 Higgs boson mass ("+pdefs.massUnits+")")
    elif (pdefs.SM == "FP"):
        first.GetXaxis().SetTitle("FP Higgs boson mass ("+pdefs.massUnits+")")
    
    pdefs.c1.SetTickx(1)    
    pdefs.c1.SetTicky(1)    
    xmin = first.GetXaxis().GetXmin()
    xmax = first.GetXaxis().GetXmax() 
    drawExclusions(xmin, xmax, min, max)
    if "p-value" in ytitle:
        #//ROOT.TBoxtwoSig(xmin, 1, xmax, ROOT::Math::normal_cdf_c(2)) twoSig.SetFillColor(220); twoSig.DrawClone();
        #//ROOT.TBoxoneSig(xmin, 1, xmax, ROOT::Math::normal_cdf_c(1)) oneSig.SetFillColor(211); oneSig.DrawClone();
        threeSig=ROOT.TLine(xmin, ROOT.Math.normal_cdf_c(3), xmax, ROOT.Math.normal_cdf_c(3)) 
        threeSig.SetLineColor(pdefs.lineAt1Color)
        threeSig.SetLineWidth(2)
        if pdefs.c1.GetGridy():
            threeSig.SetLineStyle(7)
            protrude = 1
        else: 
            threeSig.SetLineStyle(pdefs.lineAt1Style)
            protrude = 0
        latex = ROOT.TLatex() ; latex.SetTextFont(42); latex.SetTextSize(0.045); latex.SetTextColor(2);
        #protrude = pdefs.c1.GetGridy() 
        for i in xrange(9):
            if (min >  ROOT.Math.normal_cdf_c(i)) : break
            threeSig.DrawLine(xmin, ROOT.Math.normal_cdf_c(i), xmax+(xmax-xmin)*0.015*protrude, ROOT.Math.normal_cdf_c(i))
            sspam = ROOT.TString.Format("%d#sigma", i)
            latex.DrawLatex(xmax+(xmax-xmin)*0.01, ROOT.Math.normal_cdf_c(i)*1.1, sspam)
        
    
    if "CL_{S} of" in ytitle:
        xmin = first.GetXaxis().GetXmin()
        xmax = first.GetXaxis().GetXmax() 
        #//if (max > 1):
        #//    ROOT.TLineoneLine(xmin,1,xmax,1) oneLine.SetLineWidth(2);
        #//    oneLine.DrawClone()
        #//}
        lines = [ 1-0.90, 1-0.95, 1-0.99 ]
        hline = ROOT.TLine(xmin, ROOT.Math.normal_cdf_c(3), xmax, ROOT.Math.normal_cdf_c(3)) 
        hline.SetLineColor(pdefs.lineAt1Color) 
        latex = ROOT.TLatex()
        latex.SetTextFont(42); latex.SetTextSize(0.045 - 0.008*pdefs.isSquareCanvas); latex.SetTextColor(2);
        for i in xrange(3):
            hline.SetLineWidth(4) ; hline.SetLineStyle(pdefs.lineAt1Style);
            #//if (i == 1): hline.SetLineWidth(4) hline.SetLineStyle(1); }
            #//else        { hline.SetLineWidth(2) hline.SetLineStyle(2); }
            protrude = 0
            if pdefs.c1.GetGridy():
                protrude = 1 
            hline.DrawLine(xmin, lines[i], xmax+(xmax-xmin)*0.02*protrude, lines[i])
            sspam = ROOT.TString.Format("%.0f%%", 100*(1-lines[i]))
            latex.DrawLatex(xmax+(xmax-xmin)*0.01, lines[i]*1.1, sspam)
        
    




# In[7]:


def newLegend(x1, y1, x2, y2):
    pdefs.leg = ROOT.TLegend(x1,y1,x2,y2) 
    pdefs.leg.SetFillColor(0)
    pdefs.leg.SetShadowColor(0)
    pdefs.leg.SetTextFont(42)
    pdefs.leg.SetTextSize(0.05)
    return pdefs.leg


# In[8]:


def spam(text="", x1=0.17, y1=0.89, x2=0.58, y2=0.94, textAlign=-1, fill=True, fontSize=0):
   if "#splitline" in str(text): 
       if (y1 > 0.5) : y1 -= 0.065+0.04*pdefs.isTiny
       else : y2 += 0.065+0.04*pdefs.isTiny 
   elif "\n" in str(text):
      buff = str(text)
      stride = 0.0475
      if (y1 > 0.5): stride = -0.0475
      if (fontSize):
        if pdefs.isSquareCanvas: 
            stride *= fontSize/0.04
        else:
            stride *= fontSize/0.045
      lines = 0
      buff = str(text)
      #for (char *line = strtok(buff,"\n") line != 0; line = strtok((char*)0,"\n")):
      #  lines++
      lines = buff.count('\n') # lines = len(buff.split('\n')) # lines = buff(x.splitlines())

      if (stride > 0):
            lines += -1
            y1 += lines*stride
            y2 += lines*stride
            stride *= -1
      buff =str(text)
      #lines = buff.count('\n')
      lines = buff.split('\n')
      #for i in xrange(lines): #(char *line = strtok(buff,"\n") line != 0; line = strtok((char*)0,"\n")):
      for line in lines:
         if textAlign==-1 :
            spam(line,x1,y1,x2,y2,22,fill,fontSize)
         else:
            spam(line,x1,y1,x2,y2,textAlign,fill,fontSize)
         y1 += stride
         y2 += stride;
      
      return
   
   if (textAlign == -1) : textAlign=12
   #//cmsprel = new TPaveText(x1,y1,x2,y2,"brtlNDC")
   cmsprel = ROOT.TPaveText(x1,y1,x2,y2,"brtlNDC")
   if (fontSize == 0) :
    val1 = 0.045
    if pdefs.isSquareCanvas : val1 = 0.04
    val2 = 0
    if pdefs.isTiny: val2 = 0.02
    fontSize = val1 + val2
   cmsprel.SetTextSize(fontSize)
   cmsprel.SetFillColor(0)
   if (fill or text == "") :
     cmsprel.SetFillStyle(1001)
   else :
     cmsprel.SetFillStyle(0)
   cmsprel.SetLineStyle(0)
   cmsprel.SetLineColor(0)
   cmsprel.SetLineWidth(0)
   cmsprel.SetTextAlign(textAlign)
   cmsprel.SetTextFont(42)
   #//cmsprel.SetTextColor(text.Contains("PRIVATE") ? 205 : 1)
   #cmsprel.SetTextColor( TString(text ? text : pdefs.SPAM).Contains("PRIVATE") ? 205 : 1 )
   if text == "":
      cmsprel.SetTextColor(1)
   else:
      if "PRIVATE" in pdefs.SPAM:
          cmsprel.SetTextColor(205)
      else:
          cmsprel.SetTextColor(1)
   if text =="" :
     cmsprel.AddText(text)
   else:
     cmsprel.AddText(pdefs.SPAM)
   cmsprel.SetBorderSize(0)
   cmsprel.Draw("same")



# In[9]:


def finalizeNoSave(name, xmin, xmax, ymin, ymax, tspam=0, spamLow=False):
    line=ROOT.TLine(xmin,1,xmax,1) 
    line.SetLineColor(pdefs.lineAt1Color) 
    line.SetLineStyle(pdefs.lineAt1Style) 
    line.SetLineWidth(4)
    if "smcls" in name: 
        line.SetY1(0.05) ;  line.SetY2(0.05); 
        line.DrawClone()
        #//line.SetLineStyle(7) line.SetLineWidth(2);
        #//line.DrawLine(xmin,1-0.68,xmax,1-0.68)
        line.DrawLine(xmin,1-0.90,xmax,1-0.90)
        line.DrawLine(xmin,0.01,  xmax,0.01)
        #//line.DrawLine(xmin,0.003, xmax,0.003)
    elif "pval" in name:
        if "_band" in name:
            #ROOT.TLinethreeSig(xmin, ROOT::Math::normal_cdf_c(3), xmax, ROOT::Math::normal_cdf_c(3)) 
            #threeSig.SetLineColor(pdefs.lineAt1Color) threeSig.SetLineWidth(2); 
            #threeSig.SetLineStyle(pdefs.c1.GetGridy() ? 7 : pdefs.lineAt1Style); 
            threeSig=ROOT.TLine(xmin, ROOT.Math.normal_cdf_c(3), xmax, ROOT.Math.normal_cdf_c(3)) 
            threeSig.SetLineColor(pdefs.lineAt1Color)
            threeSig.SetLineWidth(2)
            if pdefs.c1.GetGridy():
               threeSig.SetLineStyle(7)
               protrude = 1
            else: 
               threeSig.SetLineStyle(pdefs.lineAt1Style)
               protrude = 0

            latex = ROOT.TLatex() ; latex.SetTextFont(42); latex.SetTextSize(0.045); latex.SetTextColor(2);
            #protrude = pdefs.c1.GetGridy() 
            for z in xrange(9):
                if (ymin >  ROOT.Math.normal_cdf_c(z)) : break
                threeSig.DrawLine(xmin, ROOT.Math.normal_cdf_c(z), xmax+(xmax-xmin)*0.015*protrude, ROOT.Math.normal_cdf_c(z))
                sspam = ROOT.TString.Format("%d#sigma", z)
                latex.DrawLatex(xmax+(xmax-xmin)*0.01, ROOT.Math.normal_cdf_c(z)*1.1, sspam)
    else:
        #print "DEBUG finalizeNoSave will line.DrawClone() "   
        #print "DEBUG line=ROOT.TLine(xmin,1,xmax,1) ",xmin," ",xmax
        # 
        line.DrawClone() # bockjoo The think red line
    #print "DEBUG finalizeNoSave will return after line.DrawClone() " 
    #return
    if (ROOT.gPad.GetLogx() and xmin <= 100 and xmax >= 600):    
        tick = ROOT.TLine() ; tick.SetLineWidth(1); tick.SetLineColor(1);
        dyh = ymax * 0.08
        dyl = ymin * 0.08 #//fabs(pdefs.c1.PixeltoY(pdefs.c1.VtoPixel(0.95)) - pdefs.c1.PixeltoY(pdefs.c1.VtoPixel(0.94)));
        if (ROOT.gPad.GetLogy() and math.log(ymax/ymin) > math.log(1000000)):
            dyh *= 2
            dyl *= 2
        if (ROOT.gPad.GetLogy() == 0):
            dyl = 0.01*(ymax-ymin)
            dyh = dyl
        if (pdefs.isTiny):
            dyh *= 2
            dyl *= 2
        for i in xrange(600): # ) (int i = 100 i < 600; i += 10)  {
            if i < 100 : continue
            if i % 10 != 0 : continue
            if (i > 400 and i % 20 == 10) : continue
            if i % 100 == 0 :
               tick.DrawLine(i, ymin, i, ymin+2*dyl)
            else:
               tick.DrawLine(i, ymin, i, ymin+dyl) 
            if i % 100 == 0 :
               tick.DrawLine(i, ymax, i, ymax+2*dyh)
            else:
               tick.DrawLine(i, ymax, i, ymax+dyh) 
        
    
    if (pdefs.leg) : pdefs.leg.Draw()
    if "pval" in name:
        print " name contains pval "
        if "_all" in name:
            print " name contains _all "
            #//spam("#splitline{Look-elsewhere effect}{not included}", 0.17 - 0.05*pdefs.isTiny, 0.21 + 0.10*pdefs.isTiny, 0.47 - 0.05*pdefs.isTiny, 0.26 + 0.10*pdefs.isTiny)
            isZoom = ("zoom" in name)
            if (pdefs.isSquareCanvas):
                print " pdefs.isSquareCanvas "
                if (pdefs.doLEESPAM):
                   #//if (pdefs.SM == "SM") spam(pdefs.LEESPAM_2L, 0.19, 0.35, 0.59, 0.40, 22, True)
                   if (pdefs.SM == "SM") : spam(pdefs.LEESPAM_3L, 0.19, 0.35, 0.59, 0.40, 22, True, 0.034)
                   if (pdefs.SM == "FP") : spam(pdefs.LEESPAM_2L, 0.19, 0.35, 0.59, 0.40, 22, True)
                   if (pdefs.SM == "SM4") : spam(pdefs.LEESPAM_2L, 0.19, 0.35, 0.59, 0.40, 22, True)
                
                #spam(tspam, .17, .15-0.00*(!pdefs.doLEESPAM), .49, .20-0.00*(!pdefs.doLEESPAM), 22)
                spam(tspam, .17, .15-0.00, .49, .20-0.00, 22)
            else:
                #//if (pdefs.doLEESPAM) spam(pdefs.LEESPAM_1L, 0.175 - 0.05*pdefs.isTiny, 0.15 + 0.10*pdefs.isTiny, 0.59 - 0.05*pdefs.isTiny, 0.20 + 0.10*pdefs.isTiny)
                #//if (pdefs.doLEESPAM) spam(pdefs.LEESPAM_1L, 0.175, 0.35 , 0.59, 0.40)
                if (pdefs.SM == "SM") : spam(pdefs.LEESPAM_3L, 0.19, 0.35, 0.51, 0.40, 22, True, 0.04)
                if (pdefs.SM == "FP") : spam(pdefs.LEESPAM_2L, 0.19, 0.35, 0.51, 0.40, 22, True)
                if (pdefs.SM == "SM4") : spam(pdefs.LEESPAM_2L, 0.19, 0.35, 0.51, 0.40, 22, True)
                #spam(tspam, .17, .15-0.000*(!pdefs.doLEESPAM), .59, .20-0.000*(!pdefs.doLEESPAM), 22)
                spam(tspam, .17, .15-0.000, .59, .20-0.000, 22)
            
            tspam = 0
        elif (tspam):
            print " tspam "
            if (pdefs.isTiny and pdefs.isSquareCanvas):
                #spam(tspam, 0.17, 0.18-0.05*(!pdefs.doLEESPAM), 0.57, 0.21-0.05*(!pdefs.doLEESPAM))
                spam(tspam, 0.17, 0.18-0.05, 0.57, 0.21-0.05)
                if (pdefs.doLEESPAM) :
                    istiny=0.0
                    if pdefs.isTiny: istiny =1.0
                    spam("Interpretation requires look-elsewhere effect correction", 0.17 + 0.05*istiny, 0.18 - 0.12*istiny, 0.94 - 0.12*istiny, 0.23 - 0.10*istiny, 22)
                print " tspam1 "
            elif "#splitline" in tspam:
                print " tspam2 "
                #spam(tspam, 0.17, 0.22-0.05*(!pdefs.doLEESPAM), 0.57, 0.25-0.05*(!pdefs.doLEESPAM))
                spam(tspam, 0.17, 0.22-0.05, 0.57, 0.25-0.05)
                if (pdefs.doLEESPAM) : spam("#splitline{Interpretation requires look-}{elsewhere effect correction}", 0.55, 0.34, 0.93, 0.38)
            else:
                print " tspam3 "
                #//bockjoo spam(tspam, 0.17 + 0.00*pdefs.isTiny, 0.23 - 0.05*pdefs.isTiny-0.05*(!pdefs.doLEESPAM), 0.94 - 0.12*pdefs.isTiny, 0.29 - 0.05*pdefs.isTiny, 22)
                istiny=0.0
                if pdefs.isTiny: istiny =1.0
                if (pdefs.doLEESPAM) : spam("Interpretation requires look-elsewhere effect correction", 0.17 + 0.00*istiny, 0.18 - 0.12*istiny, 0.94 - 0.12*istiny, 0.23 - 0.10*istiny, 22)
                if "_combe" in name:
                    spam("#splitline{Using event-by-event mass}{resolution in H #rightarrow ZZ #rightarrow 4l}", 0.55, 0.34, 0.93, 0.38)
                
            
            tspam = 0
        
    
   
    isFullComb = ("_all" in name ) or (channelFromName(name,False,True) == "Combined")
    if (tspam) :
        istiny=0.0
        if pdefs.isTiny: istiny =1.0
        spamlow = 0.0
        if spamLow : spamlow =1.0
        val1=0.0
        if "\n" in tspam: val1=1.0
        val2=0.64
        if pdefs.isSquareCanvas: val2=0.62
        val3=0.0
        if "pvala" in name or "smcls" in name: val3=-0.07
        val4=0.0
        if isFullComb: val4=1.0
        spam(tspam, 0.17 - 0.05*istiny, 0.89-0.74*spamlow + 0.04*istiny-0.005*val1,val2+val3- 0.15*val4-0.005*val1, 0.94-0.74*spamlow + 0.04*istiny)




# In[10]:


def justSave(name, xmin=0., xmax=0., ymin=0., ymax=0., tspam=0, spamLow=False):
    print "justSave saved : ",pdefs.globalPrefix+name+".eps"
    pdefs.c1.Print(pdefs.globalPrefix+name+".eps")
    #//pdefs.c1.Print(pdefs.globalPrefix+name+".png")
    #//gSystem.Exec("convert "+pdefs.globalPrefix+name+".eps "+pdefs.globalPrefix+name+".png")
    convOpt = "-q  -dBATCH -dSAFER  -dNOPAUSE  -dAlignToPixels=0 -dEPSCrop  -dPrinted -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sDEVICE=png16m"
    convCmd = Form("gs %s -sOutputFile=%s.png -q \"%s.eps\" -c showpage -c quit" % ( convOpt, pdefs.globalPrefix+name,pdefs.globalPrefix+name))
    ROOT.gSystem.Exec(convCmd)
    print "justSave saved : ",pdefs.globalPrefix+name+".png"
    if "pval_ml_" in name:
        ROOT.gSystem.Exec("epstopdf "+pdefs.globalPrefix+name+".eps --outfile="+pdefs.globalPrefix+name+".pdf")
        print "justSave saved : ",pdefs.globalPrefix+name+".pdf"
    else:
        pdefs.c1.Print(pdefs.globalPrefix+name+".pdf")
        print "justSave saved : ",pdefs.globalPrefix+name+".pdf"


def finalize(name, xmin, xmax, ymin, ymax, tspam=0, spamLow=False):
    finalizeNoSave(name,xmin, xmax, ymin, ymax, tspam, spamLow)
    justSave(name,xmin, xmax, ymin, ymax, tspam, spamLow)
    #pdefs.c1.SaveAs(pdefs.globalPrefix+name+".png") # bockjoo



# In[11]:



def slidingWindowAverage(input, slidingWindow):
    out  = input.Clone()
    isLogPlot = "smcls" in input.GetName()
    #print "DEBUG slidingWindowAverage isLogPlot " , isLogPlot
    #for (int i = 0, n = input.GetN() i < n; ++i):
    n=input.GetN()
    for i in xrange(n):
        y0 = input.GetY()[i]
        sum = 0 ; sumhi = 0 ; sumlo = 0 ; sumw = 0
        #print "DEBUG slidingWindowAverage y0 " , y0
        for j in xrange(i-slidingWindow, i+slidingWindow+1):
            #print "DEBUG slidingWindowAverage j " , j
            if (j < 0 or j >= n) : continue
            y = input.GetY()[j] ; w = 1.0 #// /(1.0 + abs(i-j));
            if (isLogPlot):
                #print "DEBUG isLogPlot next check ",abs(log(y0/y)), " > ",log(2)
                if (abs(log(y0/y)) > log(2)) : continue    
            else:
                #print "DEBUG not isLogPlot next check " , abs(y0 - y) ,  " > " , 0.1*y0
                if (abs(y0-y) > 0.1*y0) : continue    
            
            if (isLogPlot):
                sum   += w*log(y)
                sumlo += w*log(input.GetEYlow()[j])
                sumhi += w*log(input.GetEYhigh()[j])
            else:
                sum   += w*y
                sumlo += w*input.GetEYlow()[j]
                sumhi += w*input.GetEYhigh()[j]
            
            sumw  += w
        
        if (sumw == 0) : continue
        if (isLogPlot):
            out.GetY()[i] = exp(sum/sumw)
            out.GetEYlow()[i] = exp(sumlo/sumw)
            out.GetEYhigh()[i] = exp(sumhi/sumw)
        else:
            out.GetY()[i] = sum/sumw
            out.GetEYlow()[i] = sumlo/sumw
            out.GetEYhigh()[i] = sumhi/sumw
        
    
    return out


# In[12]:


def smoothSMCLs(input, slidingWindow, order):
    
    who = input.GetName()
    nonzero   = True #//TString(input.GetName()).Contains("cls");
    isLogPlot = ("smcls" in who or "pval" in who)
    absolute  = ("smcls" in who or "pval" in who)
    maxdel  = log(1.25)
    if ("smcls" in who or "pval" in who) : maxdel = 0 
    inp = input.Clone()
    out = input.Clone()
    if (absolute):
        n = input.GetN()
        for i in xrange(n):
            inp.GetEYlow() [i] = inp.GetY()[i] - inp.GetEYlow()[i]
            inp.GetEYhigh()[i] = inp.GetY()[i] + inp.GetEYhigh()[i]
            out.GetEYlow() [i] = out.GetY()[i] - out.GetEYlow()[i]
            out.GetEYhigh()[i] = out.GetY()[i] + out.GetEYhigh()[i]
        
    
    #for (int which = -1 which <= +1; ++which):
    for which in xrange(-1,2): 
        inY = inp.GetEYlow()
        if which > 0 : inY = inp.GetEYhigh()
        if which == 0 : inY = inp.GetY()
        #ouY = (which == 0 ? out.GetY() : (which > 0 ? out.GetEYhigh() : out.GetEYlow()))
        outY = out.GetEYlow()
        if which > 0 : outY = out.GetEYhigh()
        if which == 0 : outY = out.GetY()
        n=input.GetN()    
        for i in xrange(n) : #  = 0, n = input.GetN() i < n; ++i):
            y0 = inY[i] ; x0 = input.GetX()[i] ; points = 0
            if (nonzero and y0 == 0) : continue
            #for (int j = i-slidingWindow j <= i+slidingWindow; ++j):
            for j in xrange(i-slidingWindow, i+slidingWindow+1):
                if (j < 0 or j >= n) : continue
                if (nonzero and inY[j] == 0) : continue
                if (maxdel > 0):
                    delt = inp.GetY()[j]/inp.GetY()[i] 
                    if (abs(log(delt)) > maxdel) : continue
                
                xi[points] = input.GetX()[j]
                yi[points] = inY[j]
                if isLogPlot : yi[points] = log(inY[j]/y0)
                points += 1
            
            ynew = y0
            if isLogPlot: ynew = 1
            if (points > order+2):
                ynew = smoothWithPolyFit(x0, order+1, points, xi, yi)
            elif (points > 2):
                ynew = smoothWithPolyFit(x0, 1, points, xi, yi)
            elif (maxdel == 0 and nonzero and points == 1):
                ouY[i] = 0
                continue; #// kill the blip
            else : continue #// nothing to do
            ouY[i] = ynew
            if isLogPlot : ouY[i] = y0 * exp(ynew)
        
    
    if (absolute):
        n = input.GetN() # for (int i = 0, n = input.GetN() i < n; ++i):
        for i in xrange(n):
            out.GetEYlow() [i] = +out.GetY()[i] - out.GetEYlow()[i]
            out.GetEYhigh()[i] = -out.GetY()[i] + out.GetEYhigh()[i]
        
    
    #delete inp
    return out


# In[13]:


def removeGlitches(out):
    while True:
        bad = -1 ; hasgood = False
        for i in xrange(out.GetN()): # for (int i = 0 i < out.GetN(); ++i):
            if (out.GetEYlow()[i] == 0 and out.GetEYhigh()[i] == 0):
                bad = i
            else: 
                hasgood = True 
        if not hasgood: return out
        if (bad == -1) :return out
        out.RemovePoint(bad) 


# In[14]:


def loadMasses(file="masses.txt"):
    #nmasses = 0
    masses = []
    f = open(file, 'r')
    #header1 = f.readline()
    #header2 = f.readline()
    #header3 = f.readline()
    for line in f:       
       masses.append(line.strip())
       #nmasses += 1
       #print masses[nmasses-1]
    return masses


# In[15]:


def missingPoints(points):
    import math
    #if not pdefs.track_missing : print "DEBUG missingPoints track_missing is not true returns"
    #if not points : print "DEBUG missingPoints points == 0"
    if not pdefs.track_missing : return 0
    if not points : print "points missing"; return 0
    n = points.GetN() ; xj = points.GetX() ; yj = points.GetY()
    xmin = xj[0] ; xmax = xj[n-1]
    logint = True
    ret = ROOT.TGraphAsymmErrors()
    ret.SetName("missing_"+points.GetName())
    #// check if we have half-integer points
    halfint = False
    for i in xrange(n): # for (int i = 0 i < n; ++i): 
        if (points.GetX()[i] - math.floor(points.GetX()[i]) > 0.4):
            halfint = True
            break
     
    masses = loadMasses()
    nmasses = len(masses)
    nmiss = 0
    for i in xrange(nmasses): # for (int i = 0 i < nmasses; ++i):
        x = masses[i]
        if (masses[i] < xmin or masses[i] > xmax) : continue
        if ((x - math.floor(x)) > 0.4 and not halfint) : continue
        found = False ; xlo = 0 ; ylo = 0 ; xhi = 0 ; yhi = 0
        for j in xrange(n): # for (int j = 0 j < n; ++j):
            if (xj[j] < x):
                xlo = xj[j]
                ylo = yj[j]
            elif (xj[j] == x):
                found = True
                break
            elif (xhi == 0):
                xhi = xj[j] ; yhi = yj[j]
                break
        
        if (not found):
            y = ( yhi * (x-xlo) + ylo * (xhi - x) ) / (xhi - xlo)
            if (yhi > 0 and ylo > 0 and logint):
                y = exp( ( log(yhi) * (x-xlo) + log(ylo) * (xhi - x) ) / (xhi - xlo) )
            
            #//printf("Interpolated missing point %3d from (%4.0f, %7.3f) + (%4.0f, %7.3f) -. (%4.0f, %7.3f)\n", masses[i], xlo, ylo, xhi, yhi, x, y)
            nmiss += 1
            ret.Set(nmiss)
            ret.SetPoint(nmiss-1, x, y)
        
    
    ret.SetMarkerStyle(20)
    ret.SetMarkerSize(0.7)
    ret.SetMarkerColor(100)
    if (nmiss == 0): return 0
    return ret



# In[16]:


def draw2(who, fillColor68, fillColor95, lineColor, same=True, mean=False, slidingWindow=0, smoothorder=2):
    print "DEBUG draw2 1"
    mean68 = ROOT.gROOT.FindObject(who+"_median")
    if mean: mean68 = ROOT.gROOT.FindObject(who+"_mean")
    mean95 = ROOT.gROOT.FindObject(who+"_median"+"_95")
    if mean: mean95 = ROOT.gROOT.FindObject(who+"_mean"+"_95")
    if not mean68:
        if mean:
            print "MISSING ", who, "_mean"
        else:
            print "MISSING ", who, "_median"        
        return 0
    if not mean95:
        if mean:
            print "MISSING ", who, "_mean_95"
        else:
            print "MISSING ", who, "_median_95"
        return 0
    print "DEBUG draw2 2"
    mean68 = removeGlitches(mean68)
    mean95 = removeGlitches(mean95)
    if (mean68.GetN() == 1):
        mean68.SetPointError(0, 10, 10, mean68.GetErrorYlow(0), mean68.GetErrorYhigh(0))
        mean95.SetPointError(0, 10, 10, mean95.GetErrorYlow(0), mean95.GetErrorYhigh(0))
    
    meanL = mean68.Clone() # TGraphAsymmErrors *meanL = (TGraphAsymmErrors*) mean68->Clone();
    for i in xrange(meanL.GetN()): # for (int i= 0 i < meanL.GetN(); ++i): meanL.SetPointError(i, 0,0,0,0); }
        #print "DEBUG draw2 3 meanL.SetPointError ",i
        meanL.SetPointError(i, 0,0,0,0)
        
    print "DEBUG draw2 3.1 slidingWindow ",slidingWindow, " mean68.GetN() ",mean68.GetN() 
    if (slidingWindow != 0 and mean68.GetN() > 5):
        print "DEBUG draw2 4 "
        if (slidingWindow > 0):
            print "DEBUG draw2 5 "
            mean68 = slidingWindowAverage(mean68, slidingWindow)
            mean95 = slidingWindowAverage(mean95, slidingWindow)
            meanL  = slidingWindowAverage(meanL,  slidingWindow)
            print "DEBUG draw2 meanL 0",meanL.GetY()[0], " N = ",meanL.GetN()
        else:
            print "DEBUG draw2 6 "
            mean68 = smoothSMCLs(mean68, -slidingWindow, smoothorder)
            mean95 = smoothSMCLs(mean95, -slidingWindow, smoothorder)
            meanL  = smoothSMCLs(meanL,  -slidingWindow, smoothorder)        

    meanL.SetLineColor(lineColor)
    mean68.SetLineColor(lineColor)
    mean95.SetLineColor(lineColor);
    meanL.SetMarkerColor(lineColor)
    mean68.SetMarkerColor(lineColor)
    mean95.SetMarkerColor(lineColor);
    meanL.SetLineWidth(3)
    mean68.SetLineWidth(3)
    mean95.SetLineWidth(3)
    meanL.SetMarkerSize(1.6)
    mean68.SetMarkerSize(1.6)
    mean95.SetMarkerSize(1.6)
    meanL.SetLineStyle(7)
    mean68.SetLineStyle(7)
    mean95.SetLineStyle(7) 
    mean68.SetLineColor(fillColor68)
    mean95.SetLineColor(fillColor95)
    mean68.SetLineWidth(1) 
    mean95.SetLineWidth(1)
    mean68.SetFillColor(fillColor68)  
    mean95.SetFillColor(fillColor95)
    if same: print "drawing mean95.Draw(E3 SAME)"
    else: print "drawing mean95.Draw(AE3)"
    if same: mean95.Draw("E3 SAME")
    else: mean95.Draw("AE3")
    print "drawing mean68"
    mean68.Draw("E3 SAME")
    print "drawing meanL"
    meanL.Draw("LX SAME")
    return mean95


# In[17]:


def draw1(who, lineColor, option="L"):
    g = ROOT.gROOT.FindObject(who)
    if g is None:
        print "Graph ", who, " empty."
        return 0
    if g.GetN() == 0:
        print "Graph ", who, " empty."
        return 0
    g.SetLineColor(lineColor)
    g.SetMarkerColor(lineColor)
    g.SetLineWidth(3) 
    g.Draw(option+" SAME")
    return g


# In[18]:


def minMaxY(a, ymin=0.6, ymax=12, xmax=999, hardymin=-999):
    if (hardymin == -999) :
        hardymin = 0.02 
        if pdefs.SM == "SM" : hardymin = 0.08
    if not a : 
       print " a is not defined "
       return -9999.0,9999.0
    #if a is None : print " a is None ",a.GetName() ; return -9999,9999
    #ymin = 0.6
    #ymax = 12
    if hardymin < 0 : ymax = 1
    yavg = 1.; npoints = 0;
    n = a.GetN()
    for i in xrange(n): # for (int i = 0, n = a.GetN() i < n; ++i):
        yhi = a.GetY()[i] + a.GetErrorYhigh(i)
        ylo = a.GetY()[i] - a.GetErrorYlow(i)
        if (a.GetX()[i] > xmax) : continue
        npoints += 1
        if (ylo * yhi > 0) : yavg *= (ylo*yhi)
        if (yhi*3   > ymax) : ymax = yhi * 3
        if (ylo*0.6 < ymin) : ymin = ylo * 0.6
    
    if (ymin < hardymin) : ymin = hardymin
    yavg = pow(yavg, 0.5/float(npoints))
    if (ymax < 10*yavg) :ymax = 10*yavg
    if (pdefs.forceYmin != 0) : ymin = pdefs.forceYmin
    return ymin,ymax


# In[19]:


def channelFromName(who, withSpaces=False, noLumi=False):
    name = who #, space, lumi
    if "comb" in who :
        name = "Combined"
        space = ""
        lumi= "4.6-4.8 fb^{-1}"
    if "hzz" in who: 
        name = "H #rightarrow ZZ Comb." 
        if withSpaces : name =  "H #rightarrow ZZ" 
        space = "              "
        lumi= "4.7 fb^{-1}" 
    
    if "combp" in who : 
        name = "H #rightarrow ZZ + #gamma#gamma" 
        if pdefs.isSquareCanvas : name =  "ZZ + #gamma#gamma" 
        space = "     "
        lumi= "4.8 fb^{-1}" 
        if (withSpaces):
            space = ""
            lumi = ""
    
    if "combl" in who : 
        name = "H #rightarrow bb + #tau#tau + WW"  
        if pdefs.isSquareCanvas: name = "H #rightarrow bb + #tau#tau + WW"  
        space = "   "
        lumi= "4.6 fb^{-1}" 
        if (withSpaces):
            space = ""
            lumi = ""
            
    if "combs" in who : 
        name = "VBF+VH exclusive"  
        if pdefs.isSquareCanvas : name = "VBF+VH excl."  
        space = "   " ; lumi= "4.6-4.8 fb^{-1}"; 
        if (withSpaces): 
            space = "" 
            lumi = ""
    
    if "hww" in who :
        name = "H #rightarrow WW" ; space = "           "; lumi="4.6 fb^{-1}"
    if "hww2l" in who :
        name = "H #rightarrow WW #rightarrow 2l 2#nu" ; space = ""; lumi="4.6 fb^{-1}"
    if "vhww3l" in who :
        name = "WH #rightarrow 3l 3#nu" ; space = "     "; lumi="4.6 fb^{-1}"
    if "vhbb" in who :
        name = "H #rightarrow bb" ; space = "             "; lumi="4.7 fb^{-1}"
    if "hgg_vbfonly" in who :
        name = "VBF H #rightarrow #gamma#gamma" ; space = "        "; lumi="4.8 fb^{-1}"
    if "hgg" in who :
        name = "H #rightarrow #gamma#gamma" ; space = "              "; lumi="4.8 fb^{-1}"
    if "htt" in who :
        name = "H #rightarrow #tau#tau" ; space = "              "; lumi="4.6 fb^{-1}"
    if "vhtt" in who :
        name = "VH #rightarrow #tau_{h} 2l" ; space = "         "; lumi="4.6 fb^{-1}"
    if "httm" in who :
        name = "H #rightarrow #tau#tau #rightarrow #mu#mu" ; space = "    "; lumi="4.6 fb^{-1}"
    if "hzz4l" in who :
        name = "H #rightarrow ZZ #rightarrow 4l" ; space = "     "; lumi="4.7 fb^{-1}"
    if "hzz2l2q" in who :
        name = "H #rightarrow ZZ #rightarrow 2l 2q" ; space=""; lumi="4.6 fb^{-1}"
    if "hzz2l2nu" in who :
        name = "H #rightarrow ZZ #rightarrow 2l 2#nu" ; space=""; lumi="4.6 fb^{-1}"
    if "hzz2l2t" in who :
        name = "H #rightarrow ZZ #rightarrow 2l 2#tau" ; space=""; lumi="4.6 fb^{-1}"
    if (pdefs.lessSpacesInLegends):
        nsp = space.Length()-12 ; space = "";
        for i in xrange(nsp): 
            space += " " # for (int i = 0 i < nsp; ++i) space += " ";
    
    if (noLumi) : return name
    if (withSpaces) :
        if (name == "Combined" or lumi == ""): return name
        return name + space + "  ("+lumi+")"
    if (name == "Combined" and pdefs.justLumiForCombined or pdefs.channelSpamOnRightHandSide) :
        return pdefs.lumiSymbol+" = "+lumi
    return name + ", "+pdefs.lumiSymbol+" = "+lumi


# In[20]:


def colorFromName(who, dark=False):
    if ( "combp" in who ):
        return 2
        if dark: return 205
    if ( "combl" in who ):
        return 4
        if dark: return 213
    if ( "comb" in who ):
        return 1
        if dark: return 19
    if ( "hww2l" in who ):
        return 215
        if dark: return 215
    if ( "vhww" in who ):
        return 213
        if dark: return 213
    if ( "hww" in who ):
        return 4
        if dark: return 4
    if ( "hgg" in who ):
        return 209
        if dark: return 209
    if ( "vhbb" in who ):
        return 67
        if dark: return 67
    if ( "vhtt" in who ):
        return 51
        if dark: return 51
    if ( "httm" in who ):
        return 223
        if dark: return 223
    if ( "htt" in who ):
        return 221
        if dark: return 221
    if ( "hzz4l" in who ):
        return 2
        if dark: return 2
    if ( "hzz2l2q" in who ):
        return 93
        if dark: return 93
    if ( "hzz2l2nu" in who ):
        return 28
        if dark: return 28
    if ( "hzz2l2t" in who ):
        return 223
        if dark: return 223
    if ( "hzz" in who ):
        return 2
        if dark: return 205
    return 39


# In[21]:


def lineStyleFromName(who):
    if (pdefs.noLineStyles) : return 1
    # sed 's#(who.Contains(#in who #g' | sed 's#))# :#g' | awk '{print "    "$1" \""$4"\" "$2" "$3" "$5" "$6" "$7}'
    if "combzz" in who : return 1
    if "comb" in who : return 1
    if "vhbb" in who : return 3
    if "hww" in who : return 7
    if "hgg" in who : return 5
    if "htt" in who : return 2
    if "hzz4l" in who : return 1
    if "hzz2l2q" in who : return 5
    if "hzz2l2nu" in who : return 2
    if "hzz2l2t" in who : return 3
    return 39


# In[22]:


def drawOnePlot(who, what="auto"):
    #global c1
    obs = ROOT.gFile.Get(who+"_obs")
    if not obs:
        print "Missing "+who
        return
    if obs is None:
        print "Missing "+who
        return
    miss = missingPoints(obs)
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        xmax = 600.1
        if pdefs.SM == "FP" : xmax = 300.1
    col = 1 #//colorFromName(who);
    obs.SetLineWidth(2)
    obs.SetLineColor(col)
    obs.SetMarkerColor(col) 
    obs.SetMarkerStyle(21)
    obs.SetMarkerSize(0.8)
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0)
    #double ymin, ymax; minMaxY(obs, ymin, ymax);
    ymin, ymax = minMaxY(obs)
    isML = "ml_" in who or "mlz_" in who # echo 'who.Contains("ml_") or who.Contains("mlz_")' | sed 's#who.Contains(#in who #g' | sed 's#))# :#g' | awk '{print $3" "$1" "$2" "$4" "$7" "$5" "$6" "$8}' | sed 's#)##g'
    if ("ml_" in who) : ymin = 0.02
    if ("mlz_" in who): 
        ymin = -2.5 ; ymax = 5; 
        if ("comb" in who and not ("combzz" in who)) :
            ymin = -1
            ymax = 2.5
    
    if (isML):
        obs.SetFillColor(pdefs.colorFit68)
        obs.Draw("E3")
        frame0.Draw("AXIGSAME")
    obs.Draw("LPX")    
    if (miss) : miss.Draw("P")
    if (what == "auto"):
        if ("pla_" in who) : what = "P.L. Approx limit #sigma_{95%}/#sigma_{"+pdefs.SM+"}"
        if ("bayes_" in who ) : what = "Bayesian limit #sigma_{95%}/#sigma_{"+pdefs.SM+"}"
        if (isML) : what = "Best fit #sigma/#sigma_{"+pdefs.SM+"}"
    pdefs.leg = 0
    dummyBox = ROOT.TH1F("dummyBox","dummyBox",1,0.,1.)
    dummyBox.SetFillColor(pdefs.colorFit68)
    dummyBox.SetLineColor(pdefs.colorFit68)
    if (isML):
        print "isML"
        notJustComb = (who != "mlz_comb")
        issquarecanvas=0.0
        if pdefs.isSquareCanvas : issquarecanvas=1.0
        notjustcomb = 0.0
        if notJustComb : notjustcomb = 1.0
        pdefs.leg = newLegend(.65-0.09*issquarecanvas+0.07*issquarecanvas*notjustcomb,.85,.90+0.01*issquarecanvas,.92) 
        pdefs.leg.SetTextSize(0.045-0.000*issquarecanvas-0.007*issquarecanvas*notjustcomb)
        pdefs.leg.AddEntry(dummyBox, pdefs.oneSigmaFitText, "F")
        pdefs.leg.Draw()
    #//return
    #//bockjoo setCanvas(&frame0, "", ymin, ymax, what)
    #//return
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(who)
    print " myspam = ", myspam
    #//return
    finalize(who,xmin,xmax,ymin,ymax, myspam)
    #print " return "
    #return
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels()
        frame0.GetXaxis().SetNoExponent()
        finalize(who+"_logx",xmin,xmax,ymin,ymax, myspam)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom)
        frame2.Draw(); 
        obs.SetFillColor(pdefs.colorFit68)
        if (isML):
            obs.SetFillColor(pdefs.colorFit68)
            obs.Draw("E3")
            frame2.Draw("AXIGSAME")
        obs.Draw("LXP")
        if (miss) : miss.Draw("P")
        if (pdefs.leg) : pdefs.leg.Draw()
        if not ("mlz_" in who) : ymin, ymax = minMaxY(obs, ymin, ymax, pdefs.x_zoom) 
        if ("ml_" in who) : ymin = 0.02 
        setCanvas(frame2, "", ymin, ymax, what)
        finalize(who+"_zoom",xmin,pdefs.x_zoom,ymin,ymax, myspam)
        if (pdefs.doZoom2):
            frame3 = ROOT.TH1D("frame3","frame3", 1, pdefs.x_zoom2_min,pdefs.x_zoom2_max)
            frame3.Draw() 
            obs.SetFillColor(pdefs.colorFit68)
            if (isML):
                obs.SetFillColor(pdefs.colorFit68)
                obs.Draw("E3")
                frame3.Draw("AXIGSAME")
            obs.Draw("LXP")  
            if (miss) : miss.Draw("P")
            if (pdefs.leg) : pdefs.leg.Draw()
            setCanvas(frame3, "", ymin, ymax, what)
            finalize(who+"_zoom2", pdefs.x_zoom2_min,pdefs.x_zoom2_max, ymin,ymax, myspam)
    pdefs.leg = 0



# In[23]:


def minPValue(g, x_min, x_max, minPoints=0, spikeKiller=0):
    allow_halfint_pvals = True
    points = 0 ; imin = 0 ; pmin = 1.0
    n = g.GetN()
    for i in xrange(n): #    for (int i = 0, n = g.GetN() i < n; ++i):
        if (g.GetX()[i] < x_min or g.GetX()[i] > x_max) : continue
        if (not allow_halfint_pvals and float(int(g.GetX()[i]) != g.GetX()[i]) ) : continue
        if (spikeKiller > 0):
            Z  = ROOT.Math.normal_quantile_c(g.GetY()[i],1)
            ZP = 0.0
            if i > 0 : ZP = ROOT.Math.normal_quantile_c(g.GetY()[i-1],1.)
            ZN = 0.0
            if i < n-1 : ZN = ROOT.Math.normal_quantile_c(g.GetY()[i+1],1.)
            if (Z - ROOT.TMath.Min(ZP,ZN) > spikeKiller) : continue

        points += 1
        if (g.GetY()[i] < pmin):
            pmin = g.GetY()[i]
            imin = i;

    if (points < minPoints): 
        print "Plot ", g.GetName(), " has ", points, " points, less than ", minPoints 
        return -1 
    
    return pmin


# In[24]:


def drawOnePVal(who, what="Local p-value", drawToys=False):
    obs = ROOT.gFile.Get(who+"_obs")
    if not obs: print obs," not found"; return
    if obs is None : return
    exp = ROOT.gFile.Get(who+"_median")
    if not exp: exp = ROOT.gFile,Get(who+"_asimov")
    #if exp is None : exp = ROOT.gFile.Get(who+"_asimov")

    toys = 0
    toysExp = 0
    #if (drawToys and who.Index("pvala") == 0):
    if (drawToys and who.index("pvala") == 0):
        twho = who
        #twho.ReplaceAll("pvala","pval");
        twho.replace("pvala","pval");
        toys = ROOT.gFile.Get(twho+"_obs")
        if (toys == 0) : return
        toys.SetLineWidth(4)
        toys.SetMarkerStyle(20)
        toys.SetMarkerSize(1.4)
        toys.SetLineColor(62)
        toys.SetMarkerColor(62)
        toysExp = ROOT.gFile.Get(twho+"_median")
        if (toysExp != 0):
            toysExp.SetLineWidth(4)
            toysExp.SetMarkerStyle(27)
            toysExp.SetMarkerSize(1.2)
            toysExp.SetLineColor(51)
            toysExp.SetMarkerColor(51)
        
    

    miss = missingPoints(obs)
    if toys != 0 : miss = 0
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        xmax = 600.1
        if pdefs.SM == "FP" : xmax = 300.1
    col = 1 #//colorFromName(who);
    if toys == 0 :
        obs.SetLineWidth(2)
    else:
        obs.SetLineWidth(3)
    obs.SetLineColor(col)
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21)
    obs.SetMarkerSize(0.9);
    if not ( expr is None): # if exp
        smooth = smoothSMCLs(exp, 7, 2)
        if smooth is None : print "Smoothing of ", expi[i].GetName(), " returned ZERO"
        elif (pdefs.SM != "SM4") : exp = smooth
        exp.SetLineWidth(2)
        exp.SetLineStyle(1)
        exp.SetLineColor(64)
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0)
    frame0 = ROOT.TH1Dframe0("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0)
    #ymin, ymax 
    #minMaxY(obs, ymin, ymax, 999, 0); 
    ymin, ymax = minMaxY(obs, 0.0, 0.0, 999, 0); 
    ymax = 1.0
    ymin = minPValue(obs,xmin,xmax)/10
    if (ymin <= 0 or ymin > 1e-5) : ymin = 1e-5
    
    setCanvas(frame0, "", ymin, ymax, what)
    issquarecanvas = 0.0
    if pdefs.isSquareCanvas : issquarecanvas = 1.0
    frame0.GetYaxis().SetTitleOffset(1.10+0.25*issquarecanvas)
    if (toys != 0) : toys.Draw("P")
    if (toysExp !=0 ) : toysExp.Draw("P")
    obs.Draw("LXP")    
    if toys !=0 : obs.Draw("LX")    
    if not (exp is None) : exp.Draw("LX")
    if miss != 0 : miss.Draw("P")
    #TString myspam = (toys ? "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}" : pdefs.SPAM+", "+channelFromName(who)) ;
    if toys !=0 :
        myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"+", "+channelFromName(who) 
    else:
        myspam = pdefs.SPAM+", "+channelFromName(who)
    if toys !=0:
        toysexp=0.0
        if toysExp != 0 : toysexp = 1.0
        texmp=0.0
        if not (exp is None) : texmp=1.0
        pdefs.leg = newLegend(.66-0.065*issquarecanvas,.18,.93,.32+(texmp+toysexp)*0.04)
        pdefs.leg.SetTextSize(0.04-0.003*issquarecanvas);
        pdefs.leg.AddEntry(obs,  "Asymptotic Obs.", "L")
        pdefs.leg.AddEntry(toys, "Ensemble Obs.",   "P")
        if toysExp != 0 : pdefs.leg.AddEntry(toysExp, "Ensemble Exp.",   "P")
        if not (exp is None) : pdefs.leg.AddEntry(exp, "Asymptotic Exp.", "L")
    else:
        pdefs.leg = 0
    
    
    if drawToys:
        finalize(who+"_comp",xmin,xmax,ymin,ymax,myspam,True)
    else:
        finalize(who,xmin,xmax,ymin,ymax,myspam,True)
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels()
        frame0.GetXaxis().SetNoExponent()
        #finalize(who+(drawToys?"_comp_logx":"_logx"),xmin,xmax,ymin,ymax,myspam,True)
        if drawToys:
            finalize(who+"_comp_logx",xmin,xmax,ymin,ymax,myspam,True)
        else:
            finalize(who+"_logx",xmin,xmax,ymin,ymax,myspam,True)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom)
        frame2.Draw() 
        setCanvas(frame2, "", ymin, ymax, what)
        if toys != 0 : toys.Draw("P")
        if toysExp != 0 : toysExp.Draw("P")
        if toys == 0 : obs.Draw("LXP")    
        else : obs.Draw("LX")    
        if not (exp is None) : exp.Draw("LX")
        if miss != 0 : miss.Draw("P")
        #finalize(who+(drawToys?"_comp_zoom":"_zoom"),xmin,pdefs.x_zoom,ymin,ymax,myspam,True)
        if drawToys:
            finalize(who+"_comp_zoom",xmin,xmax,ymin,ymax,myspam,True)
        else:
            finalize(who+"_zoom",xmin,xmax,ymin,ymax,myspam,True)
        if (pdefs.doZoom2):
            frame2 = ROOT.TH1D("frame3","frame3", 1, pdefs.x_zoom2_min,pdefs.x_zoom2_max)
            frame3.Draw() 
            setCanvas(frame3, "", ymin, ymax, "Local p-value")
            obs.Draw("LXP")
            if not (exp is None) : exp.Draw("L")
            if miss != 0 : miss.Draw("P")
            finalize(who+"_zoom2",pdefs.x_zoom2_min,pdefs.x_zoom2_max,ymin,ymax,myspam,True)
        
    
    pdefs.leg = 0




# In[25]:


def drawOnePValBand(who, who2=""):
    obs = ROOT.gFile.Get(who+"_obs")
    m68 = ROOT.gFile.Get(who+"_median")
    m95 = ROOT.gFile.Get(who+"_median_95") 
    obs2 = ROOT.gFile.Get(who2+"_obs")
    if not obs or not m68  or not m95 : print obs, " ",m68," or ",m95," not found"; return
    if (obs is None or m68 is None or m95 is None) : return

    label1="" ;  label2=""; postfix="_band"; style2=""
    if (not (obs2 is None) and who.index("pvala") == 0 and who2.index("pval_") == 0):
        label1="Asympt. " ; label2 = "Ensemble "; postfix="_band_as"; style2="P";
        obs2.SetLineWidth(4)
        obs2.SetMarkerStyle(20)
        obs2.SetMarkerSize(1.4)
        obs2.SetLineColor(62)
        obs2.SetMarkerColor(62)
    elif (not (obs2 is None) and who.index("pval_") == 0 and who2.index("pvala_") == 0):
        label1="Ensemble " ; label2 = "Asympt. "; postfix="_band_t"; style2="L";
        obs2.SetLineColor(2)
        obs2.SetLineStyle(1)
        obs2.SetLineWidth(3)
    

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        #xmax = (pdefs.SM == "FP" ? 300.1 : 600.1); }
        xmax = 600.1
        if pdefs.SM == "FP" : xmax = 300.1
    col = 1 #//colorFromName(who);
    obs.SetLineWidth(4)
    obs.SetLineColor(1)
    obs.SetMarkerColor(1); 
    obs.SetMarkerStyle(21)
    obs.SetMarkerSize(0.8);
    m68.SetFillColor(pdefs.color68)
    m68.SetLineColor(pdefs.color50)
    m68.SetLineStyle(2)
    m68.SetLineWidth(2);
    m95.SetFillColor(pdefs.color95)
    m95.SetLineColor(pdefs.color50)
    m95.SetLineStyle(2)
    m95.SetLineWidth(2);

    
    
    
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    ymin, ymax = minMaxY(obs, 0.0, 0.0, 999, 0); 
    ymax = 1.0
    ymin = minPValue(obs,xmin,xmax)/10
    if (ymin <= 0 or ymin > 1e-5) : ymin = 1e-5
    
    setCanvas(frame0, "", ymin, ymax, "Local p-value")
    issquarecanvas = 0.0
    if pdefs.isSquareCanvas : issquarecanvas = 0.0
    frame0.GetYaxis().SetTitleOffset(1.10+0.25*issquarecanvas)
    draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, -5)
    frame0.Draw("AXIGSAME")
    if (style2 == "L"):
        obs.Draw("LXP")    
        if not (obs2 is None) :  obs2.Draw("LX")
    else:
        if not (obs2 is None) :  obs2.Draw("PX")
        obs.Draw("LX")    
    
    pdefs.leg = newLegend(.66-0.065*issquarecanvas,.18,.93,.38)
    pdefs.leg.SetTextSize(0.04-0.003*issquarecanvas);
    if obs2 is None :
       if style2 == "L" :
           pdefs.leg.AddEntry(obs, "Observed", "LP")
       else:
          pdefs.leg.AddEntry(obs, "Observed", "L")
    else :
       if style2 == "L" :
           pdefs.leg.AddEntry(obs, label1+"Obs.", "LP")
       else:
          pdefs.leg.AddEntry(obs, label1+"Obs.", "L")
    pdefs.leg.AddEntry(m68, "Expected "+pdefs.oneSigmaText,   "LF")
    pdefs.leg.AddEntry(m95, "Expected "+pdefs.twoSigmaText,   "LF")
    if not (obs2 is None) : pdefs.leg.AddEntry(obs2, label2+"Obs.",  style2)
    
    pdefs.leg.Draw()

    myspam =  "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(who)
    finalize(who+postfix,xmin,xmax,ymin,ymax,myspam,True)
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels()
        frame0.GetXaxis().SetNoExponent();
        finalize(who+postfix+"_logx",xmin,xmax,ymin,ymax,myspam,True)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom)
        frame2.Draw() 
        setCanvas(frame2, "", ymin, ymax, "Local p-value")
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, -5)
        frame0.Draw("AXIGSAME")
        if (style2 == "L"):
            obs.Draw("LXP")    
            if not (obs2 is None) : obs2.Draw("LX")
        else:
            if not (obs2 is None) : obs2.Draw("PX")
            obs.Draw("LX")    
        
        pdefs.leg.Draw()
        finalize(who+postfix+"_zoom",xmin,xmax,ymin,ymax,myspam,True)
        
    pdefs.leg = 0
    


# In[26]:


def drawOnePValFitCollage(who, whichFit, xmin, xmax, logx, postfix):
    obs = ROOT.gFile.Get("pvala_"+who+"_obs")
    exp = ROOT.gFile.Get("pvala_"+who+"_median")
    mlf = ROOT.gFile.Get(whichFit+"_"+who+"_obs")
    if not obs or not exp or not mlf : print obs, " ", exp, " or ", mlf, " not found"; return
    if obs is None or mlf is None : return

    col = 1 
    mlf.SetLineWidth(2)
    mlf.SetLineColor(col) 
    mlf.SetMarkerColor(col)
    mlf.SetMarkerStyle(21)
    mlf.SetMarkerSize(0.8)

    zero = "mlz" in whichFit
    pdefs.isTiny = True
    c1_1 = ROOT.TPad("pad1", "The pad 80% of the height",0.0,0.45,1.0,1.00,-1)
    c1_1.Draw() 
    c1_2 = ROOT.TPad("pad2", "The pad 20% of the height",0.0,0.00,1.0,0.45,-1)
    c1_2.Draw() 
    c1_1.SetLogx(logx) 
    c1_2.SetLogx(logx); 
    c1_1.SetLogy(1) 
    c1_2.SetLogy((not zero))
    c1_1.SetBottomMargin(0)
    c1_2.SetTopMargin(0.0)
    c1_2.SetBottomMargin(0.25)
    c1_1.SetLeftMargin(0.13)
    c1_2.SetLeftMargin(0.13)
    c1_1.SetRightMargin(0.04)
    c1_2.SetRightMargin(0.04)
    c1_1.cd()

    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax) 
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    ymin, ymax = minMaxY(obs, 0.0, 0.0, 999, 0)
    ymax = 1.0; ymin = 2e-5;
    setCanvas(frame0, "", ymin, ymax, "Local p-value")
    if pdefs.isSquareCanvas :
        frame0.GetYaxis().SetLabelSize(0.065)
        frame0.GetXaxis().SetLabelSize(0.065)
        frame0.GetYaxis().SetTitleSize(0.075)
        frame0.GetXaxis().SetTitleSize(0.075)
        frame0.GetYaxis().SetTitleOffset(0.68)
    else:
        frame0.GetYaxis().SetLabelSize(0.07)
        frame0.GetXaxis().SetLabelSize(0.07)
        frame0.GetYaxis().SetTitleSize(0.08)
        frame0.GetXaxis().SetTitleSize(0.08)
        frame0.GetYaxis().SetTitleOffset(0.6)
    frame0.GetXaxis().SetTitle("")
        
    #//if (exp) exp.Draw("L")
    obs.Draw("LXP")
    myspam = pdefs.SPAM+", "+channelFromName(who)
    if (pdefs.isSquareCanvas) : myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    finalizeNoSave("pvala_"+who,xmin,xmax,ymin,ymax,myspam,True)

    c1_2.cd()
    frameF = ROOT.TH1D("frame","frame", 1, xmin, xmax) 
    frameF.Draw(); 
    ROOT.gStyle.SetOptStat(0);
    yminF = 0.02 ; ymaxF = 20 
    if (zero): 
        yminF, ymaxF = minMaxY(mlf, yminF, ymaxF, 999, -1)
        ymaxF /= 2.5 ;
        yminF = -1;
    
    setCanvas(frameF, "", yminF, ymaxF, "Best fit #sigma/#sigma_{"+pdefs.SM+"} ")
    if pdefs.isSquareCanvas :
        frameF.GetYaxis().SetTitleSize(0.09)
        frameF.GetYaxis().SetLabelSize(0.07)
        frameF.GetYaxis().SetTitleOffset(0.44)
    else:
        frameF.GetYaxis().SetTitleSize(0.10)
        frameF.GetYaxis().SetLabelSize(0.08)
        frameF.GetYaxis().SetTitleOffset(0.4)

        
    frameF.GetYaxis().SetNoExponent(1)
    frameF.GetYaxis().SetNdivisions(505)
    frameF.GetXaxis().SetTitleSize(0.10)
    frameF.GetXaxis().SetLabelSize(0.08)
    frameF.GetXaxis().SetMoreLogLabels()
    frameF.GetXaxis().SetNoExponent();
    mlf.SetFillColor(pdefs.colorFit68)
    mlf.Draw("E3")    
    frameF.Draw("AXIGSAME")
    mlf.Draw("LXP")    
    issquarecanvas = 0.0
    if pdefs.isSquareCanvas : issquarecanvas = 0.0
    ZERO=0.0
    if zero: ZERO=1.0
    pdefs.leg = newLegend(.74-0.05*issquarecanvas,.27+.50*ZERO,.92+0.00*issquarecanvas,.40+.50*ZERO)
    pdefs.leg.SetTextSize(0.08);
    pdefs.leg.AddEntry(mlf, pdefs.oneSigmaFitText+" from fit", "F")
    pdefs.leg.Draw()
    finalizeNoSave(whichFit+"_"+who,xmin,xmax, yminF, ymaxF)

    pdefs.c1.cd()
    justSave("pval_"+whichFit+"_"+who+postfix,xmin,xmax,ymin,ymax)
    pdefs.leg = 0

    if (not zero):
        c1_2.cd()
        c1_2.SetLogy(0)
        c1_2.Clear() 
        frameF.Draw()
        yminF, ymaxF = minMaxY(mlf, yminF, ymaxF, 999, -1)
        ymaxF /= 3;
        setCanvas(frameF, "", 0, ymaxF, "Best fit #sigma/#sigma_{"+pdefs.SM+"} ")
        frameF.GetYaxis().SetTitle("Best fit #sigma/#sigma_{"+pdefs.SM+"} ")
        frameF.GetYaxis().SetTitleOffset(0.3+0.05*issquarecanvas)
        mlf.Draw("E3")    
        frameF.Draw("AXIGSAME")
        mlf.Draw("LXP")    
        pdefs.leg = newLegend(.74-0.05*issquarecanvas,.77,.92+0.00*issquarecanvas,.90)
        pdefs.leg.SetTextSize(0.08);
        pdefs.leg.AddEntry(mlf, pdefs.oneSigmaFitText+" from fit", "F")
        pdefs.leg.Draw()
        finalizeNoSave(whichFit+"_"+who,xmin,xmax, 0.0, ymaxF)

        pdefs.c1.cd()
        justSave("pval_"+whichFit+"_"+who+postfix+"_liny",xmin,xmax,ymin,ymax)
    

    c1_2.Delete()
    c1_1.Delete()
    pdefs.c1.Clear()
    pdefs.isTiny = False
    pdefs.leg = 0


# In[27]:


def drawOnePValFitCollage(who, whichFit="ml"):
    obs = ROOT.gFile.Get("pvala_"+who+"_obs")
    mlf = ROOT.gFile.Get(whichFit+"_"+who+"_obs")
    if not obs: print obs, " not found "; return
    if not mlf: print mlf, " not found "; return
    if ( obs is None ) : return
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200): 
        drawOnePValFitCollage(who, whichFit, 99.98, 600.1, False, "")
        drawOnePValFitCollage(who, whichFit, 99.98, 600.1, True,  "_logx")
        drawOnePValFitCollage(who, whichFit, xmin,  pdefs.x_zoom, False, "_zoom")
    else:
        drawOnePValFitCollage(who, whichFit, xmin, xmax, False, "")
        if (xmax > pdefs.x_zoom) : drawOnePValFitCollage(who, whichFit, xmin, pdefs.x_zoom, False, "_zoom")
    



# In[28]:



def drawOnePValFitCollageOld(who, what="auto"):
    obs = ROOT.gFile.Get("pvala_"+who+"_obs")
    mlf = ROOT.gFile.Get("ml_"+who+"_obs")
    if not obs or not mlf : print obs, " ",mlf, " not found "; return
    if obs is None or mlf is None : return

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        xmax = 600.1

    col = 1 
    obs.SetLineWidth(2)
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8);
    mlf.SetLineWidth(2)
    mlf.SetLineColor(col) 
    mlf.SetMarkerColor(col); 
    mlf.SetMarkerStyle(21) 
    mlf.SetMarkerSize(0.8);

    pdefs.isTiny = True
    c1_1 = ROOT.TPad("pad1", "The pad 80% of the height",0.0,0.65,1.0,1.00,-1) 
    c1_1.Draw();
    c1_2 = ROOT.TPad("pad2", "The pad 20% of the height",0.0,0.00,1.0,0.65,-1) 
    c1_2.Draw();
    if (xmax > 500 and xmin < 200): 
        c1_1.SetLogx(1) 
        c1_2.SetLogx(1)
    c1_1.SetLogy(1) 
    c1_2.SetLogy(1);
    c1_1.SetBottomMargin(0)
    c1_2.SetTopMargin(0.02)
    c1_1.SetLeftMargin(0.10)
    c1_2.SetLeftMargin(0.10)
    c1_1.SetRightMargin(0.04)
    c1_2.SetRightMargin(0.04)

    c1_1.cd()

    issquarecanvas = 0.0
    if pdefs.isSquareCanvas : issquarecanvas = 0.0
    ZERO=0.0
    if zero: ZERO=1.0

    #ROOT.TH1DframeF("frame","frame", 1, xmin, xmax) frameF.Draw(); gStyle.SetOptStat(0);
    frameF = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frameF.Draw()
    ROOT.gStyle.SetOptStat(0);
    setCanvas(frameF, "", 0.02, 20, "Best fit #sigma/#sigma_{"+pdefs.SM+"}")
    frameF.GetYaxis().SetTitleSize(0.12)
    frameF.GetYaxis().SetLabelSize(0.10)
    if pdefs.isSquareCanvas :
        frameF.GetYaxis().SetTitleOffset(0.42)
    else:
        frameF.GetYaxis().SetTitleOffset(0.4)
    frameF.GetXaxis().SetTitle("")
    frameF.GetYaxis().SetNoExponent(1)
    mlf.SetFillColor(pdefs.colorFit68)
    mlf.Draw("E3")    
    frameF.Draw("AXIGSAME")
    mlf.Draw("LXP")    
    pdefs.leg = newLegend(.74-0.05*issquarecanvas,.75,.90+0.04*issquarecanvas,.90)
    pdefs.leg.SetTextSize(0.10);
    pdefs.leg.AddEntry(mlf, pdefs.oneSigmaFitText+" from fit", "F")
    pdefs.leg.Draw()
    pdefs.leg = 0
    finalizeNoSave("ml_"+who,xmin,xmax,0.02, 20)

    c1_2.cd()
    #ROOT.TH1Dframe0("frame","frame", 1, xmin, xmax) frame0.Draw(); gStyle.SetOptStat(0);
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);

    ymin, ymax = minMaxY(obs, 0.0, 0.0, 999, 0); 
    ymax = 1.0; 
    if (ymin > 2e-6) : ymin = 2e-6;
    setCanvas(frame0, "", ymin, ymax, "Local p-value")
    if pdefs.isSquareCanvas :
        frame0.GetYaxis().SetLabelSize(0.055)
        frame0.GetXaxis().SetLabelSize(0.055)
        frame0.GetYaxis().SetTitleSize(0.065)
        frame0.GetXaxis().SetTitleSize(0.065)
        frame0.GetYaxis().SetTitleOffset(0.78)
    else:
        frame0.GetYaxis().SetLabelSize(0.06)
        frame0.GetXaxis().SetLabelSize(0.06)
        frame0.GetYaxis().SetTitleSize(0.07)
        frame0.GetXaxis().SetTitleSize(0.07)
        frame0.GetYaxis().SetTitleOffset(0.7)
        
    frame0.GetXaxis().SetMoreLogLabels() 
    frame0.GetXaxis().SetNoExponent();
    obs.Draw("LP")    
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    finalizeNoSave("pvala_"+who,xmin,xmax,ymin,ymax,myspam,True)
    pdefs.c1.cd()
    justSave("pval_fit_"+who,xmin,xmax,ymin,ymax)
    c1_2.Delete()
    c1_1.Delete()
    pdefs.c1.Clear()
    pdefs.isTiny = False
    pdefs.leg = 0



# In[29]:


def drawCompLimit(who1, who2, label="One", oldlabel="Two", what="auto"):
    obs = ROOT.gFile.Get(who1)
    old = ROOT.gFile.Get(who2)
    if not obs or not old : print obs, old, " not found "; return
    if obs is None or old is None : return
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    col = 2; dcol = 4
    obs.SetLineWidth(2)
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8);
    old.SetLineWidth(2)
    old.SetLineColor(dcol) 
    old.SetMarkerColor(dcol); 
    old.SetMarkerStyle(25) 
    old.SetMarkerSize(0.8);
    #ROOT.TH1Dframe0("frame","frame", 1, xmin, xmax) frame0.Draw(); gStyle.SetOptStat(0);
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    ymin, ymax = minMaxY(obs);
    ymin2, ymax2 = minMaxY(old);
    if (ymin2 < ymin) : ymin = ymin2 
    if (ymax2 > ymax) : ymax = ymax2
    if ("mlz_" in who1):
        ymin = -2.5
        ymax = 5
    if (( "median" in who1 and "median" in who2) or
        ("mlz" in who1 and "mlz" in who2)):
        old.SetFillStyle(1001) 
        old.SetFillColor(216);
        obs.SetFillStyle(3444) 
        obs.SetFillColor(2); 
        obs.SetLineColor(216);
        old.Draw("E3") 
        obs.Draw("E3");    
        old.Draw("LPX") 
        obs.Draw("LPX");    
    else:
        old.Draw("LP")
        obs.Draw("LP");    
    
    pdefs.leg = newLegend(.69,.75,.84,.92) 
    pdefs.leg.SetTextSize(0.05);
    pdefs.leg.AddEntry(obs,    label, "LP")
    pdefs.leg.AddEntry(old, oldlabel, "LP")
    pdefs.leg.Draw()
    setCanvas(frame0, "", ymin, ymax, "Limit")
    myspam("#splitline{PRIVATE FOR COMB.}{"+channelFromName(who1)+"}")
    finalize(who1+"_vs_"+who2,xmin,xmax,ymin,ymax,myspam)
    if (xmax >= 200 and xmin < pdefs.x_zoom):
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        #if ((who1.Contains("median") and who2.Contains("median")) or
        #    (who1.Contains("mlz") and who2.Contains("mlz"))):
        if (( "median" in who1 and "median" in who2) or
            ("mlz" in who1 and "mlz" in who2)):
            old.Draw("E3") 
            obs.Draw("E3");    
            old.Draw("LPX") 
            obs.Draw("LPX");    
        else:
            old.Draw("LP") 
            obs.Draw("LP");    
        
        pdefs.leg.Draw()
        ymin, ymax = minMaxY(obs, ymin, ymax, pdefs.x_zoom) 
        ymin2, ymax2 = minMaxY(old, ymin2, ymax2, pdefs.x_zoom);
        if (ymin2 < ymin) : ymin = ymin2 
        if (ymax2 > ymax) : ymax = ymax2;
        if ("mlz_" in who1):
            ymin = -2.5 
            ymax = 5
        setCanvas(frame2, "", ymin, ymax, what)
        finalize(who1+"_vs_"+who2+"_zoom",xmin,pdefs.x_zoom,ymin,ymax,myspam)
        if (pdefs.doZoom2):
            frame3 = ROOT.TH1D("frame3","frame3", 1,pdefs.x_zoom2_min, pdefs.x_zoom2_max) 
            frame3.Draw(); 
            #if ((who1.Contains("median") and who2.Contains("median")) or
            #    (who1.Contains("mlz") and who2.Contains("mlz"))):
            if (( "median" in who1 and "median" in who2) or
                ("mlz" in who1 and "mlz" in who2)):
                old.Draw("E3") 
                obs.Draw("E3");    
                old.Draw("LPX") 
                obs.Draw("LPX");    
            else:
                old.Draw("LP") 
                obs.Draw("LP");    
            
            pdefs.leg.Draw()
            setCanvas(frame3, "", ymin, ymax)
            finalize(who1+"_vs_"+who2+"_zoom2",pdefs.x_zoom2_min,pdefs.x_zoom2_max,ymin,ymax,myspam)
        
    
    pdefs.leg = 0


# In[30]:


def drawCompPVal(who1, who2, label="One", oldlabel="Two", what="auto", ymin=1e-4):
    obs = ROOT.gFile.Get(who1) 
    if not obs: print  "Cannot find " , who1
    if obs is None : print "Cannot find " , who1
    old = ROOT.gFile.Get(who2) 
    #if (old == 0) std::cout << "Cannot find " << who2 << std::endl;
    if old is None : print "Cannot find " , who2

    #if (obs == 0 or old == 0) return
    if obs is None or old is None : return

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    #if (xmin <= 120 and xmax > 200): xmin = 99.98 xmax = (pdefs.SM == "FP" ? 300.1 : 600.1); }
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
        
    col = 2 ; dcol = 4
    obs.SetLineWidth(2)
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8);
    old.SetLineWidth(2)
    old.SetLineColor(dcol) 
    old.SetMarkerColor(dcol); 
    old.SetMarkerStyle(25) 
    old.SetMarkerSize(0.8);
    #ROOT.TH1Dframe0("frame","frame", 1, xmin, xmax) frame0.Draw(); gStyle.SetOptStat(0);
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    ymax = 1.0
    old.Draw("LP")    
    obs.Draw("LP")    
    pdefs.leg = newLegend(.75,.14,.92,.28) 
    pdefs.leg.SetTextSize(0.05);
    pdefs.leg.AddEntry(obs,    label, "LP")
    pdefs.leg.AddEntry(old, oldlabel, "LP")
    pdefs.leg.Draw()
    setCanvas(frame0, "", ymin, ymax, "Local p-value")
    myspam("#splitline{PRIVATE FOR COMB.}{"+channelFromName(who1)+"}")
    finalize(who1+"_vs_"+who2,xmin,xmax,ymin,ymax,myspam,True)
    if (xmax >= 200 and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        finalize(who1+"_vs_"+who2+"_logx",xmin,xmax,ymin,ymax,myspam,True)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        old.Draw("LP")    
        obs.Draw("LP")  
        pdefs.leg.Draw()
        setCanvas(frame2, "", ymin, ymax, "Local p-value")
        finalize(who1+"_vs_"+who2+"_zoom",xmin,200,ymin,ymax,myspam,True)
        if (pdefs.doZoom2):
            frame3 = ROOT.TH1D("frame3","frame3", 1,pdefs.x_zoom2_min, pdefs.x_zoom2_max) 
            frame3.Draw(); 
            old.Draw("LP")    
            obs.Draw("LP")  
            pdefs.leg.Draw()
            setCanvas(frame3, "", ymin, ymax, "Local p-value")
            finalize(who1+"_vs_"+who2+"_zoom2",pdefs.x_zoom2_min,pdefs.x_zoom2_max,ymin,ymax,myspam,True)
        
    
    pdefs.leg = 0



# In[31]:


def makeAPrioriGrid(who):
    print "Asked for grid " , who #<< std::endl #//CLs_grids.ls();
    if (pdefs.CLs_grids.FindObject("grid_"+who)) :
        return pdefs.CLs_grids.FindObject("grid_"+who)
    return makeGrid(who+"_obs", who+"_median", 0.4, 2.5)



# In[32]:


pdefs.Draw_TEV=False
def drawOneCLs(who):
    obs   = ROOT.gFile.Get(who+"_obs")
    #if (obs == 0) return
    if not obs : print "drawOneCLs obs is nill" ; return
    if obs is None : return
    if pdefs.CLs_debug_apriori_grid : apriori = makeAPrioriGrid(who)
    else : apriori = 0 # None
    miss = missingPoints(obs)
    m68 = ROOT.gFile.Get(who+"_median")
    m95 = ROOT.gFile.Get(who+"_median_95") 
    if m68:
        print "missingPoints m68"
        miss68 = missingPoints(m68)
    else:
        print "m68 is 0"
        miss68 = 0
        
    if miss68:
        miss68.SetMarkerStyle(24) 
        miss68.SetMarkerSize(0.9)
        miss68.SetMarkerColor(4)
    else:
        print "DEBUG drawOneCLs miss68 is 0"
    if (pdefs.Draw_TEV):
        obsTEV =ROOT.gFile.Get("tevatron_obs")
        expTEV = ROOT.gFile.Get("tevatron_median")
        obsTEV.SetLineColor(kBlue) 
        obsTEV.SetLineWidth(5); 
        obsTEV.SetLineStyle(1);
        expTEV.SetLineColor(kBlue) 
        expTEV.SetLineWidth(5); 
        expTEV.SetLineStyle(2);
    else:
        obsTEV = 0
        expTEV = 0
        

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    col = 1; dcol = 1; smooth = 5; smoothorder = 0
    obs.SetLineWidth(2) 
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8);
    ymin, ymax =  minMaxY(obs);
    if m68 : ymin2, ymax2 = minMaxY(m68); 
    if (ymax2 > ymax) : ymax = ymax2 
    if (ymin2 < ymin) : ymin = ymin2;
    if (pdefs.Draw_TEV) : ymax *=2
        
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    frame0.GetYaxis().SetTitleOffset(1.10+0.25*pdefs.isSquareCanvas)
    what = "95% CL limit on #sigma/#sigma_{"+pdefs.SM+"}"
    if ("cls99_" in who) : what = "99% CL limit on #sigma/#sigma_{"+pdefs.SM+"}"
    if ("cls90_" in who) : what = "90% CL limit on #sigma/#sigma_{"+pdefs.SM+"}"
    #//if (who.Contains("acls_")) what = "Asymptotic "+what
    setCanvas(frame0, "", ymin, ymax, what)
    if apriori : apriori.Draw("E3")

    #draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, smooth, smoothorder) # bockjoo 68% and 95% bands and the median
    # def draw2(...)
    mean = False
    slidingWindow = smooth
    lineColor = pdefs.color50
    fillColor68 = pdefs.color68
    fillColor95 = pdefs.color95
    same = True
    mean68 = ROOT.gROOT.FindObject(who+"_median")
    if mean: mean68 = ROOT.gROOT.FindObject(who+"_mean")
    mean95 = ROOT.gROOT.FindObject(who+"_median"+"_95")
    if mean: mean95 = ROOT.gROOT.FindObject(who+"_mean"+"_95")
    if not mean68:
        if mean:
            print "MISSING ", who, "_mean"
        else:
            print "MISSING ", who, "_median"        
        return 0
    if not mean95:
        if mean:
            print "MISSING ", who, "_mean_95"
        else:
            print "MISSING ", who, "_median_95"
        return 0
    print "DEBUG draw2 2"
    mean68 = removeGlitches(mean68)
    mean95 = removeGlitches(mean95)
    if (mean68.GetN() == 1):
        mean68.SetPointError(0, 10, 10, mean68.GetErrorYlow(0), mean68.GetErrorYhigh(0))
        mean95.SetPointError(0, 10, 10, mean95.GetErrorYlow(0), mean95.GetErrorYhigh(0))
    
    meanL = mean68.Clone() # TGraphAsymmErrors *meanL = (TGraphAsymmErrors*) mean68->Clone();
    for i in xrange(meanL.GetN()): # for (int i= 0 i < meanL.GetN(); ++i): meanL.SetPointError(i, 0,0,0,0); }
        #print "DEBUG draw2 3 meanL.SetPointError ",i
        meanL.SetPointError(i, 0,0,0,0)
        
    print "DEBUG draw2 3.1 slidingWindow ",slidingWindow, " mean68.GetN() ",mean68.GetN() 
    if (slidingWindow != 0 and mean68.GetN() > 5):
        print "DEBUG draw2 4 "
        if (slidingWindow > 0):
            print "DEBUG draw2 5 "
            mean68 = slidingWindowAverage(mean68, slidingWindow)
            mean95 = slidingWindowAverage(mean95, slidingWindow)
            meanL  = slidingWindowAverage(meanL,  slidingWindow)
            print "DEBUG draw2 meanL 0",meanL.GetY()[0], " N = ",meanL.GetN()
        else:
            print "DEBUG draw2 6 "
            mean68 = smoothSMCLs(mean68, -slidingWindow, smoothorder)
            mean95 = smoothSMCLs(mean95, -slidingWindow, smoothorder)
            meanL  = smoothSMCLs(meanL,  -slidingWindow, smoothorder)        

    meanL.SetLineColor(lineColor)
    mean68.SetLineColor(lineColor)
    mean95.SetLineColor(lineColor);
    meanL.SetMarkerColor(lineColor)
    mean68.SetMarkerColor(lineColor)
    mean95.SetMarkerColor(lineColor);
    meanL.SetLineWidth(3)
    mean68.SetLineWidth(3)
    mean95.SetLineWidth(3)
    meanL.SetMarkerSize(1.6)
    mean68.SetMarkerSize(1.6)
    mean95.SetMarkerSize(1.6)
    meanL.SetLineStyle(7)
    mean68.SetLineStyle(7)
    mean95.SetLineStyle(7) 
    mean68.SetLineColor(fillColor68)
    mean95.SetLineColor(fillColor95)
    mean68.SetLineWidth(1) 
    mean95.SetLineWidth(1)
    mean68.SetFillColor(fillColor68)  
    mean95.SetFillColor(fillColor95)
    if same: print "drawing mean95.Draw(E3 SAME)"
    else: print "drawing mean95.Draw(AE3)"
    if same: mean95.Draw("E3 SAME")
    else: mean95.Draw("AE3")
    print "drawing mean68"
    mean68.Draw("E3 SAME")
    print "drawing meanL"
    meanL.Draw("LX SAME")
    # end of def draw2(...)
    #draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, smooth, smoothorder) # bockjoo 68% and 95% bands and the median
    #time.sleep(10)

    if m68: 
        m68.SetFillColor(pdefs.color68) 
        m68.SetLineColor(pdefs.color50) 
        m68.SetLineStyle(2) 
        m68.SetLineWidth(2)
    if m95: 
        m95.SetFillColor(pdefs.color95) 
        m95.SetLineColor(pdefs.color50) 
        m95.SetLineStyle(2) 
        m95.SetLineWidth(2)
    frame0.Draw("AXIGSAME")
    leg_y_hi = 0.94
    #DRAW_TEV=0.0
    #if pdefs.Draw_TEV : DRAW_TEV=1.0
    #DRAW_TEVORTEV_EXCLUDED = 0.0
    #if pdefs.Draw_TEVorpdefs.tev_excluded : DRAW_TEVORTEV_EXCLUDED = 1.0
    if pdefs.isSquareCanvas:
        leg_y_lo = 0.78- 0.05*(pdefs.lep_excluded+pdefs.tev_excluded+1.6*pdefs.Draw_TEV)
        leg_x_lo = .645-0.06*(pdefs.Draw_TEV or pdefs.tev_excluded)
    else:
        leg_y_lo = .75- 0.05*(pdefs.lep_excluded+pdefs.tev_excluded+1.6*pdefs.Draw_TEV)
        leg_x_lo = .65
        
    pdefs.leg = newLegend(leg_x_lo, leg_y_lo,.93+0.005*pdefs.isSquareCanvas, leg_y_hi) # bockjoo Legend Location
    if pdefs.isSquareCanvas or pdefs.lep_excluded : pdefs.leg.SetTextSize(0.034)
    else : pdefs.leg.SetTextSize(0.037)
    #print "DEBUG sleep 10 before drawing the observed points"
    #time.sleep(10)
    obs.Draw("LP")    # bockjoo Observed Points with a Line
    #print "DEBUG sleep 10 after drawing the observed points"
    #time.sleep(10)
    if miss68: miss68.Draw("P")
    if miss: miss.Draw("P")
    if apriori: apriori.Draw("LX")
    pdefs.leg.AddEntry(obs, "Observed", "LP")
    pdefs.leg.AddEntry(m68, "Expected "+pdefs.oneSigmaText, "LF")
    pdefs.leg.AddEntry(m95, "Expected "+pdefs.twoSigmaText, "LF")
    if (pdefs.Draw_TEV):
        pdefs.leg.AddEntry(obsTEV, "Tevatron Observed", "L")
        pdefs.leg.AddEntry(expTEV, "Tevatron Expected", "L")
    
    if (pdefs.lep_excluded) : pdefs.leg.AddEntry(pdefs.fakeLEP, "LEP excluded", "F")
    if (pdefs.tev_excluded) : pdefs.leg.AddEntry(pdefs.fakeTEV, "Tevatron excluded", "F")
    if (pdefs.cms_excluded) : pdefs.leg.AddEntry(pdefs.fakeCMS, "CMS excluded", "F")
    if (pdefs.Draw_TEV): 
        expTEV.Draw("L SAME") 
        obsTEV.Draw("L SAME")
    pdefs.leg.Draw() # bockjoo Draws legend
    #time.sleep(10)
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(who)
    if (not pdefs.c1.GetGridx() and not pdefs.c1.GetGridy()) : frame0.Draw("AXIGSAME")
    if (pdefs.channelSpamOnRightHandSide):
        spam(channelFromName(who,False,True)+" only", .65, leg_y_lo-0.01, 0.93, leg_y_lo-0.06, 22)
    print "DEBUG drawOneCLs xmin ",xmin, " xmax ",xmax," ymin ",ymin, " ymax ",ymax, " myspam ",myspam
    # Line at Y = 1
    #line=ROOT.TLine(xmin,1,xmax,1) 
    #line.SetLineColor(pdefs.lineAt1Color) 
    #line.SetLineStyle(pdefs.lineAt1Style) 
    #line.SetLineWidth(4)
    #line.DrawClone()
    #print "DEBUG drawOneCLs saving the drawings as "+pdefs.globalPrefix+who+".png"
    #pdefs.c1.SaveAs(pdefs.globalPrefix+who+".png") # bockjoo
    #time.sleep(10)
    #return
    finalize(who,xmin,xmax,ymin,ymax,myspam)
    #return
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent()
        finalize(who+"_logx",xmin,xmax,ymin,ymax,myspam)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        setCanvas(frame2, "", ymin, ymax, what)
        if not (apriori is None) : apriori.Draw("E3")
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, smooth, smoothorder)
        frame2.Draw("AXIGSAME")
        obs.Draw("LP")    
        noMoreTev = pdefs.tev_excluded and (x_zoom <= 156)
        pdefs.leg = newLegend(leg_x_lo+0.06*noMoreTev*pdefs.isSquareCanvas, leg_y_lo+0.05*noMoreTev,.93+0.005*pdefs.isSquareCanvas, leg_y_hi) 
        if pdefs.isSquareCanvas or pdefs.lep_excluded :
            pdefs.leg.SetTextSize(0.034)
        else:
            pdefs.leg.SetTextSize(0.037)
        pdefs.leg.AddEntry(obs, "Observed", "LP")
        pdefs.leg.AddEntry(m68, "Expected "+pdefs.oneSigmaText, "LF")
        pdefs.leg.AddEntry(m95, "Expected "+pdefs.twoSigmaText, "LF")
        if (pdefs.Draw_TEV):
            pdefs.leg.AddEntry(obsTEV, "Tevatron Observed", "L")
            pdefs.leg.AddEntry(expTEV, "Tevatron Expected", "L")
        
        if (pdefs.lep_excluded) : pdefs.leg.AddEntry(pdefs.fakeLEP, "LEP excluded", "F")
        if (pdefs.tev_excluded and not noMoreTev) : pdefs.leg.AddEntry(pdefs.fakeTEV, "Tevatron excluded", "F")
        if (pdefs.cms_excluded) : pdefs.leg.AddEntry(pdefs.fakeCMS, "CMS excluded", "F")

        if miss68 : miss68.Draw("P")
        if miss : miss.Draw("P")
        if (pdefs.Draw_TEV): 
            expTEV.Draw("L SAME") 
            obsTEV.Draw("L SAME")
        if apriori : apriori.Draw("LX")
        pdefs.leg.Draw()
        ymin, ymax = minMaxY(obs, ymin, ymax, pdefs.x_zoom)
        if m68 : ymin2, ymax2 = minMaxY(m68, ymin2, ymax2, pdefs.x_zoom) 
        if (ymax2 > ymax) : ymax = ymax2 
        if (ymin2 < ymin) : ymin = ymin2;
        if (pdefs.Draw_TEV) : ymax *=2
        if (not pdefs.c1.GetGridx() and not pdefs.c1.GetGridy()) : frame2.Draw("AXIGSAME")
        if (pdefs.channelSpamOnRightHandSide):
            spam(channelFromName(who,False,True)+" only", .65, leg_y_lo-0.01, 0.93, leg_y_lo-0.06, 22)
        
        finalize(who+"_zoom",xmin,pdefs.x_zoom,ymin,ymax,myspam)
        if (pdefs.doZoom2):
            frame3 = ROOT.TH1D("frame3","frame3", 1, pdefs.x_zoom2_min,pdefs.x_zoom2_max) 
            frame3.Draw(); 
            setCanvas(frame3, "", ymin, ymax, what)
            if apriori : apriori.Draw("E3")
            draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
            frame3.Draw("AXIGSAME")
            obs.Draw("LP")    
            if miss68 : miss68.Draw("P")
            if miss : miss.Draw("P")
            if apriori : apriori.Draw("LX")
            pdefs.leg.Draw()
            finalize(who+"_zoom2",pdefs.x_zoom2_min,pdefs.x_zoom2_max,ymin,ymax,myspam)
        
    
    pdefs.leg = 0



# In[33]:


def drawTwoCLs(who, who2, name2="Other", postfix="", color2=ROOT.kBlue, doObs2=True):
    obs   = ROOT.gFile.Get(who+"_obs")
    if not obs : print who+"_obs not found "; return
    if obs is None : return
    m68 = ROOT.gFile.Get(who+"_median")
    m95 = ROOT.gFile.Get(who+"_median_95") 
    obs2 = ROOT.gFile.Get(who2+"_obs")
    exp2 = ROOT.gFile.Get(who2+"_median")
    
    
    if (obs2 is None or exp2 is None): 
        print "Missing " , who2 , " to compare with " , who  
        return
    obs2.SetLineColor(color2) 
    obs2.SetLineWidth(5); 
    obs2.SetLineStyle(1);
    exp2.SetLineColor(color2) 
    exp2.SetLineWidth(5); 
    exp2.SetLineStyle(2);

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    col = 1 ; dcol = 1
    obs.SetLineWidth(2) 
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8);
    ymin, ymax  = minMaxY(obs)
    ymin2, ymax2 = minMaxY(obs2)
    if (ymax2 > ymax) : ymax = ymax2 
    if (ymin2 < ymin) : ymin = ymin2;
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    frame0.GetYaxis().SetTitleOffset(1.10+0.25*pdefs.isSquareCanvas)
    what = "95% CL limit on #sigma/#sigma_{"+pdefs.SM+"}"
    #//if (who.Contains("acls_")) what = "Asymptotic "+what
    setCanvas(frame0, "", ymin, ymax, what)
    draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
    if not (m68 is None): 
        m68.SetFillColor(pdefs.color68) 
        m68.SetLineColor(pdefs.color50) 
        m68.SetLineStyle(2) 
        m68.SetLineWidth(2)
    if not (m95 is None): 
        m95.SetFillColor(pdefs.color95) 
        m95.SetLineColor(pdefs.color50) 
        m95.SetLineStyle(2) 
        m95.SetLineWidth(2)
    frame0.Draw("AXIGSAME")
    leg_y_hi = 0.94
    if pdefs.isSquareCanvas:
        leg_y_lo = 0.78- 0.05*(1+doObs2)
        leg_x_lo = .645-0.06
    else:
        leg_y_lo = .75- 0.05*(1+doObs2)
        leg_x_lo = .65
    pdefs.leg = newLegend(leg_x_lo, leg_y_lo,.93+0.005*pdefs.isSquareCanvas, leg_y_hi) 
    if pdefs.isSquareCanvas : pdefs.leg.SetTextSize(0.034)
    else : pdefs.leg.SetTextSize(0.037)
    obs.Draw("LP")    
    exp2.Draw("LX SAME") 
    if (doObs2) : obs2.Draw("LX SAME");
    pdefs.leg.AddEntry(obs, "Observed", "LP")
    pdefs.leg.AddEntry(m68, "Expected "+pdefs.oneSigmaText, "LF")
    pdefs.leg.AddEntry(m95, "Expected "+pdefs.twoSigmaText, "LF")
    if (doObs2) : pdefs.leg.AddEntry(obs2, name2+" Obs.", "L")
    pdefs.leg.AddEntry(exp2, name2+" Exp.", "L")
    pdefs.leg.Draw()
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(who)
    finalize(who+postfix,xmin,xmax,ymin,ymax,myspam)
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        finalize(who+postfix+"_logx",xmin,xmax,ymin,ymax,myspam)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        setCanvas(frame2, "", ymin, ymax, what)
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
        frame2.Draw("AXIGSAME")
        obs.Draw("LP")    
        exp2.Draw("LX SAME") 
        if (doObs2) : obs2.Draw("LX SAME")
        pdefs.leg.Draw()
        ymin, ymax = minMaxY(obs, ymin, ymax, pdefs.x_zoom) 
        ymin, ymax = minMaxY(exp2, ymin, ymax, pdefs.x_zoom)
        if (ymax2 > ymax) : ymax = ymax2 
        if (ymin2 < ymin) : ymin = ymin2
        finalize(who+postfix+"_zoom",xmin,pdefs.x_zoom,ymin,ymax,myspam)
    
    pdefs.leg = 0




# In[34]:


def drawMethods(who, include_pla=True, bands=True):
    obs   = ROOT.gFile.Get(who+"_obs")
    if not obs : print who+"_obs not found" ; return
    if obs is None : return
    m68 = ROOT.gFile.Get(who+"_median")
    m95 = ROOT.gFile.Get(who+"_median_95")
    m68.SetFillColor(211)
    m68.SetLineColor(39)
    m95.SetFillColor(220) 
    m95.SetLineColor(39)
    isToy = not ("acls" in who)

    what = "95% CL limit on #sigma/#sigma_{"+pdefs.SM+"}"

    pla = who 
    if (isToy): 
        pla.replace("cls_","acls_")
    else:
        pla.replace("acls_","pla_")
    print "PLA is " , pla 
    bayes = who
    if isToy : bayes.replace("cls_","bayes_") 
    else : bayes.replace("acls_","bayes_") 
    pla_obs  = ROOT.gFile.Get(pla+"_obs")
    bayes_obs = ROOT.gFile.Get(bayes+"_obs")
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    obs.SetLineWidth(2) 
    obs.SetLineColor(1) 
    obs.SetMarkerColor(1) 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8)
    ymin, ymax =minMaxY(obs);
    ymin2, ymax2 =minMaxY(m68); 
    if (ymax2 > ymax) : ymax = ymax2
    if (ymin2 < ymin) : ymin = ymin2
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    setCanvas(frame0, "", ymin, ymax, what)
    pdefs.leg = newLegend(.65-0.13*pdefs.isSquareCanvas,.75 - 0.04*bands+0.01*pdefs.isSquareCanvas,.95-0.02*pdefs.isSquareCanvas,.94) 
    pdefs.leg.SetTextSize(0.037-0.003*pdefs.isSquareCanvas)
    if (bands):
        m68.SetFillColor(pdefs.color68) 
        m68.SetLineColor(pdefs.color50); 
        m68.SetLineStyle(2); 
        m68.SetLineWidth(2);
        m95.SetFillColor(pdefs.color95) 
        m95.SetLineColor(pdefs.color50); 
        m95.SetLineStyle(2); 
        m95.SetLineWidth(2);
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
        frame0.Draw("AXIGSAME")
    else:
        m68.SetLineColor(222) 
        m68.SetLineStyle(9); 
        m68.SetLineWidth(4); 
        m68.Draw("LX")
     
    obs.Draw("LP") 
    if isToy:
        pdefs.leg.AddEntry(obs, "CLs Observed", "LP")
    else:
        pdefs.leg.AddEntry(obs, "Asym. CLs Obs.", "LP")
    if (bands):
        if isToy:
            pdefs.leg.AddEntry(m68, "CLs Expected "+pdefs.oneSigmaText, "LF")
            pdefs.leg.AddEntry(m95, "CLs Expected "+pdefs.twoSigmaText, "LF")
        else:
            pdefs.leg.AddEntry(m68, "CLs Exp. "+pdefs.oneSigmaText, "LF")
            pdefs.leg.AddEntry(m95, "CLs Exp. "+pdefs.twoSigmaText, "LF")
            
    else :
        if isToy:
            pdefs.leg.AddEntry(m68, "CLs Expected", "L")
        else:
            pdefs.leg.AddEntry(m68, "Asym. CLs Exp.", "L")
    
    if (bayes_obs):
        bayes_obs.SetLineColor(215)
        bayes_obs.SetMarkerColor(215)
        bayes_obs.SetLineStyle(2)
        bayes_obs.SetMarkerStyle(24)
        bayes_obs.SetMarkerSize(0.7)
        bayes_obs.SetLineWidth(3)
        bayes_obs.Draw("LP")
        pdefs.leg.AddEntry(bayes_obs, "Bayesian Observed", "LP")
    
    if (include_pla and pla_obs != 0): 
        pla_obs.SetLineColor(2)
        pla_obs.SetLineStyle(1)
        pla_obs.SetLineWidth(2)
        pla_obs.Draw("LX")
        if isToy:
            pdefs.leg.AddEntry(pla_obs, "Asymptotic CLs Obs.", "L")
        else:
            pdefs.leg.AddEntry(pla_obs, "PL Approx. Observed", "L")
    
    pdefs.leg.Draw()
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(who)
    xname = "_comp"
    if include_pla : xname = "_comp3"
    finalize(who+xname,xmin,xmax,ymin,ymax,myspam)
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        xcut = 140
        if pdefs.SM == "SM" : xcut = pdefs.x_zoom
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent()
        finalize(who+xname+"_logx",xmin,xmax,ymin,ymax,myspam)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, xcut) 
        frame2.Draw(); 
        setCanvas(frame2, "", ymin, ymax, what)
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
        frame2.Draw("AXIGSAME")
        obs.Draw("LP")
        if (bayes_obs) : bayes_obs.Draw("LP")
        if (include_pla and pla_obs != 0) : pla_obs.Draw("LX")
        ymin, ymax = minMaxY(obs, ymin, ymax, xcut)
        ymin2, ymax2 = minMaxY(m68, ymin2, ymax2, xcut) 
        if (ymax2 > ymax) : ymax = ymax2 
        if (ymin2 < ymin) : ymin = ymin2;
        pdefs.leg.Draw()
        finalize(who+xname+"_zoom",xmin,xcut,ymin,ymax,myspam)
    
    pdefs.leg = 0



# In[35]:


def drawSMCLs(who, showAs=False):
    obs   = ROOT.gFile.Get(who+"_obs")
    m68 = ROOT.gFile.Get(who+"_median")
    m95 = ROOT.gFile.Get(who+"_median_95")
    if not obs or not m68 or not m95: print "none of ",who+"_obs", who+"_median", who+"_median_95", " is found " ; return
    if (obs is None or m68 is None or m95 is None) : return
    m68.SetFillColor(211) 
    m68.SetLineColor(39)
    m95.SetFillColor(220) 
    m95.SetLineColor(39)
    #//savGridY = pdefs.c1.GetGridy() pdefs.c1.SetGridy(0); 
    savRightMargin = pdefs.c1.GetRightMargin() 
    pdefs.c1.SetRightMargin(0.08)

    aswho = who 
    aswho.replace("smcls_", "smacls_");
    if showAs: aser = ROOT.gFile.Get(aswho+"_obs")
    else: aser = None
    if (showAs and aser is None) : print "ERROR: cannot find " , aswho+"_obs"  

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    col = colorFromName(who);
    dcol = colorFromName(who, True)
    obs.SetLineWidth(2) 
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21)
    obs.SetMarkerSize(0.8);
    if not (aser is None):
        aser.SetLineColor(215)
        aser.SetMarkerColor(215)
        aser.SetLineStyle(2)
        aser.SetMarkerStyle(24)
        aser.SetMarkerSize(0.7)
        aser.SetLineWidth(3)
    
    smooth = -5
    ymin = 5e-4; ymax = 9
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    what = "CL_{S} of "+pdefs.SM+" Higgs hypothesis"
    setCanvas(frame0, "", ymin, ymax, what)
    draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, smooth)
    frame0.Draw("AXIGSAME")
    obs.Draw("LP")
    if not (aser is None) : aser.Draw("LP")
    #//spam("#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32)
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(who)
    spam(myspam, 0.19,.88,.58-0.1*pdefs.isSquareCanvas,.93)
    #//leg = newLegend(.64,.15,.88,.31+showAs*0.04) pdefs.leg.SetTextSize(0.04);
    pdefs.leg = newLegend(.67-0.11*pdefs.isSquareCanvas,.93-.16-showAs*0.04,.91,.93) 
    pdefs.leg.SetTextSize(0.04);
    m68.SetFillColor(pdefs.color68)
    m68.SetLineColor(pdefs.color50);
    m68.SetLineStyle(2);
    m68.SetLineWidth(2);
    m95.SetFillColor(pdefs.color95)
    m95.SetLineColor(pdefs.color50); 
    m95.SetLineStyle(2);
    m95.SetLineWidth(2);
    pdefs.leg.AddEntry(obs, "Observed", "LP")
    pdefs.leg.AddEntry(m68, "Expected "+pdefs.oneSigmaText, "LF")
    pdefs.leg.AddEntry(m95, "Expected "+pdefs.twoSigmaText, "LF")
    if (showAs) : pdefs.leg.AddEntry(aser, "Asymptotic Obs.", "LP")
    frame0.Draw("AXIS SAME") #// re-draw ticks
    if showAs: finalize(who+"_comp",xmin,xmax,ymin,ymax)
    else: finalize(who,xmin,xmax,ymin,ymax)
    #//finalize(who,xmin,xmax,ymin,ymax,SPAM,True)
    if (xmax > pdefs.x_zoom and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        #//finalize(who+"_logx",xmin,xmax,ymin,ymax,SPAM,True)
        if showAs: finalize(who+"_comp_logx",xmin,xmax,ymin,ymax)
        else: finalize(who+"_logx",xmin,xmax,ymin,ymax)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        setCanvas(frame2, "", ymin, ymax, what)
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, smooth)
        frame2.Draw("AXIGSAME") #// grid
        obs.Draw("LP")
        if not (aser is None) : aser.Draw("LP")
        #//spam("#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32)
        spam(myspam, 0.19,.88,.58-0.1*pdefs.isSquareCanvas,.93)
        frame2.Draw("AXIS SAME") #// ticks
        if showAs: finalize(who+"_comp_zoom",xmin,xmax,ymin,ymax)
        else: finalize(who+"_zoom",xmin,xmax,ymin,ymax)
    
    pdefs.c1.SetRightMargin(savRightMargin) #//pdefs.c1.SetGridy(savGridY); 



# In[36]:


def drawTwoSMCLs(who, who2, name2="Other", postfix="",color2=ROOT.kBlue):
    obs   = ROOT.gFile.Get(who+"_obs")
    m68 = ROOT.gFile.Get(who+"_median")
    m95 = ROOT.gFile.Get(who+"_median_95")
    if not obs or not m68 or not m95: print "none of ",who+"_obs", who+"_median", who+"_median_95", " is found " ; return
    if (obs is None or m68 is None or m95 is None) : return
    
    m68.SetFillColor(211) 
    m68.SetLineColor(39);
    m95.SetFillColor(220) 
    m95.SetLineColor(39);
    pdefs.c1.SetRightMargin(0.08)

    obs2 = ROOT.gFile.Get(who2+"_obs")
    exp2 = ROOT.gFile.Get(who2+"_median")
    if (obs2 is None or exp2 is None): 
        print "Missing " , who2 , " to compare with " , who
        return
    obs2.SetLineColor(color2)
    obs2.SetLineWidth(5); 
    obs2.SetLineStyle(1);
    exp2.SetLineColor(color2)
    exp2.SetLineWidth(5); 
    exp2.SetLineStyle(2);

    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    col = colorFromName(who)
    dcol = colorFromName(who, True)
    obs.SetLineWidth(2) 
    obs.SetLineColor(col) 
    obs.SetMarkerColor(col); 
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(0.8);
    ymin = 8e-5
    ymax = 1
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    what = "CL_{S} of "+pdefs.SM+" Higgs hypothesis"
    setCanvas(frame0, "", ymin, ymax, what)
    draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
    frame0.Draw("AXIGSAME")
    obs.Draw("LP")
    exp2.Draw("LX SAME") 
    obs2.Draw("LX SAME");
    spam("#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32)
    pdefs.leg = newLegend(.64,.15,.88,.37) 
    pdefs.leg.SetTextSize(0.04);
    m68.SetFillColor(pdefs.color68)
    m68.SetLineColor(pdefs.color50); 
    m68.SetLineStyle(2);
    m68.SetLineWidth(2);
    m95.SetFillColor(pdefs.color95) 
    m95.SetLineColor(pdefs.color50);
    m95.SetLineStyle(2); 
    m95.SetLineWidth(2);
    pdefs.leg.AddEntry(obs, "Observed", "LP")
    pdefs.leg.AddEntry(m68, "Expected "+pdefs.oneSigmaText, "LF")
    pdefs.leg.AddEntry(m95, "Expected "+pdefs.twoSigmaText, "LF")
    pdefs.leg.AddEntry(obs2, name2+" Observed", "L")
    pdefs.leg.AddEntry(exp2, name2+" Expected", "L")
    finalize(who+postfix,xmin,xmax,ymin,ymax)
    #//finalize(who,xmin,xmax,ymin,ymax,SPAM,True)
    if (xmax >= 200 and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        #//finalize(who+"_logx",xmin,xmax,ymin,ymax,SPAM,True)
        finalize(who+postfix+"_logx",xmin,xmax,ymin,ymax)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        setCanvas(frame2, "", ymin, ymax, what)
        draw2(who, pdefs.color68, pdefs.color95, pdefs.color50, True, False, 5)
        frame2.Draw("AXIGSAME")
        obs.Draw("LP")
        exp2.Draw("LX SAME") 
        obs2.Draw("LX SAME");
        spam("#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}", 0.17,.15,.54,.20,32)
        finalize(who+postfix+"_zoom",xmin,pdefs.x_zoom,ymin,ymax)
    
    pdefs.c1.SetRightMargin(0.04)




# In[37]:


def findCrossings(who, xname, threshold, xmin, xmax, what="95% CL limit on #sigma/#sigma_{SM}"):
    obs   = ROOT.gFile.Get(who)
    if (who == 0) : return
    ihigh = -1
    ilow = -1
    n = obs.GetN()
    for i in xrange(n): # for (int i = 0, n = obs.GetN() i < n; ++i):
        if (xmin <= obs.GetX()[i] and obs.GetX()[i] <= xmax):
            if ("low" in xname):
                if (obs.GetY()[i] < threshold):
                    ilow = i
                    break
                else : 
                    ihigh = i
            else :
                if (obs.GetY()[i] < threshold):
                    ilow = i
                else : 
                    ihigh = i
                    break
    if (ilow == -1 or ihigh == -1): 
        print "didn't find points."
        return
    x1 = obs.GetX()[ilow]
    x2 = obs.GetX()[ihigh]
    y1 = obs.GetY()[ilow]
    y2 = obs.GetY()[ihigh]
    linear  = ROOT.TF1("linear", "[0] * (x - [1]) + [2]", xmin, xmax) 
    linear.SetParameters((y2-y1)/(x2-x1), x1, y1)
    linlog = ROOT.TF1("linlog", "[0] * pow([1], (x - [2])/[3])", xmin, xmax)
    linlog.SetParameters(y1, y2/y1, x1, x2-x1)
    loglog  = ROOT.TF1("linear", "[0] * pow([1], log(x/[2])/[3])", xmin, xmax)
    loglog.SetParameters(y1, y2/y1, x1, log(x2/x1))
    loglog.SetLineWidth(3) 
    loglog.SetLineColor( 63); 
    linlog.SetLineWidth(5) 
    linlog.SetLineColor(210); 
    linlog.SetLineStyle(2);
    linear.SetLineWidth(3) 
    linear.SetLineColor(  1);  
    linear.SetLineStyle(9);
    obs.SetMarkerStyle(21) 
    obs.SetMarkerSize(1.3);
    x_linear = linear.GetX(threshold, xmin, xmax)
    x_linlog = linlog.GetX(threshold, xmin, xmax)
    x_loglog = loglog.GetX(threshold, xmin, xmax)
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    linear.Draw("SAME")  
    linlog.Draw("SAME");    
    loglog.Draw("SAME");
    obs.Draw("PX")    
    ymin = threshold/2
    ymax = threshold * 3
    xleg = 0.67
    yleg = 0.18
    if ("smcls" in who): 
        ymax = 3*threshold 
        ymin = 0;  
        xleg = 0.27 
        yleg = 0.55;
    
    if ("low" in xname):
        xleg = 0.67 
        yleg = 0.66
    
    pdefs.leg = newLegend(xleg,yleg,xleg+.28,yleg+.25) 
    pdefs.leg.SetTextSize(0.04);
    pdefs.leg.AddEntry(linear, Form("Linear,  m_{H} = %.1f", x_linear), "L")
    pdefs.leg.AddEntry(linlog, Form("Lin-log, m_{H} = %.1f", x_linlog), "L")
    pdefs.leg.AddEntry(loglog, Form("Log-log, m_{H} = %.1f", x_loglog), "L")
    setCanvas(frame0, "", ymin, ymax, what)
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(who)+"}"
    finalize(who+"_crossing_"+xname,xmin,xmax,ymin,ymax, myspam)



# In[38]:


def printInfernalTable(fileName, xmin, xmax):
    pvala_obs = ROOT.gFile.Get("pvala_comb_obs")
    smcls_obs = ROOT.gFile.Get("smcls_comb_obs")
    smcls_exp = ROOT.gFile.Get("smcls_comb_median")
    cls_obs   = ROOT.gFile.Get("cls_comb_obs")
    cls_exp   = ROOT.gFile.Get("cls_comb_median")
    bayes_obs = ROOT.gFile.Get("bayes_comb_obs")
    if not cls_obs : print "cls_comb_obs not found "; return
    if (cls_obs is None) : return
    
    fout = open(fileName, "w")
    #fprintf(fout, 
    #    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n")
    #fprintf(fout, 
    #    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n")
    fout.write(    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n")
    fout.write(    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n")
    #//printf( 
    #//    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n")
    #//printf( 
    #//    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n")
    n = cls_obs.GetN()
    for i in xrange(n) : # for (int i = 0, n = cls_obs.GetN() i < n; ++i):
        mass = cls_obs.GetX()[i] 
        y_cls_obs = cls_obs.GetY()[i];
        if (mass < xmin or mass > xmax) :continue
        ipvala_obs = findBin(pvala_obs, mass) 
        if ipvala_obs != -1 : y_pvala_obs = pvala_obs.GetY()[ipvala_obs]
        else : y_pvala_obs = -1
            
        ismcls_obs = findBin(smcls_obs, mass) 
        if ismcls_obs != -1 : y_smcls_obs = smcls_obs.GetY()[ismcls_obs]
        else : y_smcls_obs = -1
        
        ismcls_exp = findBin(smcls_exp, mass) 
        if ismcls_exp != -1 : y_smcls_exp = smcls_exp.GetY()[ismcls_exp]
        else : y_smcls_exp = -1
        
        icls_exp   = findBin(cls_exp, mass)   
        if icls_exp != -1 : y_cls_exp = cls_exp.GetY()[icls_exp]
        else : y_cls_exp = -1.
        
        ibayes_obs = findBin(bayes_obs, mass) 
        if ibayes_obs != -1 : y_bayes_obs = bayes_obs.GetY()[ibayes_obs]
        else : y_bayes_obs = -1
        
    


# In[39]:


def printInfernalTable(fileName, xmin, xmax):
    pvala_obs = ROOT.gFile.Get("pvala_comb_obs")
    smcls_obs = ROOT.gFile.Get("smcls_comb_obs")
    smcls_exp = ROOT.gFile.Get("smcls_comb_median")
    cls_obs   = ROOT.gFile.Get("cls_comb_obs")
    cls_exp   = ROOT.gFile.Get("cls_comb_median")
    bayes_obs = ROOT.gFile.Get("bayes_comb_obs")
    if not cls_obs : print "cls_comb_obs not found " ; return
    if (cls_obs is None) : return
    
    fout = open(fileName, "w")
    #fprintf(fout, 
    #    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n")
    #fprintf(fout, 
    #    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n")
    fout.write(    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n")
    fout.write(    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n")
    #//printf( 
    #//    "\tHiggs boson &   observed        &   observed (expected)         &   \\multicolumn {2} {c|} {observed (expected) 95\\%% C.L. limits on $\\mu$}      \\\\\n")
    #//printf( 
    #//    "\tmass (GeV)  & $\\tilde p$-value  & $\\mathrm{CL_s}$ for $\\mu=1$   &   $\\mathrm{CL_s}$-method  &  Bayesian                  \\\\ \\hline \n")
    n = cls_obs.GetN()
    for i in xrange(n) : # for (int i = 0, n = cls_obs.GetN() i < n; ++i):
        mass = cls_obs.GetX()[i] 
        y_cls_obs = cls_obs.GetY()[i];
        if (mass < xmin or mass > xmax) :continue
        ipvala_obs = findBin(pvala_obs, mass) 
        if ipvala_obs != -1 : y_pvala_obs = pvala_obs.GetY()[ipvala_obs]
        else : y_pvala_obs = -1
            
        ismcls_obs = findBin(smcls_obs, mass) 
        if ismcls_obs != -1 : y_smcls_obs = smcls_obs.GetY()[ismcls_obs]
        else : y_smcls_obs = -1
        
        ismcls_exp = findBin(smcls_exp, mass) 
        if ismcls_exp != -1 : y_smcls_exp = smcls_exp.GetY()[ismcls_exp]
        else : y_smcls_exp = -1
        
        icls_exp   = findBin(cls_exp, mass)   
        if icls_exp != -1 : y_cls_exp = cls_exp.GetY()[icls_exp]
        else : y_cls_exp = -1.
        
        ibayes_obs = findBin(bayes_obs, mass) 
        if ibayes_obs != -1 : y_bayes_obs = bayes_obs.GetY()[ibayes_obs]
        else : y_bayes_obs = -1
        
        #fprintf(fout,
        #    "\t%3d         &    {\\bf{%5.3f}}   &     {\\bf{%5.3f}} (%5.3f)        &     {\\bf{%5.2f}} (%4.2f)    &    {\\bf{%5.2f}}       \\\\ \n",
        #    mass, y_pvala_obs, y_smcls_obs, y_smcls_exp, y_cls_obs, y_cls_exp, y_bayes_obs)
        fout.write("\t%3d         &    {\\bf{%5.3f}}   &     {\\bf{%5.3f}} (%5.3f)        &     {\\bf{%5.2f}} (%4.2f)    &    {\\bf{%5.2f}}       \\\\ \n" % (
            mass, y_pvala_obs, y_smcls_obs, y_smcls_exp, y_cls_obs, y_cls_exp, y_bayes_obs))
        #//printf(
        #//    "\t%3d         &    {\\bf{%5.3f}}   &     {\\bf{%5.3f}} (%5.3f)        &     {\\bf{%5.2f}} (%4.2f)    &    {\\bf{%5.2f}}       \\\\ \n",
        #//    mass, y_pvala_obs, y_smcls_obs, y_smcls_exp, y_cls_obs, y_cls_exp, y_bayes_obs)
        #}
    fout.close()



# In[40]:


#PLC_debug = False
#def drawCombObs(who, what="P.L. Approx limit #sigma_{95%}/#sigma_{"+pdefs.SM+"}", chann, nchann, combined, postfix=""):
def drawCombObs(who, what="P.L. Approx limit #sigma_{95%}/#sigma_{"+pdefs.SM+"}", chann=pdefs.chann, nchann=pdefs.nchann, combined=False, postfix=""):
    obs   = ROOT.gFile.Get(who+"_comb_obs")
    if not obs: print who+"_comb_obs not found " ; return
    if obs is None : return
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    isC3 = (postfix == "_c3")
    ymin, ymax =minMaxY(obs)
    if isC3: ymax *= 2
    else: ymax *= 8
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    #*obsi[pdefs.nchann_all]
    if pdefs.isSquareCanvas : 
        if pdefs.lessSpacesInLegends : legxstart = 0.60
        else : legxstart = 0.50 
    else :
        if pdefs.lessSpacesInLegends : legxstart = 0.72
        else : legxstart = 0.62
    pdefs.leg = newLegend(legxstart,.65+0.10*isC3,.95-0.02*pdefs.isSquareCanvas,.94) 
    pdefs.leg.SetTextSize(0.032);
    ntot = 0 ; nlow = 0
    for i in xrange(nchann) : # for (int i = 0 i < nchann; ++i):
        if (i == 0 and not combined) : continue
        obsi[i] = ROOT.gFile.Get(who+"_"+chann[i]+"_obs")
        if (obsi[i] == 0) : continue
        col = colorFromName(chann[i])
        #//obsi[i].SetLineWidth(2)
        if i == 0 : obsi[i].SetLineWidth(4)
        else : obsi[i].SetLineWidth(3)
        obsi[i].SetLineColor(col) 
        obsi[i].SetMarkerColor(col); 
        obsi[i].SetMarkerStyle(21) 
        obsi[i].SetMarkerSize(0.8);
        obsi[i].SetLineStyle(lineStyleFromName(chann[i]))
        if not ("_low" in str(chann[i])):
            if "cls" in who : pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "L")
            else : pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "LP")
            ntot += 1
        
        if (obsi[i].GetX()[0] < pdefs.x_zoom) : nlow += 1
    
    #for (int i = nchann-1 i >= 0; --i):
    for i in range(nchann-1,-1,-1):
        if (i == 0 and  not combined) : continue
        if (obsi[i] == 0) : continue
        if "cls" in who: obsi[i].Draw("LX")
        else : obsi[i].Draw("LPX")
    
    setCanvas(frame0, "", ymin, ymax, what)
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(chann[i])+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(chann[0])
    spam(myspam, 0.17,.89,.58,.94)
    pdefs.leg.Draw()
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax)
    if (xmax >= 200 and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        finalize(who+"_all"+postfix+"_logx",xmin,xmax,ymin,ymax)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        pdefs.leg = newLegend(legxstart,.65+0.10*isC3+(ntot-nlow)*0.032,.95-0.02*pdefs.isSquareCanvas,.94) 
        pdefs.leg.SetTextSize(0.032);
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        ymin, ymax = minMaxY(obs, ymin, ymax, 200.) 
        ymax *= 15;
        for i in range(nchann-1,-1,-1):
            if (i == 0 and not combined) : continue
            if (obsi[i] == 0) : continue
            if (obsi[i].GetX()[1] > pdefs.x_zoom) : continue
            #obsi[i].Draw(who.Contains("cls") ? "LX" : "LPX")
            if "cls" in who: obsi[i].Draw("LX")
            else : obsi[i].Draw("LPX")
    
        
        #for (int i = 0 i < nchann; ++i):
        for i in xrange(nchann) : # for (int i = 0 i < nchann; ++i):
            if (i == 0 and not combined) : continue
            if (obsi[i] == 0) : continue
            if (obsi[i].GetX()[1] > pdefs.x_zoom) : continue
            if "cls" in who: pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "L")
            else : pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "LP")
        pdefs.leg.Draw()
        #//spam(SPAM, 0.16,.85,.56,.91)
        spam(myspam, 0.17,.89,.58,.94)
        setCanvas(frame2, "", ymin, ymax, what)
        finalize(who+"_all"+postfix+"_zoom",xmin,pdefs.x_zoom,ymin,ymax)
    
    pdefs.leg = 0



# In[41]:


def drawCombBoth(who, what="P.L. Approx limit #sigma_{95%}/#sigma_{"+pdefs.SM+"}", chann=pdefs.chann, nchann=pdefs.nchann, observed=True, combined=True, postfix=""):
    obs   = ROOT.gFile.Get(who+"_comb_obs")
    if not obs: print who+"_comb_obs", "not found in ",ROOT.gFile.GetName() ; return
    if obs is None : return
    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    isC3 = (postfix == "_c3")
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    ymin, ymax = minMaxY(obs); 
    if isC3: ymax *= 2
    else: ymax *= 8
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    #*obsi[pdefs.nchann_all], *expi[pdefs.nchann_all]
    if pdefs.isSquareCanvas :
        if pdefs.lessSpacesInLegends : legxstart = 0.60
        else : legxstart = 0.50
    else :
        if pdefs.lessSpacesInLegends : legxstart = 0.72
        else : legxstart = 0.62
    legystart = .65-0.05*combined+0.14*isC3+0.03*pdefs.lessSpacesInLegends+0.06*(SM=="FP")
    pdefs.leg = newLegend(legxstart, legystart, .95-0.02*pdefs.isSquareCanvas, .94)
    pdefs.leg.SetTextSize(0.032);
    if (not observed) : pdefs.leg.SetHeader("Expected limits")
    ntot = 0 ; nlow = 0
    for i in xrange(nchann) : # for (int i = 0 i < nchann; ++i):
        if (i == 0 and not combined) : continue
        obsi[i] = ROOT.gFile.Get(who+"_"+chann[i]+"_obs")
        expi[i] = ROOT.gFile.Get(who+"_"+chann[i]+"_median")
        if (expi[i] is None or obsi[i] is None) : continue
        if ("hzz2l2q" in chann[i] or 
            (pdefs.noLineStyles == True and "comb" in chann[i]) or
            (pdefs.noLineStyles == True and "hzz" in chann[i])):
            smooth = smoothSMCLs(expi[i], 7, 2)
            if (smooth is None) : print "Smoothing of " , expi[i].GetName() , " returned ZERO" 
            else : expi[i] = smooth
        
        col = colorFromName(chann[i]) #//if (col == 96) col = 67;
        if i == 0 : obsi[i].SetLineWidth(4)
        else : obsi[i].SetLineWidth(3)
        obsi[i].SetLineColor(col)  
        obsi[i].SetMarkerColor(col); 
        obsi[i].SetMarkerStyle(21)
        obsi[i].SetMarkerSize(0.8);
        if i == 0 :
            expi[i].SetLineWidth(4)
        else : 
            if observed: 
               expi[i].SetLineWidth(3)
            else:
               expi[i].SetLineWidth(4)
        expi[i].SetLineColor(col) 
        expi[i].SetMarkerColor(col); 
        expi[i].SetMarkerStyle(21) 
        expi[i].SetMarkerSize(0.8);
        if observed or  pdefs.noLineStyles : obsi[i].SetLineStyle(1)
        else : obsi[i].SetLineStyle(2)
        if observed or not pdefs.noLineStyles : expi[i].SetLineStyle(2)
        else : expi[i].SetLineStyle(1)
        if (i == 0 and not (expi[i] is None) ):
            if (observed):
                if pdefs.lessSpacesInLegends :
                    if "cls" in who:
                        pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" obs.", "L")
                        pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" exp.", "L")
                    else:
                        pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" obs.", "LP")
                        pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" exp.", "LP")
                else:
                    if "cls" in who:
                        pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" observed", "L")
                        pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" expected", "L")
                    else:
                        pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" observed", "LP")
                        pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" expected", "LP")
            
            else:
                if "cls" in who:
                    pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True), "L")
                else:
                    pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True), "LP")
            ntot +=1
        else:
            if not ("_low" in chann[i]):
                if "cls" in who:
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "L")
                else:
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "LP")
                ntot +=1          
        
        if (obsi[i].GetX()[0] < pdefs.x_zoom) : nlow += 1
    
    #for (int i = nchann-1 i >= 0; --i):
    for i in range(nchann-1,-1,-1):
        if (expi[i] is None or obsi[i] is None) : continue
        if "cls" in who: expi[i].Draw("LX")
        else: expi[i].Draw("LPX")
        if (observed) :
            if "cls" in who: obsi[i].Draw("LX")
            else: obsi[i].Draw("LPX")    
    setCanvas(frame0, "", ymin, ymax, what)
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(chann[0])+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(chann[0])
    spam(myspam, 0.17,.89,.58,.94)
    pdefs.leg.Draw()
    if observed: finalize(who+"_all2"+postfix,xmin,xmax,ymin,ymax)
    else : finalize(who+"_allexp"+postfix,xmin,xmax,ymin,ymax)
    if (xmax >= 200 and xmin < pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        if observed: finalize(who+"_all2"+postfix+"_logx",xmin,xmax,ymin,ymax)
        else : finalize(who+"_allexp"+postfix+"_logx",xmin,xmax,ymin,ymax)
        pdefs.c1.SetLogx(0)
        xmin = xmin0
        pdefs.leg = newLegend(legxstart,legystart+(ntot-nlow)*0.032,.95-0.02*pdefs.isSquareCanvas,.94) 
        pdefs.leg.SetTextSize(0.032);
        if (not observed) : pdefs.leg.SetHeader("Expected limits")
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom) 
        frame2.Draw(); 
        ymin, ymax=minMaxY(obs, ymin, ymax, 200.) 
        if isC3: ymax *= 2
        else: ymax *= 15
        for i in range(nchann-1,-1,-1): #for (int i = nchann-1 i >= 0; --i):
            if (i == 0 and not combined) : continue
            if (expi[i] is None or obsi[i] is None) : continue
            if (obsi[i].GetX()[1] > pdefs.x_zoom) : continue
            if "cls" in who: expi[i].Draw("LX")
            else: expi[i].Draw("LPX")
            if (observed) :
               if "cls" in who: obsi[i].Draw("LX")
               else: obsi[i].Draw("LPX")    
            
        for i in xrange(nchann) : # for (int i = 0 i < nchann; ++i):
            if (i == 0 and not combined) : continue
            if (expi[i] is None or obsi[i] is None) : continue
            if (obsi[i].GetX()[1] > pdefs.x_zoom) : continue
            if (i == 0 and not ( expi[i] is None)):
                if (observed):
                    #pdefs.leg.AddEntry(obsi[i], channelFromName(pdefs.chann[i],/*withspaces=*/True)+(pdefs.lessSpacesInLegends?" obs.":" observed"), who.Contains("cls") ? "L" : "LP")
                    #pdefs.leg.AddEntry(expi[i], channelFromName(pdefs.chann[i],/*withspaces=*/True)+(pdefs.lessSpacesInLegends?" exp.":" expected"), who.Contains("cls") ? "L" : "LP")
                    if pdefs.lessSpacesInLegends :
                        if "cls" in who:
                            pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" obs.", "L")
                            pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" exp.", "L")
                        else:
                            pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" obs.", "LP")
                            pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" exp.", "LP")
                    else:
                        if "cls" in who:
                            pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" observed", "L")
                            pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" expected", "L")
                        else:
                            pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" observed", "LP")
                            pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True)+" expected", "LP")

                else:
                    #pdefs.leg.AddEntry(expi[i], channelFromName(pdefs.chann[i],/*withspaces=*/True), who.Contains("cls") ? "L" : "LP")
                    if "cls" in who:
                        pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True), "L")
                    else:
                        pdefs.leg.AddEntry(expi[i], channelFromName(chann[i],True), "LP")
                
            else:
                #pdefs.leg.AddEntry(obsi[i], channelFromName(pdefs.chann[i],/*withspaces=*/True), who.Contains("cls") ? "L" : "LP")
                if "cls" in who:
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "L")
                else:
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "LP")
            
        
        pdefs.leg.Draw()
        spam(myspam, 0.17,.89,.58,.94)
        setCanvas(frame2, "", ymin, ymax, what)
        #finalize(who+(observed?"_all2"+postfix+"_zoom":"_allexp"+postfix+"_zoom"),xmin,pdefs.x_zoom,ymin,ymax)
        if observed: finalize(who+"_all2"+postfix+"_zoom",xmin,xmax,ymin,ymax)
        else : finalize(who+"_allexp"+postfix+"_zoom",xmin,xmax,ymin,ymax)
    
    pdefs.leg = 0



# In[42]:


def drawCombPVal(who, chann, nchann, postfix="", toywho=""):
    obs = ROOT.gFile.Get(who+"_"+chann[0]+"_obs")
    if not obs : print who+"_"+chann[0]+"_obs is null in ",ROOT.gFile.GetName()
    if not obs : return
    #if obs is None : print "DEBUG "+who+"_"+chann[0]+"_obs is null"
    #if obs is None : return
    exp = ROOT.gFile.Get(who+"_"+chann[0]+"_median")
    if (exp is None) : exp =ROOT.gFile.Get(who+"_"+chann[0]+"_asimov")

    #TGraphAsymmErrors *toys = toywho != "" ? (TGraphAsymmErrors *) gFile->Get(toywho+"_obs") : 0;
    if toywho != "" : toys = ROOT.gFile.Get(toywho+"_obs")
    else : toys = 0
    if (toys):
        toys.SetLineWidth(4)
        toys.SetMarkerStyle(20)
        if toys.GetN() > 3 and pdefs.SM == "SM" : toys.SetMarkerSize(1.3)
        else : toys.SetMarkerSize(1.8)
        toys.SetLineColor(92)
        toys.SetMarkerColor(92)
    


    xmin = obs.GetX()[0]             - obs.GetErrorXlow(0)
    xmin0 = xmin
    xmax = obs.GetX()[obs.GetN()-1] + obs.GetErrorXhigh(obs.GetN()-1)
    isC3 = (postfix == "_c3")
    if (xmin <= 120 and xmax > 200):
        xmin = 99.98
        if pdefs.SM == "FP" : xmax = 300.1
        else : xmax = 600.1
    if ("_toy_" in chann[0]): 
        xmin = pdefs.x_zoom2_min 
        xmax = pdefs.x_zoom2_max
    ymin, ymax = minMaxY(obs, ymin, ymax, 999, 0)
    ymax = 1.0; ymin = 8e-7;
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    print "DEBUG setCanvas  for Local p-value"
    setCanvas(frame0, "", ymin, ymax, "Local p-value")
    #*obsi[pdefs.nchann_all]
    if pdefs.isSquareCanvas :
        if pdefs.lessSpacesInLegends : legxstart = 0.60
        else : legxstart = 0.50
    else :
        if pdefs.lessSpacesInLegends : legxstart = 0.72
        else : legxstart = 0.62
    legyend = .43-0.10*isC3 
    if (toys) : legyend += 0.04
    pdefs.leg = newLegend(legxstart,.15,.95-0.02*pdefs.isSquareCanvas,legyend)
    pdefs.leg.SetTextSize(0.032);
    for i in xrange(nchann) : #for (int i = 0 i < nchann; ++i):
        obsi[i] = ROOT.gFile.Get(who+"_"+chann[i]+"_obs")
        if (obsi[i] is None) : continue
        col = colorFromName(chann[i])
        #//obsi[i].SetLineWidth(2)
        #obsi[i].SetLineWidth(i == 0 ? 4 : 4)
        if i == 0 : obsi[i].SetLineWidth(4)
        else : obsi[i].SetLineWidth(4)
        obsi[i].SetLineColor(col) 
        obsi[i].SetMarkerColor(col); 
        obsi[i].SetMarkerStyle(21) 
        obsi[i].SetMarkerSize(0.8);
        obsi[i].SetLineStyle(lineStyleFromName(chann[i]))
        if not ("_low" in chann[i]): #if (!TString(pdefs.chann[i]).Contains("_low")):
            if (i == 0 and not (exp is None)):
                if pdefs.lessSpacesInLegends :
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" obs.", "L")
                    pdefs.leg.AddEntry(exp,     "Exp. for SM Higgs", "L")
                else:
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" observed", "L")
                    pdefs.leg.AddEntry(exp,     "Expected for SM Higgs", "L")
                if (toys) : pdefs.leg.AddEntry(toys, "Comb. ensemble", "P")
            else : 
                pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "L")
            
        
    
    if (not (exp is None)):
        exp.SetLineWidth(4) 
        exp.SetLineStyle(7); 
        exp.SetLineColor(1);
        smooth = smoothSMCLs(exp, 7, 2)
        if (smooth == 0) : print "Smoothing of " , expi[i].GetName() , " returned ZERO"
        else : exp = smooth
        exp.Draw("LX")
    
    if (toys) : toys.Draw("P")
    for i in range(nchann-1,-1,-1): # for (int i = nchann-1 i >= 0; --i):
        if (obsi[i] == 0) : continue
        obsi[i].Draw("LX")
    
    pdefs.leg.Draw()
    myspam = "#splitline{"+pdefs.SPAM+"}{"+channelFromName(chann[0])+"}"
    if (pdefs.isSquareCanvas and "Preliminary" in pdefs.SPAM) : myspam = pdefs.SPAM2L + "\n" + channelFromName(chann[0])
    finalize(who+"_all"+postfix,xmin,xmax,ymin,ymax,myspam,True)
    if (xmax >= pdefs.x_zoom and xmin <= pdefs.x_zoom):
        pdefs.c1.SetLogx(1)
        frame0.GetXaxis().SetMoreLogLabels() 
        frame0.GetXaxis().SetNoExponent();
        finalize(who+"_all"+postfix+"_logx",xmin,xmax,ymin,ymax,myspam,True)
        pdefs.c1.SetLogx(0)
        pdefs.leg = newLegend(legxstart,.153,.95-0.02*pdefs.isSquareCanvas,legyend) 
        pdefs.leg.SetTextSize(0.031);
        xmin = xmin0
        frame2 = ROOT.TH1D("frame2","frame2", 1, xmin, pdefs.x_zoom)
        frame2.Draw(); 
        setCanvas(frame2, "", ymin, ymax, "Local p-value")
        if not (exp is None) : exp.Draw("LX")
        if (toys) : toys.Draw("P")
        for i in range(nchann-1,-1,-1): # for (int i = nchann-1 i >= 0; --i):
            if (obsi[i] == 0) : continue
            if (obsi[i].GetX()[1] > pdefs.x_zoom) :continue
            obsi[i].Draw("LX")
        
        for i in xrange(nchann) : #for (int i = 0 i < nchann; ++i):
            if (obsi[i] == 0) : continue
            if (obsi[i].GetX()[1] > pdefs.x_zoom) : continue
            if (i == 0 and not ( exp is None)):
                #pdefs.leg.AddEntry(obsi[i], channelFromName(pdefs.chann[i],True)+(pdefs.lessSpacesInLegends?" obs.":" observed"), "L")
                #pdefs.leg.AddEntry(exp,     (pdefs.lessSpacesInLegends?"Exp. for SM Higgs":"Expected for SM Higgs"), "L")
                if pdefs.lessSpacesInLegends :
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" obs.", "L")
                    pdefs.leg.AddEntry(exp,     "Exp. for SM Higgs", "L")
                else:
                    pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True)+" observed", "L")
                    pdefs.leg.AddEntry(exp,     "Expected for SM Higgs", "L")
                if (toys) : pdefs.leg.AddEntry(toys, "Comb. ensemble", "P")
            else :
                pdefs.leg.AddEntry(obsi[i], channelFromName(chann[i],True), "L")
            
        
        pdefs.leg.Draw()
        finalize(who+"_all"+postfix+"_zoom",xmin,pdefs.x_zoom,ymin,ymax,myspam,True)
    
    pdefs.leg = 0



# In[43]:


def getNSigmaLine(g, sigma):
    if not g: print "g is 0" ; return 0
    if (g is None) : return 0
    if (g == 0) : return 0
    n = g.GetN()
    ret = ROOT.TGraphAsymmErrors(n)
    for i in xrange(n): # for (int i = 0 i < n; ++i):
        if sigma > 0 : ret.SetPoint(i, g.GetX()[i], g.GetY()[i] + g.GetErrorYhigh(i))
        else : ret.SetPoint(i, g.GetX()[i], g.GetY()[i] - g.GetErrorYlow(i))
        ret.SetPointError(i, 0, 0, 0, 0)
    
    if g.GetLineWidth() > 1 : ret.SetLineWidth(g.GetLineWidth() - 1)
    else : ret.SetLineWidth(1)
    ret.SetLineColor(g.GetLineColor())
    ret.SetLineStyle(2)
    return ret


# In[44]:


def drawCombMuHat(who, xmin, xmax, ymin, ymax, label="", fill=True):
    obs = ROOT.gFile.Get(who+"_comb_obs")
    if not obs : print "obs ",obs ; return
    if (obs is None) : return
    if (xmax == 600) : xmax = 600.1
    isZ = "mlz" in who
    nch = 6
    muhats = [ "vhbb", "htt",  "hzz4l", "hww", "hgg","comb" ]
    linec= [  209,     14,     215,    206,   217,  1     ]
    fillc[nch] = [   82,     17,      64,    208,    90,  1     ]
    fills[nch] = [ 1001,   1001,    1001,   1001,  1001,  3244  ]
    #*obsi[nchann], *obsi_h[nch], *obsi_l[nch]
    for i in xrange(nch): # for (int i = 0 i < nch; ++i):
        obsi[i] = ROOT.gFile.Get(who+"_"+muhats[i]+"_obs") 
        if (obsi[i] is None) : continue
        #//if (TString(muhats[i]) == "htt") continue
        obsi[i].SetLineColor(linec[i]) 
        if i == nch - 1 : 
            if fill : obsi[i].SetLineWidth(5);
            else : obsi[i].SetLineWidth(5+2);
        else :
            if fill : obsi[i].SetLineWidth(3);
            else : obsi[i].SetLineWidth(3+2);
        if (fill): 
            obsi[i].SetFillColor(fillc[i])
            obsi[i].SetFillStyle(fills[i])
        obsi_h[i] = getNSigmaLine(obsi[i], +1)
        obsi_l[i] = getNSigmaLine(obsi[i], -1)
        if (i == nch-1): 
            obsi_h[i].SetLineStyle(1) 
            obsi_l[i].SetLineStyle(1)
    
    pdefs.leg = newLegend(.62,.73,.95,.94)
    pdefs.leg.SetTextSize(0.032);
    if fill : pdefs.leg.AddEntry(obsi[nch-1], "Combined", "LF")
    else : pdefs.leg.AddEntry(obsi[nch-1], "Combined", "L")
    for i in xrange(nch-1) : #for (int i = 0 i < nch-1; ++i):
        if not (obsi[i] is None) :
            if fill : pdefs.leg.AddEntry(obsi[i], channelFromName(muhats[i],True), "LF")
            else : pdefs.leg.AddEntry(obsi[i], channelFromName(muhats[i],True), "L")
    
    frame0 = ROOT.TH1D("frame","frame", 1, xmin, xmax)
    frame0.Draw()
    ROOT.gStyle.SetOptStat(0);
    #for (int i = 0 i < nch; ++i): if (obsi[i]): if (fill) obsi[i].Draw("E3"); } }
    for i in xrange(nch) :
        if not (obsi[i] is None):
            if (fill) : obsi[i].Draw("E3")
    #for (int i = 0 i  < nch; ++i): if (obsi[i]): obsi_h[i].Draw("L"); obsi_l[i].Draw("L");}}
    for i in xrange(nch) :
        if not (obsi[i] is None):
            obsi_h[i].Draw("L")
            obsi_l[i].Draw("L")
    #for (int i = 0 i < nch; ++i): if (obsi[i]): obsi[i].Draw("LX"); } }
    for i in xrange(nch) :
        if not (obsi[i] is None):
            obsi[i].Draw("LX")
    setCanvas(frame0, "", ymin, ymax, "Best fit #sigma/#sigma_{SM}")
    frame0.Draw("AXIG SAME")
    spam(SPAM, 0.17,.89,.58,.94)
    pdefs.leg.Draw()
    finalize(who+"_all"+label,xmin,xmax,ymin,ymax)
    pdefs.leg = 0



# In[45]:


def cccPlot(who, where, mychann, nchann, postfix="", rMin=-2.0, rMax=5.0):
    obs   = ROOT.gFile.Get(who+"_"+mychann[0]+"_obs")
    if not obs : print "obs ", obs ; return
    if obs is None : return
    j = findBin(obs, where) 
    if (j == -1) : return;
    r0 = obs.GetY()[j]
    r0min = r0 - obs.GetErrorYlow(j)
    r0max = r0 + obs.GetErrorYhigh(j)
    halfint = (where - floor(where) > 0.3)
    pdefs.c1.SetLogx(0) 
    pdefs.c1.SetLogy(0) 
    pdefs.c1.SetTicky(0)
    pdefs.c1.SetTickx(0)
    #*points = new () int 
    npoints = 0;
    #names[pdefs.nchann_all]
    for i in range(nchann-1,0,-1): #for (int i = nchann-1 i > 0; --i):
       obsi   = ROOT.gFile.Get(who+"_"+mychann[i]+"_obs")
       if obsi is None : return
       ji = findBin(obsi, where) 
       if (ji != -1):
           npoints += 1
           points.Set(npoints)
           points.SetPoint(npoints-1, obsi.GetY()[ji], npoints - 0.5)
           points.SetPointError(npoints-1, obsi.GetErrorYlow(ji), obsi.GetErrorYhigh(ji), 0., 0.)
           names[npoints-1] = channelFromName(mychann[i],True,True)
       else:
           ji = findBin(obsi, floor(where))
           if (halfint and ji != -1 and findBin(obsi, ceil(where)) != -1):
               npoints += 1
               points.Set(npoints)
               points.SetPoint(npoints-1, 0.5*(obsi.GetY()[ji]+obsi.GetY()[ji+1]), npoints - 0.5)
               points.SetPointError(npoints-1, 0.5*(obsi.GetErrorYlow(ji) + obsi.GetErrorYlow(ji+1)), 
                                            0.5*(obsi.GetErrorYhigh(ji) + obsi.GetErrorYhigh(ji+1)), 
                                            0., 0.)
               names[npoints-1] = channelFromName(mychann[i],True,True)
       
    
    frame = ROOT.TH2F("frame","Best fit #sigma/#sigma_{SM};",1, rMin, rMax, npoints+2, 0., npoints+2);
    for k in xrange(npoints+1): # for (int k = 1 k <= npoints; ++k):
       frame.GetYaxis().SetBinLabel(k, names[k-1])
    
    frame.GetYaxis().SetTickLength(0)

    leftMargin = pdefs.c1.GetRightMargin() 
    pdefs.c1.SetLeftMargin(0.24+0.05*(where > 130));
    points.SetLineColor(kRed)
    points.SetLineWidth(5)
    points.SetMarkerStyle(21)
    points.SetMarkerSize(1.7)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.06)
    frame.Draw() 
    ROOT.gStyle.SetOptStat(0);
    ROOT.TBoxglobalFitBand(r0min, 0, r0max, npoints) 
    fake = ROOT.TH1F("fake","fake",1,0,1); 
    #//globalFitBand.SetFillStyle(3344) fake.SetFillStyle(3344);
    #//globalFitBand.SetFillColor(65)   fake.SetFillColor(65);
    globalFitBand.SetFillStyle(1001) 
    fake.SetFillStyle(1001);
    globalFitBand.SetFillColor(pdefs.colorFit68)  
    fake.SetFillColor(pdefs.colorFit68);
    globalFitBand.SetLineStyle(0)    
    globalFitBand.DrawClone()
    ROOT.TLineglobalFitLine(r0, 0, r0, npoints)
    globalFitLine.SetLineWidth(4)     
    fake.SetLineWidth(0);
    #//globalFitLine.SetLineColor(214) fake.SetLineColor(214); 
    globalFitLine.SetLineColor(1) 
    fake.SetLineColor(pdefs.colorFit68); 
    globalFitLine.DrawClone()
    globalFitLine.SetLineStyle(2)  
    points.Draw("P SAME")
    frame.Draw("AXIS SAME")
    xoff = 0 #//-0.1*SPAM.Contains("Preliminary");
    pdefs.leg = newLegend(.26+xoff,.74,.63+xoff,.935)
    pdefs.leg.SetTextSize(0.037); #// was TextSize 0.4, and goind up to 0.65+xoffs
    if (halfint):
        pdefs.leg.SetHeader(Form("   m_{H} = %.1f %s", where, pdefs.massUnits))
    else :
        pdefs.leg.SetHeader(Form("   m_{H} = %.0f %s", where, pdefs.massUnits))
    
    pdefs.leg.AddEntry(fake,  "Combined "+pdefs.oneSigmaFitCCC, "F") #//LEF
    #//pdefs.leg.AddEntry(points, "Single channel "+pdefs.oneSigmaFitCCC, "LP")
    pdefs.leg.AddEntry(points, "Single channel", "LP")
    if (SPAM.Contains("Preliminary")):
        spam("CMS Preliminary\n#sqrt{s} = 7 TeV\n"+channelFromName(mychann[0]), .65+xoff, .88, .94, .93, 22)
    else :
        spam("#splitline{"+pdefs.SPAM+"}{"+channelFromName(mychann[0])+"}", .65+xoff, .88, .94, .93, 22+10*("P" in pdefs.SPAM) )
    
    pdefs.leg.Draw()
    justSave(Form("%s_ccc_mH%.1f%s", who.Data(), where, postfix.Data()))
    pdefs.c1.SetLeftMargin(leftMargin) 
    pdefs.c1.SetTicky(1); pdefs.c1.SetTickx(1);



# In[46]:


def squareCanvas(gridx=1,gridy=1):
    ROOT.gStyle.SetCanvasDefW(600) #//Width of canvas
    ROOT.gStyle.SetPaperSize(20.,20.)
    pdefs.c1.Close()
    pdefs.c1 = ROOT.TCanvas("c1","c1")
    pdefs.c1.cd()
    pdefs.c1.SetWindowSize(600 + (600 - pdefs.c1.GetWw()), 600 + (600 - pdefs.c1.GetWh()))
    pdefs.c1.SetRightMargin(0.05)
    pdefs.c1.SetGridy(gridy)
    pdefs.c1.SetGridx(gridx);
    pdefs.isSquareCanvas = True


def rectangleCanvas(gridx=1,gridy=1):
    ROOT.gStyle.SetCanvasDefW(850) #//Width of canvas
    ROOT.gStyle.SetPaperSize(950./500.*20.,20.)
    pdefs.c1.Close()
    pdefs.c1 = ROOT.TCanvas("c1","c1")
    pdefs.c1.cd()
    pdefs.c1.SetWindowSize(950 + (950 - pdefs.c1.GetWw()), 500 + (500 - pdefs.c1.GetWh()))
    pdefs.c1.SetRightMargin(0.04)
    pdefs.c1.SetGridy(gridy) 
    pdefs.c1.SetGridx(gridx);
    pdefs.isSquareCanvas = False



# In[47]:


def theBand():
    pass # print 0


# In[48]:


def writeGrid(who, inputPrefix):
    g = makeAPrioriGrid("acls_%s" % str(who))
    if not g: print "g is ",g ; return
    if (g is None) : return
    g = slidingWindowAverage(g, 3)
    grid = open("grids/%s%s.txt" %  (inputPrefix, who ), "w")
    if not grid: 
        print "Cannot write to grid " , "grids/%s%s.txt" %  (inputPrefix, who ) # str(ROOT.TString.Format("grids/%s%s.txt", str(inputPrefix), str(who)))
        return
    n = g.GetN()
    for j in xrange(n) : #(int j = 0, n = g.GetN() j < n; ++j):
        x  = g.GetX()[j]
        y  = g.GetY()[j]
        yl = y-g.GetEYlow()[j]
        yh = y+g.GetEYhigh()[j]
        #print "DEBUG x yl y yh " , x , " " , yl , " " , y , " " , yh
        grid.write("%3.0f %5.2f  %5.2f  %5.2f\n" % (x, yl, y, yh))
        #fprintf(grid, "%3.0f %5.2f  %5.2f  %5.2f\n", x, yl, y, yh)
    
    grid.close()
    pdefs.CLs_grids.Add(g)

def writeGrids(inputPrefix):
    for i in xrange(pdefs.nchann_all) : # for (int i = 0 i < pdefs.nchann_all; ++i):
        writeGrid(pdefs.chann[i], inputPrefix)
        # commented out by Bockjoo
        # break
        
def readGrids(inputPrefix):
    for i in xrange(nchann) : # for (int i = 0 i < nchann; ++i):
        #FILE *fgrid = fopen(TString::Format("grids/%s%s.txt", inputPrefix.Data(), pdefs.chann[i]).Data(), "r")
        fgrid = open(ROOT.TString.Format("grids/%s%s.txt", str(inputPrefix), pdefs.chann[i]), "r")
        if not (fgrid is None) : continue
        #float x,y,yl,yh
        #*grid = new () int 
        points = 0;
        while(fscanf(fgrid,"%f %f %f %f", x, yl, y, yh) == 4):
            grid.Set(points+1)
            grid.SetPoint(points,x,y)
            grid.SetPointError(points,0,0,y-yl,yh-y)
            points +=1
        
        print "Input grid for " , pdefs.chann[i] , " contains " , points , " points"
        fgrid.close()
        grid.SetName(ROOT.TString.Format("grid_cls_%s", pdefs.chann[i]))
        grid.SetFillColor(16)
        grid.SetFillStyle(3244)
        pdefs.CLs_grids.Add(grid)
    



# In[49]:


def makeGrid(who, who2, lowfactor=0.2, highfactor=4):
    g1 = ROOT.gFile.Get(who)
    #if g1 is None: m1 = None
    if not g1 : m1 = 0    
    else: m1 = missingPoints(g1) 
    g2 = ROOT.gFile.Get(who2)
    #m2 = g2 ? missingPoints(g2) : 0
    if not g2 : m2 = 0 # g2 is None: m2 = None
    else: m2 = missingPoints(g2) 
    if not g1 and not g2: 
        print "Missing both '", who, "' and '", who2, "'"
        return False
    #TGraphAsymmErrors *ret = new TGraphAsymmErrors(nmasses); 
    masses = loadMasses()
    nmasses = len(masses)
    #print "DEBUG makeGrid nmasses == ", nmasses
    ret = ROOT.TGraphAsymmErrors(nmasses) 
    for i in xrange(nmasses): # for (int i = 0 i < nmasses; ++i):
        x = float(masses[i])
        #print "DEBUG makeGrid findBin g1 x ",x
        i1 = findBin(g1, x)
        #print "DEBUG makeGrid findBin m1 x ",x
        j1 = findBin(m1, x)
        #print "DEBUG makeGrid findBin g2 x ",x
        i2 = findBin(g2, x)
        #print "DEBUG makeGrid findBin m2 x ",x
        j2 = findBin(m2, x)
        #print "DEBUG makeGrid i1 j1 i2 j2 " , i1 , " " , j1 , " " , i2 , " " , j2

        ymin = 999.0 ; ymax = -999.0
        if (i1 != -1 and g1.GetY()[i1] != 0): 
            ymin = ROOT.TMath.Min(ymin, g1.GetY()[i1]) 
            ymax = ROOT.TMath.Max(ymax, g1.GetY()[i1])
        if (i2 != -1 and g2.GetY()[i2] != 0): 
            ymin = ROOT.TMath.Min(ymin, g2.GetY()[i2]) 
            ymax = ROOT.TMath.Max(ymax, g2.GetY()[i2])
        if (j1 != -1 and m1.GetY()[j1] != 0): 
            ymin = ROOT.TMath.Min(ymin, m1.GetY()[j1]) 
            ymax = ROOT.TMath.Max(ymax, m1.GetY()[j1])
        if (j2 != -1 and m2.GetY()[j2] != 0): 
            ymin = ROOT.TMath.Min(ymin, m2.GetY()[j2]) 
            ymax = ROOT.TMath.Max(ymax, m2.GetY()[j2])
        if (pdefs.SM != "SM4"):
            if (ymin < 0.2) : ymin = 0.2 
            if (ymin > 3) : ymin = 3;
            if ("comb" in who):
                if (ymax < 0.2) : ymax = 0.2 
                if (ymax > 8) :  ymax = 8
            else:
                if (ymax < 0.5) : ymax = 0.5 
            
        else:
            if (ymin < 0.02) : ymin = 0.02 
            if (ymin > 3) : ymin = 3
            if ("comb" in who):
                if (ymax < 0.02) : ymax = 0.02 
                if (ymax > 8) : ymax = 8
            else:
                if (ymax < 0.05) : ymax = 0.05 
            
        
        y = math.sqrt(ymin*ymax)
        ret.SetPoint(i, float(x), float(y))
        ret.SetPointError(i, 0, 0, y - lowfactor*ymin, highfactor*ymax - y)
    
    #//ret.SetLineStyle(1)
    #//ret.SetMarkerStyle(0)
    ret.SetFillColor(16)
    ret.SetFillStyle(3244)
    ret.SetLineColor(223)
    ret.SetLineWidth(4)
    ret.SetLineStyle(3)
    return ret


