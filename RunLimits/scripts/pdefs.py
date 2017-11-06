import ROOT
# globals
c1 = ROOT.TCanvas("c1","c1")
leg = 0
globalPrefix = ""

SPAM = "CMS Preliminary,  #sqrt{s} = 7 TeV"
SPAM2L = "CMS Preliminary\n#sqrt{s} = 7 TeV"
#//TString SPAM = "CMS Priv,  #sqrt{s} = 7 TeV";
#//TString SPAM = "PRIVATE,  #sqrt{s} = 7 TeV";
SM = "SM"
oneSigmaText="(68%)" #//"#pm 1#sigma";
twoSigmaText="(95%)" #//"#pm 2#sigma";
oneSigmaFitText="68% CL band"
oneSigmaFitCCC=" (68%)"
massUnits="GeV"
lumiSymbol="L" #// L_{int}
justLumiForCombined = True
doSquares = True
noLineStyles = False
channelSpamOnRightHandSide = False
doLEESPAM = False
LEESPAM_1L = "NOT USED"
LEESPAM_2L = "NOT USED"
LEESPAM_3L = "Global significance\n0.8#sigma for 110-600 GeV range\n2.1#sigma for 110-145 GeV range"

isSquareCanvas = False
isTiny = False
lessSpacesInLegends = False
track_missing = False

x_zoom = 145.
doZoom2 = False
x_zoom2_min = 160
x_zoom2_max = 300
forceYmin = 0 #//0.08;
forceYmax = 0 #//12.8;

nchann_all = 1+5+3+3+7
nchann = 14
nchann2 = 6
nchann3 = 3
nchann4 = 13
chann = [ "comb", "vhbb", "htt", "hgg", "hww", "hzz", "combp", "combs", "combl","httm", "vhtt", "vhww3l", "htt0", "hww2l", "hzz4l", "hzz2l2t", "hzz2l2q_low", "hzz2l2q", "hzz2l2nu"]
chann2 = [ "comb", "vhbb", "htt", "hgg", "hww", "hzz" ]
chann3 = [ "comb", "combp_full", "combl_full" ]
chann4 = [ "comb", "vhbb", "htt0", "httm", "vhtt", "hgg", "hww2l", "vhww3l", "hzz4l", "hzz2l2t", "hzz2l2q_low", "hzz2l2q", "hzz2l2nu" ]
color68 = 80
color95 = 90
color50 = ROOT.TColor.GetColor(20,20,20)
colorFit68 = 80
lineAt1Style = 1
lineAt1Color = 2
cms_excluded = False
lep_excluded = False
tev_excluded = False
tev_excluded_alsolow = False
fakeLEP = 0
fakeTEV = 0
fakeCMS = 0
CLs_debug_apriori_grid=True
CLs_grids = ROOT.TList()
Draw_TEV=False
PLC_debug = False

do_bands_nosyst = True
do_bands_mean = True
do_bands_median = True
do_bands_ntoys = True
do_bands_asimov = True
do_bands_cputime = False
do_bands_realtime = False
do_bands_95 = True


pdefs = "This is the global plot variables"
band_safety_crop = 0
use_precomputed_quantiles = False
precomputed_median_only = False
zero_is_valid = False
seed_is_channel = False
halfint_masses  = False #// find the halfling!
#class Enum(tuple): __getattr__ = tuple.index
#ObsAvgMode = Enum(['MeanObs', 'LogMeanObs', 'MedianObs'])    
#obs_avg_mode = ObsAvgMode.MeanObs
#BandType = Enum(['Mean', 'Median', 'Quantile', 'Observed', 'Asimov', 'CountToys', 'MeanCPUTime', 'MeanRealTime', 'AdHoc', 'ObsQuantile']) 
