# Code written by:  zaixing.mao@cern.ch && edward.laird@cern.ch from Brown U.                                                                  
#!/usr/bin/env python                                                                                                                          
import ROOT as r
import numpy
from sys import argv, exit, stdout, stderr
import math

if len(argv) < 2:
   print 'Usage:python xs_calculator_prefit.py DYCrossSection[optional]'

if len(argv)>1:
   FoundXS= numpy.array([argv[1]],dtype=float)
else:
   FoundXS=1.00

histDict = {}

#variable = ['InvariantMass_of_muon_pair_with_opposite_sign_in_one_tau_jet_region','InvariantMass_of_muon_pair_with_opposite_sign_in_two_tau_jet_region','InvariantMass_of_muon_pair_with_opposite_sign_in_zero_tau_jet_region','InvariantMass_of_muon_pair_with_same_sign_in_one_tau_jet_region','InvariantMass_of_muon_pair_with_same_sign_in_two_tau_jet_region','InvariantMass_of_muon_pair_with_same_sign_in_zero_tau_jet_region','InvariantMass_of_tau_pair_with_opposite_sign_region','InvariantMass_of_tau_pair_with_same_sign_region','InvariantMass_of_taumu_with_same_sign','InvariantMass_of_taumu_with_opposite_sign','pt_distribution_of_Leading_#mu_with_opposite_sign_in_zero_tau_jet_region','pt_distribution_of_Leading_#mu_with_same_sign_in_zero_tau_jet_region','pt_distribution_of_SubLeading_#mu_with_opposite_sign_in_zero_tau_jet_region','pt_distribution_of_SubLeading_#mu_with_same_sign_in_zero_tau_jet_region','pt_distribution_of_Leading_#tau_with_opposite_sign_region',                         
variable=['pt_distribution_of_Leading_#tau_with_same_sign_region','pt_distribution_of_SubLeading_#tau_with_opposite_sign_region','pt_distribut
ion_of_SubLeading_#tau_with_same_sign_region','InvariantMass_of_4_particle_with_totCharge==0_and_mucharge==0','InvariantMass_of_4_particle_wit
h_totCharge==0','InvariantMass_of_4_particle_with_totCharge==0_and_mucharge_not_equal_to_zero_i.e_signal_region','et_distribution_of_Leading_t
au_with_opposite_sign_region','et_distribution_of_Leading_tau_with_same_sign_region','eta_distribution_of_SubLeading_#tau_with_opposite_sign_r
egion','eta_distribution_of_SubLeading_#tau_with_same_sign_region','eta_distribution_of_Leading_#tau_with_opposite_sign_region','eta_distribut
ion_of_Leading_#tau_with_same_sign_region','eta_distribution_of_SubLeading_#mu_with_same_sign_in_zero_tau_jet_region','eta_distribution_of_Lea
ding_#mu_with_same_sign_in_zero_tau_jet_region','eta_distribution_of_SubLeading_#mu_with_opposite_sign_in_zero_tau_jet_region','eta_distributi
on_of_Leading_#mu_with_opposite_sign_in_zero_tau_jet_region']
f = r.TFile("prefit.root","recreate")
################################################                                                                                               
# Sevreal Histograms are initiated/produced here                                                                                               
defaultOrder = [('WJets',  r.TColor.GetColor(100,182,232)),
                ('TTJets', r.TColor.GetColor(155,152,204)),
                ('DY', r.TColor.GetColor(250,202,255)),
                ('WW', r.TColor.GetColor(248,206,104)),
                ('ZZ', r.TColor.GetColor(200,225,60)),
                ('HH', r.TColor.GetColor(210,70,80)),
                ('ZZTo4L', r.TColor.GetColor(30,50,120)),
                ('ST', r.TColor.GetColor(230,90,115)),
                ('WZ', r.TColor.GetColor(240,200,108))]



def buildHistDict(nbins,x_low,x_high):

    histDict = {}
    for iSample, iColor in defaultOrder:
       name = iSample
       histDict[name] = r.TH1F(name, '', nbins, x_low,x_high)
       histDict[name].SetFillColor(iColor)
       histDict[name].SetMarkerColor(iColor)
       histDict[name].SetMarkerStyle(21)
       histDict[name].SetLineColor(r.kBlack)

    name = 'bkg_'
    histDict[name] = r.TH1F(name, '', nbins, x_low,x_high)
    histDict[name].Sumw2()
    histDict[name].SetFillColor(r.kGray+2)
    histDict[name].SetLineColor(r.kGray+2)
    histDict[name].SetFillStyle(3344)
    name = 'data_'
    histDict[name] = r.TH1F(name, '', nbins, float(x_low),float(x_high))
    histDict[name].Sumw2()
    histDict[name].SetMarkerStyle(8)
    histDict[name].SetMarkerSize(0.9)
    histDict[name].SetMarkerColor(r.kBlack)
    return histDict
################################################                                                                                               

def setMyLegend(lPosition, lHistList):
    l = r.TLegend(lPosition[0], lPosition[1], lPosition[2], lPosition[3])
    l.SetFillStyle(0)
    l.SetBorderSize(0)
    for i in range(len(lHistList)):
                    l.AddEntry(lHistList[i][0], lHistList[i][1], lHistList[i][2])
    return l

def getBins(hist, x_low, x_high):
    bin_low = 0
    bin_high = 0

    for i in range(hist.GetNbinsX()):
        if hist.GetBinCenter(i+1) >= x_low and bin_low == -1:
            bin_low = i+1
        if hist.GetBinCenter(i+1) >= x_high and bin_high == -1:
            bin_high = i
        if bin_low != 0 and bin_high != 0:
            return bin_low, bin_high
    return -1,-1

def buildStackDict(histDict):
    stackDict = {}
    stackDict['s'] = r.THStack()

    for iSample, iColor in defaultOrder:
        scale = 1.0
        stackDict['s'].Add(histDict[iSample])
        histDict['bkg_'].Add(histDict[iSample])
    return stackDict

def FillHisto(input, output, weight = 1.0):
#    print 'inFillHisto',input,'ou== ', output                                                                                                 
    for i in range(input.GetNbinsX()):
       currentValue = output.GetBinContent(i+1)
       currentError = output.GetBinError(i+1)
       output.SetBinContent(i+1, currentValue+input.GetBinContent(i+1)*weight)
       output.SetBinError(i+1, math.sqrt((input.GetBinError(i+1))**2 + currentError**2))

def buildLegendDict(histDict, position):
    legendDict = {}
    histList = {'T': []}
#    histList['T'].append((histDict['data_'+sign], 'Observed', 'lep'))                                                                         
    for iSample, iColor in reversed(defaultOrder):
       histList['T'].append((histDict[iSample], iSample, 'f'))

    legendDict['T'] = setMyLegend(position, histList['T'])
    return legendDict


def xs_calculator(fileList = []):

   for i in range(len(variable)) :
      for iFileName, iFileLocation in fileList:
         ifile_ = r.TFile(iFileLocation)
         if ifile_.Get(variable[i]) :
            hist = ifile_.Get(variable[i])
            nbins = hist.GetNbinsX()
            x_low = hist.GetBinLowEdge(1)
            x_high = hist.GetBinWidth(hist.GetNbinsX())+hist.GetBinLowEdge(hist.GetNbinsX())
            histDict = buildHistDict(int(nbins),x_low,x_high)
            break
   #loop over all the samples                                                                                                                 
      for iFileName, iFileLocation in fileList:
         ifile = r.TFile(iFileLocation)
         weight = 1.
         tauWeight = 1.
         if ifile.Get(variable[i]) :
            FillHisto(ifile.Get(variable[i]), histDict[iFileName], tauWeight)
      stackDict = buildStackDict(histDict)
      legendDict = buildLegendDict(histDict, (0.6, 0.8 - 0.06*4, 0.85, 0.8))

      pdf = ''
      c = r.TCanvas("c","Test", 800, 600)
      pad1 = r.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
      pad1.SetBottomMargin(0)
      pad1.SetGridx()
      pad1.Draw()
      pad1.cd()
      max_t = 1.2*max(stackDict['s'].GetMaximum(), histDict['data_'].GetMaximum())
      stackDict['s'].Draw('hist H')
      stackDict['s'].SetTitle('%s;%s;events' %(variable[i],variable[i]))
      stackDict['s'].SetMaximum(max_t)
      stackDict['s'].GetYaxis().SetTitleOffset(1.2)
      histDict['data_'].Draw('same PE')
      histDict['bkg_'].Draw('E2 same')
      legendDict['T'].Draw('same')
      c.cd()
      pad2 = r.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
      pad2.SetTopMargin(0)
      pad2.SetBottomMargin(0.2)
      pad2.SetGridx()
      pad2.Draw()
      pad2.cd()
      h3 = histDict['data_'].Clone("h3")
      h3.SetLineColor(r.kBlack)
    h3.SetMinimum(0.)
      h3.SetMaximum(2.)
      h3.Sumw2()
      h3.SetStats(0)
      h3.Divide(histDict['bkg_'])
      h3.SetMarkerStyle(21)
      h3.Draw("ep")
      h3.GetXaxis().SetTitle("")
      h3.GetYaxis().SetTitle("ratio data/Mc ");
      h3.GetYaxis().SetNdivisions(505);
      h3.GetYaxis().SetTitleSize(20);
      h3.GetYaxis().SetTitleFont(43);
      h3.GetYaxis().SetTitleOffset(1.55);
      h3.GetYaxis().SetLabelFont(43)
      h3.GetYaxis().SetLabelSize(15);

      h3.GetXaxis().SetTitleSize(20);
      h3.GetXaxis().SetTitleFont(43);
      h3.GetXaxis().SetTitleOffset(4.);
      h3.GetXaxis().SetLabelFont(43);
      h3.GetXaxis().SetLabelSize(.08)
      pdf += variable[i]+'.gif'
      c.SaveAs('%s' %pdf)
      f.Write()
      f.Close()



dirName = '.'

fileList = [('data_', '%s/DoubleMuon_controlPlot.root' %dirName),
            ('HH', '%s/HH_controlPlot.root' %dirName),
            ('ZZTo4L', '%s/zzTo4L_controlPlot.root' %dirName),
            ('WW', '%s/WW_controlPlot.root' %dirName),
          ('TTJets', '%s/TT_controlPlot.root' %dirName),
#            ('WJets', '%s/WJets_controlPlot.root' %dirName),                                                                                  
            ('DY', '%s/DY_controlPlot.root' %dirName),
            ('WZ', '%s/WZ_controlPlot.root' %dirName),
            ('ZZ', '%s/ZZ_controlPlot.root' %dirName),
            ('ST', '%s/ST_controlPlot.root' %dirName)
           ]

xs_calculator(fileList = fileList)



