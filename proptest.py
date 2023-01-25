import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel
import sys
import argparse
import numpy as np

#python3 propdraw.py ias_study_2018B_template_2018B_15mars_eta_1_pubins_5.root 

if __name__ == "__main__":


     c = TCanvas("c1")
     c.SetBottomMargin(0.13)
     c.SetTopMargin(0.2)
     c.SetGridx()
     c.SetGridy()


     ias_intervals = 5
     ias_bins = [16.26,23.06,27.97,32.88,43.15]
     means_ias = [0.21,0.25,0.27,0.31,0.34]
     err_mean_x_base = [0.0,0.0,0.0,0.0,0.0]
     error_mean_base_cutih = [0.0021,0.0025,0.0027,0.0031,0.0034]
     print(len(ias_bins))
     p = TGraphErrors(len(ias_bins),np.array(ias_bins),np.array(means_ias),np.array(err_mean_x_base),np.array(error_mean_base_cutih))
     binx = p.GetXaxis().GetNbins()
     print(binx)
     ias_axis_names = [''] * 5
     for m in range(ias_intervals):
         print("looking for ias bin ",m)
         ias_axis_names[m] = "PU_" + str(ias_bins[m])
     

     for x in range(ias_intervals):
         p.GetXaxis().SetBinLabel(p.GetXaxis().FindBin(ias_bins[x]),ias_axis_names[x])




     p.Draw("AP")
     labellumi = TPaveLabel(0.6,0.75,0.5,0.9,"SingleMu 2018B, 13.6 TeV","brNDC")
     labellumi.SetFillColor(0)
     labellumi.SetTextColor(1)
     labellumi.SetFillStyle(0)
     labellumi.SetBorderSize(0)
     labellumi.SetTextSize(0.25)
     labellumi.SetTextAlign(12)
     labellumi.Draw("same")

      #if(saveplot):
     c.SaveAs("test1.pdf")
