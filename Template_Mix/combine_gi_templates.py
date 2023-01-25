import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math


def fillHisto(hist,hist_merged,l=1):
    for i in range(p.GetNbinsX()+1):
        for j in range(p.GetNbinsY()+1):
            for k in range(p.GetNbinsZ()+1):
                hist_merged.SetBinContent
                h.SetBinContent(i,j,k,RNG.Poisson(p.GetBinContent(i,j,k)))

#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     gi1 = "Calibration_GiTemplate_PU_1"
     gi2 = "Calibration_GiTemplate_PU_2"
     gi3 = "Calibration_GiTemplate_PU_3"
     gi4 = "Calibration_GiTemplate_PU_4"
     gi5 = "Calibration_GiTemplate_PU_5"


     gi_merged1 = "Calibration_GiTemplate_PU_1_merged" 
     gi_merged2 = "Calibration_GiTemplate_PU_2_merged" 
     gi_merged3 = "Calibration_GiTemplate_PU_3_merged" 
     gi_merged4 = "Calibration_GiTemplate_PU_4_merged" 
     gi_merged5 = "Calibration_GiTemplate_PU_5_merged"

     for i in range(999):
         filename = "Histos_" + str(i) +".root"
         file0 = ROOT.TFile.Open(filename)
         hist = file0.Get()

         
     
