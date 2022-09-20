import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3
import sys
import argparse

def poissonHisto(h,p,RNG):
    for i in range(p.GetNbinsX()+1):
        for j in range(p.GetNbinsY()+1):
            for k in range(p.GetNbinsZ()+1):
                h.SetBinContent(i,j,k,RNG.Poisson(p.GetBinContent(i,j,k)))

#python3 propdraw.py template_2018B_15mars_selection_5_bin_eta_2p1_pmin_10.root


P_Min = 1
P_Max = 16 # 1 + 14 + 1; final one is for pixel!
P_NBins = 15 #15th bin = pixel; 0 is underflow
Path_Min = 0.2
Path_Max = 1.6
Path_NBins = 42
Charge_Min = 0
Charge_Max = 5000
Charge_NBins = 500 


if __name__ == "__main__":

     file_trf = sys.argv[1]
     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", sys.argv[1])
     else:
         print("Successfully opened root file containing templates")

     AllIncTemplate = file0.Get("Charge_Vs_Path")
     if not AllIncTemplate:
         print("Could not open all inclusive template from root file")
     else:
         print("Successfully opened template from root file")

     nPE = 100

     hfile = TFile( 'PE_templates.root', 'RECREATE' )

     PE_templates = [0] * nPE

     RNG = TRandom3()
     
     for i in range(nPE):
         name = "PE_" + str(i)
         PE_templates[i] = TH3D (name,name,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates[i],AllIncTemplate,RNG)
         print("Done generating pseudo experiment #",i)         
     
     

     hfile.Write()
