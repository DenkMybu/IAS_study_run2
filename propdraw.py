import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3
import sys
import argparse

def poissonHisto(h,p,RNG):
    for i in range(p.GetNbinsX()+1):
        for j in range(p.GetNbinsY()+1):
            for k in range(p.GetNbinsZ()+1):
                h.SetBinContent(i,j,k,RNG.Poisson(p.GetBinContent(i,j,k)))


#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

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

     Tpt_PU_0_20 = file0.Get("Charge_Vs_Path_PU_between_0_20")
 
     Tpt_PU_20_25 = file0.Get("Charge_Vs_Path_PU_between_20_25")
     Tpt_PU_25_30 = file0.Get("Charge_Vs_Path_PU_between_25_30")
     Tpt_PU_30_35 = file0.Get("Charge_Vs_Path_PU_between_30_35")
     Tpt_PU_35_99 = file0.Get("Charge_Vs_Path_PU_between_35_99")



     nPE = 100

     hfile = TFile( 'PE_templates_2018B.root', 'RECREATE' )

     PE_templates = [0] * nPE
     PE_templates_0_20 = [0] * nPE
    
     PE_templates_20_25 = [0] * nPE
     PE_templates_25_30 = [0] * nPE
     PE_templates_30_35 = [0] * nPE
     PE_templates_35_99 = [0] * nPE


     RNG = TRandom3()
     
     for i in range(nPE):
         name = "PE_" + str(i)
         PE_templates[i] = TH3D (name,name,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates[i],AllIncTemplate,RNG)
         print("Done generating pseudo experiment #",i)         

         name_0_20 = "PE_" + str(i) + "_PU_0_20"
         PE_templates_0_20[i]= TH3D (name_0_20,name_0_20,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates_0_20[i],Tpt_PU_0_20,RNG)
         print("Done generating pseudo experiment #",i, " for tpt 0_20")         
          
         name_20_25 = "PE_" + str(i) + "_PU_20_25"
         PE_templates_20_25[i]= TH3D (name_20_25,name_20_25,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates_20_25[i],Tpt_PU_20_25,RNG)
         print("Done generating pseudo experiment #",i, " for tpt 20_25")         


         name_25_30 = "PE_" + str(i) + "_PU_25_30"
         PE_templates_25_30[i]= TH3D (name_25_30,name_25_30,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates_25_30[i],Tpt_PU_25_30,RNG)
         print("Done generating pseudo experiment #",i, " for tpt 25_30")         

         name_30_35 = "PE_" + str(i) + "_PU_30_35"
         PE_templates_30_35[i]= TH3D (name_30_35,name_30_35,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates_30_35[i],Tpt_PU_30_35,RNG)
         print("Done generating pseudo experiment #",i, " for tpt 30_35")         

         name_35_99 = "PE_" + str(i) + "_PU_35_99"
         PE_templates_35_99[i]= TH3D (name_35_99,name_35_99,P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max)
         poissonHisto(PE_templates_35_99[i],Tpt_PU_35_99,RNG)
         print("Done generating pseudo experiment #",i, " for tpt 35_99")         
         

     hfile.Write()
