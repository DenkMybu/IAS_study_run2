import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math



#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     gi1 = "analyzer/BaseName/Calibration_GiTemplate_PU_1"
     gi2 = "analyzer/BaseName/Calibration_GiTemplate_PU_2"
     gi3 = "analyzer/BaseName/Calibration_GiTemplate_PU_3"
     gi4 = "analyzer/BaseName/Calibration_GiTemplate_PU_4"
     gi5 = "analyzer/BaseName/Calibration_GiTemplate_PU_5"

 
     for i in range(1,30):
         #filename = "/opt/sbg/cms/ui2_data1/rhaeberl/Prod/ProdDecember2022/SingleMuon/UL2018_RunA/Histos_" + str(i) +".root"
         filename = "/opt/sbg/cms/ui3_data1/rhaeberl/HSCP_prod/CompareFramewoks/framework_run316722/Histos_" + str(i) + ".root"
         file0 = ROOT.TFile.Open(filename)
         if not file0:
             print("Could not open file ", filename)
             continue

         print("Current directory: '{}'\n".format(ROOT.gDirectory.GetName()))


        
         path = "analyzer/BaseName" 
         #file0.cd(path)
         print("Current directory: '{}'\n".format(ROOT.gDirectory.GetName()))

         hist1 = file0.Get(gi1)
         hist2 = file0.Get(gi2)
         hist3 = file0.Get(gi3)
         hist4 = file0.Get(gi4)
         hist5 = file0.Get(gi5)

         GiHist1 = hist1.Clone("Clone_Calibration_GiTemplate_PU_1")
         GiHist2 = hist2.Clone("Clone_Calibration_GiTemplate_PU_2")
         GiHist3 = hist3.Clone("Clone_Calibration_GiTemplate_PU_3")
         GiHist4 = hist4.Clone("Clone_Calibration_GiTemplate_PU_4")
         GiHist5 = hist5.Clone("Clone_Calibration_GiTemplate_PU_5")

         outname = "Run316722/GiTemplates_Histos_PostS" + str(i) +".root"
         hfile = TFile( outname, 'RECREATE' ) 
        
         hist1.Write()
         hist2.Write()
         hist3.Write()
         hist4.Write()
         hist5.Write()
 

         hfile.Write() 
     
