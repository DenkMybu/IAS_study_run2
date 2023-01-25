import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math



#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     gi1 = "Charge_Vs_Path_PU_between_0_20"
     gi2 = "Charge_Vs_Path_PU_between_20_25"
     gi3 = "Charge_Vs_Path_PU_between_25_30"
     gi4 = "Charge_Vs_Path_PU_between_30_35"
     gi5 = "Charge_Vs_Path_PU_between_35_110"
     filename = "/opt/sbg/cms/ui14_data1/rhaeberl/CMSSW_10_6_27/src/template_2018A_15mars_selection_5_bin_eta_2p1_pmin_20_pmax_48_final_HSCP_all_preselection.root"

     file0 = ROOT.TFile.Open(filename)
     if not file0:
         print("Could not open file ", filename)

 
     hist1 = file0.Get(gi1)
     hist2 = file0.Get(gi2)
     hist3 = file0.Get(gi3)
     hist4 = file0.Get(gi4)
     hist5 = file0.Get(gi5)

         #filename = "/opt/sbg/cms/ui2_data1/rhaeberl/Prod/ProdDecember2022/SingleMuon/UL2018_RunA/Histos_" + str(i) +".root"

     GiHist1 = hist1.Clone("GiTemplate_PU_1")
     GiHist2 = hist2.Clone("GiTemplate_PU_2")
     GiHist3 = hist3.Clone("GiTemplate_PU_3")
     GiHist4 = hist4.Clone("GiTemplate_PU_4")
     GiHist5 = hist5.Clone("GiTemplate_PU_5")

     outname = "Run316722/GiTemplates.root"
     hfile = TFile( outname, 'RECREATE' ) 
        
     GiHist1.Write()
     GiHist2.Write()
     GiHist3.Write()
     GiHist4.Write()
     GiHist5.Write()

     hfile.Write() 
     
