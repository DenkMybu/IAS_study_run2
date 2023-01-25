import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math



#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple
def poissonHisto(p):
    total_count = 0
    for i in range(p.GetNbinsX()+1):
        for j in range(p.GetNbinsY()+1):
            for k in range(p.GetNbinsZ()+1):
                total_count += p.GetBinContent(i,j,k)
    return total_count

if __name__ == "__main__":

     gi1 = "Charge_Vs_Path_PU_between_0_20"
     gi2 = "Charge_Vs_Path_PU_between_20_25"
     gi3 = "Charge_Vs_Path_PU_between_25_30"
     gi4 = "Charge_Vs_Path_PU_between_30_35"
     gi5 = "Charge_Vs_Path_PU_between_35_110"


     all_entries_PU =[0]*5
     filename = "/opt/sbg/cms/ui14_data1/rhaeberl/CMSSW_10_6_27/src/template_2018A_15mars_selection_5_bin_eta_2p1_pmin_20_pmax_48_final_HSCP_isotk15_minireliso002_validation.root"
     file0 = ROOT.TFile.Open(filename)
     if not file0:   
         print("Could not open file ", filename)
         #print("Current directory: '{}'\n".format(ROOT.gDirectory.GetName()))

     hist1 = file0.Get(gi1)
     hist2 = file0.Get(gi2)
     hist3 = file0.Get(gi3)
     hist4 = file0.Get(gi4)
     hist5 = file0.Get(gi5)

     all_entries_PU[0] += poissonHisto(hist1)
     all_entries_PU[1] += poissonHisto(hist2)
     all_entries_PU[2] += poissonHisto(hist3)
     all_entries_PU[3] += poissonHisto(hist4)
     all_entries_PU[4] += poissonHisto(hist5)

     for u in range(5):
         print("In total, PU bin ",u, " has nb of entries = ", all_entries_PU[u] , " on all 2018A era")
