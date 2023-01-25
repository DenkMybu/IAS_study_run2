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

     gi1 = "analyzer/BaseName/Calibration_GiTemplate_PU_1"
     gi2 = "analyzer/BaseName/Calibration_GiTemplate_PU_2"
     gi3 = "analyzer/BaseName/Calibration_GiTemplate_PU_3"
     gi4 = "analyzer/BaseName/Calibration_GiTemplate_PU_4"
     gi5 = "analyzer/BaseName/Calibration_GiTemplate_PU_5"

     all_entries_PU =[0]*5

     for i in range(3,999):
         filename = "/opt/sbg/cms/ui3_data1/rhaeberl/HSCP_prod/prodJan2023_CMSSW_10_6_30/HSCP_framework/UL2018_RunA/0000/Histos_" + str(i) +".root"

         file0 = ROOT.TFile.Open(filename)
         if not file0:
             print("Could not open file ", filename)
             continue

         #print("Current directory: '{}'\n".format(ROOT.gDirectory.GetName()))


        
         path = "analyzer/BaseName" 
         #file0.cd(path)
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
