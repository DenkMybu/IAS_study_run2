import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math




#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     #file_trf = sys.argv[1]
     DeDxMap_ = TH3D
     if not DeDxMap_:
         print("not opened") 
     '''
     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", file_trf)
     else:
         print("Successfully opened root : ",file_trf)

     #hist_num_raw = file0.Get("P_single_pt50_ias_qtl_80_90") 
     #hist_denom_raw = file0.Get("P_triple_pt50_ias_qtl_80_90")
     '''


     '''
     str_denom = "IAS_TRIPLE_PT_50_60_selection"

     
     hist_num_raw = file0.Get(str_denom) 

     hist_num = hist_num_raw.Clone("Clone_num")
     hist_num.Sumw2()

     print("BEFORE SCALING")
     print("Integral : ", hist_num.Integral())

     hist_num.Scale(1./hist_num.Integral())
     print("AFTER RENORMALIZATION")

     print("Integral : ", hist_num.Integral())
     '''
