import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math




#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     file_trf = sys.argv[1]

     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", file_trf)
     else:
         print("Successfully opened root : ",file_trf)


     '''
     split_name_num = file_trf.split('.')[0].split('_')
     split_name_denom = file_trf2.split('.')[0].split('_')

     if(len(split_name_num) == 22):
         hname = "Stddev_multiple_" + split_name_num[20] + "_" + split_name_num[21]
     '''


     hfile = TFile("test_stddev_tenth_stat2_tenth_stat1.root", 'RECREATE' )
     str_num = "IAS_1_VS_Delta_Ias"
 
     hist_num_raw = file0.Get(str_num) 




     name_xaxis = hist_num_raw.GetXaxis().GetTitle()
     print("name of x axis = ", name_xaxis)
     nb_xaxis = hist_num_raw.GetNbinsX()
     print("nb of bins X : ", nb_xaxis)
     nb_yaxis = hist_num_raw.GetNbinsY()

     std_devs = TH1D("std_devs","std_devs",nb_xaxis,0,1)
     std_devs.Sumw2()
     std_devs.SetTitle("")
     std_devs.GetYaxis().SetTitle("std_dev_y")
     std_devs.GetYaxis().SetTitleSize(15)
     std_devs.GetYaxis().SetTitleFont(43)
     std_devs.GetYaxis().SetTitleOffset(2)
     std_devs.GetYaxis().SetLabelFont(43)
     std_devs.GetYaxis().SetLabelSize(15)

     '''
     upper_lim_ratio = 0
     lower_lim_ratio = 0.005
     std_devs.SetMinimum(lower_lim_ratio)
     std_devs.SetMaximum(upper_lim_ratio)
     '''


     std_devs.GetYaxis().SetNdivisions(505)

     std_devs.GetXaxis().SetTitle("IAS")
     std_devs.GetXaxis().SetTitleSize(20)

     #std_devs.SetStats(0)


     std_devs.SetMarkerStyle(1)
     std_devs.SetMarkerSize(1)
     std_devs.SetLineColor(1)

     
     for i in range(nb_xaxis+1):
         hist_num_raw.GetXaxis().SetRange(i,i)
         print("bin #",i,"center at ", hist_num_raw.GetXaxis().GetBinCenter(i),  " std dev : ", hist_num_raw.GetStdDev(2))
         if(hist_num_raw.GetStdDev(2) != 0):
             std_devs.SetBinContent(i,hist_num_raw.GetStdDev(2))
             std_devs.SetBinError(i,hist_num_raw.GetStdDevError(2))

     #std_devs.Fit("pol3","LL","", 0,1)
     std_devs.Draw("ep")
     hfile.Write()

