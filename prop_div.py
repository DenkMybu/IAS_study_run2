import ROOT
from ROOT import TCanvas, TFile, TH2F, TH2F,TH2D,TH3D, TGraph, TGraphErrors, TGraphAsymmErrors, TString, TPaveLabel, TRandom3, TPad, TColor, TLine, TH1D, TLegend
import sys
import argparse
import math




#check if normalised ou non -> tirage inutile si oui (dans l'histo 3D ou a la fin)

#python3 prop_div.py ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root triple

#python prop_div.py triple


if __name__ == "__main__":

     
     file_trf = "ias_study_2018A_template_noiso_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root"

     file_trf_iso = "ias_study_2018A_template_iso50_2018A_15mars_eta_1_pubins_5_minp_10_final_HSCP_iso50.root"

     file0 = ROOT.TFile.Open(file_trf)
     if not file0:
         print("could not open file ", file_trf)
     else:
         print("Successfully opened root : ",file_trf)

     file1 = ROOT.TFile.Open(file_trf_iso)
     if not file1:
         print("could not open file ", file_trf_iso)
     else:
         print("Successfully opened root : ",file_trf_iso)


     ias_all_single_noiso = file0.Get("Ias_all_base_cutIH") 
     ias_all_triple_noiso = file0.Get("Ias_all_triple_cutIH")
       
     ias_s_noiso = ias_all_single_noiso.Clone("Ias_all_base_ihcut_noiso")
     ias_t_noiso = ias_all_triple_noiso.Clone("Ias_all_triple_ihcut_noiso")

     ias_all_single_iso = file1.Get("Ias_all_base_cutIH") 
     ias_all_triple_iso = file1.Get("Ias_all_triple_cutIH") 
     ias_s_iso = ias_all_single_iso.Clone("Ias_all_base_ihcut_iso")
     ias_t_iso = ias_all_triple_iso.Clone("Ias_all_triple_ihcut_iso")
    
     print(sys.argv[1])
     print(sys.argv[2])
     arg1 = sys.argv[1]
     arg2 = sys.argv[2]
     print(type(arg1))
     print(type(arg2))
      
     if arg2 == "iso":
         if arg1 == "triple":
             hfile = TFile( 'Ratio_plots_PU_triple_iso.root', 'RECREATE' )
         elif arg1 == "base":
             hfile = TFile( 'Ratio_plots_PU_base_iso.root', 'RECREATE' )
         else:
             print("problem with arg 1")
     elif arg2 == "noiso":
         if arg1 == "triple":
             hfile = TFile( 'Ratio_plots_PU_triple_noiso.root', 'RECREATE' )
         elif arg1 == "base":
             hfile = TFile( 'Ratio_plots_PU_base_noiso.root', 'RECREATE' )
         else:
             print("problem with arg 1")
     else:
         print("problem with arg2")

     ias_s_noiso.Scale(1./ias_s_noiso.Integral()) 
     ias_t_noiso.Scale(1./ias_t_noiso.Integral())
     ias_s_iso.Scale(1./ias_s_iso.Integral())
     ias_t_iso.Scale(1./ias_t_iso.Integral())

     ias_s_noiso.Write()
     ias_t_noiso.Write()
     ias_s_iso.Write()
     ias_t_iso.Write()
     
     all_iso_diff = TCanvas("iso_noiso_single","iso_noiso_single",800,800)
     #all_iso_diff = TH1D("ias_var_noiso_write_tpt","aha",50,0,1)
     ias_s_iso.SetMarkerStyle(2)
     ias_s_iso.SetMarkerColor(2)
     ias_s_iso.SetLineColor(2)
     ias_s_iso.Draw()
     ias_s_noiso.SetMarkerStyle(2)
     ias_s_noiso.SetMarkerColor(4)
     ias_s_noiso.SetLineColor(4)
     ias_s_noiso.Draw("same")
     lej_s = TLegend(0.1,0.7,0.5,0.9)
     lej_s.AddEntry(ias_s_iso,"iso","l")
     lej_s.AddEntry(ias_s_noiso,"noiso","l")
     lej_s.Draw("same")

     all_iso_diff.Write()
     #all_iso_diff.Write()
     #diff_iso_tpt.Write()
     #all_iso_diff2 = TH1D("ias_var_iso_write_tpt","aha",50,0,1)
     all_iso_diff2 = TCanvas("iso_noiso_triple","iso_noiso_triple",800,800)
     ias_t_iso.SetMarkerStyle(2)
     ias_t_iso.SetMarkerColor(2)
     ias_t_iso.SetLineColor(2)
     ias_t_iso.Draw()
     ias_t_noiso.SetMarkerStyle(2)
     ias_t_noiso.SetMarkerColor(4)
     ias_t_noiso.SetLineColor(4)
     ias_t_noiso.Draw("same")
     lej_t = TLegend(0.1,0.7,0.5,0.9)
     lej_t.AddEntry(ias_t_iso,"iso","l")
     lej_t.AddEntry(ias_t_noiso,"noiso","l")
     lej_t.Draw("same")
     all_iso_diff2.Write()

     if(sys.argv[1] == "triple"):
         if(sys.argv[2] == "noiso"):
             hden = file0.Get("Ias_when_PU_between_0_20_triple_cutIH")
             hnum = file0.Get("Ias_when_PU_between_25_30_triple_cutIH")
         else:
             hden = file1.Get("Ias_when_PU_between_0_20_triple_cutIH")
             hnum = file1.Get("Ias_when_PU_between_25_30_triple_cutIH")
       
     else:
         if(sys.argv[2] == "noiso"):
             hden = file0.Get("Ias_when_PU_between_0_20_base_cutIH")
             hnum = file0.Get("Ias_when_PU_between_25_30_base_cutIH")
         else:
             hden = file1.Get("Ias_when_PU_between_0_20_base_cutIH")
             hnum = file1.Get("Ias_when_PU_between_25_30_base_cutIH")
     

     all_triple = file0.Get("Ias_all_triple_nocut")
     all_base = file0.Get("Ias_all_base_nocut")


     low_bin = 0.4
     low_triple = all_triple.GetXaxis().FindBin(low_bin)
     up_triple = all_triple.GetXaxis().FindBin(1)

     low_base = all_base.GetXaxis().FindBin(low_bin)
     up_base = all_base.GetXaxis().FindBin(1)
     
     itg_base = all_base.Integral(low_base,up_base)

     itg_triple = all_triple.Integral(low_triple,up_triple)
     
     print("Integral",low_bin," to 1 of single tpt : ", itg_base, " and triple tpt : ", itg_triple)
     
     print("Difference ITG(simple) - ITG(triple) = ", itg_base - itg_triple)

     hname = hnum.GetName()
     hnameden = hden.GetName()
 
     hnum_norm = hnum.Clone("hnum_normalised")
     hden_norm = hden.Clone("hden normalised")
    


     hnum_norm.Scale(1./hnum_norm.Integral())
     hden_norm.Scale(1./hden_norm.Integral())
     print("hden tot intgral> ", hden_norm.Integral())
 
     print("hden tot intgral between bins 1-49> ", hden_norm.Integral(1,49))

     hnum_norm.SetTitle("")
     hnum_norm.GetYaxis().SetTitle("#events/tot")
     hnum_norm.GetYaxis().SetTitleOffset(2.8)
     hnum_norm.GetYaxis().SetTitleFont(43)
     hnum_norm.GetYaxis().SetTitleSize(20)
     
     hnum_norm.SetLineColor(4)
     hnum_norm.SetLineWidth(2)
     hnum_norm.GetYaxis().SetTitleSize(15)
     hnum_norm.GetYaxis().SetTitleFont(43)
     hnum_norm.GetYaxis().SetTitleOffset(1)
     hnum_norm.Write()
    
     hden_norm.SetLineColor(2)
     hden_norm.SetLineWidth(2)


     hden_norm.Write()
     
     xtitle = "IAS"



     c =  TCanvas("c","canvas",800,800)
     pad1 =  TPad("pad1","pad1",0,0.3,1,0.95)
     pad1.SetLogy()

     #pad1.SetBottomMargin(0)
     pad1.SetGridx()
     pad1.Draw()
     pad1.cd()

     labellumi = TPaveLabel(0.6,0.95,0.5,1,"SingleMu 2018A, 13.6 TeV","brNDC")
     labellumi.SetFillColor(0)
     labellumi.SetTextColor(1)
     labellumi.SetFillStyle(0)
     labellumi.SetBorderSize(0)
     labellumi.SetTextSize(0.25)
     labellumi.SetTextAlign(12)
     labellumi.Draw("same")

     hnum_norm.SetStats(0)
     hnum_norm.Draw()
     hden_norm.Draw("same")

     axis = hnum_norm.GetYaxis()
     axis.ChangeLabel(1,-1,-1,-1,-1,-1,"")
     axis.SetLabelFont(43)
     axis.SetLabelSize(15)
     axis.SetTitleOffset(1.5)


     '''
     c.cd()
     pad2 = TPad("pad2","pad2",0,0.3,1,0.5)
     #pad2.SetTopMargin(0)
     pad2.SetBottomMargin(0.1)
     pad2.SetLogy()
     pad2.SetGridx()
     pad2.Draw()
     pad2.cd()
 
     ratio_plot = hnum_norm.Clone("h3")

     ratio_plot.SetLineColor(1)
     ratio_plot.SetTitle("")
          
     ratio_plot.GetYaxis().SetTitle(" blue / red ")
     ratio_plot.GetYaxis().SetNdivisions(505)
     ratio_plot.GetYaxis().SetTitleSize(20)
     ratio_plot.GetYaxis().SetTitleFont(43)
     ratio_plot.GetYaxis().SetTitleOffset(1.2)
     ratio_plot.GetYaxis().SetLabelFont(43)
     ratio_plot.GetYaxis().SetLabelSize(15)

     ratio_plot.GetXaxis().SetTitleSize(20)
     ratio_plot.GetXaxis().SetTitleFont(43)
     #ratio_plot.GetXaxis().SetTitleOffset(3.5)
     ratio_plot.GetXaxis().SetLabelFont(43)
     ratio_plot.GetXaxis().SetLabelSize(15)
     ratio_plot.SetMinimum(0.3)
     ratio_plot.SetMaximum(3)
     ratio_plot.Sumw2()
     ratio_plot.SetStats(0)
     ratio_plot.Divide(hden_norm)
     ratio_plot.SetMarkerStyle(2)
     ratio_plot.SetMarkerSize(1)
     line = TLine(0,1,1,1)
     line.SetLineStyle(4)

     line_q1 = TLine(0.115,0.3,0.115,3)
     line_q1.SetLineStyle(3)
     line_q1.SetLineColor(2)
     line_q1.SetLineWidth(2)
     
     line_q2 = TLine(0.19,0.3,0.19,3)
     line_q2.SetLineStyle(3)
     line_q2.SetLineColor(2)
     line_q2.SetLineWidth(2)
     
     ratio_plot.Draw("ep")
     line_q1.Draw("same")
     line_q2.Draw("same")
     line.Draw("same")
     '''
     c.cd()
     pad3 = TPad("pad3","pad3",0,0.07,1,0.27)
     #pad3.SetTopMargin(0)
     #pad3.SetLogy()
     #pad3.SetBottomMargin(0.1)
     
     pad3.SetGridx()
     pad3.Draw()
     pad3.cd()

     ratio_itg = TH1D("ratio_itg","ratio_itg",50,0,1)
     ratio_itg.SetTitle("")
     ratio_itg.Sumw2()
     ratio_itg.GetYaxis().SetTitle("Right Integrals blue / red ")
     ratio_itg.SetMinimum(0.5)
     ratio_itg.SetMaximum(2)
     ratio_itg.SetMarkerStyle(2)
     ratio_itg.SetMarkerSize(1)
     
     print("there are ", ratio_itg.GetNbinsX(), " x bins")
     print("ITG num 50-50 = ", hnum_norm.Integral(50,50))
     print("ITG den 50-50 = ", hden_norm.Integral(50,50))
     print("ITG num  0-0 = ", hnum_norm.Integral(0,0))
     print("ITG den  0-0 = ", hden_norm.Integral(0,0))


     for i in range(0,ratio_itg.GetNbinsX()+1):
         err1 = ROOT.Double(0)
         it = hnum_norm.IntegralAndError(i,50,err1)
         #print("it num",i, " = ", it)
         err2 = ROOT.Double(0)
         it2 = hden_norm.IntegralAndError(i,50,err2)
         #print("it den",i, " = ", it2)
         if(it2 !=0):
             ratio = it/it2

             if(it != 0 and it2 != 0):
                 if ( ( (err1/it)**2 + (err2/it2)**2 )- (2*(err1/it)*(err2/it2)) > 0 ):
                      err_ratio = ratio * math.sqrt( ((err1/it)**2) + ((err2/it2)**2) - (2*(err1/it)*(err2/it2)))
                 else:
                     err_ratio = 0

             ratio_itg.SetBinContent(i,ratio)
             ratio_itg.SetBinError(i,err_ratio)
             print("bin #",i, " gets ratio = ",it, " / ", it2, " = ", it/it2)
             #SetBinError -> propagate errors each number


     ratio_itg.GetYaxis().SetNdivisions(505)
     
     ratio_itg.GetYaxis().SetLabelFont(43)
     ratio_itg.GetYaxis().SetLabelSize(15)
     ratio_itg.GetYaxis().SetTitleOffset(1.2)
     ratio_itg.GetYaxis().SetTitleSize(20)

     ratio_itg.GetXaxis().SetLabelFont(43)
     ratio_itg.GetXaxis().SetLabelSize(15)
     #ratio_itg.GetXaxis().SetTitleOffset(1)
     ratio_itg.GetXaxis().SetTitle("IAS")
     ratio_itg.GetXaxis().SetTitleSize(20)

     ratio_itg.SetStats(0)

     ratio_itg.Draw("ep")
     line2 = TLine(0,1,0.6,1)
     line2.SetLineStyle(4)
     line2.Draw("same")

     line3 = TLine(0.115,0.5,0.115,2)
     line3.SetLineStyle(3)
     line3.SetLineColor(2)
     line3.SetLineWidth(2)
     line3.Draw("same")
     
     line4 = TLine(0.19,0.5,0.19,2)
     line4.SetLineStyle(3)
     line4.SetLineColor(2)
     line4.SetLineWidth(2)
     line4.Draw("same")
     c.Write()


     c.SaveAs(hname+"_over_"+hnameden+".pdf")
     hfile.Write()

