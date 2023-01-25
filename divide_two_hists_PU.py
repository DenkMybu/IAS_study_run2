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


     #hist_num_raw = file0.Get("P_single_pt50_ias_qtl_80_90") 
     #hist_denom_raw = file0.Get("P_triple_pt50_ias_qtl_80_90")

     '''
     str_num = "P_LOWPU_pt50_cutih"
     str_denom = "P_HIGHPU_pt50_cutih"
     '''

     '''
     str_num = "P_PU_28_29_pt50_cutih"
     str_denom = "P_PU_30_31_pt50_cutih"
     '''

     str_num = "PU_after_presel"
     str_denom = "PU_before_presel"

     hist_num_raw = file0.Get(str_num) 
     hist_denom_raw = file0.Get(str_denom)

     name_xaxis = hist_num_raw.GetXaxis().GetTitle()
     print("name of x axis = ", name_xaxis)
     num_binmax = 100
     

     num_nbin = hist_num_raw.GetNbinsX()
     print("For num , there are : ", num_nbin, " bins in X axis")

     low_bin = 0
     high_bin = 100

     low_val_num_raw = hist_num_raw.GetXaxis().FindBin(low_bin)
     up_val_num_raw = hist_num_raw.GetXaxis().FindBin(high_bin)

     low_val_denom_raw = hist_denom_raw.GetXaxis().FindBin(low_bin)
     up_val_denom_raw = hist_denom_raw.GetXaxis().FindBin(high_bin)

    
     itg_num_raw = hist_num_raw.Integral(low_val_num_raw,up_val_num_raw)
     itg_denom_raw = hist_denom_raw.Integral(low_val_denom_raw,up_val_denom_raw)
     print("BEFORE SCALING")
     print("Integral ",low_bin," to ", high_bin, " from single  : ", itg_num_raw, " and triple : ", itg_denom_raw)
     
     print("Difference ITG(simple) - ITG(triple) = ", itg_num_raw - itg_denom_raw)

     hist_num = hist_num_raw.Clone("Clone_num")
     hist_denom = hist_denom_raw.Clone("Clone_denom")

     hname = hist_num_raw.GetName()
     hnameden = hist_denom_raw.GetName()

     hist_num.Sumw2()
     hist_denom.Sumw2() 
     hfile = TFile( hname+"_over_"+hnameden+".root", 'RECREATE' )

     hist_num.Scale(1./hist_num.Integral())
     hist_denom.Scale(1./hist_denom.Integral())


     title_rawdiv = hname + "/" + hnameden    
     raw_divide = TCanvas(title_rawdiv,title_rawdiv,800,800)
     #all_iso_diff = TH1D("ias_var_noiso_write_tpt","aha",50,0,1)
     hist_num.SetMarkerStyle(2)
     hist_num.SetMarkerColor(4)
     hist_num.SetLineColor(4)
     hist_num.SetLineWidth(1)
     hist_num.Draw()
     hist_denom.SetMarkerStyle(2)
     hist_denom.SetMarkerColor(2)
     hist_denom.SetLineColor(2)
     hist_denom.SetLineWidth(1)
     hist_denom.Draw("same")
     lej_denom = TLegend(0.1,0.7,0.5,0.9)
     lej_denom.AddEntry(hist_num,hname,"l")
     lej_denom.AddEntry(hist_denom,hnameden,"l")
     lej_denom.Draw("same")
     raw_divide.Write()


     low_val_num = hist_num.GetXaxis().FindBin(low_bin)
     up_val_num = hist_num.GetXaxis().FindBin(high_bin)

     low_val_denom = hist_denom.GetXaxis().FindBin(low_bin)
     up_val_denom = hist_denom.GetXaxis().FindBin(high_bin)

    
     itg_num = hist_num.Integral(low_val_num,up_val_num)
     itg_denom = hist_denom.Integral(low_val_denom,up_val_denom)
     
     print("Integral ",low_bin," to ", high_bin, " from single  : ", itg_num, " and triple : ", itg_denom)
     
     print("Difference ITG(simple) - ITG(triple) = ", itg_num - itg_denom)


     hist_num.GetYaxis().SetTitle("#events/tot")
     hist_num.GetYaxis().SetTitleOffset(2.8)
     hist_num.GetYaxis().SetTitleFont(43)
     hist_num.GetYaxis().SetTitleSize(20)
     
     hist_num.SetLineColor(4)
     hist_num.SetLineWidth(2)
     hist_num.GetXaxis().SetTitle("NPV")
     hist_num.GetXaxis().SetTitleSize(15)
     hist_num.GetXaxis().SetTitleFont(43)
     hist_num.GetXaxis().SetTitleOffset(2)
     hist_num.Write()
  

     hist_denom.GetYaxis().SetTitle("#events/tot")
     hist_denom.GetYaxis().SetTitleOffset(1)
     hist_denom.GetYaxis().SetTitleFont(43)
     hist_denom.GetYaxis().SetTitleSize(20)

     hist_denom.GetXaxis().SetTitle("NPV")
     hist_denom.GetXaxis().SetTitleSize(15)
     hist_denom.GetXaxis().SetTitleFont(43)
     hist_denom.GetXaxis().SetTitleOffset(1)

     hist_denom.SetLineColor(2)
     hist_denom.SetLineWidth(2)
     hist_denom.Write()
     

     c =  TCanvas("c","canvas",800,800)
     pad1 =  TPad("pad1","pad1",0,0.3,1,0.95)
     pad1.SetLogy()

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

     hist_num.SetStats(0)
     hist_num.Draw()
     hist_denom.Draw("same")


     lej_brut = TLegend(0.5,0.8,0.6,0.9)
     lej_brut.AddEntry(hist_num,str_num,"l")
     lej_brut.AddEntry(hist_denom,str_denom,"l")
     lej_brut.Draw("same")

     axis = hist_num.GetYaxis()
     axis.ChangeLabel(1,-1,-1,-1,-1,-1,"")
     axis.SetLabelFont(43)
     axis.SetLabelSize(15)
     axis.SetTitleOffset(1.5)

     c.cd()

     '''
     pad2 = TPad("pad2","pad2",0,0.3,1,0.5)
     #pad2.SetTopMargin(0)
     pad2.SetBottomMargin(0.1)
     #pad2.SetLogy()
     pad2.SetGridx()
     pad2.Draw()
     pad2.cd()
 
     ratio_plot = hist_num.Clone("h3")
     
     #ratio_plot.SetLineColor(1)
     ratio_plot.SetTitle("")


          
     ratio_plot.GetYaxis().SetTitle(" blue / red ")
     ratio_plot.GetYaxis().SetNdivisions(505)
     ratio_plot.GetYaxis().SetTitleSize(15)
     ratio_plot.GetYaxis().SetTitleFont(43)
     ratio_plot.GetYaxis().SetTitleOffset(2)
     ratio_plot.GetYaxis().SetLabelFont(43)
     ratio_plot.GetYaxis().SetLabelSize(15)
     ratio_plot.GetXaxis().SetTitle("")
     ratio_plot.GetXaxis().SetTitleSize(15)
     ratio_plot.GetXaxis().SetTitleFont(43)
     ratio_plot.GetXaxis().SetTitleOffset(2)
     ratio_plot.GetXaxis().SetLabelFont(43)
     ratio_plot.GetXaxis().SetLabelSize(15)
     ratio_plot.SetMinimum(0.3)
     ratio_plot.SetMaximum(3)
     ratio_plot.Sumw2()
     ratio_plot.SetStats(0)
     ratio_plot.Divide(hist_denom)
     ratio_plot.SetMarkerStyle(2)
     ratio_plot.SetMarkerSize(1)
     ratio_plot.SetLineWidth(1)
     ratio_plot.SetLineColor(15)
     ratio_plot.SetMarkerColor(15)
    
     ratio_plot.Draw("ep")

     line = TLine(0,1,1000,1)
     line.SetLineStyle(4)

     
     line_q1 = TLine(200,0.3,200,3)
     line_q1.SetLineStyle(3)
     line_q1.SetLineColor(2)
     line_q1.SetLineWidth(2)
     
     line_q2 = TLine(1000,0.3,1000,3)
     line_q2.SetLineStyle(3)
     line_q2.SetLineColor(2)
     line_q2.SetLineWidth(2)
     
     line_q1.Draw("same")
     line_q2.Draw("same")
     line.Draw("same")
     
     c.cd()
     '''
     pad3 = TPad("pad3","pad3",0,0.07,1,0.27)
     #pad3.SetTopMargin(0)
     #pad3.SetLogy()
     #pad3.SetBottomMargin(0.1)
     
     pad3.SetGridx()
     pad3.Draw()
     pad3.cd()

     ratio_itg = TH1D("ratio_itg","ratio_itg",num_nbin,0,num_binmax)
     ratio_itg.Sumw2()
     ratio_itg.SetTitle("")
     ratio_itg.GetYaxis().SetTitle("Right Integrals blue / red ")
     ratio_itg.GetYaxis().SetNdivisions(505)
     ratio_itg.GetYaxis().SetTitleSize(15)
     ratio_itg.GetYaxis().SetTitleFont(43)
     ratio_itg.GetYaxis().SetTitleOffset(2)
     ratio_itg.GetYaxis().SetLabelFont(43)
     ratio_itg.GetYaxis().SetLabelSize(15)


     upper_lim_ratio = 2.2
     lower_lim_ratio = 0.8

     ratio_itg.SetMinimum(lower_lim_ratio)
     ratio_itg.SetMaximum(upper_lim_ratio)
     ratio_itg.SetMarkerStyle(1)
     ratio_itg.SetMarkerSize(1)
     ratio_itg.SetLineColor(1) 
     print("there are ", ratio_itg.GetNbinsX(), " x bins")
     print("ITG num bins 200-200 = ", hist_num.Integral(200,200))
     print("ITG den bins 200-200 = ", hist_denom.Integral(200,200))
     print("ITG num bins 0-0 = ", hist_num.Integral(0,0))
     print("ITG den bins 0-0 = ", hist_denom.Integral(0,0))


     tot_it_num = hist_num_raw.Integral()
     tot_it_denom = hist_denom_raw.Integral()
 
     print("Before renormalization, integral num = ", tot_it_num, " and integral denom = " , tot_it_denom)   

     tot_it_num_renorm = hist_num.Integral()
     tot_it_denom_renorm = hist_denom.Integral()
     print("After renormalization, integral num = ", tot_it_num_renorm, " and integral denom = " , tot_it_denom_renorm)   

     for i in range(0,ratio_itg.GetNbinsX()+1):
      
         err1 = ROOT.Double(0)
         it = hist_num.IntegralAndError(i,ratio_itg.GetNbinsX()+1,err1)
         err1_flat = ROOT.Double(0)
         it_flat = hist_num_raw.IntegralAndError(i,ratio_itg.GetNbinsX()+1,err1_flat)
  
         #print("it num",i, " = ", it)
         err2 = ROOT.Double(0)
         it2 = hist_denom.IntegralAndError(i,ratio_itg.GetNbinsX()+1,err2)
         err2_flat = ROOT.Double(0)
         it2_flat = hist_denom_raw.IntegralAndError(i,ratio_itg.GetNbinsX()+1,err2_flat)
         #print("it den",i, " = ", it2)
         if(it2 !=0):
             print("BIN #",i, " , x center = ", ratio_itg.GetBinCenter(i))
             print("error on it num = ", err1)
             print("error on it denom = ", err2)
             ratio = it/it2
             if(it != 0 and it2 != 0):
                 if ( ( (err1/it)**2 + (err2/it2)**2 )- (2*(err1/it)*(err2/it2)) > 0 ):
                     err_ratio = ratio * math.sqrt( ((err1/it)**2) + ((err2/it2)**2))# - (2*(err1/it)*(err2/it2)))
                 else:
                     err_ratio = 0
                 ratio_itg.SetBinContent(i,ratio)
                 ratio_itg.SetBinError(i,err_ratio)
                
                 print("Before renorm ratio = ",it_flat, " / ", it2_flat ," = ", it_flat/it2_flat)
                 print("After renorm ratio = ",it, " / ", it2, " = ", it/it2, " with prop error = ", err_ratio)




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
     line2 = TLine(0,1,1000,1)
     line2.SetLineStyle(4)
     line2.Draw("same")

     line3 = TLine(200,lower_lim_ratio,200,upper_lim_ratio)
     line3.SetLineStyle(3)
     line3.SetLineColor(2)
     line3.SetLineWidth(2)
     line3.Draw("same")
     
     line4 = TLine(1000,lower_lim_ratio,1000,upper_lim_ratio)
     line4.SetLineStyle(3)
     line4.SetLineColor(2)
     line4.SetLineWidth(2)
     line4.Draw("same")
     c.Write()


     c.SaveAs(hname+"_over_"+hnameden+".pdf")
     hfile.Write()

