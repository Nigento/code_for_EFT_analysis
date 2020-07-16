#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <sstream>
#include "TStyle.h"
#include <iostream>
#include <TStyle.h>
#include <string>
#include <TGraph.h>
#include <string>
#include <TF1.h>
using namespace std;

TH1D* GetHistoWeight(TTree* t, string variable, int nbins, double xmin, double xmax, string cut, string name)
{
        string sxmin, sxmax, snbins;
        stringstream ss[3];

        ss[0] << xmin;
        ss[0] >> sxmin;
        ss[1] << xmax;
        ss[1] >> sxmax;
        ss[2] << nbins;
        ss[2] >> snbins;

        string variablenew = variable + " >> h(" + snbins + "," + sxmin + "," + sxmax + ")";

        string cutnew = "1 * (" + cut + ")";
//
        t->Draw(variablenew.c_str(), cutnew.c_str());
        TH1D *histo = (TH1D*)gDirectory->Get("h");

		if (histo->GetEntries()==0) return histo;

		double underflow = histo->GetBinContent(0);
		cout << "underflow="<<underflow<<endl;
		double val = 0;
		if (underflow>0) {
			val = histo->GetBinContent(1);
			histo->SetBinContent(1, val+underflow);
			 histo->SetBinContent(0, 0);
		}
		double overflow = histo->GetBinContent(nbins+1);
		if (overflow>0) {
		  val = histo->GetBinContent(nbins);
		  histo->SetBinContent(nbins+1, 0);
		  histo->SetBinContent(nbins, val+overflow);
		}

		cout << "Area="<<histo->Integral()<<endl;
		cout << "Nevents="<<histo->GetEntries()<<endl;
        histo->SetName(name.c_str());
        histo->SetTitle(name.c_str());

        return histo;
}

void Ratio_EFT_SM(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, string variable, string EFT, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string file_name, string lepton)
{
  //TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  //TTree* tree_file = new TTree("events","events");

  //tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m1= GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");

  Ratio_p2_SM->Divide(Histo_SM);
  Ratio_p1_SM->Divide(Histo_SM);
  Ratio_m2_SM->Divide(Histo_SM);
  Ratio_m1_SM->Divide(Histo_SM);

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);

  Ratio_m1_SM->SetLineColor(kBlack);
  Ratio_m1_SM->SetLineWidth(2);


  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());

  double max = (Ratio_p2_SM->GetMaximum()>Ratio_p1_SM->GetMaximum()) ? Ratio_p2_SM->GetMaximum() : Ratio_p1_SM->GetMaximum();

  Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->Draw("");
  Ratio_p1_SM->Draw("SAME");
  Ratio_m1_SM->Draw("SAME");
  Ratio_m2_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;
  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\Lambda^{2} = -2 (TeV^{-2})";


  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");

  legend->Draw("SAME");

  TGraph** ratio_Histo = new TGraph*[nbins];

    //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  TFile* ratio_file = new TFile(("signal_proc_"+variable+"_"+EFT+"_"+lepton+"_TF1.root").c_str(),"RECREATE");
  ratio_file->cd();
  Canvas->Print(file_name.c_str());

  for (int i = 1 ; i<=nbins ; i++)
  {
    string number_plot = to_string(i);
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_par1_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -3 , 3);

    ratio_Histo[i]->SetPoint(0,-2,Histo_EFT_m2->GetBinContent(i)/Histo_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(1,-1,Histo_EFT_m1->GetBinContent(i)/Histo_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i)/Histo_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(3,1,Histo_EFT_p1->GetBinContent(i)/Histo_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(4,2,Histo_EFT_p2->GetBinContent(i)/Histo_SM->GetBinContent(i));

    //tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());
    /*TLegend* legend2 = new TLegend(0.6, 0.7, 0.89, 0.89, "");
    legend2->SetTextSize(0.05);
    legend2->AddEntry(ratio_Histo[i]->GetName(),"0 < #phi^{*} < 1.25","l" );

    legend2->Draw("SAME");
    //Canvas->Print(file_name_eft.c_str());
*/
    //ratio_Histo[i]->Write(("bin_content_par1_"+number_plot).c_str());
    ratio_formula->Write();
  }


  //file_output->Write();
  ratio_file->Close();




}

void Compare_3Histos(TTree* t1, TTree* t2, TTree* t3, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  cout << "max="<<max<<endl;

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.7,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.5,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.5,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");


  Histo_3->SetLineColor(kOrange);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.7;
	 lx1 = 0.6;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.6;
	 ly0 = 0.6;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Histo_1->GetName(), "C_{bW}^{I}/#\Lambda^{2} = 0 (TeV^{-2})","l");
  legend->AddEntry(Histo_2->GetName(), "C_{bW}^{I}/#\Lambda^{2} = -2 (TeV^{-2})","l");
  legend->AddEntry(Histo_3->GetName(), "C_{bW}^{I}/#\Lambda^{2} = 2 (TeV^{-2})","l");

  legend->Draw("SAME");


  Canvas->Print(Name.c_str());
}

void Compare_1Histos(TTree* t1, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  Histo_1->SetStats(kFALSE);

  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.75;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }

   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");

   legend->Draw("SAME");

   Canvas->Print(Name.c_str());
}

void Compare_2Histos(TTree* t1, TTree* t2, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);

  double a = Histo_1->Integral();
  double b = Histo_2->Integral();
  Histo_1->Scale(1/a);
  Histo_2->Scale(1/b);
  cout<<legendX<<endl;
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;

}

void Compare_4Histos(TTree* t1, TTree* t2, TTree* t3, TTree* t4, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string legendEntry4, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "Histo_4");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);

  double a = Histo_1->Integral();
  double b = Histo_2->Integral();
  double c = Histo_3->Integral();
  double d = Histo_4->Integral();

  Histo_1->Scale(1/a);
  Histo_2->Scale(1/b);
  if (c>0) Histo_3->Scale(1/c);
  Histo_4->Scale(1/d);
  cout << "a="<<a<<" b="<<b<<" c="<<c<<endl;

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  if (c>0) max = (max>Histo_3->GetMaximum()) ? max : Histo_3->GetMaximum();
  if (d>0) max = (max>Histo_4->GetMaximum()) ? max : Histo_4->GetMaximum();
  cout << "max="<<max<<endl;

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  if (c>0){
    Histo_3->SetLineColor(kOrange);
    Histo_3->SetLineWidth(2);
    Histo_3->Draw("SAME");
  }

  Histo_4->SetLineColor(kGreen);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), "MadGraph EFT = 0", "l");
  legend->AddEntry(Histo_2->GetName(), "MadSpin EFT = 0", "l");
  legend->AddEntry(Histo_3->GetName(), "MadGraph EFT = 2", "l");
  legend->AddEntry(Histo_4->GetName(), "MadSpin EFT = 2", "l");
  legend->Draw("SAME");

  Canvas->Print(Name.c_str());

  cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;
  if (c>0) cout << "Histo3 mean: "<<Histo_3->GetMean()<<endl;

}

int main (){

string suffix[17];
suffix[0] = "output_madgraph_SM.root";
suffix[1] = "output_madgraph_cbwi_m1.root";
suffix[2] = "output_madgraph_cbwi_m2.root";
suffix[3] = "output_madgraph_cbwi_p1.root";
suffix[4] = "output_madgraph_cbwi_p2.root";
suffix[5] = "output_madgraph_cptbi_m1.root";
suffix[6] = "output_madgraph_cptbi_m2.root";
suffix[7] = "output_madgraph_cptbi_p1.root";
suffix[8] = "output_madgraph_cptbi_p2.root";
suffix[9] = "output_madgraph_ctw_m1.root";
suffix[10] = "output_madgraph_ctw_m2.root";
suffix[11] = "output_madgraph_ctw_p1.root";
suffix[12] = "output_madgraph_ctw_p2.root";
suffix[13] = "output_madgraph_ctwi_m1.root";
suffix[14] = "output_madgraph_ctwi_m2.root";
suffix[15] = "output_madgraph_ctwi_p1.root";
suffix[16] = "output_madgraph_ctwi_p2.root";
/*suffix[17] = "output_madspin_SM.root";
suffix[18] = "output_madspin_cbwi_m1.root";
suffix[19] = "output_madspin_cbwi_m2.root";
suffix[20] = "output_madspin_cbwi_p1.root";
suffix[21] = "output_madspin_cbwi_p2.root";
suffix[22] = "output_madspin_cptbi_m1.root";
suffix[23] = "output_madspin_cptbi_m2.root";
suffix[24] = "output_madspin_cptbi_p1.root";
suffix[25] = "output_madspin_cptbi_p.root";
suffix[26] = "output_madspin_ctw_m1.root";
suffix[27] = "output_madspin_ctw_m2.root";
suffix[28] = "output_madspin_ctw_p1.root";
suffix[29] = "output_madspin_ctw_p2.root";
suffix[30] = "output_madspin_ctwi_m1.root";
suffix[31] = "output_madspin_ctwi_m2.root";
suffix[32] = "output_madspin_ctwi_p1.root";
suffix[33] = "output_madspin_ctwi_p2.root";
suffix[34] = "/heppy/output_t_chan_MC.root";*/

	TFile* fInput[17];
	TTree* tInput[17];
	string inputName;

  for (int i=0; i<17; i++)
    {
      inputName = "data/madgraph/output/" + suffix[i];
      fInput[i] = new TFile(inputName.c_str(),"READ");
      tInput[i] = (TTree*) fInput[i]->Get("Tree");
    }
/*
  //Plot_SM
  //Compare_2Histos(tInput[0], tInput[17], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_CosTheta.pdf"  );
  Compare_2Histos(tInput[0], tInput[17], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_CosThetaStar.pdf");
  //Compare_2Histos(tInput[0], tInput[17], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_PhiStar.pdf");
  //Compare_2Histos(tInput[0], tInput[17], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_top_pt.pdf");
  Compare_2Histos(tInput[0], tInput[17], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_W_pt.pdf");
  Compare_2Histos(tInput[0], tInput[17], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[17], suffix[17], "results/dim6top_compareSM_lepton_pt.pdf");
	Compare_2Histos(tInput[0], tInput[17], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_lepton_E_Wframe.pdf");
  Compare_2Histos(tInput[0], tInput[17], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_top_mass.pdf");
  Compare_2Histos(tInput[0], tInput[17], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], "results/dim6top_compareSM_W_mass.pdf");

  //MadGraph

  //Plots ctW
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "cosTheta", 20, -1, 1, "1", "cos#theta", "a.u.", "legendUpLeft", "Operateur EFT", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_CosTheta.pdf");/*
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_CosThetaStar.pdf");
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_PhiStar.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_top_pt.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "Number of events", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_W_pt.pdf");
  /*Compare_3Histos(tInput[0], tInput[10], tInput[12], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_lepton_pt.pdf");
	Compare_3Histos(tInput[0], tInput[10], tInput[12], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "a.u.", "legendUpRight", "Operateur EFT", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_W_mass.pdf");
  Compare_3Histos(tInput[0], tInput[10], tInput[12], "W_transverse_mass",20, 50, 140, "1", "M_{T,W} (GeV)", "number of events", "legendUpRight", "Operateur de dimension 6", suffix[0], suffix[10], suffix[12], "results/madgraph_dim6top_ctW_W_transverse_mass.pdf");

	//Plots ctWI*/
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_CosTheta.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_CosThetaStar.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_PhiStar.pdf");
  /*(tInput[0], tInput[14], tInput[16], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[14], tInput[16], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[14], tInput[16], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[14], tInput[16], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0],suffix[14], suffix[16], "results/madgraph_dim6top_ctWI_lepton_E_Wframe.pdf");
  *///Compare_3Histos(tInput[0], tInput[14], tInput[16], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_ctwI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[14], tInput[16], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[14], suffix[16], "results/madgraph_dim6top_W_mass.pdf");

	//Plots cbWI
	//Compare_3Histos(tInput[0], tInput[2], tInput[4], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_CosTheta.pdf");
	*/Compare_3Histos(tInput[0], tInput[2], tInput[4], "cosThetaStar", 20, -1, 1, "1", "cos(#theta*)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_CosThetaStar.pdf");
	/*Compare_3Histos(tInput[0], tInput[2], tInput[4], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_PhiStar.pdf");
  Compare_3Histos(tInput[0], tInput[2], tInput[4], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[2], tInput[4], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[2], tInput[4], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[2], tInput[4], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbWI_lepton_E_Wframe.pdf");
  *///Compare_3Histos(tInput[0], tInput[2], tInput[4], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbwI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[2], tInput[4], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[2], suffix[4], "results/madgraph_dim6top_cbwI_W_mass.pdf");

  //Plots cptbI
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_CosTheta.pdf");
	*///Compare_3Histos(tInput[0], tInput[6], tInput[8], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_CosThetaStar.pdf");
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_PhiStar.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_top_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_W_pt.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_lepton_pt.pdf");
	//Compare_3Histos(tInput[0], tInput[6], tInput[8], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_lepton_E_Wframe.pdf");
  //Compare_3Histos(tInput[0], tInput[6], tInput[8], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_top_mass.pdf");
  /*Compare_3Histos(tInput[0], tInput[6], tInput[8], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[6], suffix[8], "results/madgraph_dim6top_cptbI_W_mass.pdf");

  //MadSpin

  //Plots ctW
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[27], tInput[29], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27],suffix[29], "results/madspin_dim6top_ctW_top_pt.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[27], tInput[29], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[27], tInput[29], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[27], tInput[29], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[27], suffix[29], "results/madspin_dim6top_ctW_W_mass.pdf");


  //Plots ctWI
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0],  suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[31], tInput[33], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_top_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[31], tInput[33], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[31], tInput[33], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[31],  suffix[33], "results/madspin_dim6top_ctWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[31], tInput[33], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[31], suffix[33], "results/madspin_dim6top_ctWI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[31], tInput[33], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[31], suffix[33], "results/madspin_dim6top_ctWI_W_mass.pdf");


	//Plots cbWI
	//Compare_3Histos(tInput[17], tInput[19], tInput[21], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_CosTheta.pdf");
	Compare_3Histos(tInput[17], tInput[19], tInput[21], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19],suffix[21], "results/madspin_dim6top_cbWI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[19], tInput[21], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_PhiStar.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_top_pt.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[19], tInput[21], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[19], tInput[21], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbWI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbwI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[19], tInput[21], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[19], suffix[21], "results/madspin_dim6top_cbwI_W_mass.pdf");

  //Plots cptbI
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_CosTheta.pdf");
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top",  suffix[0],  suffix[23], suffix[25], "results/madspin_dim6top_cptbI_CosThetaStar.pdf");
	Compare_3Histos(tInput[17], tInput[23], tInput[25], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_PhiStar.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_top_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[23],  suffix[25], "results/madspin_dim6top_cptbI_W_pt.pdf");
  //Compare_3Histos(tInput[17], tInput[23], tInput[25], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top",  suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_lepton_pt.pdf");
	//Compare_3Histos(tInput[17], tInput[23], tInput[25], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_lepton_E_Wframe.pdf");
  Compare_3Histos(tInput[17], tInput[23], tInput[25], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_top_mass.pdf");
  Compare_3Histos(tInput[17], tInput[23], tInput[25], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[23], suffix[25], "results/madspin_dim6top_cptbI_W_mass.pdf");

  //Madgraph + MadSpin

  //cbwi
  //Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[4], tInput[21], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[4], suffix[21], "results/compare_dim6top_cbwI_W_mass.pdf");

  //cptbI
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[8], tInput[25], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[8], suffix[25], "results/compare_dim6top_cptbI_W_mass.pdf");

  //ctw
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "Operateur de dimension 6", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[12], tInput[29], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[12], suffix[29], "results/compare_dim6top_ctw_W_mass.pdf");


  //ctwI
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_CosTheta.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "cosThetaStar", 20, -1, 1, "1", "cos#theta*", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_CosThetaStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "PhiStar", 20, 0, 2*TMath::Pi(), "1", "#phi* (rad)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_PhiStar.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "top_pt", 20, 0, 400, "1", "Top Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_top_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "W_pt", 20, 0, 100, "1", "W Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_W_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "lepton_pt", 20, 0, 100, "1", "Lepton Pt (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_lepton_pt.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "lepton_E_Wframe", 20, 25, 60, "1", "Lepton Energy (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_lepton_E_Wframe.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "top_mass",20, 100, 400, "1", "Top Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_top_mass.pdf");
  Compare_4Histos(tInput[0], tInput[17], tInput[16], tInput[33], "W_mass",20, 0, 140, "1", "W Mass (GeV)", "normalized", "legendUpRight", "dim6top", suffix[0], suffix[17], suffix[16], suffix[33], "results/compare_dim6top_ctwI_W_mass.pdf");
*/

  //-----------------Ratio EFT/SM-------------//

  //ctWI

  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "cosTheta","ctwi", 5, -1, 1,"1", "cos#theta", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosTheta.pdf");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "PhiStar","ctwi",20, 0, 6.2831 ,"nature_lepton == 1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","elec");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "PhiStar","ctwi",20, 0, 6.2831 ,"nature_lepton == 2", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_PhiStar.pdf","muon");
  //Ratio_EFT_SM(tInput[0], tInput[15], tInput[16], tInput[13], tInput[14], "cosThetaStar","ctwi", 3, -1, 1 ,"1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_ctwi_cosThetaStar.pdf");

/*  Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosTheta","cbwi", 5, -1, 1,"1", "cos#theta", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosTheta.pdf");
  Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "PhiStar","cbwi", 5, 0, 2*TMath::Pi(),"1", "#phi^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_PhiStar.pdf");*/
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 20, -1, 1 ,"nature_lepton==1", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "elec");
  //Ratio_EFT_SM(tInput[0], tInput[3], tInput[4], tInput[1], tInput[2], "cosThetaStar","cbwi", 20, -1, 1 ,"nature_lepton==2", "cos#theta^{*}", "ratio EFT/SM","results/ratio_madgraph/ratio_cbwi_cosThetaStar.pdf", "muon");
  /*Compare_3Histos(tInput[0], tInput[13], tInput[14], "cosTheta", 20, -1, 1, "1", "cos#theta", "normalized", "legendUpLeft", "dim6top", suffix[0], suffix[13], suffix[14], "results/madgraph_dim6top_ctWI_CosTheta_test.pdf");*/

  return 0;
}
