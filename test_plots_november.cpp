//This macro plots many different things that can be interesting to look at :)
#include "binning_histos.h"
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>

void plot_eec_distribution(){

// --- Input files --- CHANGE PATH IF NECESSARY
  TString backup_file = "/data_CMS/cms/zaidan/test_for_code_mods/backup/alldata_hist_3d_gen_aggr_n1_b_bjet_80_140.root"; // "/data_CMS/cms/zaidan/third_run/hist_3d_gen_aggr_n1_b_bjet_80_140.root";
    TString november_file ="/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_gen_aggr_n1_MC_bjet_80_140.root";    //"/data_CMS/cms/zaidan/test_for_code_mods/november/hist_3d_gen_aggr_n1_b_bjet_80_140.root";

// --- Histogram name (must exist in both files) ---
    TString hist_name = "h3D";  


// --- Open files ---
    TFile *f1 = TFile::Open(backup_file);    //I'll do it in general
    TFile *f2 = TFile::Open(november_file);

// --- Get histograms ---
    TH3D *h_3D_1 = (TH3D*) f1->Get(hist_name);
    TH3D *h_3D_2 = (TH3D*) f2->Get(hist_name+"_b1");

// --- Project in Y--- 
    TH1D *h_1D_1 = (TH1D*)h_3D_1->ProjectionY("h_1D_1", 1, bins_mb, 2,2);
    TH1D *h_1D_2 = (TH1D*)h_3D_2->ProjectionY("h_1D_2", 1, bins_mb, 2,2);

// --- Normalization ---
    h_1D_1->Scale(1.0 / h_1D_1->Integral(), "width");
    h_1D_2->Scale(1.0 / h_1D_2->Integral(), "width");

// --- Plot ---

    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_1D_1->SetStats(0);
    h_1D_1->SetTitle("EEC - comparing changed file vs backup");
    h_1D_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1D_1->GetXaxis()->CenterTitle(true);
    h_1D_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1D_1->GetYaxis()->CenterTitle(true);
    h_1D_1->GetYaxis()->SetRangeUser(0, 10.);
    
    h_1D_1->SetMarkerColor(kOrange+1);
    h_1D_1->SetLineColor(kOrange+1);
    h_1D_1->SetMarkerStyle(20);
    h_1D_1->Draw("P0 E HIST");


    h_1D_2->SetMarkerColor(kViolet);
    h_1D_2->SetLineColor(kViolet);
    h_1D_2->SetMarkerStyle(20);
    h_1D_2->Draw("P0 E HIST SAME");

    
// --- Text ---
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.65, 0.82, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.65, 0.77, "MC simulation");
    test_info_text->DrawLatex(0.65, 0.72, "single B jets");
    test_info_text->Draw("SAME");

// --- Legend ---
    TLegend *leg = new TLegend(0.15,0.7,0.25,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_1D_1, "backup");
    leg->AddEntry(h_1D_2, "modified");
    leg->Draw("SAME");

    c->Print("/grid_mnt/vol_home/llr/cms/zaidan/CMSSW_15_1_0/src/forZoe/Plotting/test_for_code_mods/eec_two_files_test.pdf");
    //c->Print("/grid_mnt/vol_home/llr/cms/zaidan/CMSSW_15_1_0/src/forZoe/Plotting/test_for_code_mods/eec_two_files_test.root");


// --- Clean ---
    //f1->Close();
    //f2->Close();
}







//______________________________________________________________________________________________
//____________________General effects on the EEC(dr) distribution_______________________________
//______________________________________________________________________________________________

//Effect of aggregation
void eec_aggr_effect(){
  //Select options
  bool show_other = true;
  bool show_gen = false;
  
  //Get your files
  TFile *file_aggr = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_aggr_BDT_n1_data_highEn_80_140.root";    //"/data_CMS/cms/meuli/testplots/wp09/hist_3d_b_bjet_80_140.root", "read");
  
  TFile *file_aggr_other_wp = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_noaggr_BDT_n1_data_highEn_80_140.root";    //"/data_CMS/cms/meuli/testplots/FN_new/hist_3d_b_bjet_80_140.root", "read");
  
  TFile *file_aggr_ideal = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_aggr_BDT_ideal_aggr_n1_data_highEn_80_140.root";    //"/data_CMS/cms/meuli/testplots/FN_ideal/hist_3d_idealaggr_b_bjet_80_140.root", "read");
  
  TFile *file_aggr_gen = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_gen_aggr_n1_MC_bjet_80_140.root";    //"/data_CMS/cms/meuli/testplots/wp09/hist_3d_gen_b_bjet_80_140.root", "read");
  
  TFile *file_noaggr = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_noaggr_n1_MC_bjet_80_140.root";    //"/data_CMS/cms/meuli/testplots/wp09/hist_3d_b_bjet_noaggr_80_140.root", "read");
    
  TFile *file_noaggr_gen = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_gen_aggr_n1_MC_bjet_80_140.root";    //"/data_CMS/cms/meuli/testplots/wp09/hist_3d_gen_b_bjet_noaggr_80_140.root", "read");
    
  //Get the histograms
  TH3D *h_3D_aggr = (TH3D*)file_aggr->Get("h3D");
  TH3D *h_3D_aggr_other_wp = (TH3D*)file_aggr_other_wp->Get("h3D");
  TH3D *h_3D_noaggr = (TH3D*)file_noaggr->Get("h3D");
  TH3D *h_3D_aggr_ideal = (TH3D*)file_aggr_ideal->Get("h3D");
  TH3D *h_3D_aggr_gen = (TH3D*)file_aggr_gen->Get("h3D");
  TH3D *h_3D_noaggr_gen = (TH3D*)file_noaggr_gen->Get("h3D");
  
  //Project into 1D EEC(dr) distributions
  TH1D *h_aggr = (TH1D*)h_3D_aggr->ProjectionY("h_aggr", 1, bins_mb, 2,2);
  TH1D *h_aggr_other_wp = (TH1D*)h_3D_aggr_other_wp->ProjectionY("h_aggr_other_wp", 1, bins_mb, 2,2);
  TH1D *h_aggr_ideal = (TH1D*)h_3D_aggr_ideal->ProjectionY("h_aggr_ideal", 1, bins_mb, 2,2);
  TH1D *h_aggr_gen = (TH1D*)h_3D_aggr_gen->ProjectionY("h_aggr_gen", 1, bins_mb, 2,2);
  TH1D *h_no_aggr = (TH1D*)h_3D_noaggr->ProjectionY("h_noaggr", 1, bins_mb, 2, 2);
  TH1D *h_no_aggr_gen = (TH1D*)h_3D_noaggr_gen->ProjectionY("h_noaggr_gen", 1, bins_mb, 2, 2);

  //Normalize
  h_aggr->Scale(1/h_aggr->Integral(), "width");
  h_aggr_other_wp->Scale(1/h_aggr_other_wp->Integral(), "width");
  h_aggr_gen->Scale(1/h_aggr_gen->Integral(), "width");
  h_aggr_ideal->Scale(1/h_aggr_ideal->Integral(), "width");
  h_no_aggr->Scale(1/h_no_aggr->Integral(), "width");
  h_no_aggr_gen->Scale(1/h_no_aggr_gen->Integral(), "width");
  
  //Plot
  TCanvas *c = new TCanvas("c", " ",170,800,700,504);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  c->SetLogx();
  
  h_aggr->SetStats(0);
  h_aggr->SetTitle("Effect of aggregation on EEC");
  h_aggr->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
  h_aggr->GetXaxis()->CenterTitle(true);
  h_aggr->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
  h_aggr->GetYaxis()->CenterTitle(true);
  h_aggr->GetYaxis()->SetRangeUser(0, 10.);
  h_aggr->SetMarkerColor(kOrange+1);
  h_aggr->SetLineColor(kOrange+1);
  h_aggr->SetMarkerStyle(20);
  h_aggr->Draw("P0 E HIST");
  
  if(show_other){
    h_aggr_other_wp->SetMarkerColor(kRed);
    h_aggr_other_wp->SetLineColor(kRed);
    h_aggr_other_wp->SetMarkerStyle(20);
    h_aggr_other_wp->Draw("P0 E HIST SAME");
  }
  
  h_no_aggr->SetMarkerColor(kViolet);
  h_no_aggr->SetLineColor(kViolet);
  h_no_aggr->SetMarkerStyle(20);
  h_no_aggr->Draw("P0 E HIST SAME");
  
    h_aggr_gen->SetMarkerColor(kAzure+4);
    h_aggr_gen->SetLineColor(kAzure+4);
    h_aggr_gen->SetMarkerStyle(20);
    if(show_gen) h_aggr_gen->Draw("P0 E HIST SAME");

    h_aggr_ideal->SetMarkerColor(kAzure+1);
    h_aggr_ideal->SetLineColor(kAzure+1);
    h_aggr_ideal->SetMarkerStyle(20);
    h_aggr_ideal->Draw("P0 E HIST SAME");

    
    h_no_aggr_gen->SetMarkerColor(kGreen+3);
    h_no_aggr_gen->SetLineColor(kGreen+3);
    h_no_aggr_gen->SetMarkerStyle(20);
    //h_no_aggr_gen->Draw("P0 E HIST SAME");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.65, 0.82, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.65, 0.77, "MC simulation");
    test_info_text->DrawLatex(0.65, 0.72, "single B jets");
    if(!show_other) test_info_text->DrawLatex(0.65, 0.67, "BDT score > -0.9");
    //test_info_text->DrawLatex(0.15, 0.67, "trkMatchSta >= 100");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.15,0.7,0.25,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    if(show_other){
        leg->AddEntry(h_aggr, "reco aggregated (BDT > -0.9)");
        leg->AddEntry(h_aggr_other_wp, "reco aggregated (FNscore > 0.3)");
    }
    else leg->AddEntry(h_aggr, "reco aggregated");
    leg->AddEntry(h_no_aggr, "reco not aggregated");
    leg->AddEntry(h_aggr_ideal, "reco ideal aggregated");
    if(show_gen) leg->AddEntry(h_aggr_gen, "gen aggregated");
    //leg->AddEntry(h_no_aggr_gen, "gen not aggregated");
    leg->Draw("same");

    c->Print("/data_CMS/cms/meuli/testplots/FN_new/eec_aggregation_effect_FNscore_forpp.pdf");

}

//Effect of applying the b-tagger
void eec_tag_effect(){

    //Get files
    TFile *file_tag = new TFile("/data_CMS/cms/zaidan/test_for_code_mods/run_with_mod_code/hist_3d_aggr_BDT_n1_MC_bjet_80_140.root", "read");
    
    TFile *file_notag = new TFile("/data_CMS/cms/zaidan/second_run/hist_3d_gen_aggr_n1_b_notag_notag_notag_notag_notag_dijet_80_140.root", "read");

    //Get histograms
    TH3D *h_3D_tag = (TH3D*)file_tag->Get("h3D");
    TH3D *h_3D_notag = (TH3D*)file_notag->Get("h3D");

    //Project to 1D EEC(dr)
    TH1D *h_tag = (TH1D*)h_3D_tag->ProjectionY("h_tag", 1, bins_mb, 2,2);
    TH1D *h_no_tag = (TH1D*)h_3D_notag->ProjectionY("h_no_tag", 1, bins_mb, 2, 2);

    //Normalize
    h_tag->Scale(1/h_tag->Integral(), "width");
    h_no_tag->Scale(1/h_no_tag->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_tag->SetStats(0);
    h_tag->SetTitle("Effect of b-tagging on EEC");
    h_tag->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_tag->GetXaxis()->CenterTitle(true);
    h_tag->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_tag->GetYaxis()->CenterTitle(true);
    h_tag->GetYaxis()->SetRangeUser(0, 10.);
    h_tag->SetMarkerColor(kBlue);
    h_tag->SetLineColor(kBlue);
    h_tag->SetMarkerStyle(20);
    h_tag->Draw("P0 E HIST");

    h_no_tag->SetMarkerColor(kGreen);
    h_no_tag->SetLineColor(kGreen);
    h_no_tag->SetMarkerStyle(20);
    h_no_tag->Draw("P0 E HIST SAME");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "gen MC");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.7,0.7,0.85,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h_tag, "tagged");
    leg->AddEntry(h_no_tag, "not tagged");
    leg->Draw("same");

    c->Print("/data_CMS/cms/zaidan/second_run/eec_tagging_effect.pdf");

}

//Effect of using a different energy weight exponent n
void eec_n_effect(){

    //Get your files
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_gen_b_bjet_80_140.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_gen_b_bjet_n2_80_140.root", "read");

    //Get histograms
    TH3D *h_3D_1 = (TH3D*)file_1->Get("h3D");
    TH3D *h_3D_2 = (TH3D*)file_2->Get("h3D");

    //Project to 1D EEC(dr) distribution
    TH1D *h_1 = (TH1D*)h_3D_1->ProjectionY("h_1", 1, bins_mb, 2,2);
    TH1D *h_2 = (TH1D*)h_3D_2->ProjectionY("h_2", 1, bins_mb, 2, 2);

    //Normalize
    h_1->Scale(1/h_1->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_1->SetStats(0);
    h_1->SetTitle("");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetRangeUser(0, 11.);
    h_1->SetMarkerColor(kMagenta);
    h_1->SetLineColor(kMagenta);
    h_1->SetMarkerStyle(20);
    h_1->Draw("P0 E HIST");

    h_2->SetMarkerColor(kBlue);
    h_2->SetLineColor(kBlue);
    h_2->SetMarkerStyle(20);
    h_2->Draw("P0 E HIST SAME");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.16, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.16, 0.75, "gen MC");
    test_info_text->DrawLatex(0.16, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.65,0.7,0.8,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h_1, "n = 1");
    leg->AddEntry(h_2, "n = 2");
    leg->Draw("same");

    //c->Print("/data_CMS/cms/meuli/wp095/biggerbins/eec_n_effect.png");
    c->Print("/data_CMS/cms/meuli/wp095/biggerbins/eec_n_effect.pdf");

}

//Ratio of 3-point EEC to 2-point EEC
void eec_ratio(){

    //Get files
    TFile *file_e2c = new TFile("/data_CMS/cms/meuli/testplots/wp095/hist_3d_gen_b_bjet_80_140.root", "read");
    
    TFile *file_e3c = new TFile("/data_CMS/cms/meuli/testplots/wp095/hist_3d_gen_b_bjet_80_140_e3c.root", "read");

    //Get histograms
    TH3D *h_e2c_3D = (TH3D*)file_e2c->Get("h3D");
    TH3D *h_e3c_3D = (TH3D*)file_e3c->Get("h3D");

    //Project to 1D EEC(dr) distributions
    TH1D *h_e2c = (TH1D*)h_e2c_3D->ProjectionY("h_e2c", 1, bins_mb, 2, 2);
    TH1D *h_e3c = (TH1D*)h_e3c_3D->ProjectionY("h_e3c", 1, bins_mb, 2, 2);

    TH1D *h_ratio = (TH1D*)h_e3c->Clone("h_ratio");
    h_ratio->Divide(h_e2c);

    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_ratio->SetStats(0);
    h_ratio->SetTitle("Ratio of EECs");
    h_ratio->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_ratio->GetXaxis()->CenterTitle(true);
    h_ratio->GetYaxis()->SetTitle("\\frac{E3C(\\Delta\\mbox{r})}{E2C(\\Delta\\mbox{r})}");
    h_ratio->GetYaxis()->CenterTitle(true);
    h_ratio->GetXaxis()->SetRangeUser(0.,0.2);
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->SetLineColor(kBlue);
    h_ratio->SetMarkerStyle(20);
    h_ratio->Draw("P0 E HIST");

    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.2, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.2, 0.75, "gen MC");
    test_info_text->DrawLatex(0.2, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");




    c->Print("/data_CMS/cms/meuli/testplots/wp095/eec_ratio.pdf");

}

//Effect of the jet pT on the EEC(dr) distribution
void pT_effect(){

    //Get file
  TFile *file = new TFile("/data_CMS/cms/zaidan/third_run/hist_3d_gen_aggr_n1_b_dijet_80_140.root", "read");    //"/data_CMS/cms/meuli/testplots/wp095/hist_3d_gen_b_bjet_80_140.root", "read");

    //Get histogram
    TH3D *h_3D = (TH3D*)file->Get("h3D");

    //Select the different pT bins and project to 1D EEC(dr) distributions
    TH1D *h_1 = (TH1D*)h_3D->ProjectionY("h_1", 1, bins_mb, 1,1);
    TH1D *h_2 = (TH1D*)h_3D->ProjectionY("h_2", 1, bins_mb, 2,2);
    TH1D *h_3 = (TH1D*)h_3D->ProjectionY("h_3", 1, bins_mb, 3,3);

    //Normalize
    h_1->Scale(1/h_1->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");
    h_3->Scale(1/h_3->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_1->SetStats(0);
    h_1->SetTitle("Effect of jet p_{T} on EEC");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetRangeUser(0, 10.);
    h_1->GetXaxis()->SetRange(2, bins_dr-1);
    h_1->SetMarkerColor(kBlue);
    h_1->SetLineColor(kBlue);
    h_1->SetMarkerStyle(20);
    h_1->Draw("P0 E HIST");

    h_2->GetXaxis()->SetRange(2, bins_dr-1);
    h_2->SetMarkerColor(kGreen);
    h_2->SetLineColor(kGreen);
    h_2->SetMarkerStyle(20);
    h_2->Draw("P0 E HIST SAME");

    
    h_3->GetXaxis()->SetRange(2, bins_dr-1);
    h_3->SetMarkerColor(kRed);
    h_3->SetLineColor(kRed);
    h_3->SetMarkerStyle(20);
    h_3->Draw("P0 E HIST SAME");
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.75, "gen MC");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.6,0.7,0.8,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);
    leg->AddEntry(h_1, "80 < p_{T} < 100 GeV");
    leg->AddEntry(h_2, "100 < p_{T} < 120 GeV");
    leg->AddEntry(h_3, "120 < p_{T} < 140 GeV");
    leg->Draw("same");

    c->Print("/data_CMS/cms/zaidan/third_run/eec_pt_effect_dijet.pdf");

}

//Effect of jet flavour on the EEC(dr) distribution
void flavour_effect(){

    //Get files
    TFile *file_b = new TFile("/data_CMS/cms/zaidan/third_run/flavour/hist_3d_gen_aggr_n1_b_notag_notag_dijet_80_140.root", "read");
    
    TFile *file_c = new TFile("/data_CMS/cms/zaidan/third_run/flavour/hist_3d_gen_aggr_n1_c_notag_notag_dijet_80_140.root", "read");
    
    TFile *file_other = new TFile("/data_CMS/cms/zaidan/third_run/flavour/hist_3d_gen_aggr_n1_light_notag_notag_dijet_80_140.root", "read");

    //Get histograms
    TH3D *h_3D_b = (TH3D*)file_b->Get("h3D");
    TH3D *h_3D_c = (TH3D*)file_c->Get("h3D");
    TH3D *h_3D_other = (TH3D*)file_other->Get("h3D");

    //If needed, also plot the inclusive
    bool inclusive = false;
    TH3D *h_3D_inclusive = (TH3D*)h_3D_other->Clone("h_3D_inclusive");

    TH1D *h_b = (TH1D*)h_3D_b->ProjectionY("h_b", 1, bins_mb, 2,2);
    TH1D *h_c = (TH1D*)h_3D_c->ProjectionY("h_c", 1, bins_mb, 2, 2);
    TH1D *h_other = (TH1D*)h_3D_other->ProjectionY("h_other", 1, bins_mb, 2, 2);
    TH1D *h_inclusive;
    
    if(inclusive){
        h_3D_inclusive->Add(h_3D_b);
        h_3D_inclusive->Add(h_3D_c);
        h_inclusive = (TH1D*)h_3D_inclusive->ProjectionY("h_inclusive", 1, bins_mb, 2, 2);
        h_inclusive->Scale(1/h_inclusive->Integral(), "width");
    }

    //Normalize
    h_b->Scale(1/h_b->Integral(), "width");
    h_c->Scale(1/h_c->Integral(), "width");
    h_other->Scale(1/h_other->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_b->SetStats(0);
    h_b->SetTitle("Effect of jet flavour on EEC");
    h_b->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_b->GetXaxis()->CenterTitle(true);
    h_b->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_b->GetYaxis()->CenterTitle(true);
    h_b->GetYaxis()->SetRangeUser(0., 12.);
    h_b->SetMarkerColor(kBlue);
    h_b->SetLineColor(kBlue);
    h_b->SetMarkerStyle(20);
    h_b->Draw("P0 E HIST");
    
    h_c->SetMarkerColor(kRed);
    h_c->SetLineColor(kRed);
    h_c->SetMarkerStyle(20);
    h_c->Draw("P0 E HIST SAME");

    h_other->SetMarkerColor(kGreen);
    h_other->SetLineColor(kGreen);
    h_other->SetMarkerStyle(20);
    h_other->Draw("P0 E HIST SAME");

    if(inclusive){
        h_inclusive->SetMarkerColor(kMagenta);
        h_inclusive->SetLineColor(kMagenta);
        h_inclusive->SetMarkerStyle(24);
        h_inclusive->Draw("P0 E HIST SAME");
    }

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "gen MC");
    test_info_text->DrawLatex(0.15, 0.7, "aggregation for b-jets");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.7,0.7,0.85,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h_b, "b jets");
    leg->AddEntry(h_c, "c jets");
    leg->AddEntry(h_other, "light jets");
    if(inclusive) leg->AddEntry(h_inclusive, "inclusive jets");
    leg->Draw("same");

    if(!inclusive) c->Print("/data_CMS/cms/zaidan/third_run/flavour/eec_flavour_effect.pdf");
    else c->Print("/data_CMS/cms/zaidan/third_run/flavour/eec_flavour_effect_inclusive.pdf");

}

//Comparison of EEC(dr) distribution for aggregated single-b jets and double-b jets
void eec_b_vs_moreb(){

    //Get files
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_b_bjet_80_140.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_moreb_bjet_80_140.root", "read");

    //Get histograms
    TH3D *h_3D_1 = (TH3D*)file_1->Get("h3D");
    TH3D *h_3D_2 = (TH3D*)file_2->Get("h3D");

    //Project into 1D EEC(dr) distribution
    TH1D *h_1 = (TH1D*)h_3D_1->ProjectionY("h_1", 1, bins_mb, 2,2);
    TH1D *h_2 = (TH1D*)h_3D_2->ProjectionY("h_2", 1, bins_mb, 2, 2);

    //Normalize
    h_1->Scale(1/h_1->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_1->SetStats(0);
    h_1->SetTitle("");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetRangeUser(0, 10.);
    h_1->SetMarkerColor(kBlue);
    h_1->SetLineColor(kBlue);
    h_1->SetMarkerStyle(20);
    h_1->Draw("P0 E HIST");

    h_2->SetMarkerColor(kRed);
    h_2->SetLineColor(kRed);
    h_2->SetMarkerStyle(20);
    h_2->Draw("P0 E HIST SAME");

    
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.17, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.17, 0.75, "reco MC");
    test_info_text->DrawLatex(0.17, 0.7, "aggregated jets");
    test_info_text->DrawLatex(0.17, 0.65, "BDT score > -0.95");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);
    leg->AddEntry(h_1, "single B");
    leg->AddEntry(h_2, "more B");
    leg->Draw("same");


    c->Print("/data_CMS/cms/meuli/wp095/biggerbins/compare_B.pdf");

}

//______________________________________________________________________________________________
//____________________Stacked distributions_____________________________________________________
//______________________________________________________________________________________________

//Makes a bigger overflow bin for the mB axis
void rebin_mass(TH3D* &h){
    Int_t bins_mb = h->GetNbinsX();
    Int_t bins_pt = h->GetNbinsZ();
    Int_t bins_y = h->GetNbinsY();

    Int_t last_mb_bin = 12;

    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_y = 1; ibin_y <= bins_y; ibin_y++){
            Float_t last_bin = 0;
            for(Int_t imb_bin = last_mb_bin; imb_bin <= bins_mb; imb_bin++){
                last_bin += h->GetBinContent(imb_bin, ibin_y, ibin_pt);
            }
            h->GetXaxis()->SetRange(1, last_mb_bin);
            h->SetBinContent(last_mb_bin, ibin_y, ibin_pt);
        }
    }
    //h->ProjectionX("h_proj", 1, bins_y, 2,2)->Draw("E0 HIST"); 
}

//Makes a bigger overflow bin for the eec weight axis
void rebin_eec(THnD* &h){

    Int_t last_bin = 5;

    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){
            for(Int_t imb_bin = 1; imb_bin <= bins_mb; imb_bin++){
                Float_t last_bin_content = 0;
                Int_t bin[4] = {imb_bin, ibin_dr, last_bin, ibin_pt};
                for(Int_t ibin_eec = last_bin; ibin_eec <= bins_eec; ibin_eec++){
                    Int_t bin_index[4] = {imb_bin, ibin_dr, ibin_eec, ibin_pt};
                    last_bin_content += h->GetBinContent(bin_index);
                }
                h->GetAxis(eec_dim)->SetRange(1, last_bin);
                h->SetBinContent(bin, last_bin_content);
            
            }
        }
    //h->ProjectionX("h_proj", 1, bins_y, 2,2)->Draw("E0 HIST"); 
    }
}

//Draws the mB stacked histogram
void draw_hist_stack(TString &dataset, TString &pT_selection){
    //Create canvas
    TCanvas *c = new TCanvas("c_stack", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);

    //Create a stack plot
    THStack *st = new THStack("st", "Mass stacked histogram (" + dataset + " sample)");

    //Get histograms
    TFile *file = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_b_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h = (TH3D*)file->Get("h3D");

    TFile *file2 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_moreb_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h2 = (TH3D*)file2->Get("h3D");

    TFile *file3 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_mc_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h3 = (TH3D*)file3->Get("h3D");

    TFile *file4 = new TFile("/data_CMS/cms/meuli/wp095/biggerbins/hist_3d_other_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h4 = (TH3D*)file4->Get("h3D");

    //Rebin with a bigger overflow bin if needed
    /*for (auto hi : {h, h2, h3, h4 }) {
         rebin_mass(hi);
         hi->Sumw2();
        }*/

    //Project into 1D EEC(dr) distributions
    bins_dr = h->GetNbinsX();
    TH1D *hx = h->ProjectionX("hx", 1, bins_dr, 2, 2);
    TH1D *h2x = h2->ProjectionX("h2x", 1, bins_dr, 2, 2);
    TH1D *h3x = h3->ProjectionX("h3x", 1, bins_dr, 2, 2);
    TH1D *h4x = h4->ProjectionX("h4x", 1, bins_dr, 2, 2);

    //Add the stacked histograms
    hx->SetTitle("mB distribution");
    hx->SetStats(0);
    hx->SetLineColor(kBlue);
    hx->SetFillColor(kBlue);
    hx->SetMarkerColor(kBlue);
    
    h2x->SetLineColor(kRed);
    h2x->SetFillColor(kRed);
    h2x->SetMarkerColor(kRed);
    st->Add(h2x);
    st->Add(hx);

    h3x->SetLineColor(kBlack);
    h3x->SetMarkerStyle(22);
    h3x->SetMarkerColor(kBlack);
    h3x->SetMarkerSize(1);
    h3x->SetFillColor(0);
    h3x->GetXaxis()->SetTitle("m_{B}");
    h3x->GetXaxis()->CenterTitle(true);
    //h3x->GetYaxis()->SetRangeUser(0, 8000);

    h4x->SetLineColor(kGreen);
    h4x->SetFillColor(kGreen);
    h4x->SetMarkerColor(kGreen);
    st->Add(h4x);

    //Plot
    st->Draw("HP");
    h3x->Draw("HP, same");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.65, 0.5, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.65, 0.45, "reco MC");
    test_info_text->DrawLatex(0.65, 0.4, "aggregated jets");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.65,0.65,0.85,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h3x, "full MC");
    leg->AddEntry(hx, "1 B");
    leg->AddEntry(h2x, "more B");
    leg->AddEntry(h4x, "other");
    leg->Draw("same");

    

    //Save
    c->Print("/data_CMS/cms/meuli/wp095/biggerbins/mB_distr_" + dataset + "_" + pT_selection + "_stack.pdf");
}

//Draw the EEC(dr) stacked histogram
void draw_hist_stack_eec(TString &dataset, TString &pT_selection){
    //Create canvas
    TCanvas *c = new TCanvas("c_stack", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);


    //Create a stack plot
    THStack *st = new THStack("st", "EEC stacked histogram (" + dataset + " sample)");

    //Get histograms
    TFile *file = new TFile("hist_3d_b_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h = (TH3D*)file->Get("h3D");
    TH1D *hx = h->ProjectionY("hx", 1, mb_bins, 2, 2);

    TFile *file2 = new TFile("hist_3d_moreb_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h2 = (TH3D*)file2->Get("h3D");
    TH1D *h2x = h2->ProjectionY("h2x", 1, mb_bins, 2, 2);

    TFile *file3 = new TFile("hist_3d_mc_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h3 = (TH3D*)file3->Get("h3D");
    TH1D *h3x = h3->ProjectionY("h3x", 1, mb_bins, 2, 2);

    TFile *file4 = new TFile("hist_3d_other_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h4 = (TH3D*)file4->Get("h3D");
    TH1D *h4x = h4->ProjectionY("h4x", 1, mb_bins, 2, 2);

    //Add histograms
    hx->SetTitle("EEC distribution");
    hx->SetStats(0);
    hx->SetLineColor(kBlue);
    hx->SetFillColor(kBlue);
    hx->SetMarkerColor(kBlue);
    
    h2x->SetLineColor(kRed);
    h2x->SetFillColor(kRed);
    h2x->SetMarkerColor(kRed);
    st->Add(h2x);
    st->Add(hx);

    h3x->SetLineColor(kBlack);
    h3x->SetMarkerStyle(22);
    h3x->SetMarkerColor(kBlack);
    h3x->SetMarkerSize(1);
    h3x->SetFillColor(0);
    //h3x->GetYaxis()->SetRangeUser(0, 8000);

    h4x->SetLineColor(kGreen);
    h4x->SetFillColor(kGreen);
    h4x->SetMarkerColor(kGreen);
    st->Add(h4x);


    //Plot
    st->Draw("HP");
    h3x->Draw("HP, same");

    //Legend
    TLegend *leg = new TLegend(0.6,0.6,0.8,0.8, "Legend"); 
    leg->AddEntry(h3x, "full MC");
    leg->AddEntry(hx, "1 B");
    leg->AddEntry(h2x, "more B");
    leg->AddEntry(h4x, "other");
    leg->Draw("same");

    

    //Save
    c->Print("eec_distr_" + dataset + "_" + pT_selection + "_stack.png");
}

//Draw the pT stacked histogram
void draw_hist_stack_pt(TString &dataset, TString &pT_selection){
    //Create canvas
    TCanvas *c = new TCanvas("c_stack", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);

    //Create a stack plot
    THStack *st = new THStack("st", "pT stacked histogram (" + dataset + " sample)");

    //Get histograms
    TFile *file = new TFile("hist_3d_b_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h = (TH3D*)file->Get("h3D");
    TH1D *hx = h->ProjectionZ("hx", 1, mb_bins, 1, bins_dr);

    TFile *file2 = new TFile("hist_3d_moreb_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h2 = (TH3D*)file2->Get("h3D");
    TH1D *h2x = h2->ProjectionZ("h2x", 1, mb_bins, 1, bins_dr);

    TFile *file3 = new TFile("hist_3d_mc_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h3 = (TH3D*)file3->Get("h3D");
    TH1D *h3x = h3->ProjectionZ("h3x", 1, mb_bins, 1, bins_dr);

    TFile *file4 = new TFile("hist_3d_other_" + dataset + "_" + pT_selection + ".root", "read");
    TH3D *h4 = (TH3D*)file4->Get("h3D");
    TH1D *h4x = h4->ProjectionZ("h4x", 1, mb_bins,1, bins_dr);

    //Add histograms
    hx->SetTitle(" pT distribution");
    hx->SetStats(0);
    hx->SetLineColor(kBlue);
    hx->SetFillColor(kBlue);
    hx->SetMarkerColor(kBlue);
    
    h2x->SetLineColor(kRed);
    h2x->SetFillColor(kRed);
    h2x->SetMarkerColor(kRed);
    st->Add(h2x);
    st->Add(hx);

    h3x->SetLineColor(kBlack);
    h3x->SetMarkerStyle(22);
    h3x->SetMarkerColor(kBlack);
    h3x->SetMarkerSize(1);
    h3x->SetFillColor(0);
    //h3x->GetYaxis()->SetRangeUser(0, 8000);

    h4x->SetLineColor(kGreen);
    h4x->SetFillColor(kGreen);
    h4x->SetMarkerColor(kGreen);
    st->Add(h4x);


    //Plot
    st->Draw("HP");
    h3x->Draw("HP, same");

    //Legend
    TLegend *leg = new TLegend(0.6,0.6,0.8,0.8, "Legend"); 
    leg->AddEntry(h3x, "full MC");
    leg->AddEntry(hx, "1 B");
    leg->AddEntry(h2x, "more B");
    leg->AddEntry(h4x, "other");
    leg->Draw("same");

    

    //Save
    c->Print("pt_distr_" + dataset + "_" + pT_selection + "_stack.png");
}

//______________________________________________________________________________________________
//____________________EEC simple distributions and filling methods______________________________
//______________________________________________________________________________________________

//Approximate a 3D distribution into the usuale 1D EEC(dr) distribution
void project_eec_1D(TH3D* &h, TH1D* &h1D, Int_t &pt_bin, const Double_t eec_last_bin = 1000){
    h->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D = (TH2D*)h->Project3D("yx");
    Double_t dr, eec;
    for(Int_t ibin_dr=1; ibin_dr <= bins_dr; ibin_dr++){
        dr = h_2D->GetXaxis()->GetBinLowEdge(ibin_dr)+dr_shiftbin;
        Float_t bin_err = 0;
        for(Int_t ibin_eec=1; ibin_eec <= bins_eec; ibin_eec++){
            eec = h_2D->GetYaxis()->GetBinLowEdge(ibin_eec)+eec_step/2;
            if(ibin_eec == bins_eec) eec = eec_last_bin;
            h1D->Fill(dr, h_2D->GetBinContent(ibin_dr, ibin_eec)*eec);
            bin_err += pow(h_2D->GetBinError(ibin_dr, ibin_eec)*eec, 2);

        }
        h1D->SetBinError(ibin_dr, sqrt(bin_err));
    }
}

//Draw a simple (approximated) EEC(dr) distribution from a 4D one in (mB, dr, energy weight, jtpt)
void draw_simple_eec(){
    //Select the pT bin
    Int_t pt_bin = 2;

    //Get File
    TFile *file = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins/hist_4d_data_HighEGJet_80_140.root", "read");

    //Get 4D histograms and set range
    THnD *h_4D = (THnD*)file->Get("h4D");    
    h_4D->GetAxis(mb_dim)->SetRange(1, bins_mb);

    //Print one random bin content (for checks)
    Int_t indx[4] = {5,2,2,2};
    std::cout << h_4D->GetBinContent(indx) << std::endl;

    //Project from 4D to 1D
    TH3D *h_3D = (TH3D*)h_4D->Projection(dr_dim, eec_dim, pt_dim, "E")->Clone("h_3D");
    h_3D->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D = (TH2D*)h_3D->Project3D("yx")->Clone("h_2D");

    TH1D* h_1 = new TH1D("h_1", "h_1", bins_dr, dr_binsVector);
    
    h_1->Sumw2();

    recover_eec_distr(h_1, h_2D);

    //Normalize
    h_1->Scale(1/h_1->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_1->SetStats(0);
    h_1->SetTitle("EEC distribution from 4D histograms");
    h_1->GetXaxis()->SetTitle("\\Delta r");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta r)");
    h_1->GetYaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetRangeUser(0, 10.);
    h_1->GetXaxis()->SetRange(2, bins_dr-1);
    h_1->SetMarkerColor(kBlue);
    h_1->SetLineColor(kBlue);
    h_1->SetMarkerStyle(20);
    h_1->Draw("P0 E HIST");

    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "reco MC");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single b jets");
    test_info_text->Draw("SAME");


    c->Print("/data_CMS/cms/meuli/wp09_4d/eec_simple.png");

}

//Draw an approximated EEC(dr) distribution from a 3D one in (dr, energy weight, jtpt) after unfolding
void draw_eec_1d(){
    //Select the pT bin
    Int_t pt_bin = 2;

    //Get the file
    TFile *file = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D.root", "read");

    //Get the histogram
    TH3D *h_3D = (TH3D*)file->Get("h_data_unfolded"); 

    //Approximate the 1D EEC(dr) distribution
    TH1D* h_1 = new TH1D("h_1", "h_1", bins_dr, dr_binsVector);
    
    h_1->Sumw2();

    project_eec_1D(h_3D, h_1, pt_bin);

    //Normalize
    h_1->Scale(1/h_1->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_1->SetStats(0);
    h_1->SetTitle("");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    //h_1->GetYaxis()->SetRangeUser(0, 10.);
    h_1->SetMarkerColor(kBlue);
    h_1->SetLineColor(kBlue);
    h_1->SetMarkerStyle(20);
    h_1->Draw("P0 E HIST");

    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.6, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.6, 0.75, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.6, 0.7, "unfolded single b jets");
    test_info_text->Draw("SAME");


    c->Print("/data_CMS/cms/meuli/wp09_4d/eec_1d_signal.pdf");


}

//Draw a 2D histogram in dr and energy weight from a 3D one after unfolding
void draw_eec_2d(){
    //Select pT bin
    Int_t pt_bin = 2;

    //Get file
    TFile *file = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D.root", "read");

    //Get and project histogram
    TH3D *h_3D = (TH3D*)file->Get("h_data_unfolded"); 

    h_3D->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D* h_1 = (TH2D*)h_3D->Project3D("yx")->Clone("h_1");
    
    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();
    c->SetLogz();

    h_1->SetStats(0);
    h_1->SetTitle("Unfolded 2D data distribution for 100 < p_{T} < 120 GeV");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("(p_{T,i}p_{T,j})");
    h_1->GetYaxis()->CenterTitle(true);
    //h_1->GetYaxis()->SetRangeUser(0, 10.);
    h_1->Draw("colz");

    /*TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.15, 0.7, "unfolded single b jets");
    test_info_text->Draw("SAME");*/


    c->Print("/data_CMS/cms/meuli/wp09_4d/eec_2d_signal.pdf");


}

//Compare direct EEC(dr) filling with the approximation from a 4D histogram
void eec_4d_vs_direct_filling(){
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_gen_b_bjet_super_big_eec_80_140.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_b_bjet_80_140.root", "read");

    Int_t pt_bin = 2;

    THnD *h_4D = (THnD*)file_1->Get("h4D");
    THnD *h_4D_simple = (THnD*)h_4D->Clone("h_4D_simple");
    TH3D *h_3D_2 = (TH3D*)file_2->Get("h3D");
    h_4D->GetAxis(mb_dim)->SetRange(1, bins_mb);
    h_4D->GetAxis(pt_dim)->SetRange(pt_bin, pt_bin);
    //rebin_eec(h_4D);
    h_4D_simple->GetAxis(mb_dim)->SetRange(1, bins_mb);
    h_4D_simple->GetAxis(pt_dim)->SetRange(pt_bin, pt_bin);

    TH2D *h_2D = (TH2D*)h_4D->Projection(eec_dim, dr_dim)->Clone("h_2D");
    TH2D *h_2D_simple = (TH2D*)h_4D_simple->Projection(eec_dim, dr_dim)->Clone("h_2D_simple");
    TH2D *h_2D_750 = (TH2D*)h_2D->Clone("h_2D_750");


    TH1D *h_1 = new TH1D("h_1", "h_1", bins_dr, dr_binsVector);
    recover_eec_distr(h_1, h_2D, 500);
    TH1D *h_1_simple = new TH1D("h_1_simple", "h_1_simple", bins_dr, dr_binsVector);
    recover_eec_distr(h_1_simple, h_2D_simple, 1000);
    TH1D *h_1_750 = new TH1D("h_1_750", "h_1_750", bins_dr, dr_binsVector);
    recover_eec_distr(h_1_750, h_2D_750, 750);
    TH1D *h_2 = (TH1D*)h_3D_2->ProjectionY("h_2", 1, bins_mb, 2, 2);

    h_1->Scale(1/h_1->Integral(), "width");
    h_1_simple->Scale(1/h_1_simple->Integral(), "width");
    h_1_750->Scale(1/h_1_750->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");

    TCanvas *c = new TCanvas("c_plot", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();


    h_1->SetStats(0);
    h_1->SetTitle("EEC for 100 < pT < 120 GeV, single b jets, dijet reco MC");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    //h_1->GetYaxis()->SetRangeUser(0, 15.);
    h_1->SetMarkerColor(kMagenta);
    h_1->SetLineColor(kMagenta);
    h_1->SetMarkerStyle(20);
    h_1->Draw("PE1 HIST");

    h_1_simple->SetMarkerColor(kGreen);
    h_1_simple->SetLineColor(kGreen);
    h_1_simple->SetMarkerStyle(20);
    h_1_simple->Draw("PE1 HIST SAME");

    h_1_750->SetMarkerColor(kOrange);
    h_1_750->SetLineColor(kOrange);
    h_1_750->SetMarkerStyle(20);
    h_1_750->Draw("PE1 HIST SAME");

    h_2->SetMarkerColor(kBlue);
    h_2->SetLineColor(kBlue);
    h_2->SetMarkerStyle(20);
    h_2->Draw("PE1 HIST SAME");

    //Legend
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9, ""); 
    leg->AddEntry(h_1, "from 4D with 500 overflow bin, filled with 500");
    leg->AddEntry(h_1_750, "from 4D with 500 overflow bin, filled with 750");
    leg->AddEntry(h_1_simple, "from 4D with 1000 overflow bin");
    leg->AddEntry(h_2, "directly");
    leg->Draw("same");


    c->Draw();


    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/compare_filling_methods.pdf");

}

//Compare the same as eec_4d_vs_direct_filling but with different energy weights of the overflow bin
void compare_filling_methods(){
    //Get data samples (before template fit and unfolding)
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_data_HighEGJet_80_140.root", "read");
    TFile *file_1_low = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_data_lowEGJet_80_140.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_data_HighEGJet_80_140.root", "read");
    TFile *file_2_low = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_data_lowEGJet_80_140.root", "read");

    //Select pT bin
    Int_t pt_bin = 2;

    //Get and project/approximate EEC(dr) histograms
    THnD *h_4D = (THnD*)file_1->Get("h4D");
    THnD *h_4D_low = (THnD*)file_1_low->Get("h4D");
    h_4D->Add(h_4D_low);

    TH3D *h_3D_2 = (TH3D*)file_2->Get("h3D");
    TH3D *h_3D_2_low = (TH3D*)file_2_low->Get("h3D");
    h_3D_2->Add(h_3D_2_low);

    h_4D->GetAxis(mb_dim)->SetRange(1, bins_mb);
    TH3D *h_3D = (TH3D*)h_4D->Projection(dr_dim, eec_dim, pt_dim)->Clone("h_3D");
    h_4D->GetAxis(pt_dim)->SetRange(pt_bin, pt_bin);

    TH2D *h_2D = (TH2D*)h_4D->Projection(eec_dim, dr_dim)->Clone("h_2D");
    
    //Approximate from 4D with different values of the last bin energy weighting
    TH1D *h_1 = new TH1D("h_1", "h_1", bins_dr, dr_binsVector);
    TH1D *h_1_750 = new TH1D("h_1_750", "h_1_750", bins_dr, dr_binsVector);
    TH1D *h_1_1000 = new TH1D("h_1_1000", "h_1_1000", bins_dr, dr_binsVector);
    TH1D *h_1_1500 = new TH1D("h_1_1500", "h_1_1500", bins_dr, dr_binsVector);

    project_eec_1D(h_3D, h_1, pt_bin, 450);
    project_eec_1D(h_3D, h_1_750, pt_bin, 750);
    project_eec_1D(h_3D, h_1_1000, pt_bin, 1000);
    project_eec_1D(h_3D, h_1_1500, pt_bin, 1500);

    TH1D *h_2 = (TH1D*)h_3D_2->ProjectionY("h_2", 1, bins_mb, 2, 2);


    h_1->Scale(1/h_1->Integral(), "width");
    h_1_750->Scale(1/h_1_750->Integral(), "width");
    h_1_1000->Scale(1/h_1_1000->Integral(), "width");
    h_1_1500->Scale(1/h_1_1500->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");

    TCanvas *c = new TCanvas("c_plot", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();


    h_1->SetStats(0);
    h_1->SetTitle("");//EEC for 100 < pT < 120 GeV, single b jets, dijet reco MC");
    h_1->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetRangeUser(0, 15.);
    h_1->SetMarkerColor(kMagenta);
    h_1->SetLineColor(kMagenta);
    h_1->SetMarkerStyle(20);
    h_1->Draw("PE1 HIST");

    h_1_750->SetMarkerColor(kRed);
    h_1_750->SetLineColor(kRed);
    h_1_750->SetMarkerStyle(20);
    h_1_750->Draw("PE1 HIST SAME");
    
    h_1_1000->SetMarkerColor(kGreen);
    h_1_1000->SetLineColor(kGreen);
    h_1_1000->SetMarkerStyle(20);
    h_1_1000->Draw("PE1 HIST SAME");
    
    h_1_1500->SetMarkerColor(kOrange);
    h_1_1500->SetLineColor(kOrange);
    h_1_1500->SetMarkerStyle(20);
    h_1_1500->Draw("PE1 HIST SAME");

    h_2->SetMarkerColor(kBlue);
    h_2->SetLineColor(kBlue);
    h_2->SetMarkerStyle(24);
    h_2->Draw("PE1 HIST SAME");

    //Legend
    TLegend *leg = new TLegend(0.45,0.68,0.65,0.83, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20); 
    leg->AddEntry(h_1, "from 2D with 500 overflow bin, filled with 450");
    leg->AddEntry(h_1_750, "from 4D with 500 overflow bin, filled with 750");
    leg->AddEntry(h_1_1000, "from 2D with 500 overflow bin, filled with 1000");
    leg->AddEntry(h_1_1500, "from 2D with 500 overflow bin, filled with 1500");
    leg->AddEntry(h_2, "directly");
    leg->Draw("same");

    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");


    c->Draw();


    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/compare_filling_methods.pdf");

}

//______________________________________________________________________________________________
//____________________Uncertainties_____________________________________________________________
//______________________________________________________________________________________________

//Get the relative uncertainty between a nominal histogram and the variation
TH1D* get_relative_uncertainty(TH1D* h_nominal, TH1D* h_variation){
    TH1D* h_uncertainty = (TH1D*)h_nominal->Clone("h_uncertainty");
    h_uncertainty->Reset();

    Int_t nbins = h_uncertainty->GetNbinsX();

    Float_t max = 0;

    for(Int_t i = 1; i <= nbins; i++){
        Float_t difference = std::abs(h_variation->GetBinContent(i) - h_nominal->GetBinContent(i));
        Float_t uncertainty = difference/h_nominal->GetBinContent(i);
        h_uncertainty->SetBinContent(i, uncertainty);
        if(max < uncertainty) max = uncertainty;
    }

    std::cout << "max relative uncertainty = " << max << std::endl;
    return h_uncertainty;
}

//Get the statistical (relative) uncertainty of a nominal histogram h_nominal
TH1D* get_statistical_uncertainty(TH1D* h_nominal){
    TH1D* h_uncertainty = (TH1D*)h_nominal->Clone("h_uncertainty");
    h_uncertainty->Reset();

    Int_t nbins = h_uncertainty->GetNbinsX();

    Float_t max = 0;

    for(Int_t i = 1; i <= nbins; i++){
        Float_t uncertainty = h_nominal->GetBinError(i)/h_nominal->GetBinContent(i);
        h_uncertainty->SetBinContent(i, uncertainty);
        if(max < uncertainty) max = uncertainty;
    }

    std::cout << "max statistical uncertainty = " << max << std::endl;
    return h_uncertainty;
}

//Plot the systematics uncertainty related to the approximation of 1D EEC(dr) from a 4D histogram
void eec_4d_filling_uncertainty(){
    //Get HighEG and LowEG files (before template fit and unfolding)
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_data_HighEGJet_80_140.root", "read");
    TFile *file_1_low = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_data_lowEGJet_80_140.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_data_HighEGJet_80_140.root", "read");
    TFile *file_2_low = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_data_lowEGJet_80_140.root", "read");

    //Select pT bin
    Int_t pt_bin = 2;

    //Get and add data histograms
    THnD *h_4D = (THnD*)file_1->Get("h4D");
    THnD *h_4D_low = (THnD*)file_1_low->Get("h4D");
    h_4D->Add(h_4D_low);


    //Project/approximate to 1D EEC(dr) histogram
    TH3D *h_3D_2 = (TH3D*)file_2->Get("h3D");
    TH3D *h_3D_2_low = (TH3D*)file_2_low->Get("h3D");
    h_3D_2->Add(h_3D_2_low);

    h_4D->GetAxis(mb_dim)->SetRange(1, bins_mb);
    TH3D *h_3D = (TH3D*)h_4D->Projection(dr_dim, eec_dim, pt_dim)->Clone("h_3D");
    h_4D->GetAxis(pt_dim)->SetRange(pt_bin, pt_bin);

    TH2D *h_2D = (TH2D*)h_4D->Projection(eec_dim, dr_dim)->Clone("h_2D");

    TH1D *h_1 = new TH1D("h_1", "h_1", bins_dr, dr_binsVector);
    project_eec_1D(h_3D, h_1, pt_bin);
    TH1D *h_2 = (TH1D*)h_3D_2->ProjectionY("h_2", 1, bins_mb, 2, 2);

    //Normalize
    h_1->Scale(1/h_1->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c_plot", " ",170,800,700,504);
    /*c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();*/
    TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
    TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
    pad1->SetTopMargin(0.01);
    pad1->SetBottomMargin(0.3);
    pad2->SetBottomMargin(0.01);
    pad1->SetLogx();
    pad2->SetLogx();


    pad2->cd();
    h_1->SetStats(0);
    h_1->SetTitle("");//"EEC for 100 < pT < 120 GeV, single b jets, dijet reco MC");
    //h_1->GetXaxis()->SetTitle("\\Delta r");
    //h_1->GetXaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_1->GetYaxis()->CenterTitle(true);
    h_1->GetYaxis()->SetRangeUser(0, 15.);
    h_1->SetMarkerColor(kMagenta);
    h_1->SetLineColor(kMagenta);
    h_1->SetMarkerStyle(20);
    h_1->Draw("PE1 HIST");

    h_2->SetMarkerColor(kBlue);
    h_2->SetLineColor(kBlue);
    h_2->SetMarkerStyle(20);
    h_2->Draw("PE1 HIST SAME");

    //Legend
    TLegend *leg = new TLegend(0.6,0.7,0.8,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);
    leg->AddEntry(h_1, "from 2D with 500 overflow bin");
    leg->AddEntry(h_2, "directly");
    leg->Draw("same");

    
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");

    //Plot relative uncertainty
    TH1D* h_uncertainty = (TH1D*)get_relative_uncertainty(h_1, h_2);

    pad1->cd();
    h_uncertainty->SetStats(0);
    h_uncertainty->SetTitle("Relative uncertainty");
    gStyle->SetTitleSize(0.08, "t");
    h_uncertainty->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_uncertainty->GetXaxis()->CenterTitle(true);
    h_uncertainty->GetXaxis()->SetTitleSize(0.08);
    /*h_uncertainty->GetYaxis()->SetTitle("Relative uncertainty");
    h_uncertainty->GetYaxis()->CenterTitle(true);*/
    h_uncertainty->GetYaxis()->SetRangeUser(0.,0.1);
    h_uncertainty->SetMarkerStyle(kFullCircle);
    h_uncertainty->SetMarkerColor(kBlack);
    h_uncertainty->SetLineColor(kBlack);
    h_uncertainty->SetLineStyle(6);
    h_uncertainty->SetMarkerSize(1);
    h_uncertainty->GetXaxis()->SetTitleOffset(1.5);
    h_uncertainty->GetXaxis()->CenterTitle(true);
    h_uncertainty->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_uncertainty->GetXaxis()->SetLabelSize(0.1);
    h_uncertainty->GetXaxis()->SetTitleSize(0.08);
    h_uncertainty->GetYaxis()->SetTitle("");
    h_uncertainty->GetYaxis()->SetTitleOffset(0.4);
    h_uncertainty->GetYaxis()->CenterTitle(true);
    h_uncertainty->GetYaxis()->SetNdivisions(8);
    h_uncertainty->Draw("PE1 L HIST SAME");


    c->cd();
    pad1->Draw();
    pad2->Draw();
    c->Draw();


    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/compare_overflow_bin.png");
    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/compare_overflow_bin.pdf");

}
//Plot the systematics uncertainty related to the model used for the template fit
void eec_template_uncertainty(){
    //Get files
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D_herwig_template.root", "read");

    //Select pT bin
    Int_t pt_bin = 2;

    //Get and project histograms
    TH3D *h_pythia = (TH3D*)file_1->Get("h_data_unfolded");
    TH3D *h_herwig = (TH3D*)file_2->Get("h_data_unfolded");

    TH1D *h_pythia_1d = new TH1D("h_pythia_1d", "h_pythia_1d", bins_dr, dr_binsVector);
    TH1D *h_herwig_1d = new TH1D("h_herwig_1d", "h_herwig_1d", bins_dr, dr_binsVector);

    project_eec_1D(h_pythia, h_pythia_1d, pt_bin);
    project_eec_1D(h_herwig, h_herwig_1d, pt_bin);

    //Normalize
    h_pythia_1d->Scale(1/h_pythia_1d->Integral(), "width");
    h_herwig_1d->Scale(1/h_herwig_1d->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c_plot", " ",170,800,700,504);
    /*c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();*/
    TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
    TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
    pad1->SetTopMargin(0.01);
    pad1->SetBottomMargin(0.3);
    pad2->SetBottomMargin(0.01);
    pad1->SetLogx();
    pad2->SetLogx();


    pad2->cd();
    h_pythia_1d->SetStats(0);
    h_pythia_1d->SetTitle("");//"EEC for 100 < pT < 120 GeV, single b jets, unfolded and b-tag efficiency corrected");
    //h_pythia_1d->GetXaxis()->SetTitle("\\Delta r");
    //h_pythia_1d->GetXaxis()->CenterTitle(true);
    h_pythia_1d->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_pythia_1d->GetYaxis()->CenterTitle(true);
    h_pythia_1d->GetYaxis()->SetRangeUser(0, 15.);
    h_pythia_1d->SetMarkerColor(kMagenta);
    h_pythia_1d->SetLineColor(kMagenta);
    h_pythia_1d->SetMarkerStyle(20);
    h_pythia_1d->Draw("PE1 HIST");

    h_herwig_1d->SetMarkerColor(kBlue);
    h_herwig_1d->SetLineColor(kBlue);
    h_herwig_1d->SetMarkerStyle(20);
    h_herwig_1d->Draw("PE1 HIST SAME");

    //Legend
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);
    leg->AddEntry(h_pythia_1d, "Pythia template");
    leg->AddEntry(h_herwig_1d, "Herwig template");
    leg->Draw("same");

    //Calculate and plot uncertainty
    TH1D* h_uncertainty = (TH1D*)get_relative_uncertainty(h_pythia_1d, h_herwig_1d);

    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single B jets");
    test_info_text->DrawLatex(0.15, 0.65, "unfolded (Pythia) and b-tag corrected");
    test_info_text->Draw("SAME");

    pad1->cd();
    h_uncertainty->SetStats(0);
    h_uncertainty->SetTitle("Relative uncertainty");
    gStyle->SetTitleSize(0.08, "t");
    h_uncertainty->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_uncertainty->GetXaxis()->CenterTitle(true);
    h_uncertainty->GetXaxis()->SetTitleSize(0.08);
    /*h_uncertainty->GetYaxis()->SetTitle("Relative uncertainty");
    h_uncertainty->GetYaxis()->CenterTitle(true);*/
    h_uncertainty->GetYaxis()->SetRangeUser(0.,0.225);
    h_uncertainty->SetMarkerStyle(kFullCircle);
    h_uncertainty->SetMarkerColor(kBlack);
    h_uncertainty->SetLineColor(kBlack);
    h_uncertainty->SetLineStyle(6);
    h_uncertainty->SetMarkerSize(1);
    h_uncertainty->GetXaxis()->SetTitleOffset(1.5);
    h_uncertainty->GetXaxis()->CenterTitle(true);
    h_uncertainty->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_uncertainty->GetXaxis()->SetLabelSize(0.1);
    h_uncertainty->GetXaxis()->SetTitleSize(0.08);
    h_uncertainty->GetYaxis()->SetTitle("");
    h_uncertainty->GetYaxis()->SetTitleOffset(0.4);
    h_uncertainty->GetYaxis()->CenterTitle(true);
    h_uncertainty->GetYaxis()->SetNdivisions(8);

    h_uncertainty->Draw("PE1 L HIST");

    
    c->cd();
    pad1->Draw();
    pad2->Draw();
    c->Draw();


    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/template_uncertainty.png");
    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/template_uncertainty.pdf");

}

//Plot the systematics uncertainty related to the model used for the response matrix
void eec_4d_unfolding_uncertainty(){
    //Get files
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D_herwig.root", "read");

    //Select pT bin
    Int_t pt_bin = 2;

    //Get and project histograms
    TH3D *h_pythia = (TH3D*)file_1->Get("h_data_unfolded");
    TH3D *h_herwig = (TH3D*)file_2->Get("h_data_unfolded");

    TH1D *h_pythia_1d = new TH1D("h_pythia_1d", "h_pythia_1d", bins_dr, dr_binsVector);
    TH1D *h_herwig_1d = new TH1D("h_herwig_1d", "h_herwig_1d", bins_dr, dr_binsVector);


    project_eec_1D(h_pythia, h_pythia_1d, pt_bin);
    project_eec_1D(h_herwig, h_herwig_1d, pt_bin);

    //Normalize
    h_pythia_1d->Scale(1/h_pythia_1d->Integral(), "width");
    h_herwig_1d->Scale(1/h_herwig_1d->Integral(), "width");

    //Plot
    TCanvas *c = new TCanvas("c_plot", " ",170,800,700,504);
    /*c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();*/
    TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
    TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
    pad1->SetTopMargin(0.01);
    pad1->SetBottomMargin(0.3);
    pad2->SetBottomMargin(0.01);
    pad1->SetLogx();
    pad2->SetLogx();


    pad2->cd();
    h_pythia_1d->SetStats(0);
    h_pythia_1d->SetTitle("");//"EEC for 100 < pT < 120 GeV, single b jets, unfolded and b-tag efficiency corrected");
    //h_pythia_1d->GetXaxis()->SetTitle("\\Delta r");
    //h_pythia_1d->GetXaxis()->CenterTitle(true);
    h_pythia_1d->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_pythia_1d->GetYaxis()->CenterTitle(true);
    h_pythia_1d->GetYaxis()->SetRangeUser(0, 15.);
    h_pythia_1d->SetMarkerColor(kMagenta);
    h_pythia_1d->SetLineColor(kMagenta);
    h_pythia_1d->SetMarkerStyle(20);
    h_pythia_1d->Draw("PE1 HIST");

    h_herwig_1d->SetMarkerColor(kBlue);
    h_herwig_1d->SetLineColor(kBlue);
    h_herwig_1d->SetMarkerStyle(20);
    h_herwig_1d->Draw("PE1 HIST SAME");

    //Legend
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);
    leg->AddEntry(h_pythia_1d, "Pythia response");
    leg->AddEntry(h_herwig_1d, "Herwig response");
    leg->Draw("same");

    //Get and plot uncertainty
    TH1D* h_uncertainty = (TH1D*)get_relative_uncertainty(h_pythia_1d, h_herwig_1d);

    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.82, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.77, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.15, 0.72, "aggregated single B jets");
    test_info_text->DrawLatex(0.15, 0.67, "unfolded and b-tag corrected");
    test_info_text->Draw("SAME");

    pad1->cd();
    h_uncertainty->SetStats(0);
    h_uncertainty->SetTitle("Relative uncertainty");
    gStyle->SetTitleSize(0.08, "t");
    h_uncertainty->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_uncertainty->GetXaxis()->CenterTitle(true);
    h_uncertainty->GetXaxis()->SetTitleSize(0.08);
    /*h_uncertainty->GetYaxis()->SetTitle("Relative uncertainty");
    h_uncertainty->GetYaxis()->CenterTitle(true);*/
    h_uncertainty->GetYaxis()->SetRangeUser(0.,0.4);
    h_uncertainty->SetMarkerStyle(kFullCircle);
    h_uncertainty->SetMarkerColor(kBlack);
    h_uncertainty->SetLineColor(kBlack);
    h_uncertainty->SetLineStyle(6);
    h_uncertainty->SetMarkerSize(1);
    h_uncertainty->GetXaxis()->SetTitleOffset(1.5);
    h_uncertainty->GetXaxis()->CenterTitle(true);
    h_uncertainty->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_uncertainty->GetXaxis()->SetLabelSize(0.1);
    h_uncertainty->GetXaxis()->SetTitleSize(0.08);
    h_uncertainty->GetYaxis()->SetTitle("");
    h_uncertainty->GetYaxis()->SetTitleOffset(0.4);
    h_uncertainty->GetYaxis()->CenterTitle(true);
    h_uncertainty->GetYaxis()->SetNdivisions(8);

    h_uncertainty->Draw("PE1 L HIST");

    
    c->cd();
    pad1->Draw();
    pad2->Draw();
    c->Draw();


    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/model_uncertainty.png");
    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/model_uncertainty.pdf");

}

//Plot a summary plot of all the uncertainties
void eec_systematic_uncertainties(){
    //Get model uncertainty
    TFile *file_1_model = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D.root", "read");
    
    TFile *file_2_model = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D_herwig.root", "read");

    Int_t pt_bin = 2;

    TH3D *h_pythia = (TH3D*)file_1_model->Get("h_data_unfolded");
    TH3D *h_herwig = (TH3D*)file_2_model->Get("h_data_unfolded");

    TH1D *h_pythia_1d = new TH1D("h_pythia_1d", "h_pythia_1d", bins_dr, dr_binsVector);
    TH1D *h_herwig_1d = new TH1D("h_herwig_1d", "h_herwig_1d", bins_dr, dr_binsVector);


    project_eec_1D(h_pythia, h_pythia_1d, pt_bin);
    project_eec_1D(h_herwig, h_herwig_1d, pt_bin);

    h_pythia_1d->Scale(1/h_pythia_1d->Integral(), "width");
    h_herwig_1d->Scale(1/h_herwig_1d->Integral(), "width");

    //Get template fit uncertainty
    TFile *file_1_template = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D.root", "read");
    
    TFile *file_2_template = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_after_unfolding_4D_herwig_template.root", "read");

    TH3D *h_pythia_template = (TH3D*)file_1_template->Get("h_data_unfolded");
    TH3D *h_herwig_template = (TH3D*)file_2_template->Get("h_data_unfolded");

    TH1D *h_pythia_template_1d = new TH1D("h_pythia_template_1d", "h_pythia_template_1d", bins_dr, dr_binsVector);
    TH1D *h_herwig_template_1d = new TH1D("h_herwig_template_1d", "h_herwig_template_1d", bins_dr, dr_binsVector);


    project_eec_1D(h_pythia_template, h_pythia_template_1d, pt_bin);
    project_eec_1D(h_herwig_template, h_herwig_template_1d, pt_bin);

    h_pythia_template_1d->Scale(1/h_pythia_template_1d->Integral(), "width");
    h_herwig_template_1d->Scale(1/h_herwig_template_1d->Integral(), "width");

    //Get filling uncertainty
    TFile *file_1 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_data_HighEGJet_80_140.root", "read");
    TFile *file_1_low = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_4d_data_lowEGJet_80_140.root", "read");
    
    TFile *file_2 = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_data_HighEGJet_80_140.root", "read");
    TFile *file_2_low = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/hist_3d_data_lowEGJet_80_140.root", "read");

    THnD *h_4D = (THnD*)file_1->Get("h4D");
    THnD *h_4D_low = (THnD*)file_1_low->Get("h4D");
    h_4D->Add(h_4D_low);

    THnD *h_4D_simple = (THnD*)h_4D->Clone("h_4D_simple");

    TH3D *h_3D_2 = (TH3D*)file_2->Get("h3D");
    TH3D *h_3D_2_low = (TH3D*)file_2_low->Get("h3D");
    h_3D_2->Add(h_3D_2_low);

    h_4D->GetAxis(mb_dim)->SetRange(1, bins_mb);
    TH3D *h_3D = (TH3D*)h_4D->Projection(dr_dim, eec_dim, pt_dim)->Clone("h_3D");
    h_4D->GetAxis(pt_dim)->SetRange(pt_bin, pt_bin);
    rebin_eec(h_4D);
    h_4D_simple->GetAxis(mb_dim)->SetRange(1, bins_mb);
    h_4D_simple->GetAxis(pt_dim)->SetRange(pt_bin, pt_bin);

    TH2D *h_2D = (TH2D*)h_4D->Projection(eec_dim, dr_dim)->Clone("h_2D");


    TH1D *h_1 = new TH1D("h_1", "h_1", bins_dr, dr_binsVector);
    project_eec_1D(h_3D, h_1, pt_bin);
    TH1D *h_2 = (TH1D*)h_3D_2->ProjectionY("h_2", 1, bins_mb, 2, 2);

    h_1->Scale(1/h_1->Integral(), "width");
    h_2->Scale(1/h_2->Integral(), "width");

    //Get uncertainties
    TH1D* h_uncertainty_model = (TH1D*)get_relative_uncertainty(h_pythia_1d, h_herwig_1d);
    TH1D* h_uncertainty_template = (TH1D*)get_relative_uncertainty(h_pythia_template_1d, h_herwig_template_1d);
    TH1D* h_uncertainty_filling = (TH1D*)get_relative_uncertainty(h_1, h_2);
    TH1D* h_statistical = (TH1D*)get_statistical_uncertainty(h_1);

    //plot
    TCanvas *c = new TCanvas("c_plot", " ",900,400,840,304);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    //Legend
    TLegend *leg = new TLegend(0.6,0.7,0.8,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);


    
    h_statistical->Scale(10); //The statistical uncertainty is scaled by a factor 10 for visualization purposes
    h_statistical->SetStats(0);
    h_statistical->SetTitle("");
    h_statistical->GetYaxis()->SetRangeUser(0.,0.4);
    h_statistical->SetMarkerStyle(kFullCircle);
    h_statistical->SetMarkerColor(kGray);
    h_statistical->SetLineColor(kGray);
    //h_statistical->SetLineStyle(6);
    h_statistical->SetMarkerSize(0);
    h_statistical->SetFillColor(kGray);
    h_statistical->GetXaxis()->CenterTitle(true);
    h_statistical->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_statistical->GetYaxis()->SetTitle("Relative uncertainty");
    h_statistical->GetYaxis()->CenterTitle(true);
    h_statistical->GetYaxis()->SetNdivisions(8);
    leg->AddEntry(h_statistical, "Statistical uncertainty x 10");
    
    h_uncertainty_model->SetMarkerStyle(kFullCircle);
    h_uncertainty_model->SetLineStyle(6);
    h_uncertainty_model->SetMarkerColor(kCyan);
    h_uncertainty_model->SetLineColor(kCyan);
    h_uncertainty_model->SetMarkerSize(1);
    leg->AddEntry(h_uncertainty_model, "Model uncertainty");

    h_uncertainty_template->SetMarkerStyle(kFullCircle);
    h_uncertainty_template->SetLineStyle(6);
    h_uncertainty_template->SetMarkerColor(kOrange);
    h_uncertainty_template->SetLineColor(kOrange);
    h_uncertainty_template->SetMarkerSize(1);
    leg->AddEntry(h_uncertainty_template, "Template fit uncertainty");

    h_uncertainty_filling->SetStats(0);
    h_uncertainty_filling->SetMarkerStyle(kFullCircle);
    h_uncertainty_filling->SetMarkerColor(kMagenta);
    h_uncertainty_filling->SetLineColor(kMagenta);
    h_uncertainty_filling->SetLineStyle(6);
    h_uncertainty_filling->SetMarkerSize(1);
    leg->AddEntry(h_uncertainty_filling, "Filling uncertainty");

    h_statistical->Draw("HIST SAME");
    h_uncertainty_model->Draw("P0 E L HIST SAME");
    h_uncertainty_template->Draw("P0 E L HIST SAME");
    h_uncertainty_filling->Draw("P0 E L HIST SAME");

    
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.8, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.75, "Data HighEGJet + LowEGJet");
    test_info_text->DrawLatex(0.15, 0.7, "aggregated single B jets");
    test_info_text->Draw("SAME");

    leg->Draw("SAME");

    
    c->Draw();


    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/systematic_uncertaintie.pdf");

}

//______________________________________________________________________________________________
//____________________Plots related to pair matching____________________________________________
//______________________________________________________________________________________________

//Compare the normalized EEC(dr) distribution for matched/non-matched pairs
void compare_matched_nonmatched_norm(){
    //Select pT bin
    Int_t pt_bin = 2;

    //Get file
    TFile *file_response = new TFile("/data_CMS/cms/meuli/wp09_4d/histos_bjet_b_response_v9.root", "read");
    
    //Select matched, non-matched or full EEC pairs distributions at reco or gen
    TH3D *h_3D_matched_reco = (TH3D*)file_response->Get("h_full_purity_numerator_eecpt")->Clone("h_3D_matched_reco");
    h_3D_matched_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_matched_reco = (TH2D*)h_3D_matched_reco->Project3D("yx")->Clone("h_2D_matched_reco");
    
    TH3D *h_3D_matched_gen = (TH3D*)file_response->Get("h_full_efficiency_numerator_eecpt")->Clone("h_3D_matched_gen");
    h_3D_matched_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_matched_gen = (TH2D*)h_3D_matched_gen->Project3D("yx")->Clone("h_2D_matched_gen");


    TH3D *h_3D_nonmatched_reco = (TH3D*)file_response->Get("h_notmatched_reco")->Clone("h_3D_nonmatched_reco");
    h_3D_nonmatched_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_nonmatched_reco = (TH2D*)h_3D_nonmatched_reco->Project3D("yx")->Clone("h_2D_nonmatched_reco");
    
    TH3D *h_3D_nonmatched_gen = (TH3D*)file_response->Get("h_notmatched_gen")->Clone("h_3D_nonmatched_gen");
    h_3D_nonmatched_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_nonmatched_gen = (TH2D*)h_3D_nonmatched_gen->Project3D("yx")->Clone("h_2D_nonmatched_gen");


    TH3D *h_3D_all_reco = (TH3D*)file_response->Get("h_full_purity_denominator_eecpt")->Clone("h_3D_all_reco");
    h_3D_all_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_all_reco = (TH2D*)h_3D_all_reco->Project3D("yx")->Clone("h_2D_all_reco");
    
    TH3D *h_3D_all_gen = (TH3D*)file_response->Get("h_full_efficiency_denominator_eecpt")->Clone("h_3D_all_gen");
    h_3D_all_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_all_gen = (TH2D*)h_3D_all_gen->Project3D("yx")->Clone("h_2D_all_gen");

    //Approximate to 1D EEC(dr) distributions
    TH1D *h_matched_reco = new TH1D("h_matched_reco", "h_matched_reco", bins_dr, dr_binsVector);
    TH1D *h_matched_gen = new TH1D("h_matched_gen", "h_matched_gen", bins_dr, dr_binsVector);
    TH1D *h_nonmatched_reco = new TH1D("h_nonmatched_reco", "h_nonmatched_reco", bins_dr, dr_binsVector);
    TH1D *h_nonmatched_gen = new TH1D("h_nonmatched_gen", "h_nonmatched_gen", bins_dr, dr_binsVector);
    TH1D *h_all_reco = new TH1D("h_all_reco", "h_all_reco", bins_dr, dr_binsVector);
    TH1D *h_all_gen = new TH1D("h_all_gen", "h_all_gen", bins_dr, dr_binsVector);


    recover_eec_distr(h_matched_reco, h_2D_matched_reco, 1000);
    recover_eec_distr(h_matched_gen, h_2D_matched_gen, 1000);
    recover_eec_distr(h_nonmatched_reco, h_2D_nonmatched_reco, 1000);
    recover_eec_distr(h_nonmatched_gen, h_2D_nonmatched_gen, 1000);
    recover_eec_distr(h_all_reco, h_2D_all_reco, 1000);
    recover_eec_distr(h_all_gen, h_2D_all_gen, 1000);

    
    //Simple projection without energy weighting
    /*TH1D *h_matched_reco = (TH1D*)h_2D_matched_reco->ProjectionX("h_matched_reco", 1,5);
    TH1D *h_matched_gen = (TH1D*)h_2D_matched_gen->ProjectionX("h_matched_gen", 1,5);
    TH1D *h_nonmatched_reco = (TH1D*)h_2D_nonmatched_reco->ProjectionX("h_nonmatched_reco", 1,5);
    TH1D *h_nonmatched_gen = (TH1D*)h_2D_nonmatched_gen->ProjectionX("h_nonmatched_gen", 1,5);
    TH1D *h_all_reco = (TH1D*)h_2D_all_reco->ProjectionX("h_all_reco", 1,5);
    TH1D *h_all_gen = (TH1D*)h_2D_all_gen->ProjectionX("h_all_gen", 1,5);*/

    //Normalize and get the maximum y value for the plot
    double ymax = 0.;
        for (auto h : {h_matched_reco, h_matched_gen, 
                       h_nonmatched_reco, h_nonmatched_gen, 
                       h_all_reco, h_all_gen
                    }) {
                        h->GetXaxis()->SetRange(1, bins_dr);
                        h->Scale(1/h->Integral(), "width");
                        h->GetXaxis()->SetRange(2, bins_dr);
                        if (h->GetMaximum()>ymax) ymax =  h->GetMaximum();
                    }

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_matched_reco->SetStats(0);
    h_matched_reco->SetTitle("Effect of pair matching on EEC");
    h_matched_reco->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_matched_reco->GetXaxis()->CenterTitle(true);
    h_matched_reco->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})\\mbox{ Normalised}");
    h_matched_reco->GetYaxis()->CenterTitle(true);
    h_matched_reco->GetYaxis()->SetRangeUser(0, 1.1*ymax);
    h_matched_reco->SetMarkerColor(kGreen);
    h_matched_reco->SetLineColor(kGreen);
    h_matched_reco->SetMarkerStyle(24);
    h_matched_reco->Draw("P0 E HIST");

    h_nonmatched_reco->SetMarkerColor(kBlue+1);
    h_nonmatched_reco->SetLineColor(kBlue+1);
    h_nonmatched_reco->SetMarkerStyle(24);
    h_nonmatched_reco->Draw("P0 E HIST SAME");
    
    h_matched_gen->SetMarkerColor(kRed);
    h_matched_gen->SetLineColor(kRed);
    h_matched_gen->SetMarkerStyle(24);
    h_matched_gen->Draw("P0 E HIST SAME");

    h_nonmatched_gen->SetMarkerColor(kMagenta+2);
    h_nonmatched_gen->SetLineColor(kMagenta+2);
    h_nonmatched_gen->SetMarkerStyle(24);
    h_nonmatched_gen->Draw("P0 E HIST SAME");

    h_all_reco->SetMarkerColor(kCyan);
    h_all_reco->SetLineColor(kCyan);
    h_all_reco->SetMarkerStyle(21);
    h_all_reco->Draw("P0 E HIST SAME");

    
    h_all_gen->SetMarkerColor(kMagenta);
    h_all_gen->SetLineColor(kMagenta);
    h_all_gen->SetMarkerStyle(21);
    h_all_gen->Draw("P0 E HIST SAME");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.7, 0.62, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.7, 0.57, "MC simulation");
    test_info_text->DrawLatex(0.7, 0.52, "single B jets");
    test_info_text->DrawLatex(0.7, 0.47, "BDT score > -0.3");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.65,0.7,0.75,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h_all_reco, "all reco");
    leg->AddEntry(h_matched_reco, "reco matched");
    leg->AddEntry(h_nonmatched_reco, "reco non-matched");
    leg->AddEntry(h_all_gen, "all gen");
    leg->AddEntry(h_matched_gen, "gen matched");
    leg->AddEntry(h_nonmatched_gen, "gen non-matched");
    leg->Draw("same");

    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins/eec_compare_matched_notmatched_norm_v8.pdf");

}

//Compare the simple (not energy weighted) dr distribution for matched/non-matched pairs and plots the
//ratio of matched over all
void compare_matched_nonmatched_ratio(){
    //Select pT bin
    Int_t pt_bin = 2;

    //Get file
    TFile *file_response = new TFile("/data_CMS/cms/meuli/wp09_4d/histos_bjet_b_response_v9.root", "read");
    
    //Get matched, non-matched and full distributions
    TH3D *h_3D_matched_reco = (TH3D*)file_response->Get("h_full_purity_numerator_eecpt")->Clone("h_3D_matched_reco");
    h_3D_matched_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_matched_reco = (TH2D*)h_3D_matched_reco->Project3D("yx")->Clone("h_2D_matched_reco");
    
    TH3D *h_3D_matched_gen = (TH3D*)file_response->Get("h_full_efficiency_numerator_eecpt")->Clone("h_3D_matched_gen");
    h_3D_matched_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_matched_gen = (TH2D*)h_3D_matched_gen->Project3D("yx")->Clone("h_2D_matched_gen");


    TH3D *h_3D_nonmatched_reco = (TH3D*)file_response->Get("h_notmatched_reco")->Clone("h_3D_nonmatched_reco");
    h_3D_nonmatched_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_nonmatched_reco = (TH2D*)h_3D_nonmatched_reco->Project3D("yx")->Clone("h_2D_nonmatched_reco");
    
    TH3D *h_3D_nonmatched_gen = (TH3D*)file_response->Get("h_notmatched_gen")->Clone("h_3D_nonmatched_gen");
    h_3D_nonmatched_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_nonmatched_gen = (TH2D*)h_3D_nonmatched_gen->Project3D("yx")->Clone("h_2D_nonmatched_gen");

    
    TH3D *h_3D_all_reco = (TH3D*)file_response->Get("h_full_purity_denominator_eecpt")->Clone("h_3D_all_reco");
    h_3D_all_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_all_reco = (TH2D*)h_3D_all_reco->Project3D("yx")->Clone("h_2D_all_reco");
    
    TH3D *h_3D_all_gen = (TH3D*)file_response->Get("h_full_efficiency_denominator_eecpt")->Clone("h_3D_all_gen");
    h_3D_all_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_all_gen = (TH2D*)h_3D_all_gen->Project3D("yx")->Clone("h_2D_all_gen");

    //Recover energy weighted EEC(dr) distributions
    /*TH1D *h_matched_reco = new TH1D("h_matched_reco", "h_matched_reco", bins_dr, dr_binsVector);
    TH1D *h_matched_gen = new TH1D("h_matched_gen", "h_matched_gen", bins_dr, dr_binsVector);
    TH1D *h_nonmatched_reco = new TH1D("h_nonmatched_reco", "h_nonmatched_reco", bins_dr, dr_binsVector);
    TH1D *h_nonmatched_gen = new TH1D("h_nonmatched_gen", "h_nonmatched_gen", bins_dr, dr_binsVector);
    TH1D *h_all_reco = new TH1D("h_all_reco", "h_all_reco", bins_dr, dr_binsVector);
    TH1D *h_all_gen = new TH1D("h_all_gen", "h_all_gen", bins_dr, dr_binsVector);


    recover_eec_distr(h_matched_reco, h_2D_matched_reco, 1000);
    recover_eec_distr(h_matched_gen, h_2D_matched_gen, 1000);
    recover_eec_distr(h_nonmatched_reco, h_2D_nonmatched_reco, 1000);
    recover_eec_distr(h_nonmatched_gen, h_2D_nonmatched_gen, 1000);
    recover_eec_distr(h_all_reco, h_2D_all_reco, 1000);
    recover_eec_distr(h_all_gen, h_2D_all_gen,1000);*/

    //Simple projection as dr distributions
    TH1D *h_matched_reco = (TH1D*)h_2D_matched_reco->ProjectionX("h_matched_reco", 1,5);
    TH1D *h_matched_gen = (TH1D*)h_2D_matched_gen->ProjectionX("h_matched_gen", 1,5);
    TH1D *h_nonmatched_reco = (TH1D*)h_2D_nonmatched_reco->ProjectionX("h_nonmatched_reco", 1,5);
    TH1D *h_nonmatched_gen = (TH1D*)h_2D_nonmatched_gen->ProjectionX("h_nonmatched_gen", 1,5);
    TH1D *h_all_reco = (TH1D*)h_2D_all_reco->ProjectionX("h_all_reco", 1,5);
    TH1D *h_all_gen = (TH1D*)h_2D_all_gen->ProjectionX("h_all_gen", 1,5);

    //(Normalize) and get the maximum y axis value for the plots
    double ymax = 0.;
    for (auto h : {h_matched_reco, h_matched_gen, h_nonmatched_reco,
                       h_nonmatched_gen, h_all_reco, h_all_gen
                    }) {
                        //h->Scale(1/h->Integral(), "width");
                    if (h->GetMaximum()>ymax) ymax =  h->GetMaximum();
    }
    
    //Plot
    TCanvas *c = new TCanvas("c", "", 800, 600);
    TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
    TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
    pad1->SetTopMargin(0.01);
    pad1->SetBottomMargin(0.3);
    pad2->SetBottomMargin(0.01);
    pad1->SetLogx();
    pad2->SetLogx();



    pad2->cd();
    h_matched_reco->SetStats(0);
    //h_matched_reco->SetTitle("Effect of pair matching on EEC");
    h_matched_reco->SetTitle("");
    h_matched_reco->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_matched_reco->GetXaxis()->CenterTitle(true);
    //h_matched_reco->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})");
    h_matched_reco->GetYaxis()->SetTitle("Entries (with MC weight)");
    h_matched_reco->GetYaxis()->CenterTitle(true);
    h_matched_reco->GetYaxis()->SetRangeUser(0, 1.1*ymax);
    h_matched_reco->SetMarkerColor(kCyan);
    h_matched_reco->SetLineColor(kCyan);
    h_matched_reco->SetMarkerStyle(24);
    h_matched_reco->Draw("P0 E1 HIST");

    /*h_nonmatched_reco->SetMarkerColor(kBlue+1);
    h_nonmatched_reco->SetLineColor(kBlue+1);
    h_nonmatched_reco->SetMarkerStyle(24);
    h_nonmatched_reco->Draw("P0 E1 HIST SAME");*/
    
    h_matched_gen->SetMarkerColor(kMagenta);
    h_matched_gen->SetLineColor(kMagenta);
    h_matched_gen->SetMarkerStyle(24);
    h_matched_gen->Draw("P0 E1 HIST SAME");

    
    /*h_nonmatched_gen->SetMarkerColor(kMagenta+2);
    h_nonmatched_gen->SetLineColor(kMagenta+2);
    h_nonmatched_gen->SetMarkerStyle(24);
    h_nonmatched_gen->Draw("P0 E1 HIST SAME");*/

    h_all_reco->SetMarkerColor(kCyan);
    h_all_reco->SetLineColor(kCyan);
    h_all_reco->SetMarkerStyle(21);
    h_all_reco->Draw("P0 E1 HIST SAME");

    
    h_all_gen->SetMarkerColor(kMagenta);
    h_all_gen->SetLineColor(kMagenta);
    h_all_gen->SetMarkerStyle(21);
    h_all_gen->Draw("P0 E1 HIST SAME");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.82, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.15, 0.77, "MC simulation");
    test_info_text->DrawLatex(0.15, 0.72, "single B jets");
    test_info_text->DrawLatex(0.15, 0.67, "BDT score > -0.9");
    //test_info_text->DrawLatex(0.15, 0.62, "3D response matrix");
    //test_info_text->DrawLatex(0.15, 0.67, "trkMatchSta >= 100");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.135,0.4,0.25,0.57, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h_all_reco, "all reco");
    leg->AddEntry(h_matched_reco, "reco matched");
    //leg->AddEntry(h_nonmatched_reco, "reco non-matched");
    leg->AddEntry(h_all_gen, "all gen");
    leg->AddEntry(h_matched_gen, "gen matched");
    //leg->AddEntry(h_nonmatched_gen, "gen non-matched");
    leg->Draw("same");


    TLegend *leg_ratio = new TLegend(0.65, 0.3, 0.9, 0.5);//0.5, 0.3, 0.85, 0.5);
    leg_ratio->SetBorderSize(1);
    
    TH1D *h_reco_ratio = (TH1D*)h_matched_reco->Clone("h_reco_ratio");
    TH1D *h_gen_ratio = (TH1D*)h_matched_reco->Clone("h_gen_ratio");

    h_reco_ratio->Divide(h_all_reco);
    h_gen_ratio->Divide(h_all_gen);

    h_reco_ratio->SetTitle("");
    h_reco_ratio->SetStats(0);
    h_reco_ratio->SetMarkerStyle(kFullCircle);
    h_reco_ratio->SetMarkerColor(kCyan);
    h_reco_ratio->SetLineColor(kCyan);
    h_reco_ratio->SetMarkerSize(1);
    h_reco_ratio->GetYaxis()->SetRangeUser(0.,1.1);
    h_reco_ratio->GetYaxis()->SetTitle("ratio");
    h_reco_ratio->GetYaxis()->SetLabelSize(0.08);
    h_reco_ratio->GetYaxis()->SetTitleSize(0.08);
    h_reco_ratio->GetXaxis()->SetTitleOffset(1.5);
    h_reco_ratio->GetXaxis()->CenterTitle(true);
    h_reco_ratio->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_reco_ratio->GetXaxis()->SetLabelSize(0.1);
    //h_reco_ratio->GetXaxis()->SetTitleSize(text_size);
    //h_reco_ratio->GetXaxis()->SetTitleOffset(3.5);
    h_reco_ratio->GetXaxis()->SetTitleSize(0.08);
    h_reco_ratio->GetYaxis()->SetTitleOffset(0.4);
    h_reco_ratio->GetYaxis()->CenterTitle(true);
    h_reco_ratio->GetYaxis()->SetNdivisions(8);
    leg_ratio->AddEntry(h_reco_ratio, "matched reco / tot reco", "pe1");

    h_gen_ratio->SetMarkerStyle(kFullCircle);
    h_gen_ratio->SetMarkerColor(kMagenta);
    h_gen_ratio->SetLineColor(kMagenta);
    h_gen_ratio->SetMarkerSize(1);
    leg_ratio->AddEntry(h_gen_ratio, "matched gen / tot gen", "pe1");

    pad1->cd();
    h_reco_ratio->Draw("pe1 same");
    h_gen_ratio->Draw("pe1 same");
    leg_ratio->Draw();

    c->cd();
    pad1->Draw();
    pad2->Draw();
    c->Draw();

    c->Print("/data_CMS/cms/meuli/wp09_4d/eec_compare_matched_notmatched_v8.pdf");

}

//Reduced version of compare_matched_nonmatched_norm
void compare_matched_nonmatched_reduced(){
    //Select pT bin
    Int_t pt_bin = 2;

    //Get file
    TFile *file_response = new TFile("/data_CMS/cms/meuli/wp09_4d/newbins2/histos_bjet_b_response_rerun_old.root", "read");
    
    //Get matched, non-matched and full distributions 
    TH3D *h_3D_matched_reco = (TH3D*)file_response->Get("h_full_purity_numerator_eecpt")->Clone("h_3D_matched_reco");
    h_3D_matched_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_matched_reco = (TH2D*)h_3D_matched_reco->Project3D("yx")->Clone("h_2D_matched_reco");
    
    TH3D *h_3D_matched_gen = (TH3D*)file_response->Get("h_full_efficiency_numerator_eecpt")->Clone("h_3D_matched_gen");
    h_3D_matched_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_matched_gen = (TH2D*)h_3D_matched_gen->Project3D("yx")->Clone("h_2D_matched_gen");

    
    TH3D *h_3D_nonmatched_reco = (TH3D*)file_response->Get("h_notmatched_reco")->Clone("h_3D_nonmatched_reco");
    h_3D_nonmatched_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_nonmatched_reco = (TH2D*)h_3D_nonmatched_reco->Project3D("yx")->Clone("h_2D_nonmatched_reco");
    
    TH3D *h_3D_nonmatched_gen = (TH3D*)file_response->Get("h_notmatched_gen")->Clone("h_3D_nonmatched_gen");
    h_3D_nonmatched_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_nonmatched_gen = (TH2D*)h_3D_nonmatched_gen->Project3D("yx")->Clone("h_2D_nonmatched_gen");

    
    TH3D *h_3D_all_reco = (TH3D*)file_response->Get("h_full_purity_denominator_eecpt")->Clone("h_3D_all_reco");
    h_3D_all_reco->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_all_reco = (TH2D*)h_3D_all_reco->Project3D("yx")->Clone("h_2D_all_reco");
    
    TH3D *h_3D_all_gen = (TH3D*)file_response->Get("h_full_efficiency_denominator_eecpt")->Clone("h_3D_all_gen");
    h_3D_all_gen->GetZaxis()->SetRange(pt_bin, pt_bin);
    TH2D *h_2D_all_gen = (TH2D*)h_3D_all_gen->Project3D("yx")->Clone("h_2D_all_gen");

    //Recover the 1D EEC(dr) distributions
    TH1D *h_matched_reco = new TH1D("h_matched_reco", "h_matched_reco", bins_dr, dr_binsVector);
    TH1D *h_matched_gen = new TH1D("h_matched_gen", "h_matched_gen", bins_dr, dr_binsVector);
    TH1D *h_nonmatched_reco = new TH1D("h_nonmatched_reco", "h_nonmatched_reco", bins_dr, dr_binsVector);
    TH1D *h_nonmatched_gen = new TH1D("h_nonmatched_gen", "h_nonmatched_gen", bins_dr, dr_binsVector);
    TH1D *h_all_reco = new TH1D("h_all_reco", "h_all_reco", bins_dr, dr_binsVector);
    TH1D *h_all_gen = new TH1D("h_all_gen", "h_all_gen", bins_dr, dr_binsVector);


    recover_eec_distr(h_matched_reco, h_2D_matched_reco, 1000);
    recover_eec_distr(h_matched_gen, h_2D_matched_gen, 1000);
    recover_eec_distr(h_nonmatched_reco, h_2D_nonmatched_reco, 1000);
    recover_eec_distr(h_nonmatched_gen, h_2D_nonmatched_gen, 1000);
    recover_eec_distr(h_all_reco, h_2D_all_reco, 1000);
    recover_eec_distr(h_all_gen, h_2D_all_gen, 1000);

    //Simple dr projection without energy weighting
    /*TH1D *h_matched_reco = (TH1D*)h_2D_matched_reco->ProjectionX("h_matched_reco", 1,5);
    TH1D *h_matched_gen = (TH1D*)h_2D_matched_gen->ProjectionX("h_matched_gen", 1,5);
    TH1D *h_nonmatched_reco = (TH1D*)h_2D_nonmatched_reco->ProjectionX("h_nonmatched_reco", 1,5);
    TH1D *h_nonmatched_gen = (TH1D*)h_2D_nonmatched_gen->ProjectionX("h_nonmatched_gen", 1,5);
    TH1D *h_all_reco = (TH1D*)h_2D_all_reco->ProjectionX("h_all_reco", 1,5);
    TH1D *h_all_gen = (TH1D*)h_2D_all_gen->ProjectionX("h_all_gen", 1,5);*/

    //Normalise and get the maximum y value for plotting
    double ymax = 0.;
        for (auto h : {h_matched_reco, h_matched_gen, 
                       h_nonmatched_reco, h_nonmatched_gen, 
                       h_all_reco, h_all_gen
                    }) {
                        h->GetXaxis()->SetRange(1, bins_dr);
                        h->Scale(1/h->Integral(), "width");
                        h->GetXaxis()->SetRange(2, bins_dr);
                        if (h->GetMaximum()>ymax) ymax =  h->GetMaximum();
                    }

    //Plot
    TCanvas *c = new TCanvas("c", " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetLogx();

    h_nonmatched_reco->SetStats(0);
    h_nonmatched_reco->SetTitle("Effect of pair matching on EEC");
    h_nonmatched_reco->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_nonmatched_reco->GetXaxis()->CenterTitle(true);
    h_nonmatched_reco->GetYaxis()->SetTitle("EEC(\\Delta\\mbox{r})\\mbox{ Normalised}");
    h_nonmatched_reco->GetYaxis()->CenterTitle(true);
    h_nonmatched_reco->GetYaxis()->SetRangeUser(0, 1.1*ymax);
    h_nonmatched_reco->SetMarkerColor(kBlue+1);
    h_nonmatched_reco->SetLineColor(kBlue+1);
    h_nonmatched_reco->SetMarkerStyle(24);
    h_nonmatched_reco->Draw("P0 E HIST");

    h_all_reco->SetMarkerColor(kCyan);
    h_all_reco->SetLineColor(kCyan);
    h_all_reco->SetMarkerStyle(21);
    h_all_reco->Draw("P0 E HIST SAME");

    
    h_all_gen->SetMarkerColor(kMagenta);
    h_all_gen->SetLineColor(kMagenta);
    h_all_gen->SetMarkerStyle(21);
    h_all_gen->Draw("P0 E HIST SAME");

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.7, 0.62, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.7, 0.57, "MC simulation");
    test_info_text->DrawLatex(0.7, 0.52, "single B jets");
    test_info_text->DrawLatex(0.7, 0.47, "BDT score > -0.9");
    test_info_text->Draw("SAME");

    //Legend
    TLegend *leg = new TLegend(0.65,0.7,0.75,0.85, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h_all_reco, "all reco");
    leg->AddEntry(h_nonmatched_reco, "reco non-matched");
    leg->AddEntry(h_all_gen, "all gen");
    leg->Draw("same");

    c->Print("/data_CMS/cms/meuli/wp09_4d/newbins2/eec_compare_matched_notmatched_norm_reduced.pdf");

}


void test_plots_november(){

    //Select a datasample and pT range for the stack histograms
    TString dataset = "bjet";
    TString pT_selection = "80_140";

    //Call the plots that you want to see

    plot_eec_distribution();
    //flavour_effect();
}
