#include "tTree.h"
#include "TFile.h"

//Rebin the mB axis to have a bigger overflow bin
void rebin_mass(TH3D* &h){
    Int_t bins_mb = h->GetNbinsX();
    Int_t bins_pt = h->GetNbinsZ();
    Int_t bins_dr = h->GetNbinsY();

    Int_t last_mb_bin = 12;

    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){
            Float_t last_bin = 0;
            for(Int_t imb_bin = last_mb_bin; imb_bin <= bins_mb; imb_bin++){
                last_bin += h->GetBinContent(imb_bin, ibin_dr, ibin_pt);
            }
            h->GetXaxis()->SetRange(1, last_mb_bin);
            h->SetBinContent(last_mb_bin, ibin_dr, ibin_pt);
        }
    } 
}

//Rebin the dr axis to make the bins start later
void rebin_dr(TH3D* &h){

    Int_t first_bin = 3;
    Int_t bins_dr = h->GetNbinsY();
    Int_t bins_mb = h->GetNbinsX();
    Int_t bins_pt = h->GetNbinsZ();

    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t imb_bin = 1; imb_bin <= bins_mb; imb_bin++){
            Float_t first_bin_content = 0;
            for(Int_t ibin_dr = 1; ibin_dr <= first_bin; ibin_dr++){
                first_bin_content += h->GetBinContent(imb_bin, ibin_dr, ibin_pt);
            }
            h->GetYaxis()->SetRange(first_bin, bins_dr);
            h->SetBinContent(imb_bin, first_bin, ibin_pt, first_bin_content);
            
            
        }
    }
}

//Performs 3D template fit of mB for every dr and jtpt bin
//do_hist(TString &filename, TString folder, int &dataType, bool &isMC, Float_t &pT_low, Float_t &pT_high, Int_t &n, bool &btag, Int_t &lowEn, Int_t &highEn, bool aggregated = true,  bool ideal_aggr = false){                                                                                                    
 

void do_template_fit(TString pT_selection, TString folder, TString &fout_name, bool &also_bjet){


     //Open dataset
    TFile *file_data = new TFile(folder + "hist_3d_aggr_BDT_n1_data_highEn" + "_" + pT_selection + ".root", "read");
    TH3D *h_data = (TH3D*)file_data->Get("h3D_all")->Clone("h_data");

    TFile *file_data_low = new TFile(folder + "hist_3d_aggr_BDT_n1_data_lowEn" + "_" + pT_selection + ".root", "read");
    TH3D *h_data_lowEG = (TH3D*)file_data_low->Get("h3D_all")->Clone("h_data_lowEG");                                                                                   
    h_data->Add(h_data_lowEG);
 
    //Open file and get histograms from bjet
    TFile *file_bjet = new TFile(folder + "hist_3d_aggr_BDT_n1_MC_bjet" + "_" + pT_selection + ".root", "read");
    TH3D *h_bjet_b1 = (TH3D*)file_bjet->Get("h3D_b1")->Clone("h_bjet_b1");
    TH3D *h_bjet_b2 = (TH3D*)file_bjet->Get("h3D_b2")->Clone("h_bjet_b2");
    TH3D *h_bjet_c = (TH3D*)file_bjet->Get("h3D_c")->Clone("h_bjet_c");
    TH3D *h_bjet_light = (TH3D*)file_bjet->Get("h3D_light")->Clone("h_bjet_light");
    TH3D *h_bjet_all = (TH3D*)file_bjet->Get("h3D_all")->Clone("h_bjet_all");

    TH3D *h_bjet_other = (TH3D*)h_bjet_c->Clone("h_bjet_other");
    h_bjet_other->Add(h_bjet_light); 
    
    //Open file and get histograms from dijet
    TFile *file_dijet = new TFile(folder + "hist_3d_aggr_BDT_n1_MC_dijet" + "_" + pT_selection + ".root", "read");
    TH3D *h_dijet_b1 = (TH3D*)file_dijet->Get("h3D_b1")->Clone("h_dijet_b1");
    TH3D *h_dijet_b2 = (TH3D*)file_dijet->Get("h3D_b2")->Clone("h_dijet_b2");
    TH3D *h_dijet_c = (TH3D*)file_dijet->Get("h3D_c")->Clone("h_dijet_c");
    TH3D *h_dijet_light = (TH3D*)file_dijet->Get("h3D_light")->Clone("h_dijet_light");
    TH3D *h_dijet_all = (TH3D*)file_dijet->Get("h3D_all")->Clone("h_dijet_all");  

    TH3D *h_dijet_other = (TH3D*)h_dijet_c->Clone("h_dijet_other");
    h_dijet_other->Add(h_dijet_light); 

    // Fix statistical uncertainties on all histograms
    for (auto h : {h_bjet_b1, h_bjet_b1, h_bjet_other, h_bjet_all, h_dijet_b1, h_dijet_b1, h_dijet_other, h_dijet_all, h_data}){
      h->Sumw2();                                                                                                                                                        
    }           


    /*
    // Fix statistical uncertainties on all histograms
    for (auto h : {hb_bjet, hmoreb_bjet, hrest_bjet, hmc_bjet,
                   hb_dijet, hmoreb_dijet, hrest_dijet, hmc_dijet, h_data
                }) {
         //rebin_mass(h);
         //rebin_dr(h);
         h->Sumw2();
        }
    */
    

    //Get number of bin entries
    Int_t bins_pt = h_data->GetNbinsZ();
    Int_t bins_dr = h_dijet_b1->GetNbinsY();
    Int_t mb_bins = h_data->GetNbinsX();

    
    
    //Create the histograms storing the fit parameters (the dr are stored in the X, the jtpt in the Y)
    TH2D *h_sig_fraction = (TH2D *) h_data->Project3D("zy")->Clone("h_sig_fraction");
    h_sig_fraction->SetTitle("x=observable, y=jtpt, fraction / error");
    h_sig_fraction->Reset();

    TH2D *h_bkg_bb_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_bb_fraction");
    TH2D *h_bkg_rest_fraction = (TH2D *) h_sig_fraction->Clone("h_bkg_rest_fraction");

    TH2D *h_sig_fraction_error = (TH2D *) h_sig_fraction->Clone("h_sig_fraction_error");
    TH2D *h_bkg_bb_fraction_error = (TH2D *) h_sig_fraction->Clone("h_bkg_bb_fraction_error");
    TH2D *h_bkg_rest_fraction_error = (TH2D *) h_sig_fraction->Clone("h_bkg_rest_fraction_error"); 

    //Store the true parameters and errors
    TH2D *h_sig_frac_true = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true");
    TH2D *h_sig_frac_true_error = (TH2D *) h_sig_fraction->Clone("h_sig_frac_true_error");

    TH2D *h_bkg_bb_true = (TH2D *) h_sig_fraction->Clone("h_bkg_bb_true");
    TH2D *h_bkg_bb_true_error = (TH2D *) h_sig_fraction->Clone("h_bkg_bb_true_error");


    //Vector to test the convergence
    std::vector<std::pair<int, int>> non_converge_bins;


    // Fitting - loop over dr and jtpt entries
    for(Int_t ibin_pt = 1; ibin_pt <= bins_pt; ibin_pt++){
        for(Int_t ibin_dr = 1; ibin_dr <= bins_dr; ibin_dr++){

            // Make slices for data
            TH1D *h_data_mb = (TH1D *) h_data->ProjectionX(Form("h_data_mb_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);

            // Make slices for bjet
            TH1D *h_mc_bjet = (TH1D *) h_bjet_all->ProjectionX(Form("h_mc_bjet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            TH1D *h_b_bjet = (TH1D *) h_bjet_b1->ProjectionX(Form("h_b_bjet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            TH1D *h_moreb_bjet = (TH1D *) h_bjet_b2->ProjectionX(Form("h_moreb_bjet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            TH1D *h_rest_bjet = (TH1D *) h_bjet_other->ProjectionX(Form("h_rest_bjet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);

            // Make slices for dijet
            TH1D *h_mc_dijet = (TH1D *) h_dijet_all->ProjectionX(Form("h_mc_dijet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            TH1D *h_b_dijet = (TH1D *) h_dijet_b1->ProjectionX(Form("h_b_dijet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            TH1D *h_moreb_dijet = (TH1D *) h_dijet_b2->ProjectionX(Form("h_moreb_dijet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);
            TH1D *h_rest_dijet = (TH1D *) h_dijet_other->ProjectionX(Form("h_rest_dijet_%d_%d", ibin_dr, ibin_pt), ibin_dr, ibin_dr, ibin_pt, ibin_pt);

            // Calculate true fractions to be used as initial values for the fit (the true fractions are the dijet qcd ones)
            double int0 = h_b_dijet->Integral();
            double int1 = h_moreb_dijet->Integral();
            double int2 = h_rest_dijet->Integral();

            std::cout << "int0=" << int0 << std::endl;
            std::cout << "int1=" << int1 << std::endl;
            std::cout << "int2=" << int2 << std::endl;


            double sig_fraction_true = int0 / (int0 + int1 + int2); 
            double bkg_bb_fraction_true = int1 / (int0 + int1 + int2);
            std::cout << "sig_fraction_true=" << sig_fraction_true << std::endl;
            std::cout << "bkg_bb_fraction_true=" << bkg_bb_fraction_true << std::endl;

            //Add signal
            TH1D *h_sig = (TH1D*)h_b_dijet->Clone(Form("h_sig_%d_%d", ibin_dr, ibin_pt));
            Int_t h_sig_bins = h_sig->GetNbinsX();
            if(also_bjet) h_sig->Add(h_b_bjet, 1);
            h_sig->Scale(1/h_sig->Integral(1, h_sig_bins));

            //Add background from bb jets
            TH1D *h_bkg = (TH1D*)h_moreb_dijet->Clone(Form("h_bkg_%d_%d", ibin_dr, ibin_pt));
            Int_t h_bkg_bins = h_bkg->GetNbinsX();
            if (also_bjet) h_bkg->Add(h_moreb_bjet, 1);
            h_bkg->Scale(1/h_bkg->Integral(1, h_bkg_bins));

            //Add background from other non-B hadrons
            TH1D *h_rest = (TH1D*)h_rest_dijet->Clone(Form("h_rest_%d_%d", ibin_dr, ibin_pt));
            Int_t h_rest_bins = h_rest->GetNbinsX();
            if (also_bjet) h_rest->Add(h_rest_bjet, 1);
            h_rest->Scale(1/h_rest->Integral(1, h_rest_bins));


            //Add signal and rest with the correct scaling (from Lida)
            // add together light+c to sig with nominal MC ratio
            // a sig + b bb + c light = n, a+b+c=1, a+c=1-b
            // c'=c/(a+c)=(1-a-b)/(1-b)
            // a'=a/(a+c)=a/(1-b)=1-c'
            double cl_frac = (1 - sig_fraction_true - bkg_bb_fraction_true) / (1 - bkg_bb_fraction_true); // nominal
            // cl_frac = 2*cl_frac; // systematic
            // cl_frac = 0.; // systematic
            double sig_frac = 1 - cl_frac;
            h_sig->Add(h_sig, h_rest, sig_frac, 1-sig_frac);
            

            //Fitting

            // Create the observable
            Double_t min_mb = h_data_mb->GetXaxis()->GetBinLowEdge(1);
            Double_t max_mb = h_data_mb->GetXaxis()->GetBinUpEdge(mb_bins);
            RooRealVar mb(Form("mb_%d_%d", ibin_dr, ibin_pt), "mb", min_mb, max_mb); //this sets a variable able to float in the range, the initial value is set in the middle of the range
            mb.setBins(mb_bins); //Create a uniform binning under name 'name' for this variable.
     

            // Create the RooDataHist object for the observed data + templates
            RooDataHist *dh_data_mb = new RooDataHist(Form("dh_data_mb_%d_%d", ibin_dr, ibin_pt), "dh_data_mb", mb, RooFit::Import(*h_data_mb));
            RooDataHist *dh_sig_mb = new RooDataHist(Form("dh_sig_mb_%d_%d", ibin_dr, ibin_pt), "dh_sig_mb", mb, RooFit::Import(*h_sig));
            RooDataHist *dh_bkg_mb = new RooDataHist(Form("h_bkg_mb_%d_%d", ibin_dr, ibin_pt), "dh_bkg_mb", mb, RooFit::Import(*h_bkg));

            // Create the RooHistPdf objects for the template PDFs
            RooHistPdf sig_template(Form("sig_template_%d_%d", ibin_dr, ibin_pt), "sig_template", mb, *dh_sig_mb);
            RooHistPdf bkg_template(Form("bkg_template_%d_%d", ibin_dr, ibin_pt), "bkg_template", mb, *dh_bkg_mb);
            // Create list of templates
            RooArgList template_list(sig_template, bkg_template, "template_list");

            // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
            // debug 
            // sig_fraction_true = 0.5;
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d_%d", ibin_dr, ibin_pt), "sig_fraction_val",1-bkg_bb_fraction_true,0.,1.);//

            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("model0_%d_%d", ibin_dr, ibin_pt), "model0", template_list, sig_fraction_val, true); //fit the scaling when adding the single b with the more b distribution
            RooFitResult* result = model0.fitTo(*dh_data_mb, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name
                                                                                                                                                                                              //instead of sign histogram here we would put the data histogram
            Int_t status = result->status();
            /*result->Print();

            std::cout << "covariance matrix:" << std::endl;
            (result->covarianceMatrix().Print());*/

            //Check if it converged for a dr and jtpt bin
            if (status != 0) {
                std::cout << "\n\n\n\n!!!Fitting for ipt = " << ibin_pt 
                          << ", ix = " << ibin_dr 
                          << " did not converge\n\n\n\n" << std::endl;
                non_converge_bins.push_back(std::pair<int, int>(ibin_pt, ibin_dr));
                continue;
            }
    
            // Get the fitted parameter values
            double a = sig_fraction_val.getValV();
            double da = sig_fraction_val.getError();

            //Print some check
            std::cout << "a = " << a << " da = " << da << std::endl;

            Double_t p0, p1, p2, errP0, errP1, errP2;
            p0 = a*sig_frac;
            p1 = 1-a;
            p2 = a*(1-sig_frac);

            errP0 = da*sig_frac;
            errP1 = da;
            errP2 = da*(1-sig_frac);
            // std::cout << "errP0 = " << errP0 << std::endl;
            std::cout << "a'/a = " << p0/sig_fraction_true << std::endl;
            std::cout << "c'/c = " << p2/(1-sig_fraction_true-bkg_bb_fraction_true) << std::endl;
            std::cout << "c=" << (1-sig_fraction_true-bkg_bb_fraction_true) << ", c'=" << p2 << std::endl;
            std::cout << "a=" << sig_fraction_true << ", a'=" << p0 << std::endl;
            std::cout << "b=" << bkg_bb_fraction_true << ", b'=" << p1 << std::endl;

        
            //save the fit, rescaling for the signal only (not rest background)
            h_sig_fraction->SetBinContent(ibin_dr, ibin_pt, p0);
            h_sig_fraction->SetBinError(ibin_dr, ibin_pt, errP0);
            h_sig_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP0);
            h_bkg_bb_fraction->SetBinContent(ibin_dr, ibin_pt, p1);
            h_bkg_bb_fraction->SetBinError(ibin_dr, ibin_pt, errP1);
            h_bkg_bb_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP1);
            h_bkg_rest_fraction->SetBinContent(ibin_dr, ibin_pt, p2);
            h_bkg_rest_fraction->SetBinError(ibin_dr, ibin_pt, errP2);
            h_bkg_rest_fraction_error->SetBinContent(ibin_dr, ibin_pt, errP2);


            //save the true fraction
            h_sig_frac_true->SetBinContent(ibin_dr, ibin_pt, sig_fraction_true);
            h_bkg_bb_true->SetBinContent(ibin_dr, ibin_pt, bkg_bb_fraction_true);

            }
        

    }

    // Save histograms
    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");
    for (auto h : {h_data,
                   h_bjet_b1, h_bjet_b2,
                   h_bjet_other, h_bjet_all,
                   h_dijet_b1, h_dijet_b2,
                   h_dijet_other, h_dijet_all,
                  }) {
                   h->Write();
        }
    for (auto h : {h_sig_fraction, h_sig_fraction_error,
                   h_bkg_bb_fraction, h_bkg_bb_fraction_error,
                   h_bkg_rest_fraction, h_bkg_rest_fraction_error,
                   h_sig_frac_true, h_sig_frac_true_error,
                   h_bkg_bb_true
                   }) {
                    h->Write();
    }

    //See if some bins did not converge
    for (auto p : non_converge_bins) {
        std::cout << "Fit did not converge for (" << p.first << ", " << p.second << ")" << std::endl;
    }
        
    fout->Close();

}

//Draws the result of the template fit
void draw_template_fit_result(TString fout_name, TString &folder, bool &isMC, TString &pT_selection, TString &pT_selection_label, bool &also_bjet, Int_t &pt_bin){

  TString label;

  if (!isMC){
    label = "data";
  }
  else{ label = "MC dijet"; 
      }

    //Define the canvas
    TCanvas *c = new TCanvas("c_" + label, " ",170,800,700,504);
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetTitle("Fitted parameter distribution");
    c->SetLogx();

    //Get fractions for the jtpt bin pt_bin
    TFile *file = new TFile(fout_name, "read");
    //signal fraction
    TH2D *h_2D = (TH2D*)file->Get("h_sig_fraction");
    TH1D *h = (TH1D*)h_2D->ProjectionX("h", pt_bin,pt_bin);
    TH2D *htrue_2D = (TH2D*)file->Get("h_sig_frac_true");
    TH1D *htrue = (TH1D*)htrue_2D->ProjectionX("htrue", pt_bin,pt_bin);
    //background fraction
    TH2D *hbkg_bb_2D = (TH2D*)file->Get("h_bkg_bb_fraction");
    TH1D *hbkg_bb = (TH1D*)hbkg_bb_2D->ProjectionX("hbkg_bb", pt_bin,pt_bin);
    TH2D *hbkg_bb_true_2D = (TH2D*)file->Get("h_bkg_bb_true");
    TH1D *hbkg_bb_true = (TH1D*)hbkg_bb_true_2D->ProjectionX("hbkg_bb_true",pt_bin,pt_bin);


    //Draw results
    h->SetStats(0);
    h->SetTitle("");
    if(!isMC) label = "data";
    h->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h->GetYaxis()->SetRangeUser(0,1);
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->SetTitle("fitted parameter");
    h->GetYaxis()->CenterTitle(true);
    h->SetMarkerColor(kBlue);
    h->SetLineColor(kBlue);
    
    
    h->Draw("EH");
    htrue->SetMarkerColor(kRed);
    htrue->SetLineColor(kRed);
    htrue->SetLineStyle(9);
    htrue->Draw("same");

    hbkg_bb->SetStats(0);
    hbkg_bb->SetMarkerColor(kCyan);
    hbkg_bb->SetLineColor(kCyan);
    hbkg_bb->Draw("EH SAME");
    hbkg_bb_true->SetMarkerColor(kMagenta);
    hbkg_bb_true->SetLineColor(kMagenta);
    hbkg_bb_true->SetLineStyle(9);
    hbkg_bb_true->Draw("same");

    

    //Text
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.15, 0.5, "100 < p_{T} < 120 GeV");
    if(label == "data") test_info_text->DrawLatex(0.15, 0.45, "Data HighEGJet + LowEGJet");
    else test_info_text->DrawLatex(0.15, 0.45, label);
    
    if(also_bjet) test_info_text->DrawLatex(0.15, 0.4, "Template from MC dijet + bjet");
    else test_info_text->DrawLatex(0.15, 0.4, "Template from MC dijet");
    test_info_text->Draw("SAME");


    TLegend *leg = new TLegend(0.406,0.39,0.7,0.55, ""); //get position from legend.C file (lower left corner is 0,0)
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.50);
    leg->AddEntry(h, "fitted signal");
    leg->AddEntry(htrue, "MC signal");
    leg->AddEntry(hbkg_bb, "fitted more B background");
    leg->AddEntry(hbkg_bb_true, "MC more B background");
    leg->Draw("same");

    c->Print(folder + "sign_frac_result_" + label + "_" + pT_selection + "_fromtree.pdf");

}

//Draw a preliminary version of the EEC(dr) after template fit for checking
void draw_eec(TString fout_name, TString &folder, bool &isMC, TString &pT_selection, TString &pT_selection_label, bool &norm, bool &all, bool &eff_corr, Int_t &pt_bin){
    //Define the canvas
    TCanvas *c = new TCanvas("c", " ",170,800,800,504);
    c->SetLogx();
    //c->SetLogy();
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetFrameBorderMode(0);
    c->SetFrameBorderMode(0);
    c->SetTitle("EEC corrected for single-b signal fraction, all mB bins");


    // Get signal and backgrounds fractions histograms
    TFile *file = new TFile(fout_name, "read");
    TH2D *h_sigfrac = (TH2D*)file->Get("h_sig_fraction");
    TH2D *h_bkg_bb_fraction = (TH2D*)file->Get("h_bkg_bb_fraction");
    TH2D *h_bkg_rest_fraction = (TH2D*)file->Get("h_bkg_rest_fraction");

    // Get data
    TH3D *h3D = (TH3D*)file->Get("h_data");

    //Get full MC samples
    TH3D *h3D_dijet_all = (TH3D*)file->Get("h_dijet_all");
    TH2D *h2D_dijet_all = (TH2D*)h3D_dijet_all->Project3D("zy")->Clone("h2D_dijet_all");
    TH1D *h1D_dijet_all = (TH1D*)h2D_dijet_all->ProjectionX("h1D_dijet_all", pt_bin, pt_bin);

    TH3D *h3D_bjet_all = (TH3D*)file->Get("h_bjet_all");
    TH2D *h2D_bjet_all = (TH2D*)h3D_bjet_all->Project3D("zy")->Clone("h2D_bjet_all");
    TH1D *h1D_bjet_all = (TH1D*)h2D_bjet_all->ProjectionX("h1D_bjet_all", pt_bin, pt_bin);

    TH1D *h1D_all = (TH1D*)h1D_dijet_all->Clone("h1D_all");
    h1D_all->Add(h1D_bjet_all);

    
    //Get 1 b samples
    TH3D *h3D_dijet_b1 = (TH3D*)file->Get("h_dijet_b1");
    TH2D *h2D_dijet_b1 = (TH2D*)h3D_dijet_b1->Project3D("zy")->Clone("h2D_dijet_b1");
    TH1D *h1D_dijet_b1 = (TH1D*)h2D_dijet_b1->ProjectionX("h1D_dijet_b1", pt_bin, pt_bin);

    TH3D *h3D_bjet_b1 = (TH3D*)file->Get("h_bjet_b1");
    TH2D *h2D_bjet_b1 = (TH2D*)h3D_bjet_b1->Project3D("zy")->Clone("h2D_bjet_b1");
    TH1D *h1D_bjet_b1 = (TH1D*)h2D_bjet_b1->ProjectionX("h1D_bjet_b1", pt_bin, pt_bin);

    TH1D *h1D_b1 = (TH1D*)h1D_dijet_b1->Clone("h1D_b1");
    h1D_b1->Add(h1D_bjet_b1);
    
    //Get more b samples
    TH3D *h3D_dijet_b2 = (TH3D*)file->Get("h_dijet_b2");
    TH2D *h2D_dijet_b2 = (TH2D*)h3D_dijet_b2->Project3D("zy")->Clone("h2D_dijet_b2");
    TH1D *h1D_dijet_b2 = (TH1D*)h2D_dijet_b2->ProjectionX("h1D_dijet_b2", pt_bin, pt_bin);

    TH3D *h3D_bjet_b2 = (TH3D*)file->Get("h_bjet_b2");
    TH2D *h2D_bjet_b2 = (TH2D*)h3D_bjet_b2->Project3D("zy")->Clone("h2D_bjet_b2");
    TH1D *h1D_bjet_b2 = (TH1D*)h2D_bjet_b2->ProjectionX("h1D_bjet_b2", pt_bin, pt_bin);

    TH1D *h1D_b2 = (TH1D*)h1D_dijet_b2->Clone("h1D_b2");
    h1D_b2->Add(h1D_bjet_b2);
	
    //Get other background samples
    TH3D *h3D_dijet_other = (TH3D*)file->Get("h_dijet_other");
    TH2D *h2D_dijet_other = (TH2D*)h3D_dijet_other->Project3D("zy")->Clone("h2D_dijet_other");
    TH1D *h1D_dijet_other = (TH1D*)h2D_dijet_other->ProjectionX("h1D_dijet_other", pt_bin, pt_bin);

    TH3D *h3D_bjet_other = (TH3D*)file->Get("h_bjet_other");
    TH2D *h2D_bjet_other = (TH2D*)h3D_bjet_other->Project3D("zy")->Clone("h2D_bjet_other");
    TH1D *h1D_bjet_other = (TH1D*)h2D_bjet_other->ProjectionX("h1D_bjet_other", pt_bin, pt_bin);

    TH1D *h1D_other = (TH1D*)h1D_dijet_other->Clone("h1D_other");
    h1D_other->Add(h1D_bjet_other);
    
    //Get the efficiency correction factor data
    //    TFile *file3D_1b_notag = new TFile(folder + "hist_3d_gen_aggr_n1_b_dijet_notag_" + pT_selection + ".root", "read"); //ATTENTION
    TFile *file3D_1b_notag = new TFile(folder + "hist_3d_gen_aggr_n1_MC_dijet_notag_" + pT_selection + ".root", "read");
    //"hist_3d_gen_aggr_n1_b_notag_dijet_" + pT_selection + ".root", "read"); //ATTENTION
    TH3D *h3D_1b_notag = (TH3D*)file3D_1b_notag->Get("h3D_b1");
    TH2D *h1b_notag = (TH2D*)h3D_1b_notag->Project3D("zy")->Clone("h1b_notag");
 
    TFile *file3D_1b_tag = new TFile(folder + "hist_3d_gen_aggr_n1_MC_dijet_" + pT_selection + ".root", "read");
    TH3D *h3D_1b_tag = (TH3D*)file3D_1b_tag->Get("h3D_b1");
    TH2D *h1b_tag = (TH2D*)h3D_1b_tag->Project3D("zy")->Clone("h1b_tag");

    //Prepare scaled EEC(dr) histograms for data
    TH2D *h2D_scaled = (TH2D*)h3D->Project3D("zy")->Clone("h2D_scaled");
    TH2D *h2D_bkg_bb_scaled = (TH2D*)h2D_scaled->Clone("h2D_bkg_bb_scaled");
    TH2D *h2D_bkg_rest_scaled = (TH2D*)h2D_scaled->Clone("h2D_bkg_rest_scaled");

    //save efficiency correction
    TH2D *h_eff = (TH2D*)h_sigfrac->Clone("h_eff");
    h_eff->Reset();

    //Get number of bins
    int pt_bins = h3D->GetNbinsZ();
    int dr_bins = h3D->GetNbinsY();
    int mB_bins = h3D->GetNbinsX();


    //Get the tagging efficiency
    for(int bin_pt = 1; bin_pt <= pt_bins; bin_pt++){
        for(int bin_dr = 1; bin_dr <= dr_bins; bin_dr++){

            Float_t tag_eff = h1b_tag->GetBinContent(bin_dr, bin_pt)/h1b_notag->GetBinContent(bin_dr, bin_pt);
            h_eff->SetBinContent(bin_dr, bin_pt, tag_eff);
        
        }
    }

    //Multiply data by signal and background fraction and efficiency correct
    h2D_scaled->Multiply(h_sigfrac);
    h2D_bkg_bb_scaled->Multiply(h_bkg_bb_fraction);
    h2D_bkg_rest_scaled->Multiply(h_bkg_rest_fraction);

    //Correct already for b-tagging efficiency
    if(eff_corr){
        h2D_scaled->Divide(h_eff);
        h2D_bkg_bb_scaled->Divide(h_eff);
        h2D_bkg_rest_scaled->Divide(h_eff);
    }

    //Create slices in jtpt and scale them
    TH1D *heec = (TH1D*)h2D_scaled->ProjectionX("heec", pt_bin, pt_bin);

    TH1D *heec_bkg_bb = (TH1D*)h2D_bkg_bb_scaled->ProjectionX("heec_bkg_bb", pt_bin, pt_bin);

    TH1D *heec_bkg_rest = (TH1D*)h2D_bkg_rest_scaled->ProjectionX("heec_bkg_rest", pt_bin, pt_bin);

    //Normalise
    if(norm){
        heec->Scale(1/heec->Integral(), "width");
        heec_bkg_bb->Scale(1/heec_bkg_bb->Integral(), "width");
        heec_bkg_rest->Scale(1/heec_bkg_rest->Integral(), "width");
        h1D_b1->Scale(1/h1D_b1->Integral(), "width");
        h1D_b2->Scale(1/h1D_b2->Integral(), "width");
        h1D_other->Scale(1/h1D_other->Integral(), "width");
        h1D_all->Scale(1/h1D_all->Integral(), "width");
    }

    //Plot
    
    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9, "Legend");

    TString label;
    
    //Fitted histograms
    heec->SetStats(0);
    if(!isMC) label = "data";
    if(!eff_corr) heec->SetTitle("EEC distribution integrated, " + label + " (not efficiency corrected) "+ pT_selection_label + " (MC from dijet), pt bin = " + pt_bin);
    else heec->SetTitle("EEC distribution integrated, " + label + ", after efficiency correction, "+ pT_selection_label + " (MC from dijet), pt_bin = " + pt_bin);
    heec->GetXaxis()->SetTitle("\\Delta r");
    heec->GetXaxis()->CenterTitle(true);
    if(norm) heec->GetYaxis()->SetTitle("eec (norm)");
    else heec->GetYaxis()->SetTitle("eec");
    heec->GetYaxis()->CenterTitle(true);
    //if (!norm) heec->GetYaxis()->SetRangeUser(10000, 7000000);
    if(norm) heec->GetYaxis()->SetRangeUser(0, 20);
    //heec->GetXaxis()->SetRangeUser(0.001, 1);
    heec->SetMarkerStyle(20);
    heec->SetMarkerColor(kRed);
    heec->SetLineColor(kRed);
    leg->AddEntry(heec, "signal fitted (data)");
    heec->Draw("P0EHIST");
    
    //Plot all the data and MC distributions
    if(all){
        heec_bkg_bb->SetMarkerStyle(20);
        heec_bkg_bb->SetMarkerColor(kBlue);
        heec_bkg_bb->SetLineColor(kBlue);
        leg->AddEntry(heec_bkg_bb, "more B background fitted (data)");
        heec_bkg_bb->Draw("P0EHIST SAME");
        heec_bkg_rest->SetMarkerStyle(20);
        heec_bkg_rest->SetMarkerColor(kYellow);
        heec_bkg_rest->SetLineColor(kYellow);
        leg->AddEntry(heec_bkg_rest, "other background fitted (data)");
        heec_bkg_rest->Draw("P0EHIST SAME");

        //Comparisons
        h1D_b1->Draw("P0EHIST SAME");
        h1D_b1->SetMarkerStyle(24);
        h1D_b1->SetMarkerColor(kRed+1);
        h1D_b1->SetLineColor(kRed+1);
        h1D_b2->Draw("P0EHIST SAME");
        h1D_b2->SetMarkerStyle(24);
        h1D_b2->SetMarkerColor(kBlue+1);
        h1D_b2->SetLineColor(kBlue+1);
        h1D_other->Draw("P0EHIST SAME");
        h1D_other->SetMarkerStyle(24);
        h1D_other->SetMarkerColor(kOrange);
        h1D_other->SetLineColor(kOrange);
        h1D_all->Draw("P0EHIST SAME");
        h1D_all->SetMarkerStyle(3);
        h1D_all->SetMarkerColor(kMagenta);
        h1D_all->SetLineColor(kMagenta);
    

        leg->AddEntry(h1D_all, "MC not fitted");
        leg->AddEntry(h1D_b1, "1 B");
        leg->AddEntry(h1D_b2, "more B");
        leg->AddEntry(h1D_other, "other");

        leg->Draw("SAME");
    }

    TString plot_name = "eec_template_fit_check_" + label + "_" + pT_selection;
    if(eff_corr) plot_name += "_effcorr_";
    if(norm) plot_name += "norm.pdf";
    c->Print(folder + plot_name);

    TFile *file_efficiency = new TFile(folder + "file_efficiency_" + label + "_" + pT_selection + ".root", "recreate");
    h_eff->Write();
    file_efficiency->Close();
}

void do_template_fit_3d(bool isMC = false){

  TString pT_selection = "80_140";
  TString pT_selection_label = "80 < pT < 140 GeV";
    
  TString folder = "/data_CMS/cms/zaidan/test_for_code_mods/smaller_bins/";
  TString fout_name = folder + "histos_3d_from_templ_data_" + pT_selection + ".root";
  
    //Template fit using also the bjet distribution
    bool also_bjet = true;

    //Normalise the eec plot
    bool norm = true;

    //Plot also the mc eec
    bool all = true;

    //Correct for tagging efficiency
    bool eff_corr = false;

    //Select jtpt bin for plotting
    Int_t pt_bin = 2;

    //Calculate and plot eec
    do_template_fit(pT_selection, folder, fout_name, also_bjet);
    draw_template_fit_result(fout_name, folder, isMC,  pT_selection, pT_selection_label, also_bjet, pt_bin);
    draw_eec(fout_name, folder, isMC, pT_selection, pT_selection_label, norm, all, eff_corr, pt_bin);
    


} 
