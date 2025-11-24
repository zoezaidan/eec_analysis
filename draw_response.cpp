//Draws the response matrices for any dimension in the unfolding obtained from create_responses.c
#include "binning_histos.h"

//A variant of recover_eec_distribution
void project_eec_2D(TH3D* &h, TH2D* &h2D, Double_t eec_max_fill = 1000){

    for(Int_t pt_bin = 1; pt_bin <= jtpt_bins; pt_bin++){
        TH3D *h_3D = (TH3D*)h->Clone("h_3D");
        h_3D->GetZaxis()->SetRange(pt_bin, pt_bin);
        Double_t eec;
        TH2D *h_2D = (TH2D*)h_3D->Project3D("yx");
        for(Int_t ibin_dr=1; ibin_dr <= bins_dr; ibin_dr++){
            Float_t bin_err = 0;
            for(Int_t ibin_eec=1; ibin_eec <= bins_eec; ibin_eec++){
                eec = h_2D->GetYaxis()->GetBinLowEdge(ibin_eec)+eec_step/2;
                if(ibin_eec == bins_eec) eec = eec_max_fill;
                h2D->SetBinContent(ibin_dr, pt_bin, h_2D->GetBinContent(ibin_dr, ibin_eec)*eec);
                bin_err += pow(h_2D->GetBinError(ibin_dr, ibin_eec)*eec, 2);
    
            }
            h2D->SetBinError(ibin_dr, pt_bin, sqrt(bin_err));
            //std::cout << bin_err << std::endl;
        }
    }
}
//Normalise the rows so that the intregral for each column is one 
void normalise_rows(TH2D* &h){
    Int_t rows = h->GetNbinsX();
    Int_t cols = h->GetNbinsY();
    for(Int_t i = 1; i <= cols; i++){
        Float_t integral = 0;
        for(Int_t j = 1; j <= rows; j++){
            integral += h->GetBinContent(j, i);
        }
        
        for(Int_t j = 1; j <= rows; j++){
            Float_t bin_cont_before = h->GetBinContent(j, i);
            h->SetBinContent(j, i, bin_cont_before/integral);
        }
    }

    //h->Draw("colz");
}

//Draws a 2D response matrix (both for 2D and 1D unfolding where the response matrix is still filled as a "fake" 2D)
void draw_response_2D(TString &sample, TString &label, TString &folder, bool &btag, Int_t &n, bool matched){

    Float_t font_size = 26.;
    gStyle->SetPalette(57);
    gStyle->SetPaintTextFormat(".2f"); 


    TString observable = "eec";
    TString xlabel = "#Deltar";

    TString fin_name;
    if(matched) fin_name = "histos_response_2D_";
    else fin_name = "histos_response_1D_";

    if(!btag) label += "_notag"; 

    fin_name += TString(Form("n%i_", n)) + sample + "_" + label + ".root";

    std::cout << "File in: " << folder+fin_name << std::endl;
    TFile *fin = new TFile(folder+fin_name);
    TH2D *h_purity;
    TH2D *h_efficiency;

    //Get purity and efficiency
    h_purity = (TH2D *) fin->Get("h_full_purity_" + observable + "pt");
    h_efficiency = (TH2D *) fin->Get("h_full_efficiency_" + observable + "pt");
    
    RooUnfoldResponse *response = (RooUnfoldResponse *) fin->Get("response_full_" + observable + "pt");

    //Draw purity and efficiency
    TCanvas *c_purity = new TCanvas("c_purity", "purity", 800, 600);
    c_purity->SetLogx();
    h_purity->SetTitle("");
    h_purity->SetStats(0);
    h_purity->GetZaxis()->SetTitle("Reconstruction purity");
    h_purity->GetXaxis()->SetTitle("\\Delta\\mbox{r}");
    h_purity->GetXaxis()->CenterTitle(true);
    h_purity->GetYaxis()->SetTitle("p_{T}");
    h_purity->Draw("colz texte");
    //c_purity->Draw();
    c_purity->Print(folder + "purity_"+sample+"_"+label+"_purity_"+observable+".pdf");

    //TCanvas *c_efficiency = new TCanvas("c_efficiency", "efficiency", 800, 600);
    // h_efficiency->GetZaxis()->SetTitle("Reconstruction efficiency");
    // h_efficiency->Draw("colz texte");
    // c_efficiency->Draw();
    // // c_efficiency->Print(folder+sample+"_"+label+"_efficiency_"+observable+".png");
    

    // Draw the response matrix
    Int_t nbins_x = h_purity->GetNbinsX();
    Int_t nbins_pt = h_purity->GetNbinsY(); 

    TMatrixD response_matrix = response->Mresponse();
    TH2D *response_histogram = new TH2D(response_matrix);

    if (true) {
        TCanvas *c_response = new TCanvas("c_response", "response", 800, 600);
        c_response->SetRightMargin(0.2);
        // c_response->SetLogz();
        //response_histogram->SetTitle(sample + " response matrix");
        response_histogram->SetStats(0);
        response_histogram->GetXaxis()->SetTitle("Detector level weighted " + xlabel + " * p_{T}^{jet} bins");
        response_histogram->GetYaxis()->SetTitle("Particle level weighted " + xlabel + " * p_{T}^{jet} bins");
        response_histogram->GetZaxis()->SetTitle("Migration probability");
        response_histogram->GetZaxis()->SetTitleOffset(1);
        response_histogram->Draw("colz");

        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x;
            TLine *vline = new TLine(coord, 0, coord, nbins_x*nbins_pt);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_x*nbins_pt, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }

        if(matched) c_response->Print(folder + "plot_"+sample+"_"+label+"_response_2D.pdf");
        else c_response->Print(folder + "plot_"+sample+"_"+label+"_response_1D.pdf");
    }

    // DRAW SPECIFIC PT BIN (from Lida)
    if (false) {
        int ibin_pt = 2;
        TCanvas *c_response_i = new TCanvas(Form("c_response_%d",ibin_pt), "response", 800, 600);
        response_histogram->GetXaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
        response_histogram->GetYaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
        response_histogram->GetXaxis()->SetTitle("Detector level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
        response_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
        for (int i=1;i<=nbins_x+1;i++) {
            TString label = Form("%.2f",h_purity->GetXaxis()->GetBinLowEdge(i));
            response_histogram->GetYaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
            response_histogram->GetXaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
        }
        response_histogram->GetZaxis()->SetRangeUser(0, 0.7);

        response_histogram->Draw("colz");
        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x;
            TLine *vline = new TLine(coord, 0, coord, nbins_x*nbins_pt);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_x*nbins_pt, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }

    }

    // Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    for (int i=0; i < singular_values.GetNrows(); i++) {
        double val = singular_values[i];
        // std::cout << val << std::endl;
    }
    std::cout << "Largest value = " << singular_values.Max() 
              << "\nSmallest value = " << singular_values.Min()
              << "\nCondition number = " << cond_number
              << std::endl;


    //Print values from the diagonal
    /*Int_t tot_bins = nbins_pt*nbins_x*nbins_y;
    for(Int_t i = 1; i <= tot_bins; i++){
        std::cout << "diagonal entry = " << i << ", value = " << response_histogram->GetBinContent(i,i) << std::endl;
    }*/

}

//Draws the 2D response matrix from 1D unfolding as it should be understood just with nine bins
void draw_response_1D(TString &sample, TString &label, TString &folder, bool &btag, Int_t &n){

    Float_t font_size = 26.;
    gStyle->SetPalette(57);
    gStyle->SetPaintTextFormat(".2f"); 


    TString observable = "eec";
    TString xlabel = "#Deltar";
    
    TString fin_name = "histos_response_1D_";

    if(!btag) label += "_notag"; 

    fin_name += TString(Form("n%i_", n)) + sample + "_" + label + ".root";

    std::cout << "File in: " << folder+fin_name << std::endl;
    TFile *fin = new TFile(folder+fin_name);

    //Get purity/efficiency and the response matrix
    TH2D *h_purity;
    TH2D *h_efficiency;

    
    h_purity = (TH2D *) fin->Get("h_full_purity_" + observable + "pt");
    h_efficiency = (TH2D *) fin->Get("h_full_efficiency_" + observable + "pt");
    

    
    RooUnfoldResponse *response = (RooUnfoldResponse *) fin->Get("response_full_" + observable + "pt");

    // Draw the response matrix
    Int_t nbins_x = h_purity->GetNbinsX();
    Int_t nbins_pt = h_purity->GetNbinsY(); 

    TMatrixD response_matrix = response->Mresponse();
    TH2D *response_histogram_first = (TH2D*)response->Hresponse()->Clone("response_histogram_first");//new TH2D(response_matrix);//
    TH2D *response_histogram = new TH2D("response_histogram", "response_histogram", 3, 80, 140, 3, 80, 140);
    

    //Adds the bin content of the small bins to 3 big bins in jtpt
    for(Int_t pt_bin_col = 0; pt_bin_col < nbins_pt; pt_bin_col++){
        for(Int_t pt_bin_row = 0; pt_bin_row < nbins_pt; pt_bin_row++){
            Float_t entry = 0;
            for(Int_t dr_bin = 1; dr_bin <= nbins_x; dr_bin++){
                Int_t dr_shift_row = nbins_x * pt_bin_row;
                Int_t dr_shift_col = nbins_x * pt_bin_col;
                entry += response_histogram_first->GetBinContent(dr_shift_row + dr_bin, dr_shift_col + dr_bin);

            }
            response_histogram->SetBinContent(pt_bin_row+1, pt_bin_col+1, entry);
        }
    }

    //Normalise as a response matrix would be normalised
    normalise_rows(response_histogram);

    //Draw response matrix
    if (true) {
        TCanvas *c_response = new TCanvas("c_response", "response", 800, 600);
        c_response->SetRightMargin(0.2);
        // c_response->SetLogz();
        //response_histogram->SetTitle(sample + " response matrix");
        //response_histogram->Scale(1/response_histogram->Integral());
        response_histogram->SetTitle("");
        response_histogram->SetStats(0);
        response_histogram->GetXaxis()->SetTitle("Detector level  p_{T}^{jet} (GeV)");
        response_histogram->GetYaxis()->SetTitle("Particle level p_{T}^{jet} (GeV)");
        response_histogram->GetZaxis()->SetTitle("Migration probability");
        response_histogram->GetZaxis()->SetTitleOffset(1);
        response_histogram->Draw("colz");

        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x;
            TLine *vline = new TLine(coord, 0, coord, nbins_pt);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_pt, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }
        c_response->Print(folder+"plot_"+sample+"_"+label+"_reduced_response"+".pdf");
    }

    // DRAW SPECIFIC PT BIN (from Lida)
    if (false) {
        int ibin_pt = 2;
        TCanvas *c_response_i = new TCanvas(Form("c_response_%d",ibin_pt), "response", 800, 600);
        response_histogram->GetXaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
        response_histogram->GetYaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
        response_histogram->GetXaxis()->SetTitle("Detector level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
        response_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
        for (int i=1;i<=nbins_x+1;i++) {
            TString label = Form("%.2f",h_purity->GetXaxis()->GetBinLowEdge(i));
            response_histogram->GetYaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
            response_histogram->GetXaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
        }
        response_histogram->GetZaxis()->SetRangeUser(0, 0.7);

        response_histogram->Draw("colz");
        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x;
            TLine *vline = new TLine(coord, 0, coord, nbins_pt);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_pt, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }
    }

    // Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    for (int i=0; i < singular_values.GetNrows(); i++) {
        double val = singular_values[i];
        // std::cout << val << std::endl;
    }
    std::cout << "Largest value = " << singular_values.Max() 
              << "\nSmallest value = " << singular_values.Min()
              << "\nCondition number = " << cond_number
              << std::endl;


    //Print values from the diagonal
    /*Int_t tot_bins = nbins_pt*nbins_x*nbins_y;
    for(Int_t i = 1; i <= tot_bins; i++){
        std::cout << "diagonal entry = " << i << ", value = " << response_histogram->GetBinContent(i,i) << std::endl;
    }*/

}

//Draws the 3D response matrix 
void draw_response_3D(TString observable, TString &sample, TString &label, TString &folder, bool btag, Int_t n){

    Float_t font_size = 26.;
    gStyle->SetPalette(57);
    gStyle->SetPaintTextFormat(".2f"); 

    TString xlabel;

    if (observable=="eec") xlabel = "eec";
    else if (observable=="mb") xlabel = "mB";

    TString fin_name = "histos_response_3D_";

    if(!btag) label += "_notag"; 

    fin_name += TString(Form("n%i_", n)) + sample + "_" + label + ".root";

    std::cout << "File in: " << folder+fin_name << std::endl;
    TFile *fin = new TFile(folder+fin_name);

    //Define purity and efficiency histograms
    TH3D *h_purity;
    TH3D *h_efficiency;
    TH3D *h_purity_numerator;
    TH3D* h_purity_denominator;

    h_purity = (TH3D *) fin->Get("h_full_purity_" + observable + "pt");
    h_efficiency = (TH3D *) fin->Get("h_full_efficiency_" + observable + "pt");
    
    h_purity_numerator = (TH3D *) fin->Get("h_full_purity_numerator_" + observable + "pt");
    h_purity_denominator = (TH3D *) fin->Get("h_full_purity_denominator_" + observable + "pt");

    TH3D *h_efficiency_numerator = (TH3D *) fin->Get("h_full_efficiency_numerator_" + observable + "pt");
    TH3D *h_efficiency_denominator = (TH3D *) fin->Get("h_full_efficiency_denominator_" + observable + "pt");

    TH2D *h_purity_numerator_2d = (TH2D*)h_purity->Project3D("zx")->Clone("h_purity_numerator_2d");
    h_purity_numerator_2d->Reset();

    TH2D *h_purity_denominator_2d = (TH2D*)h_purity_numerator_2d->Clone("h_purity_denominator_2d");

    project_eec_2D(h_purity_numerator, h_purity_numerator_2d);
    project_eec_2D(h_purity_denominator, h_purity_denominator_2d);

    h_purity_numerator_2d->Divide(h_purity_numerator_2d, h_purity_denominator_2d, 1., 1., "b");

    TH2D *h_purity_2d = (TH2D*)h_purity_denominator_2d->Clone("h_purity_2d");


    //Get the response matrix
    RooUnfoldResponse *response = (RooUnfoldResponse *) fin->Get("response_full_" + observable + "pt");

    
    //Draw purity for the eec bin eec_bin or for all eec_bins
    Int_t eec_bin = 1;
    TCanvas *c_purity = new TCanvas("c_purity", "purity", 900, 1000, 800, 600);
    c_purity->SetLogx();
    /*c_purity->SetLogx();
    h_purity_2d->SetTitle("");
    h_purity_2d->SetStats(0);
    h_purity_2d->GetZaxis()->SetTitle("Reconstruction purity");
    h_purity_2d->GetXaxis()->SetTitle("\\Delta r");
    h_purity_2d->GetXaxis()->CenterTitle(true);
    h_purity_2d->GetYaxis()->SetTitle("p_{T}");
    h_purity_2d->Draw("colz texte");*/

    
    //Legend
    TLegend *leg = new TLegend(0.6,0.7,0.8,0.8, ""); 
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetMargin(0.20);

    
    //Projection over all eec bins
    TH1D* h_purity_numerator_1d = (TH1D *)h_purity_numerator->ProjectionX("h_purity_numerator_1d", 1, bins_eec,2,2);
    TH1D* h_purity_denominator_1d = (TH1D *)h_purity_denominator->ProjectionX("h_purity_denominator_1d", 1, bins_eec ,2,2);
    TH1D *h_purity_1d = (TH1D*)h_purity_numerator_1d->Clone("h_purity_1d");
    h_purity_1d->Divide(h_purity_denominator_1d);
    h_purity_1d->SetTitle("");
    h_purity_1d->SetStats(0);
    h_purity_1d->GetXaxis()->SetTitle("\\Delta r");
    h_purity_1d->GetXaxis()->CenterTitle(true);
    h_purity_1d->GetYaxis()->SetTitle("Purity/Efficiency");
    h_purity_1d->GetYaxis()->SetRangeUser(0., 1.);
    h_purity_1d->SetMarkerColor(kBlue);
    h_purity_1d->SetLineColor(kBlue);
    h_purity_1d->SetMarkerStyle(20);
    leg->AddEntry(h_purity_1d, "purity");
    h_purity_1d->Draw("pe1 hist same");
    //c_purity->Draw();
    //c_purity->Print(folder+"purity_"+sample+"_"+label+"_purity_"+observable+"_dr.pdf");

    //Draw efficiency in the same or a separate canvas
    //TCanvas *c_efficiency = new TCanvas("c_efficiency", "efficiency", 800, 600);TH1D* h_purity_numerator_1d = (TH1D *)h_purity_numerator->ProjectionX("h_purity_numerator_1d", 2,2,2,2);
    TH1D* h_efficiency_denominator_1d = (TH1D *)h_efficiency_denominator->ProjectionX("h_efficiency_denominator_1d", 1, bins_eec,2,2);
    TH1D* h_efficiency_numerator_1d = (TH1D *)h_efficiency_numerator->ProjectionX("h_efficiency_numerator_1d", 1, bins_eec,2,2);
    TH1D *h_efficiency_1d = (TH1D*)h_efficiency_numerator_1d->Clone("h_efficiency_1d");
    h_efficiency_1d->Divide(h_efficiency_denominator_1d);
    /*h_efficiency_1d->GetXaxis()->SetTitle("\\Delta r");
    h_efficiency_1d->GetXaxis()->CenterTitle(true);
    h_efficiency_1d->GetYaxis()->SetTitle("Efficiency");*/
    h_efficiency_1d->SetMarkerColor(kRed);
    h_efficiency_1d->SetLineColor(kRed);
    h_efficiency_1d->SetMarkerStyle(20);
    leg->AddEntry(h_efficiency_1d, "efficiency");
    h_efficiency_1d->Draw("pe1 hist same");
    // h_efficiency->GetZaxis()->SetTitle("Reconstruction efficiency");
    // h_efficiency->Draw("colz texte");
    // if (sample.Contains("herwig")) drawHeaderHerwig();
    // else drawHeaderSimulation();
    // c_efficiency->Draw();
    leg->Draw("same");
    TLatex *test_info_text = new TLatex;
    test_info_text->SetNDC();
    test_info_text->SetTextSize(0.03);
    test_info_text->DrawLatex(0.62, 0.65, "100 < p_{T} < 120 GeV");
    test_info_text->DrawLatex(0.62, 0.6, "integrated over eec");
    test_info_text->DrawLatex(0.62, 0.55, "aggregated single B jets");
    test_info_text->Draw("same");
    
    c_purity->Print(folder + "purity_efficiency_"+sample+"_"+label+"_purity_"+observable+"_dr.pdf");

    // // c_efficiency->Print("plots_an/"+sample+"_"+label+"_efficiency_"+observable+".png");

    // Draw the response matrix
    Int_t nbins_x = h_purity->GetNbinsX(); 
    Int_t nbins_y = h_purity->GetNbinsY();
    Int_t nbins_pt = h_purity->GetNbinsZ(); 

    TMatrixD response_matrix = response->Mresponse();
    TH2D *response_histogram = new TH2D(response_matrix);

    if (true) {
        TCanvas *c_response = new TCanvas("c_response", "response", 800, 600);
        c_response->SetRightMargin(0.2);
        // c_response->SetLogz();
        response_histogram->SetStats(0);
        response_histogram->GetXaxis()->SetTitle("Detector level " + xlabel + " * #Deltar * p_{T}^{jet} bins");
        response_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + " * #Deltar * p_{T}^{jet} bins");
        response_histogram->GetZaxis()->SetTitle("Migration probability");
        response_histogram->GetZaxis()->SetTitleOffset(1.);
        response_histogram->Draw("colz");

        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x * nbins_y;
            TLine *vline = new TLine(coord, 0, coord, nbins_x*nbins_pt*nbins_y);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_x*nbins_pt*nbins_y, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }

        c_response->Print(folder + "plot_"+sample+"_"+label+"_response_3D.pdf");
    }

    // DRAW SPECIFIC PT BIN (from Lida)
    if (false) {
        int ibin_pt = 2;
        TCanvas *c_response_i = new TCanvas(Form("c_response_%d",ibin_pt), "response", 800, 600);
        response_histogram->GetXaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
        response_histogram->GetYaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
        response_histogram->GetXaxis()->SetTitle("Detector level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
        response_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
        for (int i=1;i<=nbins_x+1;i++) {
            TString label = Form("%.2f",h_purity->GetXaxis()->GetBinLowEdge(i));
            response_histogram->GetYaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
            response_histogram->GetXaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
        }
        response_histogram->GetZaxis()->SetRangeUser(0, 0.7);

        response_histogram->Draw("colz");
        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x;
            TLine *vline = new TLine(coord, 0, coord, nbins_x*nbins_pt);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_x*nbins_pt, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }

    }

    // Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    for (int i=0; i < singular_values.GetNrows(); i++) {
        double val = singular_values[i];
        // std::cout << val << std::endl;
    }
    std::cout << "Largest value = " << singular_values.Max() 
              << "\nSmallest value = " << singular_values.Min()
              << "\nCondition number = " << cond_number
              << std::endl;


    //Print values from the diagonal
    /*Int_t tot_bins = nbins_pt*nbins_x*nbins_y;
    for(Int_t i = 1; i <= tot_bins; i++){
        std::cout << "diagonal entry = " << i << ", value = " << response_histogram->GetBinContent(i,i) << std::endl;
    }*/

}

void draw_response(){
  std::vector<TString> datasets{"bjet"};    //"bjet"};//, "dijet"};
    

    //Create labels
    std::vector<TString> labels_vec{"b"};//, "moreb", "other","mc"};

    TString folder = "/data_CMS/cms/zaidan/bigrerun/";

    bool btag = true;

    bool matched = true;

    Int_t n = 1; 

    TString observable = "eec";

    for(Int_t i = 0; i < datasets.size(); i++){
        for(Int_t j = 0; j < labels_vec.size(); j++){
            draw_response_3D(observable,  datasets.at(i), labels_vec.at(j), folder, btag, n);
        }

    }
}
