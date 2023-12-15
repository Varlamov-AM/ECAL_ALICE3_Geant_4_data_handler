#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"


int         main();
void        Differential_hist_create(TH1D*);
TH1D*       Convert_spectrum_with_optic_photons(TH1D*);
void        Photon_data_handler(TFile*, TChain*);
void        Electron_data_handler(TFile*, TChain*);
void        Muon_data_handler(TFile*, TChain*);
void        Pion_data_handler(TFile*, TChain*);
void        Proton_data_handler(TFile*, TChain*);
void        Kaon_data_handler(TFile*, TChain*);
void        Cluster_geometry_research(TFile*, TChain*);
void        Cluster_with_magnet(TFile*, TChain*);
void        Chi_c_data_from_Geant4_wom_magnet(TFile*, TChain*);
void        Chi_c_data_from_Geant4_wim_magnet(TFile*, TChain*);
Double_t    CB(Double_t*, Double_t*);

int Simulation_of_a_calorimeter_cell_data_handler(){
    int exeption_code  = main();
    return exeption_code;
}

int main(){

    TFile* final_result = 
        new TFile("Handled_calorimeter_energy_edeption_data.root", 
        "UPDATE");
    //     "RECREATE");

    // final_result->cd();

    // TDirectory* Electron_data = final_result->mkdir("Electron_data");
    // Electron_data->cd();

    // TDirectory* Electrons_Prime_energy_spectrum = Electron_data->mkdir("Prime_energy_spectrum");
    // TDirectory* Electrons_Smear_energy_spectrum = Electron_data->mkdir("Smear_energy_spectrum");
    // TDirectory* Electrons_Prime_optic_photons_spectrum = Electron_data->mkdir("Prime_optic_photons_spectrum");
    // TDirectory* Electrons_Smear_optic_photons_spectrum = Electron_data->mkdir("Smear_optic_photons_spectrum");
    // TDirectory* Electrons_Parametrization_of_resolution = Electron_data->mkdir("Parametrization_of_resolution");

    // final_result->cd();

    // TDirectory* Photon_data = final_result->mkdir("Photon_data");
    // Photon_data->cd();

    // TDirectory* Photons_Prime_energy_spectrum = Photon_data->mkdir("Prime_energy_spectrum");
    // TDirectory* Photons_Smear_energy_spectrum = Photon_data->mkdir("Smear_energy_spectrum");
    // TDirectory* Photons_Prime_optic_photons_spectrum = Photon_data->mkdir("Prime_optic_photons_spectrum");
    // TDirectory* Photons_Smear_optic_photons_spectrum = Photon_data->mkdir("Smear_optic_photons_spectrum");
    // TDirectory* Photons_Parametrization_of_resolution = Photon_data->mkdir("Parametrization_of_resolution");

    // final_result->cd();

    // TDirectory* Muon_data = final_result->mkdir("Muon_data");
    // Muon_data->cd();

    // TDirectory* Muons_Prime_energy_spectrum = Muon_data->mkdir("Prime_energy_spectrum");
    // TDirectory* Muons_Smear_energy_spectrum = Muon_data->mkdir("Smear_energy_spectrum");
    // TDirectory* Muons_Prime_optic_photons_spectrum = Muon_data->mkdir("Prime_optic_photons_spectrum");
    // TDirectory* Muons_Smear_optic_photons_spectrum = Muon_data->mkdir("Smear_optic_photons_spectrum");
    // TDirectory* Muons_MIP_peak_fitting = Muon_data->mkdir("MIP_peak_fitting");

    // TDirectory* Pion_data = final_result->mkdir("Pion_data");
    // Pion_data->cd();

    // TDirectory* Pions_Prime_energy_spectrum = Pion_data->mkdir("Prime_energy_spectrum");
    // TDirectory* Pions_Smear_energy_spectrum = Pion_data->mkdir("Smear_energy_spectrum");
    // TDirectory* Pions_Prime_optic_photons_spectrum = Pion_data->mkdir("Prime_optic_photons_spectrum");
    // TDirectory* Pions_Smear_optic_photons_spectrum = Pion_data->mkdir("Smear_optic_photons_spectrum");
    // TDirectory* Pions_Spectrum_parametrization = Pion_data->mkdir("Spectrum_parametrization");

    // TDirectory* Proton_data = final_result->mkdir("Proton_data");
    // Proton_data->cd();

    // TDirectory* Protons_Prime_energy_spectrum = Proton_data->mkdir("Prime_energy_spectrum");
    // TDirectory* Protons_Smear_energy_spectrum = Proton_data->mkdir("Smear_energy_spectrum");
    // TDirectory* Protons_Prime_optic_photons_spectrum = Proton_data->mkdir("Prime_optic_photons_spectrum");
    // TDirectory* Protons_Smear_optic_photons_spectrum = Proton_data->mkdir("Smear_optic_photons_spectrum");
    // TDirectory* Protons_Spectrum_parametrization = Proton_data->mkdir("Spectrum_parametrization");

    // TDirectory* Kaon_data = final_result->mkdir("Kaon_data");
    // Kaon_data->cd();

    // TDirectory* Kaons_Prime_energy_spectrum = Kaon_data->mkdir("Prime_energy_spectrum");
    // TDirectory* Kaons_Smear_energy_spectrum = Kaon_data->mkdir("Smear_energy_spectrum");
    // TDirectory* Kaons_Prime_optic_photons_spectrum = Kaon_data->mkdir("Prime_optic_photons_spectrum");
    // TDirectory* Kaons_Smear_optic_photons_spectrum = Kaon_data->mkdir("Smear_optic_photons_spectrum");
    // TDirectory* Kaons_Spectrum_parametrization = Kaon_data->mkdir("Spectrum_parametrization");

    // TDirectory* Cluster_data = final_result->mkdir("Cluster_data_for_separated_calorimeter");
    // Cluster_data->cd();

    // TDirectory* Energy_distribution_in_cells = Cluster_data->mkdir("Cluster_energy_distribution_in_cells");
    // TDirectory* Center_of_cluster_properties = Cluster_data->mkdir("Center_of_cluster_properties");
    // TDirectory* Edepted_energy_parametrization = Energy_distribution_in_cells->mkdir("Edepted_energy_parametrization");

    // TDirectory* Cluster_data_wm = final_result->mkdir("Cluster_data_for_calorimeter_with_magnet");
    // Cluster_data_wm->cd();

    // TDirectory* Energy_distribution_in_cells_wm   = Cluster_data_wm->mkdir("Cluster_energy_distribution_in_cells");
    // TDirectory* Center_of_cluster_properties_wm   = Cluster_data_wm->mkdir("Center_of_cluster_properties");
    // TDirectory* Edepted_energy_parametrization_wm = Energy_distribution_in_cells_wm->mkdir("Edepted_energy_parametrization");

    // TDirectory* Chi_c_data_from_Geant4_wom_dir = final_result->mkdir("Chi_c_data_from_Geant4_without_magnet");

    // TDirectory* Chi_c_data_from_Geant4_wim_dir = final_result->mkdir("Chi_c_data_from_Geant4_with_magnet");

    TChain* photon_edepted_energy_tree_chain = new TChain("energy_tree");
    TString photon_tree_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_fit_parameterization_28_02_2023";

    for (int i = 0; i < 16; ++i){
        photon_edepted_energy_tree_chain->
            Add(Form("%s/Job%02d/*_tree.root", photon_tree_dir_path.Data(), i + 1));
    }

    TChain* electron_edepted_energy_tree_chain = new TChain("energy_tree");
    TString electron_tree_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_electron_data_store_01_03_2023";

    for (int i = 0; i < 16; ++i){
        electron_edepted_energy_tree_chain->
            Add(Form("%s/Job%02d/*_tree.root", electron_tree_dir_path.Data(), i + 1));
    }

    TChain* muon_edepted_energy_tree_chain = new TChain("energy_tree");
    TString muon_tree_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_muon_data_store_01_03_2023";

    for (int i = 0; i < 16; ++i){
        muon_edepted_energy_tree_chain->
            Add(Form("%s/Job%02d/*_tree.root", muon_tree_dir_path.Data(), i + 1));
    }

    TChain* pion_edepted_energy_tree_chain = new TChain("energy_tree");
    TString pion_tree_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_pion_data_store_02_03_2023";

    for (int i = 0; i < 16; ++i){
        pion_edepted_energy_tree_chain->
            Add(Form("%s/Job%02d/*_tree.root", pion_tree_dir_path.Data(), i + 1));
    }

    TChain* proton_edepted_energy_tree_chain = new TChain("energy_tree");
    TString proton_tree_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_proton_data_store_04_10_2023";

    for (int i = 0; i < 16; ++i){
        proton_edepted_energy_tree_chain->
            Add(Form("%s/Job%02d/*_tree.root", proton_tree_dir_path.Data(), i + 1));
    }

    TChain* kaon_edepted_energy_tree_chain = new TChain("energy_tree");
    TString kaon_tree_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_kaon_data_store_04_10_2023";

    for (int i = 0; i < 16; ++i){
        kaon_edepted_energy_tree_chain->
            Add(Form("%s/Job%02d/*_tree.root", kaon_tree_dir_path.Data(), i + 1));
    }

    TChain* Cluster_data_for_separated_calorimeter_chain = new TChain("calorimeter_data");
    TString Cluster_data_for_separated_calorimeter_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_cluster_size_investignation_13_07_2023";
    for (int i = 0; i < 26; ++i){
        Cluster_data_for_separated_calorimeter_chain->
            Add(Form("%s/Job%02d/*.root", Cluster_data_for_separated_calorimeter_dir_path.Data(), i + 1));
    }

    TChain* Cluster_data_for_magnet_and_calorimeter_chain = new TChain("calorimeter_data");
    TString Cluster_data_for_magnet_and_calorimeter_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_batch_cluster_with_magnet_10_08_2023";
    for (int i = 0; i < 26; ++i){
        Cluster_data_for_magnet_and_calorimeter_chain->
            Add(Form("%s/Job%02d/*.root", Cluster_data_for_magnet_and_calorimeter_dir_path.Data(), i + 1));
    }

    TChain* Chi_c_data_from_Geant4_with_0dot0X0 = new TChain("calorimeter_data");
    TString Chi_c_data_from_Geant4_with_0dot0X0_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_test_chic_decays";
    for (int i = 0; i < 1; ++i){
        Chi_c_data_from_Geant4_with_0dot0X0->
            Add(Form("%s/%s", Chi_c_data_from_Geant4_with_0dot0X0_dir_path.Data(), "Calorimeter_energy_edeption_data_0dot0X0.root"));
    }

    TChain* Chi_c_data_from_Geant4_with_1dot3X0 = new TChain("calorimeter_data");
    TString Chi_c_data_from_Geant4_with_1dot3X0_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_test_chic_decays";
    for (int i = 0; i < 1; ++i){
        Chi_c_data_from_Geant4_with_1dot3X0->
            Add(Form("%s/%s", Chi_c_data_from_Geant4_with_1dot3X0_dir_path.Data(), "Calorimeter_energy_edeption_data_1dot3X0.root"));
    }

    // TString Cluster_data_for_magnet_and_calorimeter_dir_path = "/home/antony/Documents/Projects/Simulation_of_a_calorimeter_cell/build/Local_test_non_dot_beam_13_08_2023";
    // Cluster_data_for_magnet_and_calorimeter_chain->Add(Form("%s/*.root", Cluster_data_for_magnet_and_calorimeter_dir_path.Data()));
    
    // Photon_data_handler(final_result, photon_edepted_energy_tree_chain);
    // Electron_data_handler(final_result, electron_edepted_energy_tree_chain);
    // Muon_data_handler(final_result, muon_edepted_energy_tree_chain);
    // Pion_data_handler(final_result, pion_edepted_energy_tree_chain);
    // Proton_data_handler(final_result, proton_edepted_energy_tree_chain);
    // Kaon_data_handler(final_result, kaon_edepted_energy_tree_chain);
    // Cluster_geometry_research(final_result, Cluster_data_for_separated_calorimeter_chain);
    // Cluster_with_magnet(final_result, Cluster_data_for_magnet_and_calorimeter_chain);
    // Chi_c_data_from_Geant4_wom_magnet(final_result, Chi_c_data_from_Geant4_with_0dot0X0);
    Chi_c_data_from_Geant4_wim_magnet(final_result, Chi_c_data_from_Geant4_with_1dot3X0);
    

    return 0;
}

Double_t CB(Double_t * x, Double_t * par)
    {
    // Crystal Ball + pol1
    // par[0:2] - Gaussian
    // par[3:4] - Peak asymmetry parameters (Crystal Ball), usyally fixed
    Double_t m=par[1] ;
    Double_t s=par[2] ;
    Double_t n=par[3] ;
    Double_t a=par[4] ;
    Double_t dx=(x[0]-m)/s ;
    if(dx>-a)
        return par[0]*exp(-dx*dx/2.);
    else{
        Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2) ;
        Double_t B=n/TMath::Abs(a)-TMath::Abs(a) ;
        return par[0]*A*TMath::Power((B-dx),-n);
    }
}

Double_t Energy_edeption_parametrization(Double_t* x, Double_t* p){
    
    Double_t result = CB(x, p) + p[5]*TMath::Landau(x[0], p[6], p[7]);
    return result;
}


void Differential_hist_create(TH1D* input_hist){

    const int nbinsx = input_hist->GetNbinsX();

    for (int i = 0; i < nbinsx; ++i){
        input_hist->SetBinContent(i, 
            input_hist->GetBinContent(i)/(input_hist->GetBinWidth(i)));
    } 

    return;
}

void Photon_data_handler(TFile* output_data, TChain* input){

    /*
    
        With photon data we should do next steps:
    
    */
    output_data->cd();

    TDirectory* Photon_data = (TDirectory*)output_data->Get("Photon_data");
    Photon_data->cd();

    TDirectory* Prime_energy_spectrum = (TDirectory*)Photon_data->Get("Prime_energy_spectrum");
    TDirectory* Smear_energy_spectrum = (TDirectory*)Photon_data->Get("Smear_energy_spectrum");
    TDirectory* Prime_optic_photons_spectrum = (TDirectory*)Photon_data->Get("Prime_optic_photons_spectrum");
    TDirectory* Smear_optic_photons_spectrum = (TDirectory*)Photon_data->Get("Smear_optic_photons_spectrum");
    TDirectory* Parametrization_of_resolution = (TDirectory*)Photon_data->Get("Parametrization_of_resolution");

    const int n_entries = input->GetEntries();
    const int n_hist    = 16;
    const int n_params  = 20; 
    const int n_events  = 100000;

    TH1D* prime_energy_spectrum[n_hist][n_params];
    TH1D* smear_energy_spectrum[n_hist][n_params];
    TH1D* prime_optical_photons[n_hist][n_params];
    TH1D* smear_optical_photons[n_hist][n_params];

    TH1D* parametrization[n_params];

    double min_par = 0.6;
    double max_par = 1.0;

    // double opt_electron_rate[n_params] = {0.9};
    double opt_electron_rate[n_params]; // = 
    //     {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    for (int m = 0; m < n_params; ++m){
        opt_electron_rate[m] = min_par + (max_par - min_par) * m/(n_params - 1);
    }

    TH1D* parameter_depends = 
        new TH1D("Stochastic_parameter_depends", 
                 "stochastic parameter depends from opt photons rate", 
                 n_params, 
                 opt_electron_rate[0] -
                    (opt_electron_rate[n_params - 1] - opt_electron_rate[0])/(2 * n_params),
                 opt_electron_rate[n_params - 1] +
                    (opt_electron_rate[n_params - 1] - opt_electron_rate[0])/(2 * n_params));

    double edepted_energy;

    input->SetBranchAddress("Edepted_energy", &edepted_energy);

    for (int i = 0; i < n_params; ++i){

        parametrization[i] = 
            new TH1D(Form("parametrization_histogram_par%d", i), 
                     "#sigma/E", 16, 20./16 - 20./32,  20. +20./32);

        for (int j = 0; j < n_hist; ++j){
            prime_energy_spectrum[j][i] = new 
                TH1D(Form("prime_energy_spectrum_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1), 
                    Form("Prime edepted energy optical photon rate = %3.2f", opt_electron_rate[i]), 
                    2200, 0., 21.);
            smear_energy_spectrum[j][i] = new 
                TH1D(Form("smear_energy_spectrum_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1),
                    Form("Smear edepted energy optical photon rate = %3.2f", opt_electron_rate[i]), 
                    2200, 0., 21.);
            prime_optical_photons[j][i] = new 
                TH1D(Form("prime_optical_photons_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1), 
                    Form("Prime optical photon  spectrm with rate = %3.2f", opt_electron_rate[i]), 
                    22000, 0., 22000.);
            smear_optical_photons[j][i] = new 
                TH1D(Form("smear_optical_photons_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1), 
                    Form("Smear optical photon  spectrm with rate = %3.2f", opt_electron_rate[i]),
                    22000, 0., 22000.);
        }
    }

    for (int i = 0; i < n_params; ++i){
        for (int j = 0; j < n_hist; ++j){
            for (int k = n_events * j; k < n_events * (j + 1); ++k){
                input->GetEntry(k);
                prime_energy_spectrum[j][i]->Fill(edepted_energy/1000.);
                int n_pr_opt_photons = 
                    (int)edepted_energy * opt_electron_rate[i];
                prime_optical_photons[j][i]->Fill(n_pr_opt_photons);
                int n_sm_opt_photons = 
                    (int)gRandom->Poisson(n_pr_opt_photons);
                n_sm_opt_photons = 
                    (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
                smear_optical_photons[j][i]->Fill(n_sm_opt_photons);
                smear_energy_spectrum[j][i]->Fill(n_sm_opt_photons/
                                                (1000. * opt_electron_rate[i]));
                if (k % 100000 == 0){
                    std::cout << 
                        Form("i = %d in %d, j = %d in %d, k = %07d in %d\n",
                             i, n_params, j, n_hist, k, n_events * 16);
                }
            }
            smear_energy_spectrum[j][i]->Fit("gaus", "QE");
            TF1* gauss = smear_energy_spectrum[j][i]->GetFunction("gaus");
            parametrization[i]->
                SetBinContent(j + 1, gauss->GetParameter(2)/gauss->GetParameter(1));
            parametrization[i]->
                SetBinError(j + 1, 
                    TMath::Sqrt((gauss->GetParError(1)*gauss->GetParameter(2)/
                                (gauss->GetParameter(1) * gauss->GetParameter(1))) * 
                                (gauss->GetParError(1)*gauss->GetParameter(2)/
                                (gauss->GetParameter(1) * gauss->GetParameter(1))) + 
                                (gauss->GetParError(2)/
                                 gauss->GetParameter(1)) * 
                                (gauss->GetParError(2)/
                                 gauss->GetParameter(1))));
            Prime_energy_spectrum->cd();
            prime_energy_spectrum[j][i]->Write();
            Smear_energy_spectrum->cd();
            smear_energy_spectrum[j][i]->Write();
            Prime_optic_photons_spectrum->cd();
            prime_optical_photons[j][i]->Write();
            Smear_optic_photons_spectrum->cd();
            smear_optical_photons[j][i]->Write();
        }

        TF1* parametr = new TF1("parametr", "sqrt([0] * [0] / x + [1] * [1])",
                                0.001, 22.);
        parametr->SetParName(0, "b");
        parametr->SetParName(1, "c");
        parametr->SetParameter("b", 0.04);
        parametr->SetParameter("b", 0.01);
        Parametrization_of_resolution->cd();
        parametrization[i]->Fit(parametr, "E", "pe", 1., 22.);
        parametrization[i]->Write();
        parameter_depends->SetBinContent(i + 1, parametr->GetParameter("b"));
    }

    parameter_depends->Write();
    return;
}

void Electron_data_handler(TFile* output_data, TChain* input){
    /*
    
        With electron data we should do next steps:
    
    */

    output_data->cd();

    TDirectory* Electron_data = (TDirectory*)output_data->Get("Electron_data");
    Electron_data->cd();

    TDirectory* Prime_energy_spectrum = (TDirectory*)Electron_data->Get("Prime_energy_spectrum");
    TDirectory* Smear_energy_spectrum = (TDirectory*)Electron_data->Get("Smear_energy_spectrum");
    TDirectory* Prime_optic_photons_spectrum = (TDirectory*)Electron_data->Get("Prime_optic_photons_spectrum");
    TDirectory* Smear_optic_photons_spectrum = (TDirectory*)Electron_data->Get("Smear_optic_photons_spectrum");
    TDirectory* Parametrization_of_resolution = (TDirectory*)Electron_data->Get("Parametrization_of_resolution");

    const int n_entries = input->GetEntries();
    const int n_hist    = 16;
    const int n_params  = 20; 
    const int n_events  = 100000;

    TH1D* prime_energy_spectrum[n_hist][n_params];
    TH1D* smear_energy_spectrum[n_hist][n_params];
    TH1D* prime_optical_photons[n_hist][n_params];
    TH1D* smear_optical_photons[n_hist][n_params];

    TH1D* parametrization[n_params];

    double min_par = 0.6;
    double max_par = 1.0;

    // double opt_electron_rate[n_params] = {0.9};
    double opt_electron_rate[n_params]; // = 
    //     {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    for (int m = 0; m < n_params; ++m){
        opt_electron_rate[m] = min_par + (max_par - min_par) * m/(n_params - 1);
    }

    TH1D* parameter_depends = 
        new TH1D("Stochastic_parameter_depends", 
                 "stochastic parameter depends from opt photons rate", 
                 n_params, 
                 opt_electron_rate[0] -
                    (opt_electron_rate[n_params - 1] - opt_electron_rate[0])/(2 * n_params),
                 opt_electron_rate[n_params - 1] +
                    (opt_electron_rate[n_params - 1] - opt_electron_rate[0])/(2 * n_params));

    double edepted_energy;

    input->SetBranchAddress("Edepted_energy", &edepted_energy);

    for (int i = 0; i < n_params; ++i){

        parametrization[i] = 
            new TH1D(Form("parametrization_histogram_par%d", i), 
                     "#sigma/E", 16, 20./16 - 20./32,  20. +20./32);

        for (int j = 0; j < n_hist; ++j){
            prime_energy_spectrum[j][i] = new 
                TH1D(Form("prime_energy_spectrum_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1), 
                    Form("Prime edepted energy optical photon rate = %3.2f", opt_electron_rate[i]), 
                    2200, 0., 21.);
            smear_energy_spectrum[j][i] = new 
                TH1D(Form("smear_energy_spectrum_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1),
                    Form("Smear edepted energy optical photon rate = %3.2f", opt_electron_rate[i]), 
                    2200, 0., 21.);
            prime_optical_photons[j][i] = new 
                TH1D(Form("prime_optical_photons_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1), 
                    Form("Prime optical photon  spectrm with rate = %3.2f", opt_electron_rate[i]), 
                    22000, 0., 22000.);
            smear_optical_photons[j][i] = new 
                TH1D(Form("smear_optical_photons_energy_%d_opt_phot_rate_%d",
                    j + 1, 
                    i + 1), 
                    Form("Smear optical photon  spectrm with rate = %3.2f", opt_electron_rate[i]),
                    22000, 0., 22000.);
        }
    }

    for (int i = 0; i < n_params; ++i){
        for (int j = 0; j < n_hist; ++j){
            for (int k = n_events * j; k < n_events * (j + 1); ++k){
                input->GetEntry(k);
                prime_energy_spectrum[j][i]->Fill(edepted_energy/1000.);
                int n_pr_opt_photons = 
                    (int)edepted_energy * opt_electron_rate[i];
                prime_optical_photons[j][i]->Fill(n_pr_opt_photons);
                int n_sm_opt_photons = 
                    (int)gRandom->Poisson(n_pr_opt_photons);
                n_sm_opt_photons = 
                    (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
                smear_optical_photons[j][i]->Fill(n_sm_opt_photons);
                smear_energy_spectrum[j][i]->Fill(n_sm_opt_photons/
                                                (1000. * opt_electron_rate[i]));
                if (k % 100000 == 0){
                    std::cout << 
                        Form("i = %d in %d, j = %d in %d, k = %07d in %d\n",
                             i, n_params, j, n_hist, k, n_events * 16);
                }
            }
            smear_energy_spectrum[j][i]->Fit("gaus", "QE");
            TF1* gauss = smear_energy_spectrum[j][i]->GetFunction("gaus");
            parametrization[i]->
                SetBinContent(j + 1, gauss->GetParameter(2)/gauss->GetParameter(1));
            parametrization[i]->
                SetBinError(j + 1, 
                    TMath::Sqrt((gauss->GetParError(1)*gauss->GetParameter(2)/
                                (gauss->GetParameter(1) * gauss->GetParameter(1))) * 
                                (gauss->GetParError(1)*gauss->GetParameter(2)/
                                (gauss->GetParameter(1) * gauss->GetParameter(1))) + 
                                (gauss->GetParError(2)/
                                 gauss->GetParameter(1)) * 
                                (gauss->GetParError(2)/
                                 gauss->GetParameter(1))));
            Prime_energy_spectrum->cd();
            prime_energy_spectrum[j][i]->Write();
            Smear_energy_spectrum->cd();
            smear_energy_spectrum[j][i]->Write();
            Prime_optic_photons_spectrum->cd();
            prime_optical_photons[j][i]->Write();
            Smear_optic_photons_spectrum->cd();
            smear_optical_photons[j][i]->Write();
        }

        TF1* parametr = new TF1("parametr", "sqrt([0] * [0] / x + [1] * [1])",
                                0.001, 22.);
        parametr->SetParName(0, "b");
        parametr->SetParName(1, "c");
        parametr->SetParameter("b", 0.04);
        parametr->SetParameter("b", 0.01);
        // parametrization[i]->Draw("")
        parametrization[i]->Fit(parametr, "E", "pe", 1., 22.);
        Parametrization_of_resolution->cd();
        parametrization[i]->Write();
        parameter_depends->SetBinContent(i + 1, parametr->GetParameter("b"));
    }

    parameter_depends->Write();
    return;
}

void Muon_data_handler(TFile* output_data, TChain* input){

    output_data->cd();

    TDirectory* Muon_data = (TDirectory*)output_data->Get("Muon_data");
    Muon_data->cd();

    TDirectory* Prime_energy_spectrum = (TDirectory*)Muon_data->Get("Prime_energy_spectrum");
    TDirectory* Smear_energy_spectrum = (TDirectory*)Muon_data->Get("Smear_energy_spectrum");
    TDirectory* Prime_optic_photons_spectrum = (TDirectory*)Muon_data->Get("Prime_optic_photons_spectrum");
    TDirectory* Smear_optic_photons_spectrum = (TDirectory*)Muon_data->Get("Smear_optic_photons_spectrum");
    TDirectory* Muons_MIP_peak_fitting = (TDirectory*)Muon_data->Get("MIP_peak_fitting");

    const int n_entries = input->GetEntries();
    const int n_hist    = 16;
    const int n_events  = 100000;

    TH1D* prime_energy_spectrum[n_hist];
    TH1D* smear_energy_spectrum[n_hist];
    TH1D* prime_optical_photons[n_hist];
    TH1D* smear_optical_photons[n_hist];

    TH1D* landau_MPV_dep = 
        new TH1D("MIP_peak_MPV_depends_from_energy", 
                 "MIP peak MPV depends from energy", 
                 n_hist, 0.625, 1.25*n_hist + 0.625);
    
    TH1D* landau_sigma_dep = 
        new TH1D("MIP_peak_sigma_depends_from_energy", 
                 "MIP peak sigma depends from energy", 
                 n_hist, 0.625, 1.25*n_hist + 0.625);

    const double opt_electron_rate = 0.9;

    double edepted_energy;

    input->SetBranchAddress("Edepted_energy", &edepted_energy);


    for (int j = 0; j < n_hist; ++j){
        prime_energy_spectrum[j] = new 
            TH1D(Form("prime_energy_spectrum_energy_%d", j + 1), 
                Form("Prime edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        smear_energy_spectrum[j] = new 
            TH1D(Form("smear_energy_spectrum_energy_%d", j + 1),
                Form("Smear edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        prime_optical_photons[j] = new 
            TH1D(Form("prime_optical_photons_energy_%d", j + 1), 
                Form("Prime optical photon  spectrm with rate = %3.2f", opt_electron_rate), 
                22000, 0., 22000.);
        smear_optical_photons[j] = new 
            TH1D(Form("smear_optical_photons_energy_%d", j + 1), 
                Form("Smear optical photon  spectrm with rate = %3.2f", opt_electron_rate),
                22000, 0., 22000.);
    }

    
    for (int j = 0; j < n_hist; ++j){
        for (int k = n_events * j; k < n_events * (j + 1); ++k){
            input->GetEntry(k);
            prime_energy_spectrum[j]->Fill(edepted_energy/1000.);
            int n_pr_opt_photons = 
                (int)edepted_energy * opt_electron_rate;
            prime_optical_photons[j]->Fill(n_pr_opt_photons);
            int n_sm_opt_photons = 
                (int)gRandom->Poisson(n_pr_opt_photons);
            n_sm_opt_photons = 
                (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
            smear_optical_photons[j]->Fill(n_sm_opt_photons);
            smear_energy_spectrum[j]->Fill(n_sm_opt_photons/
                                            (1000. * opt_electron_rate));
            if (k % 100000 == 0){
                std::cout << 
                    Form("j = %d in %d, k = %07d in %d\n",
                          j, n_hist, k, n_events * 16);
            }
        }

        smear_energy_spectrum[j]->Fit("landau", "QE", "", 0.15, 1.);
        TF1* landau = smear_energy_spectrum[j]->GetFunction("landau");
        landau_MPV_dep->SetBinContent(j + 1, landau->GetParameter("MPV"));
        landau_sigma_dep->SetBinContent(j + 1, landau->GetParameter("Sigma"));
        landau_MPV_dep->SetBinError(j + 1, landau->GetParError(1));
        landau_sigma_dep->SetBinError(j + 1, landau->GetParError(2));
        Prime_energy_spectrum->cd();
        prime_energy_spectrum[j]->Write();
        Smear_energy_spectrum->cd();
        smear_energy_spectrum[j]->Write();
        Prime_optic_photons_spectrum->cd();
        prime_optical_photons[j]->Write();
        Smear_optic_photons_spectrum->cd();
        smear_optical_photons[j]->Write();
    }

    Muons_MIP_peak_fitting->cd();
    landau_MPV_dep->Write();
    landau_sigma_dep->Write();

    return;
}

void Pion_data_handler(TFile* output_data, TChain* input){

    output_data->cd();

    TDirectory* Pion_data = (TDirectory*)output_data->Get("Pion_data");
    Pion_data->cd();

    TDirectory* Prime_energy_spectrum = (TDirectory*)Pion_data->Get("Prime_energy_spectrum");
    TDirectory* Smear_energy_spectrum = (TDirectory*)Pion_data->Get("Smear_energy_spectrum");
    TDirectory* Prime_optic_photons_spectrum = (TDirectory*)Pion_data->Get("Prime_optic_photons_spectrum");
    TDirectory* Smear_optic_photons_spectrum = (TDirectory*)Pion_data->Get("Smear_optic_photons_spectrum");
    TDirectory* Spectrum_parametrization = (TDirectory*)Pion_data->Get("Spectrum_parametrization");

    const int n_entries = input->GetEntries();
    const int n_hist    = 16;
    const int n_events  = 100000;

    TH1D* prime_energy_spectrum[n_hist];
    TH1D* smear_energy_spectrum[n_hist];
    TH1D* prime_optical_photons[n_hist];
    TH1D* smear_optical_photons[n_hist];

    TH1D* spectrum_params[6];
    TH1D* Energy_edeption_source_ratio = new TH1D("Energy_edeption_source_ratio", 
                                                  "MIP and strong interaction ratio", n_hist, 1.25/2, 1.25 * n_hist + 1.25/2);


    for (int i = 0; i < 6; ++i){
        spectrum_params[i] = new TH1D(Form("spectrum_params_%d", i), 
                                      Form("spectrum_params_%d", i), 
                                      n_hist, 1.25/2, 1.25 * n_hist + 1.25/2);
    }


    const double opt_electron_rate = 0.9;

    double edepted_energy;

    input->SetBranchAddress("Edepted_energy", &edepted_energy);


    for (int j = 0; j < n_hist; ++j){
        prime_energy_spectrum[j] = new 
            TH1D(Form("prime_energy_spectrum_energy_%02d", j + 1), 
                Form("Prime edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        smear_energy_spectrum[j] = new 
            TH1D(Form("smear_energy_spectrum_energy_%02d", j + 1),
                Form("Smear edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        prime_optical_photons[j] = new 
            TH1D(Form("prime_optical_photons_energy_%02d", j + 1), 
                Form("Prime optical photon  spectrm with rate = %3.2f", opt_electron_rate), 
                22000, 0., 22000.);
        smear_optical_photons[j] = new 
            TH1D(Form("smear_optical_photons_energy_%02d", j + 1), 
                Form("Smear optical photon  spectrm with rate = %3.2f", opt_electron_rate),
                22000, 0., 22000.);
    }

    
    for (int j = 0; j < n_hist; ++j){
        for (int k = n_events * j; k < n_events * (j + 1); ++k){
            input->GetEntry(k);
            prime_energy_spectrum[j]->Fill(edepted_energy/1000.);
            int n_pr_opt_photons = 
                (int)edepted_energy * opt_electron_rate;
            prime_optical_photons[j]->Fill(n_pr_opt_photons);
            int n_sm_opt_photons = 
                (int)gRandom->Poisson(n_pr_opt_photons);
            n_sm_opt_photons = 
                (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
            smear_optical_photons[j]->Fill(n_sm_opt_photons);
            smear_energy_spectrum[j]->Fill(n_sm_opt_photons/
                                            (1000. * opt_electron_rate));
            if (k % 100000 == 0){
                std::cout << 
                    Form("j = %d in %d, k = %07d in %d\n",
                          j, n_hist, k, n_events * n_hist);
            }
        }

        int nbin_max = smear_energy_spectrum[j]->GetMaximumBin();
        double landau_mpv = smear_energy_spectrum[j]->GetBinCenter(nbin_max);
        double final_fit_function_par[6];

        smear_energy_spectrum[j]->Fit("landau", "QE", "", 0.5*landau_mpv, 3*landau_mpv);
        final_fit_function_par[0] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(0);
        final_fit_function_par[1] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(1);
        final_fit_function_par[2] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(2);
        smear_energy_spectrum[j]->Fit("gaus", "QE", "", 1.25 * (j + 1)/2  - 1.25 * (j + 1)/3, 1.25 * (j + 1)/2  + 1.25 * (j + 1)/2);
        final_fit_function_par[3] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(0);
        final_fit_function_par[4] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(1);
        final_fit_function_par[5] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(2);


        TF1* final_fit_function = new TF1("gaus_and_landau", "[0]*TMath::Landau(x, [1], [2]) + [3]* TMath::Gaus(x, [4], [5])", 0., 21.);
        final_fit_function->SetParameters(final_fit_function_par);
        smear_energy_spectrum[j]->Fit(final_fit_function, "QE", "", final_fit_function_par[1] - 3 * final_fit_function_par[2], 
                                                                    final_fit_function_par[4] + 3 * final_fit_function_par[5]);

        double* params = (double*)final_fit_function->GetParameters();
        double* parers = (double*)final_fit_function->GetParErrors();
        

        for (int m = 0; m < 6; ++m){
            spectrum_params[m]->SetBinContent(j + 1, params[m]);
            spectrum_params[m]->SetBinError(j + 1, parers[m]);
        }

        TF1* landau_final_part = new TF1("landau_final_part", "landau", 0., 21.);
        TF1* gaus_final_part = new TF1("gaus_final_part", "gaus", 0., 21.);

        landau_final_part->SetParameter(0, final_fit_function_par[0]);
        landau_final_part->SetParameter(1, final_fit_function_par[1]);
        landau_final_part->SetParameter(2, final_fit_function_par[2]);
        gaus_final_part->SetParameter(0, final_fit_function_par[3]);
        gaus_final_part->SetParameter(1, final_fit_function_par[4]);
        gaus_final_part->SetParameter(2, final_fit_function_par[5]);

        Energy_edeption_source_ratio->SetBinContent(j + 1, (landau_final_part->Integral(0., 21.))/(gaus_final_part->Integral(0., 21.)));

        Prime_energy_spectrum->cd();
        prime_energy_spectrum[j]->Write();
        Smear_energy_spectrum->cd();
        smear_energy_spectrum[j]->Write();
        Prime_optic_photons_spectrum->cd();
        prime_optical_photons[j]->Write();
        Smear_optic_photons_spectrum->cd();
        smear_optical_photons[j]->Write();
        Spectrum_parametrization->cd();
    }

    Energy_edeption_source_ratio->Write();

    for (int i = 0; i < 6; ++i){
        TF1* func = new TF1("fit_func", "[0] * pow(x, [1])", 1.25, 1.25 * n_hist);
        func->SetParameter(0, 1);
        func->SetParameter(1, 1);
        if (i != 1 or i != 2){
            spectrum_params[i]->Fit(func, "QE");
        } else {
            spectrum_params[i]->Fit("expo", "QE");
        }
        spectrum_params[i]->Write();
    }


    return;
}

void Proton_data_handler(TFile* output_data, TChain* input){

    output_data->cd();

    TDirectory* Proton_data = (TDirectory*)output_data->Get("Proton_data");
    Proton_data->cd();

    TDirectory* Prime_energy_spectrum = (TDirectory*)Proton_data->Get("Prime_energy_spectrum");
    TDirectory* Smear_energy_spectrum = (TDirectory*)Proton_data->Get("Smear_energy_spectrum");
    TDirectory* Prime_optic_photons_spectrum = (TDirectory*)Proton_data->Get("Prime_optic_photons_spectrum");
    TDirectory* Smear_optic_photons_spectrum = (TDirectory*)Proton_data->Get("Smear_optic_photons_spectrum");
    TDirectory* Spectrum_parametrization = (TDirectory*)Proton_data->Get("Spectrum_parametrization");

    const int n_entries = input->GetEntries();
    const int n_hist    = 16;
    const int n_events  = 100000;

    TH1D* prime_energy_spectrum[n_hist];
    TH1D* smear_energy_spectrum[n_hist];
    TH1D* prime_optical_photons[n_hist];
    TH1D* smear_optical_photons[n_hist];

    TH1D* landau_MPV_dep = 
        new TH1D("MIP_peak_MPV_depends_from_energy", 
                 "MIP peak MPV depends from energy", 
                 n_hist, 0.625, 1.25*n_hist + 0.625);
    
    TH1D* landau_sigma_dep = 
        new TH1D("MIP_peak_sigma_depends_from_energy", 
                 "MIP peak sigma depends from energy", 
                 n_hist, 0.625, 1.25*n_hist + 0.625);

    TH1D* spectrum_params[6];
    TH1D* Energy_edeption_source_ratio = new TH1D("Energy_edeption_source_ratio", 
                                                  "MIP and strong interaction ratio", n_hist, 1.25/2, 1.25 * n_hist + 1.25/2);


    for (int i = 0; i < 6; ++i){
        spectrum_params[i] = new TH1D(Form("spectrum_params_%d", i), 
                                      Form("spectrum_params_%d", i), 
                                      n_hist, 1.25/2, 1.25 * n_hist + 1.25/2);
    }


    const double opt_electron_rate = 0.9;

    double edepted_energy;

    input->SetBranchAddress("Edepted_energy", &edepted_energy);


    TCanvas* canv = new TCanvas("canv", "canv");
    canv->Divide(4, 4);


    for (int j = 0; j < n_hist; ++j){
        prime_energy_spectrum[j] = new 
            TH1D(Form("prime_energy_spectrum_energy_%02d", j + 1), 
                Form("Prime edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        smear_energy_spectrum[j] = new 
            TH1D(Form("smear_energy_spectrum_energy_%02d", j + 1),
                Form("Smear edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        prime_optical_photons[j] = new 
            TH1D(Form("prime_optical_photons_energy_%02d", j + 1), 
                Form("Prime optical photon  spectrm with rate = %3.2f", opt_electron_rate), 
                22000, 0., 22000.);
        smear_optical_photons[j] = new 
            TH1D(Form("smear_optical_photons_energy_%02d", j + 1), 
                Form("Smear optical photon  spectrm with rate = %3.2f", opt_electron_rate),
                22000, 0., 22000.);
    }

    double gaus_0[16] = {2.08456e+03, 1.10980e+03, 7.40413e+02, 5.70843e+02, 4.47007e+02, 3.54518e+02, 2.92552e+02, 2.46729e+02, 2.13662e+02, 1.86349e+02, 1.69530e+02, 1.52971e+02, 1.38974e+02, 1.24725e+02, 1.18133e+02, 1.10373e+02};
    double gaus_1[16] = {5.71469e-01, 1.05207e+00, 1.56235e+00, 2.13703e+00, 2.77405e+00, 3.27427e+00, 3.80524e+00, 4.31606e+00, 4.79937e+00, 5.35925e+00, 5.79658e+00, 6.34762e+00, 6.72503e+00, 7.31838e+00, 7.62334e+00, 8.03603e+00};
    double gaus_2[16] = {2.00289e-01, 3.73283e-01, 5.32339e-01, 6.01226e-01, 6.40397e-01, 7.93481e-01, 9.34556e-01, 1.09606e+00, 1.26166e+00, 1.38984e+00, 1.57959e+00, 1.69693e+00, 1.91486e+00, 2.04055e+00, 2.28153e+00, 2.47653e+00};
    
    for (int j = 0; j < n_hist; ++j){
        canv->cd(j + 1);
        for (int k = n_events * j; k < n_events * (j + 1); ++k){
            input->GetEntry(k);
            prime_energy_spectrum[j]->Fill(edepted_energy/1000.);
            int n_pr_opt_photons = 
                (int)edepted_energy * opt_electron_rate;
            prime_optical_photons[j]->Fill(n_pr_opt_photons);
            int n_sm_opt_photons = 
                (int)gRandom->Poisson(n_pr_opt_photons);
            n_sm_opt_photons = 
                (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
            smear_optical_photons[j]->Fill(n_sm_opt_photons);
            smear_energy_spectrum[j]->Fill(n_sm_opt_photons/
                                            (1000. * opt_electron_rate));
            if (k % 100000 == 0){
                std::cout << 
                    Form("j = %d in %d, k = %07d in %d\n",
                          j, n_hist, k, n_events * n_hist);
            }
        }

        int nbin_max = smear_energy_spectrum[j]->GetMaximumBin();
        double landau_mpv = smear_energy_spectrum[j]->GetBinCenter(nbin_max);
        double final_fit_function_par[9];

        // smear_energy_spectrum[j]->Fit("gaus", "QE", "", 1.25 * (j + 1)/2  - 1.25 * (j + 1)/3, 1.25 * (j + 1)/2  + 1.25 * (j + 1)/2);
        // final_fit_function_par[0] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(0);
        // final_fit_function_par[1] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(1);
        // final_fit_function_par[2] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(2);

        smear_energy_spectrum[j]->Fit("landau", "QE", "", 0.5*landau_mpv, 3*landau_mpv);
        final_fit_function_par[5] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(0);
        final_fit_function_par[6] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(1);
        final_fit_function_par[7] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(2);

        final_fit_function_par[8] = 0.5;

        // TF1* final_fit_function = new TF1("gaus_and_landau", "[0]*TMath::Landau(x, [1], [2]) + [3]* TMath::Gaus(x, [4], [5])", 0., 21.);
        // final_fit_function->SetParameters(final_fit_function_par);
        // smear_energy_spectrum[j]->Fit(final_fit_function, "QE", "", final_fit_function_par[1] - 3 * final_fit_function_par[2], 
        //                                                             final_fit_function_par[4] + 3 * final_fit_function_par[5]);


        TF1* final_fit_function = new TF1("gaus_and_landau", Energy_edeption_parametrization, 0., 21., 8);
        final_fit_function->FixParameter(0, gaus_0[j]);
        final_fit_function->FixParameter(1, gaus_1[j]);
        final_fit_function->FixParameter(2, gaus_2[j]);
        final_fit_function->FixParameter(3, 100);
        final_fit_function->FixParameter(4, 0.5);
        // final_fit_function->FixParameter(5, 0);
        final_fit_function->SetParameter(5, final_fit_function_par[5]);
        final_fit_function->SetParameter(6, final_fit_function_par[6]);
        final_fit_function->SetParameter(7, final_fit_function_par[7]);
        // final_fit_function->FixParameter(8, final_fit_function_par[8]);

        smear_energy_spectrum[j]->Fit(final_fit_function, "E", "",  final_fit_function_par[6] - 4 * final_fit_function_par[7],
                                                                    final_fit_function_par[6] + 20 * final_fit_function_par[7]);

        // TF1* landau_ff = new TF1("landau", "[0] * TMath::Landau(x, [1], [2])", 0., 21.);
        TF1* landau_ff = new TF1("landau", "landau", 0., 21.);
        landau_ff->SetParameter(0, (1 - final_fit_function->GetParameter(8)) * final_fit_function->GetParameter(5));
        landau_ff->SetParameter(1, final_fit_function->GetParameter(6));
        landau_ff->SetParameter(2, final_fit_function->GetParameter(7));
        // landau_ff->Draw("same");
        landau_ff->SetNpx(smear_energy_spectrum[j]->GetNbinsX());
        TH1D* landau_ff_hist = (TH1D*)landau_ff->GetHistogram();
        TH1D* clone_hist = (TH1D*)smear_energy_spectrum[j]->Clone(Form("clone_%02d", j));
        // clone_hist->Add(landau_ff_hist, -1);
        clone_hist->SetLineColor(kGreen);
        landau_ff_hist->SetLineColor(kMagenta);
        smear_energy_spectrum[j]->Draw();
        final_fit_function->Draw("same");
        // landau_ff_hist->Draw("same");
        clone_hist->Add(landau_ff_hist, -1);
        // clone_hist->Draw("same");

        
        gPad->SetLogy();
        gStyle->SetOptFit(111);
        double* params = (double*)final_fit_function->GetParameters();
        double* parers = (double*)final_fit_function->GetParErrors();
        

        for (int m = 0; m < 8; ++m){
            spectrum_params[m]->SetBinContent(j + 1, params[m]);
            spectrum_params[m]->SetBinError(j + 1, parers[m]);
        }

        TF1* landau_final_part = new TF1("landau_final_part", "landau", 0., 21.);
        TF1* gaus_final_part = new TF1("gaus_final_part", "gaus", 0., 21.);

        landau_final_part->SetParameter(0, final_fit_function_par[5]);
        landau_final_part->SetParameter(1, final_fit_function_par[6]);
        landau_final_part->SetParameter(2, final_fit_function_par[7]);

        gaus_final_part->SetParameter(0, final_fit_function_par[0]);
        gaus_final_part->SetParameter(1, final_fit_function_par[1]);
        gaus_final_part->SetParameter(2, final_fit_function_par[2]);

        Energy_edeption_source_ratio->SetBinContent(j + 1, (landau_final_part->Integral(0., 21.))/(gaus_final_part->Integral(0., 21.)));

        std::cout << "Landau pars: " << final_fit_function_par[5] << " " << final_fit_function_par[6] << " " << final_fit_function_par[7] << "\n";
        std::cout << "Landau pars: " << params[5] << " " << params[6] << " " << params[7] << "\n";

        landau_MPV_dep->SetBinContent(j + 1, params[6]);
        landau_sigma_dep->SetBinContent(j + 1, params[7]);
        landau_MPV_dep->SetBinError(j + 1, parers[6]);
        landau_sigma_dep->SetBinError(j + 1, parers[7]);

        Prime_energy_spectrum->cd();
        prime_energy_spectrum[j]->Write();
        Smear_energy_spectrum->cd();
        smear_energy_spectrum[j]->Write();
        Prime_optic_photons_spectrum->cd();
        prime_optical_photons[j]->Write();
        Smear_optic_photons_spectrum->cd();
        smear_optical_photons[j]->Write();
        Spectrum_parametrization->cd();
    }

    Spectrum_parametrization->cd();
    Energy_edeption_source_ratio->Write();
    landau_MPV_dep->Write();
    landau_sigma_dep->Write();

    // for (int i = 0; i < 6; ++i){
    //     TF1* func = new TF1("fit_func", "[0] * pow(x, [1])", 1.25, 1.25 * n_hist);
    //     func->SetParameter(0, 1);
    //     func->SetParameter(1, 1);
    //     if (i != 1 or i != 2){
    //         spectrum_params[i]->Fit(func, "QE");
    //     } else {
    //         spectrum_params[i]->Fit("expo", "QE");
    //     }
    //     spectrum_params[i]->Write();
    // }


    return;
}

void Kaon_data_handler(TFile* output_data, TChain* input){

    output_data->cd();

    TDirectory* Kaon_data = (TDirectory*)output_data->Get("Kaon_data");
    Kaon_data->cd();

    TDirectory* Prime_energy_spectrum = (TDirectory*)Kaon_data->Get("Prime_energy_spectrum");
    TDirectory* Smear_energy_spectrum = (TDirectory*)Kaon_data->Get("Smear_energy_spectrum");
    TDirectory* Prime_optic_photons_spectrum = (TDirectory*)Kaon_data->Get("Prime_optic_photons_spectrum");
    TDirectory* Smear_optic_photons_spectrum = (TDirectory*)Kaon_data->Get("Smear_optic_photons_spectrum");
    TDirectory* Spectrum_parametrization = (TDirectory*)Kaon_data->Get("Spectrum_parametrization");

    const int n_entries = input->GetEntries();
    const int n_hist    = 16;
    const int n_events  = 100000;

    TH1D* prime_energy_spectrum[n_hist];
    TH1D* smear_energy_spectrum[n_hist];
    TH1D* prime_optical_photons[n_hist];
    TH1D* smear_optical_photons[n_hist];

    TH1D* spectrum_params[6];
    TH1D* Energy_edeption_source_ratio = new TH1D("Energy_edeption_source_ratio", 
                                                  "MIP and strong interaction ratio", n_hist, 1.25/2, 1.25 * n_hist + 1.25/2);


    for (int i = 0; i < 6; ++i){
        spectrum_params[i] = new TH1D(Form("spectrum_params_%d", i), 
                                      Form("spectrum_params_%d", i), 
                                      n_hist, 1.25/2, 1.25 * n_hist + 1.25/2);
    }


    const double opt_electron_rate = 0.9;

    double edepted_energy;

    input->SetBranchAddress("Edepted_energy", &edepted_energy);


    for (int j = 0; j < n_hist; ++j){
        prime_energy_spectrum[j] = new 
            TH1D(Form("prime_energy_spectrum_energy_%02d", j + 1), 
                Form("Prime edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        smear_energy_spectrum[j] = new 
            TH1D(Form("smear_energy_spectrum_energy_%02d", j + 1),
                Form("Smear edepted energy optical photon rate = %3.2f", opt_electron_rate), 
                2200, 0., 21.);
        prime_optical_photons[j] = new 
            TH1D(Form("prime_optical_photons_energy_%02d", j + 1), 
                Form("Prime optical photon  spectrm with rate = %3.2f", opt_electron_rate), 
                22000, 0., 22000.);
        smear_optical_photons[j] = new 
            TH1D(Form("smear_optical_photons_energy_%02d", j + 1), 
                Form("Smear optical photon  spectrm with rate = %3.2f", opt_electron_rate),
                22000, 0., 22000.);
    }

    
    for (int j = 0; j < n_hist; ++j){
        for (int k = n_events * j; k < n_events * (j + 1); ++k){
            input->GetEntry(k);
            prime_energy_spectrum[j]->Fill(edepted_energy/1000.);
            int n_pr_opt_photons = 
                (int)edepted_energy * opt_electron_rate;
            prime_optical_photons[j]->Fill(n_pr_opt_photons);
            int n_sm_opt_photons = 
                (int)gRandom->Poisson(n_pr_opt_photons);
            n_sm_opt_photons = 
                (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
            smear_optical_photons[j]->Fill(n_sm_opt_photons);
            smear_energy_spectrum[j]->Fill(n_sm_opt_photons/
                                            (1000. * opt_electron_rate));
            if (k % 100000 == 0){
                std::cout << 
                    Form("j = %d in %d, k = %07d in %d\n",
                          j, n_hist, k, n_events * n_hist);
            }
        }

        int nbin_max = smear_energy_spectrum[j]->GetMaximumBin();
        double landau_mpv = smear_energy_spectrum[j]->GetBinCenter(nbin_max);
        double final_fit_function_par[6];

        smear_energy_spectrum[j]->Fit("landau", "QE", "", 0.5*landau_mpv, 3*landau_mpv);
        final_fit_function_par[0] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(0);
        final_fit_function_par[1] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(1);
        final_fit_function_par[2] = smear_energy_spectrum[j]->GetFunction("landau")->GetParameter(2);
        smear_energy_spectrum[j]->Fit("gaus", "QE", "", 1.25 * (j + 1)/2  - 1.25 * (j + 1)/3, 1.25 * (j + 1)/2  + 1.25 * (j + 1)/2);
        final_fit_function_par[3] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(0);
        final_fit_function_par[4] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(1);
        final_fit_function_par[5] = smear_energy_spectrum[j]->GetFunction("gaus")->GetParameter(2);


        TF1* final_fit_function = new TF1("gaus_and_landau", "[0]*TMath::Landau(x, [1], [2]) + [3]* TMath::Gaus(x, [4], [5])", 0., 21.);
        final_fit_function->SetParameters(final_fit_function_par);
        smear_energy_spectrum[j]->Fit(final_fit_function, "QE", "", final_fit_function_par[1] - 3 * final_fit_function_par[2], 
                                                                    final_fit_function_par[4] + 3 * final_fit_function_par[5]);

        double* params = (double*)final_fit_function->GetParameters();
        double* parers = (double*)final_fit_function->GetParErrors();
        

        // for (int m = 0; m < 6; ++m){
        //     spectrum_params[m]->SetBinContent(j + 1, params[m]);
        //     spectrum_params[m]->SetBinError(j + 1, parers[m]);
        // }

        TF1* landau_final_part = new TF1("landau_final_part", "landau", 0., 21.);
        TF1* gaus_final_part = new TF1("gaus_final_part", "gaus", 0., 21.);

        landau_final_part->SetParameter(0, final_fit_function_par[0]);
        landau_final_part->SetParameter(1, final_fit_function_par[1]);
        landau_final_part->SetParameter(2, final_fit_function_par[2]);
        gaus_final_part->SetParameter(0, final_fit_function_par[3]);
        gaus_final_part->SetParameter(1, final_fit_function_par[4]);
        gaus_final_part->SetParameter(2, final_fit_function_par[5]);

        Energy_edeption_source_ratio->SetBinContent(j + 1, (landau_final_part->Integral(0., 21.))/(gaus_final_part->Integral(0., 21.)));

        Prime_energy_spectrum->cd();
        prime_energy_spectrum[j]->Write();
        Smear_energy_spectrum->cd();
        smear_energy_spectrum[j]->Write();
        Prime_optic_photons_spectrum->cd();
        prime_optical_photons[j]->Write();
        Smear_optic_photons_spectrum->cd();
        smear_optical_photons[j]->Write();
        Spectrum_parametrization->cd();
    }

    Energy_edeption_source_ratio->Write();

    for (int i = 0; i < 6; ++i){
        TF1* func = new TF1("fit_func", "[0] * pow(x, [1])", 1.25, 1.25 * n_hist);
        func->SetParameter(0, 1);
        func->SetParameter(1, 1);
        if (i != 1 or i != 2){
            spectrum_params[i]->Fit(func, "QE");
        } else {
            spectrum_params[i]->Fit("expo", "QE");
        }
        spectrum_params[i]->Write();
    }


    return;
}

void Cluster_geometry_research(TFile* output_data, TChain* input){

    const int n_entries = input->GetEntries();
    const int n_trees   = 26;
    const int n_events  = 10000;

    TH2D* hist[n_events];

    std::cout << n_entries << "\n";
    
    const double opt_electron_rate = 1;
    // const double opt_electron_rate = 7; -- energy resolution sigma_E/E \prop 0.02/sqrt(E)
    // const double opt_electron_rate = 1.1; -- energy resolution sigma_E/E \sim to PHOS


    std::vector<double>* cells_energy = 0;
    double true_energy = 0;

    int max_cluster_size = 11;

    input->SetBranchAddress("cell_energy", &cells_energy);
    input->SetBranchAddress("initial_particle_energy", &true_energy);
    TDirectory* Cluster_data = (TDirectory*)output_data->Get("Cluster_data_for_separated_calorimeter");
    Cluster_data->cd();

    TDirectory* Energy_distribution_in_cells   = 
        (TDirectory*)Cluster_data->
            Get("Cluster_energy_distribution_in_cells");
    TDirectory* Center_of_cluster_properties   = 
        (TDirectory*)Cluster_data->
            Get("Center_of_cluster_properties");
    TDirectory* Edepted_energy_parametrization = 
        (TDirectory*)Energy_distribution_in_cells->
            Get("Edepted_energy_parametrization");

    double bin_edge[27] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 
                           1.9, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10 - 0.5, 
                           11.25 - 0.625, 12.5 - 0.625, 13.75 - 0.625, 
                           15 - 0.625, 16.25 - 0.625, 17.5 - 0.625, 
                           18.75 - 0.625, 20 - 0.625, 20. + 0.625};



    TH2D* global_cluster_capacity = 
        new TH2D("Cluster_global_capacity_for_different_cluster_size", 
                 "E_{cluster}/E_{event}(cluster size)", 
                 (max_cluster_size - 1)/2, 0, max_cluster_size - 1, 
                 26, bin_edge);

    global_cluster_capacity->SetXTitle("E_{true}, GeV");
    global_cluster_capacity->SetYTitle("Cluster size");

    TH2D* local_cluster_capacity = 
        new TH2D("Cluster_local_capacity_for_different_cluster_size", 
                 "E_{cluster}/E_{all cells}(cluster size)", 
                 (max_cluster_size - 1)/2, 0, max_cluster_size - 1, 
                 26, bin_edge);

    local_cluster_capacity->SetXTitle("E_{true}, GeV");
    global_cluster_capacity->SetYTitle("Cluster size");

    TH1D* cluster_center_of_mass_coord_sigma = 
        new TH1D("Cluster_center_of_mass_coord_sigma", "#sigma, mm", 
                 26, bin_edge);

    cluster_center_of_mass_coord_sigma->SetXTitle("E_{true}, GeV");

    TH1D* filled_cells_in_events = 
        new TH1D("filled_cells_in_events", "average number of non-empty cells in event", 
                 26, bin_edge);

    filled_cells_in_events->SetXTitle("E_{true}, GeV");

    TH2D* edepted_energy = 
        new TH2D("Edepted_energy_depends_of_true_energy", 
                 "E_{edep}(E_{true}), cluster = 3x3", 
                 26, bin_edge,
                 2000, 0., 20.);

    edepted_energy->SetXTitle("E_{true}, GeV");
    edepted_energy->SetYTitle("E_{edep}, GeV");


    TH1D* parametrization = 
        new TH1D("energy_resolution_parametrization", 
                 "photon energy resolution with 0X_{0}", 26, bin_edge);

    parametrization->SetXTitle("E_{true}, GeV");

    for (int j = 0; j < n_trees; ++j){
        
        double energy_in_hist_title = 0;
        
            if (j < 10){
                energy_in_hist_title = 0.2 * (j + 1);
            } else if (j < 18) {
                energy_in_hist_title = 1. * (j - 7);
            } else {
                energy_in_hist_title = 10 + 1.25 * (j - 17);
            }

        double average_calorimeter_energy = 0;

        std::vector<double> clusters_energy((max_cluster_size - 1)/2, 0.);
        std::vector<double> clusters_energy_error((max_cluster_size - 1)/2, 0.);
        
        long long int filled_cells = 0;

        TH1D* clusters_center = 
            new TH1D(Form("Cluster_center_x_distribution_%02d", j),
                          "Distribution of cluster center mass",
                          200, -1, 1);


        for (int k = n_events * j; k < n_events * (j + 1); ++k){
            input->GetEntry(k);
            int cal_size = static_cast<int>(sqrt(cells_energy->size()));
            hist[k] = new TH2D(Form("event%06d", k + 1), "event_calorimeter_map", 
                               cal_size, -cal_size/2, cal_size/2, 
                               cal_size, -cal_size/2, cal_size/2);
            for (int m = 0; m < cells_energy->size(); ++m){
                int n_pr_opt_photons = 
                    (int)(((*cells_energy)[m]) * opt_electron_rate);
                int n_sm_opt_photons = 
                    (int)gRandom->Poisson(n_pr_opt_photons);
                n_sm_opt_photons = 
                    (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
                (*cells_energy)[m] = n_sm_opt_photons / opt_electron_rate;
                if ((*cells_energy)[m] <= 4){
                    (*cells_energy)[m] = 0;
                } else {
                    filled_cells++;
                }      
                hist[k]->SetBinContent((m)/ cal_size + 1, 
                                       (m) % cal_size + 1, 
                                       (*cells_energy)[m]);
            }

            int max_bin_index = hist[k]->GetMaximumBin();
            int cluster_center_x = max_bin_index / cal_size - 1;
            int cluster_center_y = max_bin_index % cal_size - 1;

            double full_event_energy = hist[k]->Integral();
            double integral_error = 0;
            average_calorimeter_energy += hist[k]->Integral();

            for (int i = 0; i < 441; ++i){
                integral_error += hist[k]->GetBinError(i);
            }

            for (int p = 1; p < max_cluster_size; p+=2){
                for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                    for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                        // (q, r) -- relative_coords_for_cluster_cell;
                        clusters_energy_error[p / 2] += hist[k]->
                            GetBinError(cluster_center_x + q, 
                                        cluster_center_y + r)/(n_entries);
                        clusters_energy[p / 2] += hist[k]->
                            GetBinContent(cluster_center_x + q, 
                                          cluster_center_y + r);
                    }
                }
            }


            for (int p = 1; p < max_cluster_size; p+=2){
                double x_cm = 0;
                double weight_sum = 0;
                double cluster_energy = 0;
                for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                    for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                        // (q, r) -- relative_coords_for_cluster_cell;
                        // w -- cell weight
                        cluster_energy+=hist[k]->
                            GetBinContent(cluster_center_x + q, 
                            cluster_center_y + r);
                    }
                }
                for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                    for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                        // (q, r) -- relative_coords_for_cluster_cell;
                        // w -- cell weight
                        double w = 
                            TMath::Log(hist[k]->
                                GetBinContent(cluster_center_x + q, 
                                              cluster_center_y + r) 
                                              / 
                                              cluster_energy + 5);
                        if (w < 0) {w = 0;  }
                        x_cm += r * 22 * w;
                        weight_sum += w;
                    }
                }
                if (p == 3){
                    clusters_center->Fill(x_cm/weight_sum);
                    edepted_energy->
                        Fill(energy_in_hist_title + 0.01, cluster_energy/1000);
                }
            }

            if (k % 10000 == 0){
                std::cout << 
                    Form("j = %d in %d, k = %06d in %d\n",
                          j, n_trees, k, n_events * n_trees);
                Energy_distribution_in_cells->cd();
                hist[k]->Write();
            }
        }

        Center_of_cluster_properties->cd();
        clusters_center->Fit("gaus", "QEO");
        clusters_center->Write();
        cluster_center_of_mass_coord_sigma->
            SetBinContent(j + 1, clusters_center->
                                    GetFunction("gaus")->GetParameter(2));
        cluster_center_of_mass_coord_sigma->
            SetBinError(j + 1, clusters_center->
                                    GetFunction("gaus")->GetParError(2));


        for (int p = 1; p < max_cluster_size; p+=2){
            global_cluster_capacity->
                SetBinContent(p / 2 + 1, j + 1,  clusters_energy[p / 2]/(n_events * energy_in_hist_title * 1000));;
            global_cluster_capacity->
                SetBinError(p / 2 + 1, j + 1, clusters_energy_error[p / 2]/(1000));
            local_cluster_capacity->
                SetBinContent(p / 2 + 1, j + 1,  clusters_energy[p / 2]/(average_calorimeter_energy));;
            local_cluster_capacity->
                SetBinError(p / 2 + 1, j + 1, clusters_energy_error[p / 2]/(1000));
        }

        filled_cells_in_events->SetBinContent(j + 1, (1. * filled_cells)/n_events);
    }

    TH1D* sliced_edepted_energy[n_trees];
    Edepted_energy_parametrization->cd();

    for (int k = 0; k < n_trees; ++k){
        sliced_edepted_energy[k] = (TH1D*)edepted_energy->ProjectionY(Form("Sliced_edepted_energy_%02d", k), k + 1, k + 1);
        sliced_edepted_energy[k]->Write();
        double max_energy_in_hist = sliced_edepted_energy[k]->GetBinCenter(sliced_edepted_energy[k]->GetMaximumBin());
        sliced_edepted_energy[k]->Fit("gaus", "", "", max_energy_in_hist - 0.1, max_energy_in_hist + 0.3);
        TF1* func_test_gaus = (TF1*)sliced_edepted_energy[k]->GetFunction("gaus");
        sliced_edepted_energy[k]->Fit("gaus", "", "", func_test_gaus->GetParameter(1) - func_test_gaus->GetParameter(2), func_test_gaus->GetParameter(1) + 5 * func_test_gaus->GetParameter(2));
        TF1* func_final_gaus = (TF1*)sliced_edepted_energy[k]->GetFunction("gaus");
        parametrization->
            SetBinContent(k + 1, func_final_gaus->GetParameter(2)/func_final_gaus->GetParameter(1));
        parametrization->
            SetBinError(k + 1, 
                TMath::Sqrt((func_final_gaus->GetParError(1)*func_final_gaus->GetParameter(2)/
                            (func_final_gaus->GetParameter(1) * func_final_gaus->GetParameter(1))) * 
                            (func_final_gaus->GetParError(1)*func_final_gaus->GetParameter(2)/
                            (func_final_gaus->GetParameter(1) * func_final_gaus->GetParameter(1))) + 
                            (func_final_gaus->GetParError(2)/
                                func_final_gaus->GetParameter(1)) * 
                            (func_final_gaus->GetParError(2)/
                                func_final_gaus->GetParameter(1))));

        double energy_in_hist_title = 0;
        
            if (k < 10){
                energy_in_hist_title = 0.2 * (k + 1);
            } else if (k < 18) {
                energy_in_hist_title = 1. * (k - 7);
            } else {
                energy_in_hist_title = 10 + 1.25 * (k - 17);
            }
        // std::cout << func_final_gaus->GetParameter(1) << " " << energy_in_hist_title << " " <<  func_final_gaus->GetParameter(1)/energy_in_hist_title << "\n";
    }

    TF1* resolution = new TF1("resolution", "TMath::Sqrt([0] * [0] / x + [1] * [1])", 0., 20.);

    parametrization->Fit("resolution");
    parametrization->Write();

    Center_of_cluster_properties->cd();
    cluster_center_of_mass_coord_sigma->Write();
    
    Energy_distribution_in_cells->cd();
    filled_cells_in_events->Write();
    global_cluster_capacity->Write();
    local_cluster_capacity->Write();
    edepted_energy->Write();

    return;
}

void Cluster_with_magnet(TFile* output_data, TChain* input){

    const int n_entries = input->GetEntries();
    const int n_trees   = 26;
    const int n_events  = 10000;

    TH2D* hist[n_events];

    std::cout << n_entries << "\n";
    
    const double opt_electron_rate = 7;
    // const double opt_electron_rate = 4;

    std::vector<double>* cells_energy = 0;
    double true_energy = 0;

    int max_cluster_size = 11;

    input->SetBranchAddress("cell_energy", &cells_energy);
    input->SetBranchAddress("initial_particle_energy", &true_energy);
    TDirectory* Cluster_data = (TDirectory*)output_data->Get("Cluster_data_for_calorimeter_with_magnet");
    Cluster_data->cd();

    TDirectory* Energy_distribution_in_cells   = 
        (TDirectory*)Cluster_data->
            Get("Cluster_energy_distribution_in_cells");
    TDirectory* Center_of_cluster_properties   = 
        (TDirectory*)Cluster_data->
            Get("Center_of_cluster_properties");
    TDirectory* Edepted_energy_parametrization = 
        (TDirectory*)Energy_distribution_in_cells->
            Get("Edepted_energy_parametrization");

    double bin_edge[27] = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 
                           1.9, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10 - 0.5, 
                           11.25 - 0.625, 12.5 - 0.625, 13.75 - 0.625, 
                           15 - 0.625, 16.25 - 0.625, 17.5 - 0.625, 
                           18.75 - 0.625, 20 - 0.625, 20. + 0.625};



    TH2D* global_cluster_capacity = 
        new TH2D("Cluster_global_capacity_for_different_cluster_size", 
                 "E_{cluster}/E_{event}(cluster size)", 
                 (max_cluster_size - 1)/2, 0, max_cluster_size - 1, 
                 26, bin_edge);

    global_cluster_capacity->SetXTitle("E_{true}, GeV");
    global_cluster_capacity->SetYTitle("Cluster size");

    TH2D* local_cluster_capacity = 
        new TH2D("Cluster_local_capacity_for_different_cluster_size", 
                 "E_{cluster}/E_{all cells}(cluster size)", 
                 (max_cluster_size - 1)/2, 0, max_cluster_size - 1, 
                 26, bin_edge);

    local_cluster_capacity->SetXTitle("E_{true}, GeV");
    global_cluster_capacity->SetYTitle("Cluster size");

    TH1D* cluster_center_of_mass_coord_sigma = 
        new TH1D("Cluster_center_of_mass_coord_sigma", "#sigma, mm", 
                 26, bin_edge);

    cluster_center_of_mass_coord_sigma->SetXTitle("E_{true}, GeV");

    TH1D* filled_cells_in_events = 
        new TH1D("filled_cells_in_events", "average number of non-empty cells in event", 
                 26, bin_edge);

    filled_cells_in_events->SetXTitle("E_{true}, GeV");

    TH2D* edepted_energy = 
        new TH2D("Edepted_energy_depends_of_true_energy", 
                 "E_{edep}(E_{true}), cluster = 3x3", 
                 26, bin_edge,
                 2000, 0., 20.);

    edepted_energy->SetXTitle("E_{true}, GeV");
    edepted_energy->SetYTitle("E_{edep}, GeV");


    TH1D* parametrization = 
        new TH1D("energy_resolution_parametrization", 
                 "photon energy resolution with 1.2X_{0}", 26, bin_edge);

    parametrization->SetXTitle("E_{true}, GeV");


    for (int j = 0; j < n_trees; ++j){
        
        double energy_in_hist_title = 0;
        
            if (j < 10){
                energy_in_hist_title = 0.2 * (j + 1);
            } else if (j < 18) {
                energy_in_hist_title = 1. * (j - 7);
            } else {
                energy_in_hist_title = 10 + 1.25 * (j - 17);
            }

        double average_calorimeter_energy = 0;

        std::vector<double> clusters_energy((max_cluster_size - 1)/2, 0.);
        std::vector<double> clusters_energy_error((max_cluster_size - 1)/2, 0.);

        long long int filled_cells = 0;
        
        TH1D* clusters_center = 
            new TH1D(Form("Cluster_center_x_distribution_%02d", j),
                          "Distribution of cluster center mass",
                          200, -1, 1);


        for (int k = n_events * j; k < n_events * (j + 1); ++k){
            input->GetEntry(k);
            int cal_size = static_cast<int>(sqrt(cells_energy->size()));
            hist[k] = new TH2D(Form("event%06d", k + 1), "event_calorimeter_map", 
                               cal_size, -cal_size/2, cal_size/2, 
                               cal_size, -cal_size/2, cal_size/2);
            for (int m = 0; m < cells_energy->size(); ++m){
                int n_pr_opt_photons = 
                    (int)(((*cells_energy)[m]) * opt_electron_rate);
                int n_sm_opt_photons = 
                    (int)gRandom->Poisson(n_pr_opt_photons);
                n_sm_opt_photons = 
                    (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
                (*cells_energy)[m] = n_sm_opt_photons / opt_electron_rate;
                if ((*cells_energy)[m] <= 4){
                    (*cells_energy)[m] = 0;
                } else {
                    filled_cells++;
                }
                hist[k]->SetBinContent((m)/ cal_size + 1, 
                                       (m) % cal_size + 1, 
                                       (*cells_energy)[m]);
            }

            int max_bin_index = hist[k]->GetMaximumBin();
            int cluster_center_x = max_bin_index / cal_size - 1;
            int cluster_center_y = max_bin_index % cal_size - 1;

            double full_event_energy = hist[k]->Integral();
            double integral_error = 0;
            average_calorimeter_energy += hist[k]->Integral();

            for (int i = 0; i < 441; ++i){
                integral_error += hist[k]->GetBinError(i);
            }

            for (int p = 1; p < max_cluster_size; p+=2){
                for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                    for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                        // (q, r) -- relative_coords_for_cluster_cell;
                        clusters_energy_error[p / 2] += hist[k]->
                            GetBinError(cluster_center_x + q, 
                                        cluster_center_y + r)/(n_entries);
                        clusters_energy[p / 2] += hist[k]->
                            GetBinContent(cluster_center_x + q, 
                                          cluster_center_y + r);
                    }
                }
            }

            for (int p = 1; p < max_cluster_size; p+=2){
                double x_cm = 0;
                double weight_sum = 0;
                double cluster_energy = 0;
                for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                    for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                        // (q, r) -- relative_coords_for_cluster_cell;
                        // w -- cell weight
                        cluster_energy+=hist[k]->
                            GetBinContent(cluster_center_x + q, 
                            cluster_center_y + r);
                    }
                }
                for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                    for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                        // (q, r) -- relative_coords_for_cluster_cell;
                        // w -- cell weight
                        double w = 
                            TMath::Log(hist[k]->
                                GetBinContent(cluster_center_x + q, 
                                              cluster_center_y + r) 
                                              / 
                                              cluster_energy + 5);
                        if (w < 0) {w = 0;  }
                        x_cm += r * 22 * w;
                        weight_sum += w;
                    }
                }
                if (p == 3){
                    clusters_center->Fill(x_cm/weight_sum);
                    edepted_energy->
                        Fill(energy_in_hist_title + 0.01, cluster_energy/1000);
                }
            }

            if (k % 10000 == 0){
                std::cout << 
                    Form("j = %d in %d, k = %06d in %d\n",
                          j, n_trees, k, n_events * n_trees);
                Energy_distribution_in_cells->cd();
                hist[k]->Write();
            }
        }

        Center_of_cluster_properties->cd();
        clusters_center->Fit("gaus", "QEO");
        clusters_center->Write();
        cluster_center_of_mass_coord_sigma->
            SetBinContent(j + 1, clusters_center->
                                    GetFunction("gaus")->GetParameter(2));
        cluster_center_of_mass_coord_sigma->
            SetBinError(j + 1, clusters_center->
                                    GetFunction("gaus")->GetParError(2));


        for (int p = 1; p < max_cluster_size; p+=2){
            global_cluster_capacity->
                SetBinContent(p / 2 + 1, j + 1,  clusters_energy[p / 2]/(n_events * energy_in_hist_title * 1000));;
            global_cluster_capacity->
                SetBinError(p / 2 + 1, j + 1, clusters_energy_error[p / 2]/(1000));
            local_cluster_capacity->
                SetBinContent(p / 2 + 1, j + 1,  clusters_energy[p / 2]/(average_calorimeter_energy));;
            local_cluster_capacity->
                SetBinError(p / 2 + 1, j + 1, clusters_energy_error[p / 2]/(1000));
        }

        filled_cells_in_events->SetBinContent(j + 1, (1. * filled_cells)/n_events);

    }

    TH1D* sliced_edepted_energy[n_trees];
    Edepted_energy_parametrization->cd();

    for (int k = 0; k < n_trees; ++k){
        sliced_edepted_energy[k] = (TH1D*)edepted_energy->ProjectionY(Form("Sliced_edepted_energy_%02d", k), k + 1, k + 1);
        sliced_edepted_energy[k]->Write();
        double max_energy_in_hist = sliced_edepted_energy[k]->GetBinCenter(sliced_edepted_energy[k]->GetMaximumBin());
        sliced_edepted_energy[k]->Fit("gaus", "", "", max_energy_in_hist - 0.1, max_energy_in_hist + 0.3);
        TF1* func_test_gaus = (TF1*)sliced_edepted_energy[k]->GetFunction("gaus");
        sliced_edepted_energy[k]->Fit("gaus", "", "", func_test_gaus->GetParameter(1) - func_test_gaus->GetParameter(2), func_test_gaus->GetParameter(1) + 5 * func_test_gaus->GetParameter(2));
        TF1* func_final_gaus = (TF1*)sliced_edepted_energy[k]->GetFunction("gaus");
        parametrization->
            SetBinContent(k + 1, func_final_gaus->GetParameter(2)/func_final_gaus->GetParameter(1));
        parametrization->
            SetBinError(k + 1, 
                TMath::Sqrt((func_final_gaus->GetParError(1)*func_final_gaus->GetParameter(2)/
                            (func_final_gaus->GetParameter(1) * func_final_gaus->GetParameter(1))) * 
                            (func_final_gaus->GetParError(1)*func_final_gaus->GetParameter(2)/
                            (func_final_gaus->GetParameter(1) * func_final_gaus->GetParameter(1))) + 
                            (func_final_gaus->GetParError(2)/
                                func_final_gaus->GetParameter(1)) * 
                            (func_final_gaus->GetParError(2)/
                                func_final_gaus->GetParameter(1))));

        double energy_in_hist_title = 0;
        
            if (k < 10){
                energy_in_hist_title = 0.2 * (k + 1);
            } else if (k < 18) {
                energy_in_hist_title = 1. * (k - 7);
            } else {
                energy_in_hist_title = 10 + 1.25 * (k - 17);
            }
        // std::cout << func_final_gaus->GetParameter(1) << " " << energy_in_hist_title << " " <<  func_final_gaus->GetParameter(1)/energy_in_hist_title << "\n";
    }

    TF1* resolution = new TF1("resolution", "TMath::Sqrt([0] * [0] / x + [1] * [1])", 0., 20.);

    parametrization->Fit("resolution", "", "", 0., 20.);
    parametrization->Write();

    Center_of_cluster_properties->cd();
    cluster_center_of_mass_coord_sigma->Write();
    
    Energy_distribution_in_cells->cd();
    filled_cells_in_events->Write();
    global_cluster_capacity->Write();
    local_cluster_capacity->Write();
    edepted_energy->Write();

    return;
}

void Chi_c_data_from_Geant4_wom_magnet(TFile* output_data, TChain* input){

    double chic0_branching = 0.00083594;
    double chic1_branching = 0.0204805;
    double chic2_branching = 0.0113449;

    std::vector<double>* cells_energy = 0;
    double true_energy = 0;

    double electron_px = 0;
    double electron_py = 0;
    double electron_pz = 0;
    double electron_p0 = 0;

    double positron_px = 0;
    double positron_py = 0;
    double positron_pz = 0;
    double positron_p0 = 0;

    double gamma_px = 0;
    double gamma_py = 0;
    double gamma_pz = 0;
    double gamma_p0 = 0;

    double chi_c_mass_before_rotating = 0;
    double chi_c_pt_before_rotating = 0;

    double branchings = 0;

    int max_cluster_size = 11;

    input->SetBranchAddress("cell_energy", &cells_energy);
    input->SetBranchAddress("electron_px", &electron_px);
    input->SetBranchAddress("electron_py", &electron_py);
    input->SetBranchAddress("electron_pz", &electron_pz);
    input->SetBranchAddress("electron_p0", &electron_p0);
    input->SetBranchAddress("positron_px", &positron_px);
    input->SetBranchAddress("positron_py", &positron_py);
    input->SetBranchAddress("positron_pz", &positron_pz);
    input->SetBranchAddress("positron_p0", &positron_p0);
    input->SetBranchAddress("gamma_px", &gamma_px);
    input->SetBranchAddress("gamma_py", &gamma_py);
    input->SetBranchAddress("gamma_pz", &gamma_pz);
    input->SetBranchAddress("gamma_p0", &gamma_p0);
    input->SetBranchAddress("chi_c_mass_before_rotating", &chi_c_mass_before_rotating);
    input->SetBranchAddress("chi_c_pt_before_rotating", &chi_c_pt_before_rotating);
    input->SetBranchAddress("event_branching", &branchings);
    input->SetBranchAddress("initial_particle_energy", &true_energy);

    output_data->cd();

    TDirectory* Chi_c_data_from_Geant4_wim_dir = (TDirectory*)output_data->GetDirectory("Chi_c_data_from_Geant4_without_magnet");
    Chi_c_data_from_Geant4_wim_dir->cd();

    Int_t min_opt_electron_rate = 7;
    Int_t max_opt_electron_rate = 7;

    // const double opt_electron_rate = 7; -- energy resolution sigma_E/E \prop 0.02/sqrt(E)
    // const double opt_electron_rate = 1.1; -- energy resolution sigma_E/E \sim to PHOS

    TDirectory* opt_electron_rate_dir[max_opt_electron_rate - min_opt_electron_rate + 1];

    TDirectory* Inv_mass_spectrum_before_simulation[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Inv_mass_spectrum_after_simulation[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Inv_mass_sep_chic_after_simulation[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Energy_distribution_in_cells_chic[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Chic_parameters[max_opt_electron_rate - min_opt_electron_rate + 1];

    for (int opt_electron_rate = min_opt_electron_rate; opt_electron_rate <= max_opt_electron_rate; ++opt_electron_rate){

        opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate] = (TDirectory*)Chi_c_data_from_Geant4_wim_dir->mkdir(Form("Opt_elctron_rate_%1d_per_MeV", opt_electron_rate));        
        opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->cd();

        Inv_mass_spectrum_before_simulation[opt_electron_rate - min_opt_electron_rate] = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Inv_mass_spectrum_before_simulation");
        Inv_mass_spectrum_after_simulation[opt_electron_rate - min_opt_electron_rate]  = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Inv_mass_spectrum_after_simulation");
        Inv_mass_sep_chic_after_simulation[opt_electron_rate - min_opt_electron_rate]  = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Inv_mass_sep_chic_after_simulation");
        Energy_distribution_in_cells_chic[opt_electron_rate - min_opt_electron_rate]   = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Cluster_energy_distribution_in_cells");
        Chic_parameters[opt_electron_rate - min_opt_electron_rate]                     = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Chic_parameters_pt_depends");

        const int n_entries = input->GetEntries();

        TH2D** hist = new TH2D*[n_entries];

        std::cout << n_entries << "\n";

        std::vector<double> clusters_energy((max_cluster_size - 1)/2, 0.);
        std::vector<double> clusters_energy_error((max_cluster_size - 1)/2, 0.);

        long long int filled_cells = 0;

        int n_empty = 0;

        TH2D* edepted_energy = new TH2D("edepted_energy", "edepted_energy", 500, 0., 5., 500, 0., 5.);

        TH2D* Jpsi_mass = new TH2D("Jpmin_opt_electron_ratesi_mass", "Jpsi_mass", 100, 2.8, 3.4, 10, 0., 10.);
        TH2D* Jpsi_mass_with_elec_posi_resolution = new TH2D("Jpsi_mass_with_elec_posi_resolution", "Jpsi_mass_with_elec_posi_resolution", 100, 2.8, 3.4, 10, 0., 10.);

        TH2D* chic_mass_before = new TH2D("chi_c_mass_before", "chi_c_mass_before", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_before_slice[chic_mass_before->GetNbinsY()];

        TH2D* chic_mass_with_gamma_resolution = new TH2D("chic_mass_with_gamma_resolution", "chic_mass_with_gamma_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_with_gamma_resolution_slice[chic_mass_with_gamma_resolution->GetNbinsY()];

        TH2D* chic_mass_with_gamma_elec_posi_resolution = new TH2D("chic_mass_with_gamma_elec_posi_resolution", "chic_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_with_gamma_elec_posi_resolution_slice[chic_mass_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH2D* chic_mass_diff_with_gamma_elec_posi_resolution = new TH2D("chic_mass_diff_with_gamma_elec_posi_resolution", "chic_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_diff_with_gamma_elec_posi_resolution_slice[chic_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH2D* chic0_mass_with_gamma_elec_posi_resolution = new TH2D("separate_chic0_mass_with_gamma_elec_posi_resolution", "separate_chic0_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic1_mass_with_gamma_elec_posi_resolution = new TH2D("separate_chic1_mass_with_gamma_elec_posi_resolution", "separate_chic1_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic2_mass_with_gamma_elec_posi_resolution = new TH2D("separate_chic2_mass_with_gamma_elec_posi_resolution", "separate_chic2_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic0_mass_with_gamma_elec_posi_resolution_slice[chic0_mass_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic1_mass_with_gamma_elec_posi_resolution_slice[chic1_mass_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic2_mass_with_gamma_elec_posi_resolution_slice[chic2_mass_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH2D* chic0_mass_diff_with_gamma_elec_posi_resolution = new TH2D("separate_chic0_mass_diff_with_gamma_elec_posi_resolution", "separate_chic0_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic1_mass_diff_with_gamma_elec_posi_resolution = new TH2D("separate_chic1_mass_diff_with_gamma_elec_posi_resolution", "separate_chic1_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic2_mass_diff_with_gamma_elec_posi_resolution = new TH2D("separate_chic2_mass_diff_with_gamma_elec_posi_resolution", "separate_chic2_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic0_mass_diff_with_gamma_elec_posi_resolution_slice[chic0_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic1_mass_diff_with_gamma_elec_posi_resolution_slice[chic1_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic2_mass_diff_with_gamma_elec_posi_resolution_slice[chic2_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH1D* mass_chic0_pt_depends = new TH1D("mass_chic0_pt_depends", "mass_chic0_pt_depends", 10, -0.5, 10.5);
        TH1D* mass_chic1_pt_depends = new TH1D("mass_chic1_pt_depends", "mass_chic1_pt_depends", 10, -0.5, 10.5);
        TH1D* mass_chic2_pt_depends = new TH1D("mass_chic2_pt_depends", "mass_chic2_pt_depends", 10, -0.5, 10.5);

        TH1D* sigma_chic0_pt_depends = new TH1D("sigma_chic0_pt_depends", "sigma_chic0_pt_depends", 10, -0.5, 10.5);
        TH1D* sigma_chic1_pt_depends = new TH1D("sigma_chic1_pt_depends", "sigma_chic1_pt_depends", 10, -0.5, 10.5);
        TH1D* sigma_chic2_pt_depends = new TH1D("sigma_chic2_pt_depends", "sigma_chic2_pt_depends", 10, -0.5, 10.5);

        TH1D* efficiency_chic0_pt_depends = new TH1D("efficiency_chic0_pt_depends", "efficiency_chic0_pt_depends", 10, -0.5, 10.5);
        TH1D* efficiency_chic1_pt_depends = new TH1D("efficiency_chic1_pt_depends", "efficiency_chic1_pt_depends", 10, -0.5, 10.5);
        TH1D* efficiency_chic2_pt_depends = new TH1D("efficiency_chic2_pt_depends", "efficiency_chic2_pt_depends", 10, -0.5, 10.5);

        for (int k = 0; k < n_entries; ++k){
            input->GetEntry(k);

            if (k % 10000 == 0){
                std::cout << 
                    Form("k = %06d in %d\n", k, n_entries);
            }

            int cal_size = static_cast<int>(sqrt(cells_energy->size()));
            hist[k] = new TH2D(Form("event%06d", k + 1), "event_calorimeter_map", 
                                cal_size, -cal_size/2, cal_size/2, 
                                cal_size, -cal_size/2, cal_size/2);
            for (int m = 0; m < cells_energy->size(); ++m){
                int n_pr_opt_photons = 
                    (int)(((*cells_energy)[m]) * opt_electron_rate);
                int n_sm_opt_photons = 
                    (int)gRandom->Poisson(n_pr_opt_photons);
                n_sm_opt_photons = 
                    (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
                (*cells_energy)[m] = n_sm_opt_photons / opt_electron_rate;
                if ((*cells_energy)[m] <= 4){
                    (*cells_energy)[m] = 0;
                } else {
                    filled_cells++;
                }
                hist[k]->SetBinContent((m)/ cal_size + 1, 
                                        (m) % cal_size + 1, 
                                        (*cells_energy)[m]);
            }

            int max_bin_index = hist[k]->GetMaximumBin();
            int cluster_center_x = max_bin_index / cal_size - 1;
            int cluster_center_y = max_bin_index % cal_size - 1;

            double full_event_energy = hist[k]->Integral();
            if (full_event_energy <= 0.1){
                ++n_empty;
                continue;
            }

            double integral_error = 0;

            for (int i = 0; i < 441; ++i){
                integral_error += hist[k]->GetBinError(i);
            }

            int p = 3;
            double x_cm = 0;
            double weight_sum = 0;
            double cluster_energy = 0;
            for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                    // (q, r) -- relative_coords_for_cluster_cell;
                    // w -- cell weight
                    cluster_energy+=hist[k]->
                        GetBinContent(cluster_center_x + q, 
                        cluster_center_y + r);
                }
            }
               
            edepted_energy->
                Fill(true_energy/1000., cluster_energy/1000);
            
            TLorentzVector* elec4m = new TLorentzVector(electron_px, electron_py, electron_pz, electron_p0);
            TLorentzVector* posi4m = new TLorentzVector(positron_px, positron_py, positron_pz, positron_p0);
            TLorentzVector* gamm4m = new TLorentzVector(gamma_px, gamma_py, gamma_pz, gamma_p0);

            Jpsi_mass->Fill((*elec4m + *posi4m).M(), (*elec4m + *posi4m).Pt());
            chic_mass_before->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            
            Double_t p_gamma = gamm4m->P();
            
            gamm4m->SetPx(gamm4m->Px() * cluster_energy/(1000 * gamm4m->E()));
            gamm4m->SetPy(gamm4m->Py() * cluster_energy/(1000 * gamm4m->E()));
            gamm4m->SetPz(gamm4m->Pz() * cluster_energy/(1000 * gamm4m->E()));
            gamm4m->SetE(cluster_energy/1000);

            chic_mass_with_gamma_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);

            Double_t p_elec = elec4m->P();
            Double_t ep_mass = elec4m->M();

            Double_t sigma_p_elec = p_elec * 0.015;

            Double_t smeared_p_elec = gRandom->Gaus(p_elec, sigma_p_elec);

            elec4m->SetPx(elec4m->Px() * smeared_p_elec/p_elec);
            elec4m->SetPy(elec4m->Py() * smeared_p_elec/p_elec);
            elec4m->SetPz(elec4m->Pz() * smeared_p_elec/p_elec);

            elec4m->SetE(TMath::Sqrt(ep_mass * ep_mass + smeared_p_elec * smeared_p_elec));


            Double_t p_posi = posi4m->P();
            ep_mass = posi4m->M();

            Double_t sigma_p_posi = p_posi * 0.015;
            
            Double_t smeared_p_posi = gRandom->Gaus(p_posi, sigma_p_posi);

            posi4m->SetPx(posi4m->Px() * smeared_p_posi/p_posi);
            posi4m->SetPy(posi4m->Py() * smeared_p_posi/p_posi);
            posi4m->SetPz(posi4m->Pz() * smeared_p_posi/p_posi);

            posi4m->SetE(TMath::Sqrt(ep_mass * ep_mass + smeared_p_posi * smeared_p_posi));

            Jpsi_mass_with_elec_posi_resolution->Fill((*elec4m + *posi4m).M(), (*elec4m + *posi4m).Pt());
            chic_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            chic_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);

            if (fabs(branchings - chic0_branching) < 0.00001){
                chic0_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
                chic0_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            } else if (fabs(branchings - chic1_branching) < 0.00001){
                chic1_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
                chic1_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            } else if (fabs(branchings - chic2_branching) < 0.00001) {
                chic2_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
                chic2_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            }

            delete elec4m;
            delete posi4m;
            delete gamm4m;
        }

        delete[] hist;

        for (int i = 0; i < chic0_mass_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic0_mass_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic0_mass_with_gamma_elec_posi_resolution->ProjectionX(Form("chic0_mass_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(0, chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic0_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic0_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
        }

        for (int i = 0; i < chic0_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic0_mass_diff_with_gamma_elec_posi_resolution->ProjectionX(Form("chic0_mass_diff_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(0, chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            mass_chic0_pt_depends->SetBinContent(i + 1, gaus->GetParameter(1));
            mass_chic0_pt_depends->SetBinError(i + 1, gaus->GetParError(1));
            sigma_chic0_pt_depends->SetBinContent(i + 1, gaus->GetParameter(2));
            sigma_chic0_pt_depends->SetBinError(i + 1, gaus->GetParError(2));
            efficiency_chic0_pt_depends->SetBinContent(i + 1, (chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral(chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) - 5 * gaus->GetParameter(2)), chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) + 5 * gaus->GetParameter(2)))/(chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral())));
        }

        for (int i = 0; i < chic1_mass_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic1_mass_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic1_mass_with_gamma_elec_posi_resolution->ProjectionX(Form("chic1_mass_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic1_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic1_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
        }

        for (int i = 0; i < chic1_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic1_mass_diff_with_gamma_elec_posi_resolution->ProjectionX(Form("chic1_mass_diff_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            mass_chic1_pt_depends->SetBinContent(i + 1, gaus->GetParameter(1));
            mass_chic1_pt_depends->SetBinError(i + 1, gaus->GetParError(1));
            sigma_chic1_pt_depends->SetBinContent(i + 1, gaus->GetParameter(2));
            sigma_chic1_pt_depends->SetBinError(i + 1, gaus->GetParError(2));
            efficiency_chic1_pt_depends->SetBinContent(i + 1, (chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral(chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) - 5 * gaus->GetParameter(2)), chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) + 5 * gaus->GetParameter(2)))/(chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral())));
        }

        for (int i = 0; i < chic2_mass_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic2_mass_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic2_mass_with_gamma_elec_posi_resolution->ProjectionX(Form("chic2_mass_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic2_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic2_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            
        }

        for (int i = 0; i < chic2_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic2_mass_diff_with_gamma_elec_posi_resolution->ProjectionX(Form("chic2_mass_diff_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            mass_chic2_pt_depends->SetBinContent(i + 1, gaus->GetParameter(1));
            mass_chic2_pt_depends->SetBinError(i + 1, gaus->GetParError(1));
            sigma_chic2_pt_depends->SetBinContent(i + 1, gaus->GetParameter(2));
            sigma_chic2_pt_depends->SetBinError(i + 1, gaus->GetParError(2));
            // efficiency_chic2_pt_depends->SetBinContent(i + 1, gaus->Integral(gaus->GetParameter(1) - 5 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2))/(chic2_branching * chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral()));
            // efficiency_chic2_pt_depends->SetBinError(i + 1, gaus->IntegralError(gaus->GetParameter(1) - 5 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2))/(chic2_branching * chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral()));
            efficiency_chic2_pt_depends->SetBinContent(i + 1, (chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral(chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) - 5 * gaus->GetParameter(2)), chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) + 5 * gaus->GetParameter(2)))/(chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral())));
        }

        output_data->cd(); 

        
        opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->cd();

        Inv_mass_spectrum_before_simulation[opt_electron_rate - min_opt_electron_rate]->cd();
        Jpsi_mass->Write();
        chic_mass_before->Write();

        Inv_mass_spectrum_after_simulation[opt_electron_rate - min_opt_electron_rate]->cd();
        chic_mass_with_gamma_resolution->Write();
        chic_mass_with_gamma_elec_posi_resolution->Write();
        chic_mass_diff_with_gamma_elec_posi_resolution->Write();

        Inv_mass_sep_chic_after_simulation[opt_electron_rate - min_opt_electron_rate]->cd();
        chic0_mass_with_gamma_elec_posi_resolution->Write();
        chic1_mass_with_gamma_elec_posi_resolution->Write();
        chic2_mass_with_gamma_elec_posi_resolution->Write();
        chic0_mass_diff_with_gamma_elec_posi_resolution->Write();
        chic1_mass_diff_with_gamma_elec_posi_resolution->Write();
        chic2_mass_diff_with_gamma_elec_posi_resolution->Write();


        for (auto slice_chic0_mass_diff : chic0_mass_diff_with_gamma_elec_posi_resolution_slice){
            slice_chic0_mass_diff->Write();
        }
        for (auto slice_chic1_mass_diff : chic1_mass_diff_with_gamma_elec_posi_resolution_slice){
            slice_chic1_mass_diff->Write();
        }
        for (auto slice_chic2_mass_diff : chic2_mass_diff_with_gamma_elec_posi_resolution_slice){
            slice_chic2_mass_diff->Write();
        }

        Chic_parameters[opt_electron_rate - min_opt_electron_rate]->cd();
        mass_chic0_pt_depends->Write();
        sigma_chic0_pt_depends->Write();
        efficiency_chic0_pt_depends->Write();
        mass_chic1_pt_depends->Write();
        sigma_chic1_pt_depends->Write();
        efficiency_chic1_pt_depends->Write();
        mass_chic2_pt_depends->Write();
        sigma_chic2_pt_depends->Write();
        efficiency_chic2_pt_depends->Write();

        
        Energy_distribution_in_cells_chic[opt_electron_rate - min_opt_electron_rate]->cd();
        edepted_energy->Write();
        
    }

    return;
}

void Chi_c_data_from_Geant4_wim_magnet(TFile* output_data, TChain* input){

    double chic0_branching = 0.00083594;
    double chic1_branching = 0.0204805;
    double chic2_branching = 0.0113449;

    std::vector<double>* cells_energy = 0;
    double true_energy = 0;

    double electron_px = 0;
    double electron_py = 0;
    double electron_pz = 0;
    double electron_p0 = 0;

    double positron_px = 0;
    double positron_py = 0;
    double positron_pz = 0;
    double positron_p0 = 0;

    double gamma_px = 0;
    double gamma_py = 0;
    double gamma_pz = 0;
    double gamma_p0 = 0;

    double chi_c_mass_before_rotating = 0;
    double chi_c_pt_before_rotating = 0;

    double branchings = 0;

    int max_cluster_size = 11;

    input->SetBranchAddress("cell_energy", &cells_energy);
    input->SetBranchAddress("electron_px", &electron_px);
    input->SetBranchAddress("electron_py", &electron_py);
    input->SetBranchAddress("electron_pz", &electron_pz);
    input->SetBranchAddress("electron_p0", &electron_p0);
    input->SetBranchAddress("positron_px", &positron_px);
    input->SetBranchAddress("positron_py", &positron_py);
    input->SetBranchAddress("positron_pz", &positron_pz);
    input->SetBranchAddress("positron_p0", &positron_p0);
    input->SetBranchAddress("gamma_px", &gamma_px);
    input->SetBranchAddress("gamma_py", &gamma_py);
    input->SetBranchAddress("gamma_pz", &gamma_pz);
    input->SetBranchAddress("gamma_p0", &gamma_p0);
    input->SetBranchAddress("chi_c_mass_before_rotating", &chi_c_mass_before_rotating);
    input->SetBranchAddress("chi_c_pt_before_rotating", &chi_c_pt_before_rotating);
    input->SetBranchAddress("event_branching", &branchings);
    input->SetBranchAddress("initial_particle_energy", &true_energy);

    output_data->cd();

    TDirectory* Chi_c_data_from_Geant4_wim_dir = (TDirectory*)output_data->GetDirectory("Chi_c_data_from_Geant4_with_magnet");
    Chi_c_data_from_Geant4_wim_dir->cd();

    Int_t min_opt_electron_rate = 7;
    Int_t max_opt_electron_rate = 7;

    // const double opt_electron_rate = 7; -- energy resolution sigma_E/E \prop 0.02/sqrt(E)
    // const double opt_electron_rate = 1.1; -- energy resolution sigma_E/E \sim to PHOS

    TDirectory* opt_electron_rate_dir[max_opt_electron_rate - min_opt_electron_rate + 1];

    TDirectory* Inv_mass_spectrum_before_simulation[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Inv_mass_spectrum_after_simulation[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Inv_mass_sep_chic_after_simulation[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Energy_distribution_in_cells_chic[max_opt_electron_rate - min_opt_electron_rate + 1];
    TDirectory* Chic_parameters[max_opt_electron_rate - min_opt_electron_rate + 1];

    for (int opt_electron_rate = min_opt_electron_rate; opt_electron_rate <= max_opt_electron_rate; ++opt_electron_rate){

        opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate] = (TDirectory*)Chi_c_data_from_Geant4_wim_dir->mkdir(Form("Opt_elctron_rate_%1d_per_MeV", opt_electron_rate));        
        opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->cd();

        Inv_mass_spectrum_before_simulation[opt_electron_rate - min_opt_electron_rate] = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Inv_mass_spectrum_before_simulation");
        Inv_mass_spectrum_after_simulation[opt_electron_rate - min_opt_electron_rate]  = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Inv_mass_spectrum_after_simulation");
        Inv_mass_sep_chic_after_simulation[opt_electron_rate - min_opt_electron_rate]  = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Inv_mass_sep_chic_after_simulation");
        Energy_distribution_in_cells_chic[opt_electron_rate - min_opt_electron_rate]   = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Cluster_energy_distribution_in_cells");
        Chic_parameters[opt_electron_rate - min_opt_electron_rate]                     = (TDirectory*)opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->mkdir("Chic_parameters_pt_depends");

        const int n_entries = input->GetEntries();

        TH2D** hist = new TH2D*[n_entries];

        std::cout << n_entries << "\n";

        std::vector<double> clusters_energy((max_cluster_size - 1)/2, 0.);
        std::vector<double> clusters_energy_error((max_cluster_size - 1)/2, 0.);

        long long int filled_cells = 0;

        int n_empty = 0;

        TH2D* edepted_energy = new TH2D("edepted_energy", "edepted_energy", 500, 0., 5., 500, 0., 5.);

        TH2D* Jpsi_mass = new TH2D("Jpmin_opt_electron_ratesi_mass", "Jpsi_mass", 100, 2.8, 3.4, 10, 0., 10.);
        TH2D* Jpsi_mass_with_elec_posi_resolution = new TH2D("Jpsi_mass_with_elec_posi_resolution", "Jpsi_mass_with_elec_posi_resolution", 100, 2.8, 3.4, 10, 0., 10.);

        TH2D* chic_mass_before = new TH2D("chi_c_mass_before", "chi_c_mass_before", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_before_slice[chic_mass_before->GetNbinsY()];

        TH2D* chic_mass_with_gamma_resolution = new TH2D("chic_mass_with_gamma_resolution", "chic_mass_with_gamma_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_with_gamma_resolution_slice[chic_mass_with_gamma_resolution->GetNbinsY()];

        TH2D* chic_mass_with_gamma_elec_posi_resolution = new TH2D("chic_mass_with_gamma_elec_posi_resolution", "chic_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_with_gamma_elec_posi_resolution_slice[chic_mass_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH2D* chic_mass_diff_with_gamma_elec_posi_resolution = new TH2D("chic_mass_diff_with_gamma_elec_posi_resolution", "chic_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic_mass_diff_with_gamma_elec_posi_resolution_slice[chic_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH2D* chic0_mass_with_gamma_elec_posi_resolution = new TH2D("separate_chic0_mass_with_gamma_elec_posi_resolution", "separate_chic0_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic1_mass_with_gamma_elec_posi_resolution = new TH2D("separate_chic1_mass_with_gamma_elec_posi_resolution", "separate_chic1_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic2_mass_with_gamma_elec_posi_resolution = new TH2D("separate_chic2_mass_with_gamma_elec_posi_resolution", "separate_chic2_mass_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic0_mass_with_gamma_elec_posi_resolution_slice[chic0_mass_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic1_mass_with_gamma_elec_posi_resolution_slice[chic1_mass_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic2_mass_with_gamma_elec_posi_resolution_slice[chic2_mass_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH2D* chic0_mass_diff_with_gamma_elec_posi_resolution = new TH2D("separate_chic0_mass_diff_with_gamma_elec_posi_resolution", "separate_chic0_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic1_mass_diff_with_gamma_elec_posi_resolution = new TH2D("separate_chic1_mass_diff_with_gamma_elec_posi_resolution", "separate_chic1_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH2D* chic2_mass_diff_with_gamma_elec_posi_resolution = new TH2D("separate_chic2_mass_diff_with_gamma_elec_posi_resolution", "separate_chic2_mass_diff_with_gamma_elec_posi_resolution", 100, 3.3, 3.6, 10, 0., 10.);
        TH1D* chic0_mass_diff_with_gamma_elec_posi_resolution_slice[chic0_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic1_mass_diff_with_gamma_elec_posi_resolution_slice[chic1_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];
        TH1D* chic2_mass_diff_with_gamma_elec_posi_resolution_slice[chic2_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY()];

        TH1D* mass_chic0_pt_depends = new TH1D("mass_chic0_pt_depends", "mass_chic0_pt_depends", 10, -0.5, 10.5);
        TH1D* mass_chic1_pt_depends = new TH1D("mass_chic1_pt_depends", "mass_chic1_pt_depends", 10, -0.5, 10.5);
        TH1D* mass_chic2_pt_depends = new TH1D("mass_chic2_pt_depends", "mass_chic2_pt_depends", 10, -0.5, 10.5);

        TH1D* sigma_chic0_pt_depends = new TH1D("sigma_chic0_pt_depends", "sigma_chic0_pt_depends", 10, -0.5, 10.5);
        TH1D* sigma_chic1_pt_depends = new TH1D("sigma_chic1_pt_depends", "sigma_chic1_pt_depends", 10, -0.5, 10.5);
        TH1D* sigma_chic2_pt_depends = new TH1D("sigma_chic2_pt_depends", "sigma_chic2_pt_depends", 10, -0.5, 10.5);

        TH1D* efficiency_chic0_pt_depends = new TH1D("efficiency_chic0_pt_depends", "efficiency_chic0_pt_depends", 10, -0.5, 10.5);
        TH1D* efficiency_chic1_pt_depends = new TH1D("efficiency_chic1_pt_depends", "efficiency_chic1_pt_depends", 10, -0.5, 10.5);
        TH1D* efficiency_chic2_pt_depends = new TH1D("efficiency_chic2_pt_depends", "efficiency_chic2_pt_depends", 10, -0.5, 10.5);

        for (int k = 0; k < n_entries; ++k){
            input->GetEntry(k);

            if (k % 10000 == 0){
                std::cout << 
                    Form("k = %06d in %d\n", k, n_entries);
            }

            int cal_size = static_cast<int>(sqrt(cells_energy->size()));
            hist[k] = new TH2D(Form("event%06d", k + 1), "event_calorimeter_map", 
                                cal_size, -cal_size/2, cal_size/2, 
                                cal_size, -cal_size/2, cal_size/2);
            for (int m = 0; m < cells_energy->size(); ++m){
                int n_pr_opt_photons = 
                    (int)(((*cells_energy)[m]) * opt_electron_rate);
                int n_sm_opt_photons = 
                    (int)gRandom->Poisson(n_pr_opt_photons);
                n_sm_opt_photons = 
                    (int)gRandom->Gaus(n_sm_opt_photons, 0.007*n_sm_opt_photons);
                (*cells_energy)[m] = n_sm_opt_photons / opt_electron_rate;
                if ((*cells_energy)[m] <= 4){
                    (*cells_energy)[m] = 0;
                } else {
                    filled_cells++;
                }
                hist[k]->SetBinContent((m)/ cal_size + 1, 
                                        (m) % cal_size + 1, 
                                        (*cells_energy)[m]);
            }

            int max_bin_index = hist[k]->GetMaximumBin();
            int cluster_center_x = max_bin_index / cal_size - 1;
            int cluster_center_y = max_bin_index % cal_size - 1;

            double full_event_energy = hist[k]->Integral();
            if (full_event_energy <= 0.1){
                ++n_empty;
                continue;
            }

            double integral_error = 0;

            for (int i = 0; i < 441; ++i){
                integral_error += hist[k]->GetBinError(i);
            }

            int p = 3;
            double x_cm = 0;
            double weight_sum = 0;
            double cluster_energy = 0;
            for (int q = - (p - 1)/2; q < (p + 1)/2; ++q){
                for (int r = - (p - 1)/2; r < (p + 1)/2; ++r){
                    // (q, r) -- relative_coords_for_cluster_cell;
                    // w -- cell weight
                    cluster_energy+=hist[k]->
                        GetBinContent(cluster_center_x + q, 
                        cluster_center_y + r);
                }
            }
               
            edepted_energy->
                Fill(true_energy/1000., cluster_energy/1000);
            
            TLorentzVector* elec4m = new TLorentzVector(electron_px, electron_py, electron_pz, electron_p0);
            TLorentzVector* posi4m = new TLorentzVector(positron_px, positron_py, positron_pz, positron_p0);
            TLorentzVector* gamm4m = new TLorentzVector(gamma_px, gamma_py, gamma_pz, gamma_p0);

            Jpsi_mass->Fill((*elec4m + *posi4m).M(), (*elec4m + *posi4m).Pt());
            chic_mass_before->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            
            Double_t p_gamma = gamm4m->P();
            
            gamm4m->SetPx(gamm4m->Px() * cluster_energy/(1000 * gamm4m->E()));
            gamm4m->SetPy(gamm4m->Py() * cluster_energy/(1000 * gamm4m->E()));
            gamm4m->SetPz(gamm4m->Pz() * cluster_energy/(1000 * gamm4m->E()));
            gamm4m->SetE(cluster_energy/1000);

            chic_mass_with_gamma_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);

            Double_t p_elec = elec4m->P();
            Double_t ep_mass = elec4m->M();

            Double_t sigma_p_elec = p_elec * 0.015;

            Double_t smeared_p_elec = gRandom->Gaus(p_elec, sigma_p_elec);

            elec4m->SetPx(elec4m->Px() * smeared_p_elec/p_elec);
            elec4m->SetPy(elec4m->Py() * smeared_p_elec/p_elec);
            elec4m->SetPz(elec4m->Pz() * smeared_p_elec/p_elec);

            elec4m->SetE(TMath::Sqrt(ep_mass * ep_mass + smeared_p_elec * smeared_p_elec));


            Double_t p_posi = posi4m->P();
            ep_mass = posi4m->M();

            Double_t sigma_p_posi = p_posi * 0.015;
            
            Double_t smeared_p_posi = gRandom->Gaus(p_posi, sigma_p_posi);

            posi4m->SetPx(posi4m->Px() * smeared_p_posi/p_posi);
            posi4m->SetPy(posi4m->Py() * smeared_p_posi/p_posi);
            posi4m->SetPz(posi4m->Pz() * smeared_p_posi/p_posi);

            posi4m->SetE(TMath::Sqrt(ep_mass * ep_mass + smeared_p_posi * smeared_p_posi));

            Jpsi_mass_with_elec_posi_resolution->Fill((*elec4m + *posi4m).M(), (*elec4m + *posi4m).Pt());
            chic_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            chic_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);

            if (fabs(branchings - chic0_branching) < 0.00001){
                chic0_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
                chic0_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            } else if (fabs(branchings - chic1_branching) < 0.00001){
                chic1_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
                chic1_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            } else if (fabs(branchings - chic2_branching) < 0.00001) {
                chic2_mass_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M(), (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
                chic2_mass_diff_with_gamma_elec_posi_resolution->Fill((*elec4m + *posi4m + *gamm4m).M() - (*elec4m + *posi4m).M() + 3.096, (*elec4m + *posi4m + *gamm4m).Pt(), branchings);
            }

            delete elec4m;
            delete posi4m;
            delete gamm4m;
        }

        delete[] hist;

        for (int i = 0; i < chic0_mass_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic0_mass_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic0_mass_with_gamma_elec_posi_resolution->ProjectionX(Form("chic0_mass_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(0, chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic0_mass_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic0_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic0_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
        }

        for (int i = 0; i < chic0_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic0_mass_diff_with_gamma_elec_posi_resolution->ProjectionX(Form("chic0_mass_diff_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(0, chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            mass_chic0_pt_depends->SetBinContent(i + 1, gaus->GetParameter(1));
            mass_chic0_pt_depends->SetBinError(i + 1, gaus->GetParError(1));
            sigma_chic0_pt_depends->SetBinContent(i + 1, gaus->GetParameter(2));
            sigma_chic0_pt_depends->SetBinError(i + 1, gaus->GetParError(2));
            efficiency_chic0_pt_depends->SetBinContent(i + 1, (chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral(chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) - 5 * gaus->GetParameter(2)), chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) + 5 * gaus->GetParameter(2)))/(chic0_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral())));
        }

        for (int i = 0; i < chic1_mass_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic1_mass_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic1_mass_with_gamma_elec_posi_resolution->ProjectionX(Form("chic1_mass_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic1_mass_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic1_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic1_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
        }

        for (int i = 0; i < chic1_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic1_mass_diff_with_gamma_elec_posi_resolution->ProjectionX(Form("chic1_mass_diff_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            mass_chic1_pt_depends->SetBinContent(i + 1, gaus->GetParameter(1));
            mass_chic1_pt_depends->SetBinError(i + 1, gaus->GetParError(1));
            sigma_chic1_pt_depends->SetBinContent(i + 1, gaus->GetParameter(2));
            sigma_chic1_pt_depends->SetBinError(i + 1, gaus->GetParError(2));
            efficiency_chic1_pt_depends->SetBinContent(i + 1, (chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral(chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) - 5 * gaus->GetParameter(2)), chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) + 5 * gaus->GetParameter(2)))/(chic1_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral())));
        }

        for (int i = 0; i < chic2_mass_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic2_mass_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic2_mass_with_gamma_elec_posi_resolution->ProjectionX(Form("chic2_mass_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic2_mass_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic2_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic2_mass_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            
        }

        for (int i = 0; i < chic2_mass_diff_with_gamma_elec_posi_resolution->GetNbinsY(); ++i){
            chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i] = (TH1D*)chic2_mass_diff_with_gamma_elec_posi_resolution->ProjectionX(Form("chic2_mass_diff_with_gamma_elec_posi_resolution_slice_%02d", i + 1), i + 1, i + 1);
            TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", 0., 4.5);
            gaus->SetParameter(1, chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximum());
            gaus->SetParameter(1, chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetBinCenter(chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetMaximumBin()));
            gaus->SetParameter(2, chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->GetStdDev());
            gaus->SetParLimits(2, 0., 1.);
            chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus);
            chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Fit(gaus, "", "", gaus->GetParameter(1) - 1 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2));
            mass_chic2_pt_depends->SetBinContent(i + 1, gaus->GetParameter(1));
            mass_chic2_pt_depends->SetBinError(i + 1, gaus->GetParError(1));
            sigma_chic2_pt_depends->SetBinContent(i + 1, gaus->GetParameter(2));
            sigma_chic2_pt_depends->SetBinError(i + 1, gaus->GetParError(2));
            // efficiency_chic2_pt_depends->SetBinContent(i + 1, gaus->Integral(gaus->GetParameter(1) - 5 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2))/(chic2_branching * chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral()));
            // efficiency_chic2_pt_depends->SetBinError(i + 1, gaus->IntegralError(gaus->GetParameter(1) - 5 * gaus->GetParameter(2), gaus->GetParameter(1) + 5 * gaus->GetParameter(2))/(chic2_branching * chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral()));
            efficiency_chic2_pt_depends->SetBinContent(i + 1, (chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral(chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) - 5 * gaus->GetParameter(2)), chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->FindBin(gaus->GetParameter(1) + 5 * gaus->GetParameter(2)))/(chic2_mass_diff_with_gamma_elec_posi_resolution_slice[i]->Integral())));
        }

        output_data->cd(); 

        
        opt_electron_rate_dir[opt_electron_rate - min_opt_electron_rate]->cd();

        Inv_mass_spectrum_before_simulation[opt_electron_rate - min_opt_electron_rate]->cd();
        Jpsi_mass->Write();
        chic_mass_before->Write();

        Inv_mass_spectrum_after_simulation[opt_electron_rate - min_opt_electron_rate]->cd();
        chic_mass_with_gamma_resolution->Write();
        chic_mass_with_gamma_elec_posi_resolution->Write();
        chic_mass_diff_with_gamma_elec_posi_resolution->Write();

        Inv_mass_sep_chic_after_simulation[opt_electron_rate - min_opt_electron_rate]->cd();
        chic0_mass_with_gamma_elec_posi_resolution->Write();
        chic1_mass_with_gamma_elec_posi_resolution->Write();
        chic2_mass_with_gamma_elec_posi_resolution->Write();
        chic0_mass_diff_with_gamma_elec_posi_resolution->Write();
        chic1_mass_diff_with_gamma_elec_posi_resolution->Write();
        chic2_mass_diff_with_gamma_elec_posi_resolution->Write();


        for (auto slice_chic0_mass_diff : chic0_mass_diff_with_gamma_elec_posi_resolution_slice){
            slice_chic0_mass_diff->Write();
        }
        for (auto slice_chic1_mass_diff : chic1_mass_diff_with_gamma_elec_posi_resolution_slice){
            slice_chic1_mass_diff->Write();
        }
        for (auto slice_chic2_mass_diff : chic2_mass_diff_with_gamma_elec_posi_resolution_slice){
            slice_chic2_mass_diff->Write();
        }

        Chic_parameters[opt_electron_rate - min_opt_electron_rate]->cd();
        mass_chic0_pt_depends->Write();
        sigma_chic0_pt_depends->Write();
        efficiency_chic0_pt_depends->Write();
        mass_chic1_pt_depends->Write();
        sigma_chic1_pt_depends->Write();
        efficiency_chic1_pt_depends->Write();
        mass_chic2_pt_depends->Write();
        sigma_chic2_pt_depends->Write();
        efficiency_chic2_pt_depends->Write();

        
        Energy_distribution_in_cells_chic[opt_electron_rate - min_opt_electron_rate]->cd();
        edepted_energy->Write();
        
    }

    return;
}
