#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <utility>
#include <tuple>
#include <string>
#include "Include/GraphingUtils.h"
#include "Include/HistVariables.h"
#include "Include/DataFilePaths.h"
#include "Include/EfficiencyResults.h"
#include "Include/EfficiencyCounting.h"

using namespace std;

//Function Declarations

void EfficiencyPurityDataWrite_withtrkpt(
    float xsec_SD, 
    float xsec_DD,
    float xsec_had, 
    float xsec_alphaO,
    int N, 
    float xMax,
    int trkptcut,
    string outfolder = "DataFilesExample/",
    string outFileAdd = "");

void EfficiencyPurityDataWrite_withtrkpt_asymmetric(
    float xsec_SD, 
    float xsec_DD,
    float xsec_had, 
    int trkptcut,
    vector<pair<float, float>> HFEcutvector,
    string outfolder = "DataFilesExample/",
    string outFileAdd = "");

void EfficiencyPurityDataWrite_withtrkpt_root(
    const char* variableNamePlus = "HFEMaxPlus",
    const char* variableNameMinus = "HFEMaxMinus",
    int N = 11, 
    float xMax = 10,
    int trkptcut = -1, // Default to -1, can be changed if needed
    string outfolder = "DataFilesExample/",
    string outFileAdd = "");

void run(){

    int N = 51;
    float xMax = 50;
    float trkptcut = -1;

    vector<pair<float, float>> HFEcutvector;
    for (int i = 0; i < 5; ++i) {
        HFEcutvector.push_back({static_cast<float>(i), static_cast<float>(i)});
    }
    for (int i = 0; i < HFEcutvector.size(); ++i) {
        cout << "HFEcutvector: " << HFEcutvector[i].first << ", " << HFEcutvector[i].second << endl;
    }

    EfficiencyPurityDataWrite_withtrkpt_root("HFEMaxPlus","HFEMaxMinus", 
            N, xMax, 3, "Data/DataFiles_Final/", "Offline_allValues");
    EfficiencyPurityDataWrite_withtrkpt_root("mMaxL1HFAdcPlus","mMaxL1HFAdcMinus", 
            N, xMax, 3, "Data/DataFiles_Final/", "Online_allValues");

 /*/   for (int i = 0; i <= 6; i+=2){
        EfficiencyPurityDataWrite_withtrkpt_root("mMaxL1HFAdcPlus","mMaxL1HFAdcMinus", 
            N, xMax, i, "Data/DataFiles_Final/", "Online_allValues");
        EfficiencyPurityDataWrite_withtrkpt_root("HFEMaxPlus","HFEMaxMinus", 
            N, xMax, i, "Data/DataFiles_Final/", "Offline_allValues");
    }*/


}

void EfficiencyPurityDataWrite_withtrkpt_root(
    const char* varNamePlus = "",
    const char* varNameMinus = "",
    int N = 11, 
    float xMax = 10,
    int trkptcut = -1, // Default to -1, can be changed if needed
    string outfolder = "DataFilesExample/",
    string outFileAdd = "") 
{
    const char* inFileNameSignal = HIJINGFile; // Default signal file 
    const char* inFileNameBKGSD = StarlightSD; // Default background file for SD
    const char* inFileNameBKGDD = StarlightDD; // Default background file for DD
    const char* inFileNameaO = alphaOFile;     // Default background file for alphaO

    string rootfilename = outfolder + Form("EfficiencyPurityData_N%d_trkpt%d_", N-1, trkptcut) + outFileAdd + ".root";
    TFile *rootfile = new TFile(rootfilename.c_str(), "RECREATE");
    TTree *tree = new TTree("effTree", "Efficiency and Purity Data");

    // Variables for branches
    double cut;
    double Eff_HijingAND, Eff_HijingOR;
    double Eff_SDAND, Eff_SDOR;
    double Eff_DDAND, Eff_DDOR;
    double Eff_AlphaOAND, Eff_AlphaOOR;
    double PurityAND, PurityOR;
    int nTracks;

    // Create branches
    tree->Branch("Cut", &cut, "Cut/D");
    tree->Branch("Eff_HijingAND", &Eff_HijingAND, "Eff_HijingAND/D");
    tree->Branch("Eff_HijingOR", &Eff_HijingOR, "Eff_HijingOR/D");
    tree->Branch("Eff_SDAND", &Eff_SDAND, "Eff_SDAND/D");
    tree->Branch("Eff_SDOR", &Eff_SDOR, "Eff_SDOR/D");
    tree->Branch("Eff_DDAND", &Eff_DDAND, "Eff_DDAND/D");
    tree->Branch("Eff_DDOR", &Eff_DDOR, "Eff_DDOR/D");
    tree->Branch("Eff_AlphaOAND", &Eff_AlphaOAND, "Eff_AlphaOAND/D");
    tree->Branch("Eff_AlphaOOR", &Eff_AlphaOOR, "Eff_AlphaOOR/D");
    tree->Branch("PurityAND", &PurityAND, "PurityAND/D");
    tree->Branch("PurityOR", &PurityOR, "PurityOR/D");
    tree->Branch("nTracks", &nTracks, "nTracks/I");

    cout << "------- Writing Efficiency and Purity Data to ROOT file -------" << endl;
    cout << "Output file: " << rootfilename << endl;
    cout << "Variable Names: " << varNamePlus << ", " << varNameMinus << endl;

    for (int i = 0; i < N; ++i) {
        cut = xMax * i / (N - 1); // Uniform spacing from 0 to xMax
        Eff_HijingAND = countingTrkptVariable(inFileNameSignal, varNamePlus, varNameMinus, cut, true, true, trkptcut).first;
        Eff_HijingOR  = countingTrkptVariable(inFileNameSignal, varNamePlus, varNameMinus, cut, false, true, trkptcut).first;
        Eff_SDAND     = countingTrkptVariable(inFileNameBKGSD, varNamePlus, varNameMinus, cut, true, false, trkptcut).first;
        Eff_SDOR      = countingTrkptVariable(inFileNameBKGSD, varNamePlus, varNameMinus, cut, false, false, trkptcut).first;
        Eff_DDAND     = countingTrkptVariable(inFileNameBKGDD, varNamePlus, varNameMinus, cut, true, false, trkptcut).first;
        Eff_DDOR      = countingTrkptVariable(inFileNameBKGDD, varNamePlus, varNameMinus, cut, false, false, trkptcut).first;
        Eff_AlphaOAND = countingTrkptVariable(inFileNameaO, varNamePlus, varNameMinus, cut, true, false, trkptcut).first;
        Eff_AlphaOOR  = countingTrkptVariable(inFileNameaO,varNamePlus, varNameMinus, cut, false, false, trkptcut).first;
        
        tree->Fill();
    }

    rootfile->Write();
    rootfile->Close();
    delete rootfile;
}


void EfficiencyPurityDataWrite_withtrkpt(
    float xsec_SD, 
    float xsec_DD,
    float xsec_had, 
    float xsec_alphaO,
    int N, 
    float xMax,
    int trkptcut,
    string outfolder = "DataFilesExample/",
    string outFileAdd = "") {

    const char* inFileNameSignal = HIJINGFile; // Default signal file 
    const char* inFileNameBKGSD = StarlightSD; // Default background file for SD
    const char* inFileNameBKGDD = StarlightDD; // Default background file for DD
    const char* inFileNameaO = alphaOFile; // Default background file for DD

        
    string filename = outfolder + Form("EfficiencyPurityData_N%d_trkpt%d_", N-1, trkptcut) + outFileAdd + ".txt";
    cout << "------- Writing Efficiency and Purity Data to File -------" << endl;
    ofstream outfile(filename.c_str());
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return;
    }
    cout << "Output file: " << filename << endl;
    double x[N], y[N], y2[N];
    double EfficiencyDDAND[N], EfficiencyDDOR[N], EfficiencyaOAND[N], EfficiencyaOOR[N];
    double EfficiencyHijingAND[N], EfficiencyHijingOR[N],EfficiencySDAND[N], EfficiencySDOR[N],PurityAND[N], PurityOR[N];
    float leadingPtCut = 0.0;
    outfile << "SignalFile: " << inFileNameSignal << "\n";
    outfile << "BackgroundFileSD: " << inFileNameBKGSD << "\n";
    outfile << "BackgroundFileDD: " << inFileNameBKGDD << "\n";
    outfile << "BackgroundFileAlphaO: " << inFileNameaO << "\n";
    outfile << "xsec_SD: " << xsec_SD << " "
            << "xsec_DD: " << xsec_DD << " "
            << "xsec_had: " << xsec_had << " "
            << "xsec_alphaO: " << xsec_alphaO << " "
            << "N: " << N << " "
            << "xMax: " << xMax << " "
            << "trkptcut: " << trkptcut << "\n";
    cout << "xsec_SD: " << xsec_SD << ", xsec_DD: " << xsec_DD << ", xsec_had: " << xsec_had << ", xsec_alphaO: " << xsec_alphaO << endl;
    outfile << "HFEMaxCut "
            << "Eff_Hijing(AND) Eff_Hijing(OR) "
            << "Eff_SD(AND) Eff_SD(OR) "
            << "Eff_DD(AND) Eff_DD(OR) "
            << "Eff_AlphaO(AND) Eff_AlphaO(OR) "
            << "Purity(AND) Purity(OR)\n";
    for (int i = 0; i < N; ++i) {
        x[i] = xMax * i / (N - 1); // Uniform spacing from 0 to xMax
        EfficiencyHijingAND[i] = countingTrkptNPart(inFileNameSignal, x[i], true, true, trkptcut).first;
        EfficiencyHijingOR[i] = countingTrkptNPart(inFileNameSignal, x[i], false, true, trkptcut).first;
        EfficiencySDAND[i] = countingTrkptNPart(inFileNameBKGSD, x[i], true,false, trkptcut).first;
        EfficiencySDOR[i] = countingTrkptNPart(inFileNameBKGSD, x[i], false,false, trkptcut).first;
        EfficiencyDDAND[i] = countingTrkptNPart(inFileNameBKGDD, x[i], true,false, trkptcut).first;
        EfficiencyDDOR[i] = countingTrkptNPart(inFileNameBKGDD, x[i], false,false, trkptcut).first;
        EfficiencyaOAND[i] = countingTrkptNPart(inFileNameaO, x[i], true,false, trkptcut).first;
        EfficiencyaOOR[i] = countingTrkptNPart(inFileNameaO, x[i], false,false, trkptcut).first;

        PurityAND[i] = (xsec_had * EfficiencyHijingAND[i]) / (xsec_SD * EfficiencySDAND[i] + xsec_DD * EfficiencyDDAND[i] + xsec_had * EfficiencyHijingAND[i] + xsec_alphaO * EfficiencyaOAND[i]);
        PurityOR[i] = (xsec_had * EfficiencyHijingOR[i]) / (xsec_SD * EfficiencySDOR[i] + xsec_DD * EfficiencyDDOR[i] + xsec_had * EfficiencyHijingOR[i] + xsec_alphaO * EfficiencyaOOR[i]);
        cout << "EfficiencyAND[" << i << "] = " << EfficiencyHijingAND[i] << endl;
        cout << "EfficiencySDAND[" << i << "] = " << EfficiencySDAND[i] << endl;
        cout << "EfficiencyDDAND[" << i << "] = " << EfficiencyDDAND[i] << endl;
        cout << "EfficiencyAlphaOAND[" << i << "] = " << EfficiencyaOAND[i] << endl;
        cout << "EfficiencyOR[" << i << "] = " << EfficiencyHijingOR[i] << endl;

        cout << "PurityAND[" << i << "] = " << PurityAND[i] << endl;
        cout << "PurityOR[" << i << "] = " << PurityOR[i] << endl;
        outfile << x[i] << " "
                << EfficiencyHijingAND[i] << " "
                << EfficiencyHijingOR[i] << " "
                << EfficiencySDAND[i] << " "
                << EfficiencySDOR[i] << " "
                << EfficiencyDDAND[i] << " "
                << EfficiencyDDOR[i] << " "
                << EfficiencyaOAND[i] << " "
                << EfficiencyaOOR[i] << " "
                << PurityAND[i] << " "
                << PurityOR[i] << "\n";
    }
    outfile.close();
}

void EfficiencyPurityDataWrite_withtrkpt_asymmetric(
    float xsec_SD, 
    float xsec_DD,
    float xsec_had, 
    float xsec_alphaO,
    int trkptcut,
    vector<pair<float, float>> HFEcutvector,
    string outfolder = "DataFilesExample/",
    string outFileAdd = "") {

    const char* inFileNameSignal = HIJINGFile; // Default signal file 
    const char* inFileNameBKGSD = StarlightSD; // Default background file for SD
    const char* inFileNameBKGDD = StarlightDD; // Default background file for DD
    const char* inFileNameBKGaO = alphaOFile; // Default background file for DD

    int N = HFEcutvector.size();

    string filename = outfolder + Form("EfficiencyPurityData_N%d_trkpt%d_", N-1, trkptcut) + outFileAdd + ".txt";
    cout << "------- Writing Efficiency and Purity Data to File -------" << endl;
    ofstream outfile(filename.c_str());
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file " << filename << " for writing." << endl;
        return;
    }
    cout << "Output file: " << filename << endl;
    double xplus[N], xminus[N];
    double EfficiencyDDAND[N], EfficiencyDDOR[N];
    double EfficiencyHijingAND[N], EfficiencyHijingOR[N],EfficiencySDAND[N], EfficiencySDOR[N],BKGRejectionAND[N], BKGRejectionOR[N],PurityAND[N], PurityOR[N];
    outfile << "SignalFile: " << inFileNameSignal << "\n";
    outfile << "BackgroundFileSD: " << inFileNameBKGSD << "\n";
    outfile << "BackgroundFileDD: " << inFileNameBKGDD << "\n";
    outfile << "xsec_SD: " << xsec_SD << " "
            << "xsec_DD: " << xsec_DD << " "
            << "xsec_had: " << xsec_had << " "
            << "trkptcut: " << trkptcut << "\n";
    outfile << "HFECutHigh HFECutLow Eff_Hijing(AND) Eff_Hijing(OR) Eff_SD(AND) Eff_SD(OR) Eff_DD(AND) Eff_DD(OR) BKGRej(AND) BKGRej(OR) Purity(AND) Purity(OR)\n";
    for (int i = 0; i < N; ++i) {
        xplus[i] = HFEcutvector[i].first; // Use the first element of the pair
        xminus[i] = HFEcutvector[i].second; // Use the second element of the pair
        EfficiencyHijingAND[i] = countingTrkptNPartAsymmetric(inFileNameSignal, xplus[i], xminus[i], true, true, trkptcut).first;
        EfficiencyHijingOR[i] = countingTrkptNPartAsymmetric(inFileNameSignal, xplus[i], xminus[i], false, true, trkptcut).first;
        EfficiencySDAND[i] = countingTrkptNPartAsymmetric(inFileNameBKGSD, xplus[i], xminus[i], true, false, trkptcut).first;
        EfficiencySDOR[i] = countingTrkptNPartAsymmetric(inFileNameBKGSD, xplus[i], xminus[i], false, false, trkptcut).first;
        EfficiencyDDAND[i] = countingTrkptNPartAsymmetric(inFileNameBKGDD, xplus[i], xminus[i], true, false, trkptcut).first;
        EfficiencyDDOR[i] = countingTrkptNPartAsymmetric(inFileNameBKGDD, xplus[i], xminus[i], false, false, trkptcut).first;
        BKGRejectionAND[i] = 1.0 - EfficiencySDAND[i];
        BKGRejectionOR[i] = 1.0 - EfficiencySDOR[i];
        PurityAND[i] = 1 - (xsec_SD * EfficiencySDAND[i] + xsec_DD * EfficiencyDDAND[i]) / (xsec_SD * EfficiencySDAND[i] + xsec_DD * EfficiencyDDAND[i] + xsec_had * EfficiencyHijingAND[i]);
        PurityOR[i] = 1 - (xsec_SD * EfficiencySDOR[i] + xsec_DD * EfficiencyDDOR[i]) / (xsec_SD * EfficiencySDOR[i] + xsec_DD * EfficiencyDDOR[i] + xsec_had * EfficiencyHijingOR[i]);
        cout << "EfficiencyAND[" << i << "] = " << EfficiencyHijingAND[i] << endl;
        cout << "EfficiencyOR[" << i << "] = " << EfficiencyHijingOR[i] << endl;   
        cout << "EfficiencySDAND[" << i << "] = " << EfficiencySDAND[i] << endl;
        cout << "EfficiencySDOR[" << i << "] = " << EfficiencySDOR[i] << endl;
        cout << "PurityAND[" << i << "] = " << PurityAND[i] << endl;
        cout << "PurityOR[" << i << "] = " << PurityOR[i] << endl;
        outfile << xplus[i] << " "
                << xminus[i] << " "
                << EfficiencyHijingAND[i] << " "
                << EfficiencyHijingOR[i] << " "
                << EfficiencySDAND[i] << " "
                << EfficiencySDOR[i] << " "
                << EfficiencyDDAND[i] << " "
                << EfficiencyDDOR[i] << " "
                << BKGRejectionAND[i] << " "
                << BKGRejectionOR[i] << " "
                << PurityAND[i] << " "
                << PurityOR[i] << "\n";
    }
    outfile.close();
}