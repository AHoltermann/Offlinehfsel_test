#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <utility>
#include <tuple>
#include <string>
#include "GraphingUtils.h"
#include "HistVariables.h"
#include "EfficiencyResults.h"
#include "/home/xirong/MITHIGAnalysis2024/CommonCode/include/Messenger.h"


using namespace std;

//Function declarations
pair<double, int> countingTrkptNPart(const char* inFileName, 
    float HFEcut = 0, 
    bool HFEcutBool = false, 
    bool HijingBool = false,
    int trkptcut = 0);
pair<double, int> countingLeadingptcuts(const char* inFileName, float cut = 0, bool cutBool = false, int leadingPtSelection = 0);
pair<double, int> countingTrkptcuts(const char* inFileName, float cut = 0, bool cutBool = false, int trkptcut = 0);
pair<double, int> countingCoinc(const char* inFileName, float cut = 0, bool cutBool = false, int coincidence = 0);
pair<double, int> countingTrkptNPartAsymmetric(const char* inFileName, 
    float HFEcutplus = 0, 
    float HFEcutminus = 0,
    bool HFEcutBool = false, 
    bool HijingBool = false,
    int trkptcut = 0);
std::pair<double, int> countingTrkptVariable(const char* inFileName, 
    const char* variableNamePlus,
    const char* variableNameMinus,
    float cut = 0, 
    bool cutBool = false, 
    bool HijingBool = false,
    int trkptcut = 0);


//Write Efficiency, Purity, and BKG Rejection to a file
pair<double, int> countingLeadingptcuts(const char* inFileName, float cut = 0, bool cutBool = false, int leadingPtSelection = 0) {
    cout << "File Name: " << inFileName << endl;
    cout << Form("------- Counting for HFEMax cut %f and ptcut %d -------",cut,leadingPtSelection) << endl;

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return make_pair(0.0, 0);
    }

    TTree* tree = (TTree*)inFile->Get("Tree");
    if (!tree) {
        cout << "Error: Could not find the tree!" << endl;
        return make_pair(0.0, 0);
    }
    
    int CCFilter, PVFilter, nVtx;
    float HFEMaxPlus, HFEMaxMinus;
    float VZ, leadingPtEta1p0_sel;
    const Long64_t nEvt = tree->GetEntries();
    int eventCounter = 0, subeventCounter = 0;
    int denominator = 0;

    tree->SetBranchAddress("ClusterCompatibilityFilter", &CCFilter);
    tree->SetBranchAddress("PVFilter", &PVFilter);
    tree->SetBranchAddress("nVtx", &nVtx);
    tree->SetBranchAddress("VZ", &VZ);
    tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
    tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
    tree->SetBranchAddress("leadingPtEta1p0_sel",&leadingPtEta1p0_sel);

    for (int i = 0; i < nEvt; i++) {
        tree->GetEntry(i);
        if (leadingPtEta1p0_sel > leadingPtSelection) {
            denominator++;
            if (CCFilter == 1 && PVFilter == 1 && nVtx > 0 && abs(VZ) < 15) {
                if (cutBool == false){
                    if (HFEMaxPlus > cut || HFEMaxMinus > cut){
                        eventCounter++;
                    }
                }
                else if (cutBool == true){
                    if (HFEMaxPlus > cut && HFEMaxMinus > cut){
                        eventCounter++;
                    }
                }
            }
        }
    }

    double ratio = static_cast<double>(eventCounter) / denominator;
    cout << "Ratio between Original and Filtered: " << ratio << endl;
    cout << "Total number of events: " << nEvt << endl;
    cout << "Denominator after leadingPtSelection: " << leadingPtSelection << ": " << denominator << endl;
    cout << "Total number of events after HFEcut "<< cut << ": " <<eventCounter << endl;
    cout << endl;
    cout << "------- Count Complete -------" << endl;
    inFile->Close();
    delete inFile;
    return make_pair(ratio, eventCounter);
}
pair<double, int> countingTrkptcuts(const char* inFileName, 
    float cut = 0, 
    bool cutBool = false, 
    int trkptcut = 0) {
    cout << Form("------- Counting for HFEMax cut %f and ptcut %d -------",cut,trkptcut) << endl;

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return make_pair(0.0, 0);
    }

    TTree* tree = (TTree*)inFile->Get("Tree");
    if (!tree) {
        cout << "Error: Could not find the tree!" << endl;
        return make_pair(0.0, 0);
    }
    
    int CCFilter, PVFilter, nVtx;
    float HFEMaxPlus, HFEMaxMinus;
    float VZ, leadingPtEta1p0_sel;
    vector<float>* trkPt = nullptr;
    const Long64_t nEvt = tree->GetEntries();
    int eventCounter = 0, subeventCounter = 0;
    int denominator = 0;
    bool trkPtBool = false;

    tree->SetBranchAddress("ClusterCompatibilityFilter", &CCFilter);
    tree->SetBranchAddress("PVFilter", &PVFilter);
    tree->SetBranchAddress("nVtx", &nVtx);
    tree->SetBranchAddress("VZ", &VZ);
    tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
    tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
    tree->SetBranchAddress("leadingPtEta1p0_sel",&leadingPtEta1p0_sel);
    tree->SetBranchAddress("trkPt",&trkPt); 

    for (int i = 0; i < nEvt; i++) {
        tree->GetEntry(i);
        trkPtBool = false;
        if (trkptcut == -1){
            trkPtBool = true;
        }
        for (int j = 0; j < trkPt->size(); j++) {
            if (trkPt->at(j) > trkptcut) {
                trkPtBool = true;
            }
        }
        if (trkPtBool) {
            denominator++;
            if (CCFilter == 1 && PVFilter == 1 && nVtx > 0 && abs(VZ) < 15) {
                if (cutBool == false){
                    if (HFEMaxPlus > cut || HFEMaxMinus > cut){
                        eventCounter++;
                    }
                }
                else if (cutBool == true){
                    if (HFEMaxPlus > cut && HFEMaxMinus > cut){
                        eventCounter++;
                    }
                }
            }
        }
    }

    double ratio = static_cast<double>(eventCounter) / denominator;
    cout << "Ratio between Original and Filtered: " << ratio << endl;
    cout << "Total number of events: " << nEvt << endl;
    cout << "Denominator after trkptcuts: " << denominator << ": " << denominator << endl;
    cout << "Total number of events after HFEcut "<< cut << ": " <<eventCounter << endl;
    cout << endl;
    cout << "------- Count Complete -------" << endl;
    inFile->Close();
    delete inFile;

    return make_pair(ratio, eventCounter);
}
pair<double, int> countingCoinc(const char* inFileName,
    float cut = 0,
    bool cutBool = false,
    int coincidence = 0) {
    cout << "------- Counting (Cuts Same as Jing) -------" << endl;

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return make_pair(0.0, 0);
    }

    TTree* tree = (TTree*)inFile->Get("Tree");
    if (!tree) {
        cout << "Error: Could not find the tree!" << endl;
        return make_pair(0.0, 0);
    }

    int CCFilter, PVFilter, nVtx;
    float HFEMaxPlus, HFEMaxMinus;
    float HFEMaxPlus2, HFEMaxMinus2, HFEMaxPlus3, HFEMaxMinus3;
    float VZ, leadingPtEta1p0_sel;

    tree->SetBranchAddress("ClusterCompatibilityFilter", &CCFilter);
    tree->SetBranchAddress("PVFilter", &PVFilter);
    tree->SetBranchAddress("nVtx", &nVtx);
    tree->SetBranchAddress("VZ", &VZ);
    tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
    tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);

    tree->SetBranchAddress("HFEMaxPlus2", &HFEMaxPlus2);
    tree->SetBranchAddress("HFEMaxMinus2", &HFEMaxMinus2);
    tree->SetBranchAddress("HFEMaxPlus3", &HFEMaxPlus3);
    tree->SetBranchAddress("HFEMaxMinus3", &HFEMaxMinus3);
    tree->SetBranchAddress("leadingPtEta1p0_sel",&leadingPtEta1p0_sel);


    const Long64_t nEvt = tree->GetEntries();

    int eventCounter = 0, subeventCounter = 0;
    int eventCounterTotal = 0;
    int CC0Counter = 0, CC1Counter = 0, PV0Counter = 0, PV1Counter = 0, nVtxCounter = 0, VZCounter = 0;

    for (int i = 0; i < nEvt; i++) {
        tree->GetEntry(i);
        if (CCFilter == 1 && PVFilter == 1 && abs(VZ) < 15) {
            eventCounterTotal ++;
            if (cutBool == false) { // OR logic
                if (coincidence == 1) {
                    if (HFEMaxPlus > cut || HFEMaxMinus > cut) {
                        eventCounter++;
                    }
                } else if (coincidence == 2) {
                    if ((HFEMaxPlus > cut || HFEMaxMinus > cut) &&
                        (HFEMaxPlus2 > cut || HFEMaxMinus2 > cut)) {
                        eventCounter++;
                    }
                }
            } else if (cutBool == true) { // AND logic
                if (coincidence == 1) {
                    if (HFEMaxPlus > cut && HFEMaxMinus > cut) {
                        eventCounter++;
                    }
                } else if (coincidence == 2) {
                    if ((HFEMaxPlus > cut && HFEMaxMinus > cut) &&
                        (HFEMaxPlus2 > cut && HFEMaxMinus2 > cut)) {
                        eventCounter++;
                    }
                }
            }
        }
    }
    double ratio = static_cast<double>(eventCounter) / nEvt;
    cout << "Ratio between Original and Filtered: " << ratio << endl;
    cout << "Total number of events: " << nEvt << endl;
    cout << "Total number of events after HFEcut "<< cut << ": " <<eventCounter << endl;
    cout << endl;

    cout << "------- Count Complete -------" << endl;
    inFile->Close();
    delete inFile;
    return make_pair(ratio, eventCounter);
}
pair<double, int> countingTrkptNPart(const char* inFileName, 
    float HFEcut = 0, 
    bool HFEcutBool = false, 
    bool HijingBool = false,
    int trkptcut = 0) {
    cout << Form("------- Counting for HFEMax cut %f and ptcut %d -------",HFEcut,trkptcut) << endl;

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return make_pair(0.0, 0);
    }

    TTree* tree = (TTree*)inFile->Get("Tree");
    if (!tree) {
        cout << "Error: Could not find the tree!" << endl;
        return make_pair(0.0, 0);
    }
    
    int CCFilter, PVFilter, nVtx;
    float HFEMaxPlus, HFEMaxMinus;
    float VZ, leadingPtEta1p0_sel, Npart, Ncoll;
    vector<float>* trkPt = nullptr;
    const Long64_t nEvt = tree->GetEntries();
    int eventCounter = 0, subeventCounter = 0;
    int denominator = 0;
    bool trkPtBool = false;

    tree->SetBranchAddress("ClusterCompatibilityFilter", &CCFilter);
    tree->SetBranchAddress("PVFilter", &PVFilter);
    tree->SetBranchAddress("nVtx", &nVtx);
    tree->SetBranchAddress("VZ", &VZ);
    tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
    tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
    tree->SetBranchAddress("leadingPtEta1p0_sel",&leadingPtEta1p0_sel);
    tree->SetBranchAddress("trkPt",&trkPt); 
    tree->SetBranchAddress("Npart", &Npart);
    tree->SetBranchAddress("Ncoll", &Ncoll);


    for (int i = 0; i < nEvt; i++) {
        tree->GetEntry(i);
        trkPtBool = false;
        if (trkptcut == -1){
            trkPtBool = true;
        }
        for (int j = 0; j < trkPt->size(); j++) {
            if (trkPt->at(j) > trkptcut) {
                trkPtBool = true;
            }
        }
        if (trkPtBool) {
            denominator++;
            if (CCFilter == 1 && PVFilter == 1 && nVtx > 0 && abs(VZ) < 15  && (!HijingBool || Npart > 1)) {
                if (HFEcutBool == false){
                    if (HFEMaxPlus > HFEcut || HFEMaxMinus > HFEcut){
                        eventCounter++;
                    }
                }
                else if (HFEcutBool == true){
                    if (HFEMaxPlus > HFEcut && HFEMaxMinus > HFEcut){
                        eventCounter++;
                    }
                }
            }
        }
    }

    double ratio = static_cast<double>(eventCounter) / denominator;
    cout << "Ratio between Original and Filtered: " << ratio << endl;
    cout << "Total number of events: " << nEvt << endl;
    cout << "Denominator after trkptcuts: " << denominator << ": " << denominator << endl;
    cout << "Total number of events after HFEcut "<< HFEcut << ": " << eventCounter << endl;
    cout << endl;
    cout << "------- Count Complete -------" << endl;
    inFile->Close();
    delete inFile;

    return make_pair(ratio, eventCounter);
}
pair<double, int> countingTrkptNPartAsymmetric(const char* inFileName, 
    float HFEcutlarge = 0, 
    float HFEcutsmall = 0,
    bool HFEcutBool = false, 
    bool HijingBool = false,
    int trkptcut = 0) {
    cout << Form("------- Counting for HFEMax+ %f and HFEMax- %f and ptcut %d -------",HFEcutlarge,HFEcutsmall,trkptcut) << endl;

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return make_pair(0.0, 0);
    }

    TTree* tree = (TTree*)inFile->Get("Tree");
    if (!tree) {
        cout << "Error: Could not find the tree!" << endl;
        return make_pair(0.0, 0);
    }
    
    int CCFilter, PVFilter, nVtx;
    float HFEMaxPlus, HFEMaxMinus;
    float VZ, Npart, Ncoll;
    vector<float>* trkPt = nullptr;
    const Long64_t nEvt = tree->GetEntries();
    int eventCounter = 0, subeventCounter = 0;
    int denominator = 0;
    bool trkPtBool = false;

    tree->SetBranchAddress("ClusterCompatibilityFilter", &CCFilter);
    tree->SetBranchAddress("PVFilter", &PVFilter);
    tree->SetBranchAddress("nVtx", &nVtx);
    tree->SetBranchAddress("VZ", &VZ);
    tree->SetBranchAddress("HFEMaxPlus", &HFEMaxPlus);
    tree->SetBranchAddress("HFEMaxMinus", &HFEMaxMinus);
    tree->SetBranchAddress("trkPt",&trkPt); 
    tree->SetBranchAddress("Npart", &Npart);
    tree->SetBranchAddress("Ncoll", &Ncoll);


    for (int i = 0; i < nEvt; i++) {
        tree->GetEntry(i);
        trkPtBool = false;
        if (trkptcut == -1){
            trkPtBool = true;
        }
        for (int j = 0; j < trkPt->size(); j++) {
            if (trkPt->at(j) > trkptcut) {
                trkPtBool = true;
            }
        }
        if (trkPtBool) {
            denominator++;
            if (CCFilter == 1 && PVFilter == 1 && nVtx > 0 && abs(VZ) < 15 && (!HijingBool || Npart > 1)){
                if (HFEcutBool == false){
                    if ((HFEMaxPlus > HFEcutlarge || HFEMaxMinus > HFEcutsmall)|| (HFEMaxPlus > HFEcutsmall || HFEMaxMinus > HFEcutlarge)){
                        eventCounter++;
                    }
                }
                else if (HFEcutBool == true){
                    if ((HFEMaxPlus > HFEcutlarge && HFEMaxMinus > HFEcutsmall)|| (HFEMaxPlus > HFEcutsmall && HFEMaxMinus > HFEcutlarge)){
                        eventCounter++;
                    }
                }
            }
        }
    }

    double ratio = static_cast<double>(eventCounter) / denominator;
    cout << "Ratio between Original and Filtered: " << ratio << endl;
    cout << "Total number of events: " << nEvt << endl;
    cout << "Denominator after trkptcuts: " << denominator << ": " << denominator << endl;
    cout << "Total number of events after HFEcut : " << eventCounter << endl;
    cout << endl;
    cout << "------- Count Complete -------" << endl;
    inFile->Close();
    delete inFile;

    return make_pair(ratio, eventCounter);
}

pair<double, int> countingTrkptVariable(
    const char*     inFileName, 
    const char*     varNamePlus,
    const char*     varNameMinus,
    float           cut = 0, 
    bool            cutBool = false, 
    bool            HijingBool = false,
    int             trkptcut = 0) {

    cout << Form("------- Counting for %s and %s > %f and ptcut %d -------",
        varNamePlus, varNameMinus, cut,trkptcut) << endl;
    cout << "File Name: " << inFileName << endl;

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return make_pair(0.0, 0);
    }

    TTree* tree = (TTree*)inFile->Get("Tree");
    if (!tree) {
        cout << "Error: Could not find the tree!" << endl;
        return make_pair(0.0, 0);
    }
    
    int CCFilter, PVFilter, nVtx;
    Float_t varPlusF, varMinusF;
    Int_t varPlusI, varMinusI;

    float VZ, leadingPtEta1p0_sel, Npart, Ncoll;
    vector<float>* trkPt = nullptr;
    const Long64_t nEvt = tree->GetEntries();
    int eventCounter = 0, subeventCounter = 0, denominator = 0;
    bool trkPtBool = false;

    tree->SetBranchAddress("ClusterCompatibilityFilter", &CCFilter);
    tree->SetBranchAddress("PVFilter", &PVFilter);
    tree->SetBranchAddress("nVtx", &nVtx);
    tree->SetBranchAddress("VZ", &VZ);
    tree->SetBranchAddress("leadingPtEta1p0_sel",&leadingPtEta1p0_sel);
    tree->SetBranchAddress("trkPt",&trkPt); 
    tree->SetBranchAddress("Npart", &Npart);
    tree->SetBranchAddress("Ncoll", &Ncoll);
    
    TLeaf *leafPlus  = tree->GetLeaf(varNamePlus);   // e.g. "HFEMaxPlus"
    TLeaf *leafMinus = tree->GetLeaf(varNameMinus);  // e.g. "HFEMaxMinus"
    double varPlus, varMinus;
    for (int i = 0; i < nEvt; i++) {
        tree->GetEntry(i);

        varPlus  = leafPlus ->GetValue();   // converted to Double_t automatically
        varMinus = leafMinus->GetValue();

        trkPtBool = false;
        if (trkptcut == -1){
            trkPtBool = true;
        }
        for (int j = 0; j < trkPt->size(); j++) {
            if (trkPt->at(j) > trkptcut) {
                trkPtBool = true;
            }
        }
        if (trkPtBool) {
            denominator++;
            if (CCFilter == 1 && PVFilter == 1 && nVtx > 0 && abs(VZ) < 15  && (!HijingBool || Npart > 1)) {
                if (cutBool == false){
                    if (varPlus > cut || varMinus > cut){
                        eventCounter++;
                    }
                }
                else if (cutBool == true){
                    if (varPlus > cut && varMinus > cut){
                        eventCounter++;
                    }
                }
            }
        }
    }

    double ratio = static_cast<double>(eventCounter) / denominator;
    cout << "Selection Efficiency: " << ratio << endl;
    cout << "Number of Event: " << nEvt << endl;
    cout << "Number of Event after ptCuts: " << denominator << ": " << denominator << endl;
    cout << "Number of event after cuts "<< cut << ": " << eventCounter << endl;
    cout << endl;
    cout << "------- Count Complete -------" << endl;
    inFile->Close();
    delete inFile;
    return make_pair(ratio, eventCounter);
}