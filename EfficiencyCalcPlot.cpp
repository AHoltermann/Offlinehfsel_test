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
#include "Include/GraphingUtils.h"
#include "Include/HistVariables.h"
#include "Include/EfficiencyResults.h"
#include "Include/EfficiencyCounting.h"
#include "Include/DataFilePaths.h"
#include "Include/Color.h"

using namespace std;

EfficiencyResults ReturnEfficiency(const char* inFilename);
EfficiencyResults ReturnEfficiency1bkg(const char* inFilename);
EfficiencyResults ReturnEfficiency_asymmetric(const char* inFilename);
EfficiencyResults ReturnEfficiency_aO(const char* inFilename);
EfficiencyResults ReturnEfficiency_root(const char* inFilename);
void PlotLeadingPtCut(const char* inFileName, const char* inFileshort);
void drawROCsame(
    const vector<vector<double>>& xVecs,
    const vector<vector<double>>& yVecs,
    EfficiencyResults results,
    HistVar1D hvar,
    const vector<std::string>& legends// legends[0], legends[1]
);
double CalculatePurity(vector<double> xsec, 
                    vector<double> Eff,
                    vector<double> backgroundfraction);
void drawROC(const char* fileNameShort,vector<double> x1, vector<double> y1, 
    EfficiencyResults results,
    HistVar1D hvar);
void drawROCmulti(
    const vector<vector<double>>& xVecs,
    const vector<vector<double>>& yVecs,
    const vector<string>& labels,
    const map<string,float>& xsec,
    const HistVar1D& hvar,
    bool npartcut = false);
void drawROCmultiShapes(
    const vector<EfficiencyResults>& resultsVec,
    const vector<vector<double>>& effVecs,
    const vector<vector<double>>& purityVecs,
    const vector<string>& labels,
    const map<string,float>& xsec,
    const HistVar1D& hvar,
    const char* varName,
    string andor = "and",
    bool npartcut = false);
void drawROCasymm(const EfficiencyResults& res,
    HistVar1D hvar,
    string andor = "and");
EfficiencyResults SelectEfficiencyResultsByX(const EfficiencyResults& input, 
    const std::vector<double>& wanted_x);
    
void run(){

    //std::vector<double> wanted_x = {0,6,8,10,12,14,16,18};
   // std::vector<double> wanted_x = {0,6,10,14,16,18,20,22};
    std::vector<double> wanted_x = {0,4,6,8,10,12};

    EfficiencyResults res_Offline_allbkg= ReturnEfficiency_root(Offline_allbkg);
    EfficiencyResults res_Offline_allbkg_bestvalues = SelectEfficiencyResultsByX(res_Offline_allbkg, wanted_x);
    EfficiencyResults res_Online_allbkg= ReturnEfficiency_root(Online_allbkg);
    EfficiencyResults res_Online_allbkg_bestvalues = SelectEfficiencyResultsByX(res_Online_allbkg, wanted_x);

    EfficiencyResults res_trkpt0_Offline = ReturnEfficiency_root(trkpt0_Offline);
    EfficiencyResults res_trkpt0_Offline_bestvalues = SelectEfficiencyResultsByX(res_trkpt0_Offline, wanted_x);
    EfficiencyResults res_trkpt2_Offline = ReturnEfficiency_root(trkpt2_Offline);
    EfficiencyResults res_trkpt2_Offline_bestvalues = SelectEfficiencyResultsByX(res_trkpt2_Offline, wanted_x);
    EfficiencyResults res_trkpt4_Offline = ReturnEfficiency_root(trkpt4_Offline);
    EfficiencyResults res_trkpt4_Offline_bestvalues = SelectEfficiencyResultsByX(res_trkpt4_Offline, wanted_x);
    EfficiencyResults res_trkpt6_Offline = ReturnEfficiency_root(trkpt6_Offline);
    EfficiencyResults res_trkpt6_Offline_bestvalues = SelectEfficiencyResultsByX(res_trkpt6_Offline, wanted_x);

    EfficiencyResults res_trkpt0_Online = ReturnEfficiency_root(trkpt0_Online);
    EfficiencyResults res_trkpt0_Online_bestvalues = SelectEfficiencyResultsByX(res_trkpt0_Online, wanted_x);
    EfficiencyResults res_trkpt2_Online = ReturnEfficiency_root(trkpt2_Online);
    EfficiencyResults res_trkpt2_Online_bestvalues = SelectEfficiencyResultsByX(res_trkpt2_Online, wanted_x);
    EfficiencyResults res_trkpt4_Online = ReturnEfficiency_root(trkpt4_Online);
    EfficiencyResults res_trkpt4_Online_bestvalues = SelectEfficiencyResultsByX(res_trkpt4_Online, wanted_x);
    EfficiencyResults res_trkpt6_Online = ReturnEfficiency_root(trkpt6_Online);
    EfficiencyResults res_trkpt6_Online_bestvalues = SelectEfficiencyResultsByX(res_trkpt6_Online, wanted_x);


    map<string,float> xsec;
    xsec["Hijing"] = 1.3;
    xsec["Starlight_SD"] = 0.6;
    xsec["Starlight_DD"] = 0.0003;
    xsec["AlphaO"] = 1.6;

    vector<EfficiencyResults> res_vector_Offline = {
        res_Offline_allbkg,
        res_trkpt0_Offline,
        res_trkpt2_Offline,
        res_trkpt4_Offline
    };

    vector<EfficiencyResults> res_vector_Online = {
        res_Online_allbkg,
        res_trkpt0_Online,
        res_trkpt2_Online,
        res_trkpt4_Online,
        res_trkpt6_Online
    };

    vector<EfficiencyResults> res_vector_bestvalues_Offline = {
        res_Offline_allbkg_bestvalues,
        res_trkpt0_Offline_bestvalues,
        res_trkpt2_Offline_bestvalues,
        res_trkpt4_Offline_bestvalues
    };

    vector<EfficiencyResults> res_vector_bestvalues_Online = {
        res_Online_allbkg_bestvalues,
        res_trkpt0_Online_bestvalues,
        res_trkpt2_Online_bestvalues,
        res_trkpt4_Online_bestvalues,
        res_trkpt6_Online_bestvalues
    };

    //
    const auto& res_vector = res_vector_bestvalues_Offline;

    // Build effVecs_OR and effVecs_AND from the selected res_vector_bestvalues
    vector<vector<double>> effVecs_OR, effVecs_AND;
    for (const auto& res : res_vector) {
        effVecs_OR.push_back(res.EfficiencyHijingOR);
        effVecs_AND.push_back(res.EfficiencyHijingAND);
    }

    vector<double> vPurity_and, vPurity_or;
    vector<vector<double>> vPurityAND, vPurityOR;
    double Purity_and = 0;
    double Purity_or = 0;

    // Loop over each EfficiencyResults in res_vector
    for (size_t idx = 0; idx < res_vector.size(); ++idx) {
        auto& res = res_vector[idx];
        vPurity_and.clear();
        vPurity_or.clear();
        // Only print detailed output for the first res in the vector
        for (int i = 0; i < res.x.size(); i++) {
            cout << "Trkptcut: " << idx*2 << endl;
            cout << "Energy Cut: " << res.x[i] << endl;
            cout << "Efficiency Hijing AND: " << res.EfficiencyHijingAND[i] << endl;
            cout << "Efficiency SD AND: " << res.EfficiencySDAND[i] << endl;
            cout << "Efficiency DD AND: " << res.EfficiencyDDAND[i] << endl;
            cout << "Efficiency Hijing OR: " << res.EfficiencyHijingOR[i] << endl;
            cout << "Efficiency SD OR: " << res.EfficiencySDOR[i] << endl;
            cout << "Efficiency DD OR: " << res.EfficiencyDDOR[i] << endl;

            
            Purity_and = CalculatePurity({xsec["Hijing"], xsec["Starlight_SD"], xsec["Starlight_DD"], xsec["AlphaO"]},
                        {res.EfficiencyHijingAND[i], res.EfficiencySDAND[i],res.EfficiencyDDAND[i],res.EfficiencyAOAND[i]},
                        {1,1,1, 0});
            Purity_or = CalculatePurity({xsec["Hijing"], xsec["Starlight_SD"], xsec["Starlight_DD"], xsec["AlphaO"]},
                        {res.EfficiencyHijingOR[i], res.EfficiencySDOR[i],res.EfficiencyDDOR[i],res.EfficiencyAOOR[i]},
                        {1,1,1, 0});

            vPurity_and.push_back(Purity_and);
            vPurity_or.push_back(Purity_or);

            cout << "Purity AND for x = " << res.x[i] << ": " << Purity_and << endl;
            cout << "Purity OR for x = " << res.x[i] << ": " << Purity_or << endl;
        }
        vPurityAND.push_back(vPurity_and);
        vPurityOR.push_back(vPurity_or);
    }
    //xsec["AlphaO"] = 0; // Set AlphaO to 0 for purity calculations
     


   /* vector<double> vPurity_percent_and,vPurity_percent_or;
    vector<vector<double>> vPurity_AND,vPurity_OR;
    vector<vector<double>> vPurity_sel;
    double Purity_and = 0;
    double Purity_or = 0;

    EfficiencyResults res = res_Offline_allbkg_bestvalues;

    for (float percent = 0; percent <= 1.1; percent += 0.1){
        vPurity_percent_and.clear();
        vPurity_percent_or.clear();

        cout << "Calculating purity for percent: " << percent << endl;
        for (int i = 0; i < res.x.size(); i++) {
            cout << endl;
            cout << "----------------------------" << endl;
            cout << "Energy Cut: " << res.x[i] << endl;
            cout << "Percent Contamination: " << percent << endl;
            cout << endl;
            Purity_and = CalculatePurity({xsec["Hijing"], xsec["Starlight_SD"], xsec["Starlight_DD"], xsec["AlphaO"]},
                                    {res.EfficiencyHijingAND[i], res.EfficiencySDAND[i],res.EfficiencyDDAND[i],res.EfficiencyAOAND[i]},
                                    {1,1,1, percent});
            Purity_or = CalculatePurity({xsec["Hijing"], xsec["Starlight_SD"], xsec["Starlight_DD"], xsec["AlphaO"]},
                                    {res.EfficiencyHijingOR[i], res.EfficiencySDOR[i],res.EfficiencyDDOR[i],res.EfficiencyAOOR[i]},
                                    {1,1,1, percent});

            vPurity_percent_and.push_back(Purity_and);
            vPurity_percent_or.push_back(Purity_or);
        }
        vPurity_AND.push_back(vPurity_percent_and);
        vPurity_OR.push_back(vPurity_percent_or);
    }
*/

    HistVar1D hvar;
    //hvar.histTitle = "ROC Curve Online mMaxL1HFAdc Selection No trkpt selection No AlphaO";

    hvar.xLabel = "Efficiency";
    hvar.yLabel = "Purity";
    hvar.xmin = 0.8;
    hvar.xmax = 1;
    hvar.ymin = 0.2;
    hvar.ymax = 1.2;
    hvar.legendTitle = "";
    hvar.outFolderName = "Plots/FilterEfficiencyFinalPlots/";
   
    hvar.histTitle = "ROC Curve Offline HFAND Selection";
    hvar.outFileName = "OfflineHFEMax_trkpt_AND_bestvalues";

  drawROCmultiShapes(res_vector,
        effVecs_AND,
        vPurityAND,
        {"no trkpt selection", "track pt > 0 GeV", "track pt > 2 GeV", "track pt > 4 GeV"},
        xsec,
        hvar,
        "HFEmax",
        "and",
        true);

}

EfficiencyResults ReturnEfficiency_root(const char* inFilename){
    TFile* file = TFile::Open(inFilename);
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file!" << endl;
        return EfficiencyResults(); // Return empty results
    }
    // Get the tree
    TTree* tree = (TTree*)file->Get("effTree");
    if (!tree) {
        cerr << "Error: Cannot find tree!" << endl;
        return EfficiencyResults(); // Return esmpty results
    }

    // Variables to read5
    double Cut, Eff_HijingAND, Eff_HijingOR, Eff_SDAND, Eff_SDOR;
    double Eff_DDAND, Eff_DDOR, Eff_AlphaOAND, Eff_AlphaOOR, PurityAND, PurityOR;
    EfficiencyResults results;
    // Set branch addresses only if the branch exists
    #define SET_BRANCH_IF_EXISTS(branch, ptr) \
        if (tree->GetBranch(branch)) tree->SetBranchAddress(branch, ptr);

    SET_BRANCH_IF_EXISTS("Cut", &Cut);
    SET_BRANCH_IF_EXISTS("Eff_HijingAND", &Eff_HijingAND);
    SET_BRANCH_IF_EXISTS("Eff_HijingOR", &Eff_HijingOR);
    SET_BRANCH_IF_EXISTS("Eff_SDAND", &Eff_SDAND);
    SET_BRANCH_IF_EXISTS("Eff_SDOR", &Eff_SDOR);
    SET_BRANCH_IF_EXISTS("Eff_DDAND", &Eff_DDAND);
    SET_BRANCH_IF_EXISTS("Eff_DDOR", &Eff_DDOR);
    SET_BRANCH_IF_EXISTS("Eff_AlphaOAND", &Eff_AlphaOAND);
    SET_BRANCH_IF_EXISTS("Eff_AlphaOOR", &Eff_AlphaOOR);
    SET_BRANCH_IF_EXISTS("PurityAND", &PurityAND);
    SET_BRANCH_IF_EXISTS("PurityOR", &PurityOR);

    #undef SET_BRANCH_IF_EXISTS

    // Loop over entries
    Long64_t nentries = tree->GetEntries();
    cout << "Number of entries in the tree: " << nentries << endl;
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        results.x.push_back(Cut);
        results.EfficiencyHijingAND.push_back(Eff_HijingAND);
        results.EfficiencyHijingOR.push_back(Eff_HijingOR);
        results.EfficiencySDAND.push_back(Eff_SDAND);
        results.EfficiencySDOR.push_back(Eff_SDOR);
        results.EfficiencyDDAND.push_back(Eff_DDAND);
        results.EfficiencyDDOR.push_back(Eff_DDOR);
        results.EfficiencyAOAND.push_back(Eff_AlphaOAND);
        results.EfficiencyAOOR.push_back(Eff_AlphaOOR);
        cout << Cut << " ";
    }

    file->Close();
    delete file;
    return results;

}

EfficiencyResults SelectEfficiencyResultsByX(const EfficiencyResults& input, 
    const std::vector<double>& wanted_x) {
    EfficiencyResults selected;
    for (size_t i = 0; i < input.x.size(); ++i) {
        double xval = input.x[i];
        if (std::find(wanted_x.begin(), wanted_x.end(), static_cast<int>(xval)) != wanted_x.end()) {
            selected.x.push_back(xval);
            selected.EfficiencyHijingAND.push_back(input.EfficiencyHijingAND[i]);
            selected.EfficiencyHijingOR.push_back(input.EfficiencyHijingOR[i]);
            selected.EfficiencySDAND.push_back(input.EfficiencySDAND[i]);
            selected.EfficiencySDOR.push_back(input.EfficiencySDOR[i]);
            selected.EfficiencyDDAND.push_back(input.EfficiencyDDAND[i]);
            selected.EfficiencyDDOR.push_back(input.EfficiencyDDOR[i]);
            selected.EfficiencyAOAND.push_back(input.EfficiencyAOAND[i]);
            selected.EfficiencyAOOR.push_back(input.EfficiencyAOOR[i]);
        }
    }
    return selected;
}

EfficiencyResults ReturnEfficiency(const char* inFilename){

    float xsec_SD, xsec_DD, xsec_had, xMax;
    int N,trkptcut = 0; // Default to 0, can be changed if needed
    EfficiencyResults results;

    vector<double> x, EfficiencyHijingAND, EfficiencyHijingOR;
    vector<double> EfficiencySDAND, EfficiencySDOR;
    vector<double> BKGRejectionAND, BKGRejectionOR;
    vector<double> PurityAND, PurityOR;
    string signalFile, backgroundFile;

    ifstream infile(inFilename);
    string line, dummy;
    if (!infile.is_open()) {
        cerr << "Error: Could not open file " << inFilename << endl;
        return results; // Return empty results
    }

    getline(infile, line);
    istringstream iss1(line);
    iss1 >> dummy >> signalFile;

    getline(infile, line);
    istringstream iss2(line);
    iss2 >> dummy >> backgroundFile;

    getline(infile, line);
    istringstream iss3(line);
    iss3 >> dummy >> backgroundFile;

    getline(infile, line);
    istringstream iss(line);
    iss >> dummy >> xsec_SD >> dummy >> xsec_DD >> dummy >> xsec_had >> dummy >> N >> dummy >> xMax >> dummy >> trkptcut;
    results.xsec_SD = xsec_SD;
    results.xsec_DD = xsec_DD;
    results.xsec_had = xsec_had;
    results.N = N;
    results.xMax = xMax;
    cout << "xsec_SD: " << xsec_SD << ", xsec_DD: " << xsec_DD << ", xsec_had: " << xsec_had << ", N: " << N << ", xMax: " << xMax << ", trkptcut: " << trkptcut << endl;
    getline(infile, line);

    while (getline(infile, line)) {
        istringstream iss(line);
        double tx, tEffHijingAND, tEffHijingOR, tEffSDAND, tEffSDOR, tEffDDAND, tEffDDOR, tBKGRejAND, tBKGRejOR, tPurityAND, tPurityOR;
        iss >> tx
            >> tEffHijingAND
            >> tEffHijingOR
            >> tEffSDAND
            >> tEffSDOR
            >> tEffDDAND
            >> tEffDDOR
            >> tBKGRejAND
            >> tBKGRejOR
            >> tPurityAND
            >> tPurityOR;

        results.x.push_back(tx);
        results.EfficiencyHijingAND.push_back(tEffHijingAND);
        results.EfficiencyHijingOR.push_back(tEffHijingOR);
        results.EfficiencySDAND.push_back(tEffSDAND);
        results.EfficiencySDOR.push_back(tEffSDOR);
        results.EfficiencyDDAND.push_back(tEffDDAND);
        results.EfficiencyDDOR.push_back(tEffDDOR);
        results.BKGRejectionAND.push_back(tBKGRejAND);
        results.BKGRejectionOR.push_back(tBKGRejOR);
        results.PurityAND.push_back(tPurityAND);
        results.PurityOR.push_back(tPurityOR);
    }
    infile.close();
    return results;
}

EfficiencyResults ReturnEfficiency_aO(const char* inFilename){

    float xsec_SD, xsec_DD, xsec_had, xsec_alphaO, xMax;
    int N,trkptcut = 0; // Default to 0, can be changed if needed
    EfficiencyResults results;

    vector<double> x, EfficiencyHijingAND, EfficiencyHijingOR;
    vector<double> EfficiencySDAND, EfficiencySDOR;
    vector<double> BKGRejectionAND, BKGRejectionOR;
    vector<double> PurityAND, PurityOR;
    string signalFile, backgroundFile;

    ifstream infile(inFilename);
    string line, dummy;
    if (!infile.is_open()) {
        cerr << "Error: Could not open file " << inFilename << endl;
        return results; // Return empty results
    }

    getline(infile, line);
    istringstream iss1(line);
    iss1 >> dummy >> signalFile;

    getline(infile, line);
    iss1 >> dummy >> backgroundFile;

    getline(infile, line);
    iss1 >> dummy >> backgroundFile;

    getline(infile, line);
    iss1 >> dummy >> backgroundFile;

    getline(infile, line);
    istringstream iss(line);
    iss >> dummy >> xsec_SD >> dummy >> xsec_DD >> dummy >> xsec_had >> dummy >> xsec_alphaO >> dummy >> N >> dummy >> xMax >> dummy >> trkptcut;
    results.xsec_SD = xsec_SD;
    results.xsec_DD = xsec_DD;
    results.xsec_had = xsec_had;
    results.xsec_alphaO = xsec_alphaO;
    results.N = N;
    results.xMax = xMax;
    cout << "xsec_SD: " << xsec_SD << ", xsec_DD: " << xsec_DD << ", xsec_had: " << xsec_had 
    << ", xsec_alphaO: " << xsec_alphaO
    << ", N: " << N << ", xMax: " << xMax << ", trkptcut: " << trkptcut << endl;
    getline(infile, line);

    while (getline(infile, line)) {
        istringstream iss(line);
        double tx, tEffHijingAND, tEffHijingOR, tEffSDAND, tEffSDOR, tEffDDAND, tEffDDOR, tEffaOAND, tEffaOOR, tPurityAND, tPurityOR;
        iss >> tx
            >> tEffHijingAND
            >> tEffHijingOR
            >> tEffSDAND
            >> tEffSDOR
            >> tEffDDAND
            >> tEffDDOR
            >> tEffaOAND
            >> tEffaOOR
            >> tPurityAND
            >> tPurityOR;

        results.x.push_back(tx);
        results.EfficiencyHijingAND.push_back(tEffHijingAND);
        results.EfficiencyHijingOR.push_back(tEffHijingOR);
        results.EfficiencySDAND.push_back(tEffSDAND);
        results.EfficiencySDOR.push_back(tEffSDOR);
        results.EfficiencyDDAND.push_back(tEffDDAND);
        results.EfficiencyDDOR.push_back(tEffDDOR);
        results.EfficiencyAOAND.push_back(tEffaOAND);
        results.EfficiencyAOOR.push_back(tEffaOOR);
        results.PurityAND.push_back(tPurityAND);
        results.PurityOR.push_back(tPurityOR);

    }
    infile.close();
    return results;
}


EfficiencyResults ReturnEfficiency_asymmetric(const char* inFilename){

    float xsec_SD, xsec_DD, xsec_had;
    int N,trkptcut = 0; // Default to 0, can be changed if needed
    EfficiencyResults results;

    vector<double> xplus, xminus, EfficiencyHijingAND, EfficiencyHijingOR;
    vector<double> EfficiencySDAND, EfficiencySDOR;
    vector<double> BKGRejectionAND, BKGRejectionOR;
    vector<double> PurityAND, PurityOR;
    string signalFile, backgroundFile;

    ifstream infile(inFilename);
    string line, dummy;
    if (!infile.is_open()) {
        cerr << "Error: Could not open file " << inFilename << endl;
        return results; // Return empty results
    }


    getline(infile, line);
    istringstream iss1(line);
    iss1 >> dummy >> signalFile;

    getline(infile, line);
    istringstream iss2(line);
    iss2 >> dummy >> backgroundFile;

    getline(infile, line);
    istringstream iss3(line);
    iss3 >> dummy >> backgroundFile;

    getline(infile, line);
    istringstream iss(line);
    iss >> dummy >> xsec_SD >> dummy >> xsec_DD >> dummy >> xsec_had >> dummy >> trkptcut;
    results.xsec_SD = xsec_SD;
    results.xsec_DD = xsec_DD;
    results.xsec_had = xsec_had;

    cout << "xsec_SD: " << xsec_SD << ", xsec_DD: " << xsec_DD << ", xsec_had: " << xsec_had << ", trkptcut: " << trkptcut << endl;
    getline(infile, line);

    while (getline(infile, line)) {
        istringstream iss(line);
        double txplus, txminus, tEffHijingAND, tEffHijingOR, tEffSDAND, tEffSDOR, tEffDDAND, tEffDDOR, tBKGRejAND, tBKGRejOR, tPurityAND, tPurityOR;
        iss >> txplus
            >> txminus
            >> tEffHijingAND
            >> tEffHijingOR
            >> tEffSDAND
            >> tEffSDOR
            >> tEffDDAND
            >> tEffDDOR
            >> tBKGRejAND
            >> tBKGRejOR
            >> tPurityAND
            >> tPurityOR;
        cout << "txplus: " << txplus << ", txminus: " << txminus << endl;
        results.xplus.push_back(txplus);
        results.xminus.push_back(txminus);
        results.EfficiencyHijingAND.push_back(tEffHijingAND);
        results.EfficiencyHijingOR.push_back(tEffHijingOR);
        results.EfficiencySDAND.push_back(tEffSDAND);
        results.EfficiencySDOR.push_back(tEffSDOR);
        results.BKGRejectionAND.push_back(tBKGRejAND);
        results.BKGRejectionOR.push_back(tBKGRejOR);
        results.PurityAND.push_back(tPurityAND);
        results.PurityOR.push_back(tPurityOR);
    }
    infile.close();
    return results;
}


EfficiencyResults ReturnEfficiency1bkg(const char* inFilename){

    ifstream infile(inFilename);
    string line, dummy;
    float xsec_SD, xsec_had, xMax;
    int N,coincidence = 0; // Default to 0, can be changed if needed

    EfficiencyResults results;

    vector<double> x, EfficiencyHijingAND, EfficiencyHijingOR;
    vector<double> EfficiencySDAND, EfficiencySDOR;
    vector<double> BKGRejectionAND, BKGRejectionOR;
    vector<double> PurityAND, PurityOR;
    string signalFile, backgroundFile;

    getline(infile, line);
    istringstream iss1(line);
    iss1 >> dummy >> signalFile;
    cout << "Signal File: " << signalFile << endl;

    getline(infile, line);
    istringstream iss2(line);
    iss2 >> dummy >> backgroundFile;

    cout << "Background File SD: " << backgroundFile << endl;

    getline(infile, line);
    istringstream iss(line);
    iss >> dummy >> xsec_SD >> dummy >> xsec_had >> dummy >> N >> dummy >> xMax >> dummy >> coincidence;
    results.xsec_SD = xsec_SD;
    results.xsec_had = xsec_had;
    results.N = N;
    results.xMax = xMax;
    results.coincidence = coincidence;
    cout << "xsec_EM: " << xsec_SD << ", xsec_had: " << xsec_had << ", N: " << N << ", xMax: " << xMax << ", coincidence: " << coincidence << endl;
    getline(infile, line);

    for (int i = 0; i < N && getline(infile, line); ++i) {
        istringstream iss(line);
        double tx, tEffHijingAND, tEffHijingOR, tEffSDAND, tEffSDOR, tBKGRejAND, tBKGRejOR, tPurityAND, tPurityOR;
        iss >> tx
            >> tEffHijingAND
            >> tEffHijingOR
            >> tEffSDAND
            >> tEffSDOR
            >> tBKGRejAND
            >> tBKGRejOR
            >> tPurityAND
            >> tPurityOR;
    
        results.x.push_back(tx);
        results.EfficiencyHijingAND.push_back(tEffHijingAND);
        results.EfficiencyHijingOR.push_back(tEffHijingOR);
        results.EfficiencySDAND.push_back(tEffSDAND);
        results.EfficiencySDOR.push_back(tEffSDOR);
        results.BKGRejectionAND.push_back(tBKGRejAND);
        results.BKGRejectionOR.push_back(tBKGRejOR);
        results.PurityAND.push_back(tPurityAND);
        results.PurityOR.push_back(tPurityOR);
    }
    infile.close();
    return results;
}


void PlotLeadingPtCut(const char* inFileName, const char* inFileshort) {
    
    int HFEmaxcut = 4;
    int xmax = 30;

    TH1F* h = new TH1F("h", Form("Efficiency of Leading Track p_{T} cuts (%s);Track p_{T} Cut;Efficiency", inFileshort), 30, 1, 31);
    for (int i = 0; i <= xmax; ++i) {
        double eff = countingTrkptcuts(inFileName, HFEmaxcut, true, i).first; // Assuming countingJingcuts returns a pair with efficiency as the first element
        h->SetBinContent(i, eff);
    }
    
    TCanvas* c = new TCanvas("cEff", "Efficiency vs Leading p_{T} Cut", 800, 600);
    h->Draw("hist");
    h->GetXaxis()->SetRangeUser(0, xmax);
    h->GetYaxis()->SetRangeUser(0.0, 1.1);
    gStyle->SetOptStat(0);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.6, 0.85, Form("HFE_max+ and - > %i", HFEmaxcut));
    latex.DrawLatex(0.23, 0.43, "Cuts on Num only:");
    latex.DrawLatex(0.23, 0.40, "PVFilter == 1");
    latex.DrawLatex(0.23, 0.36, "ClusterCompatibilityFilter == 1");
    latex.DrawLatex(0.23, 0.32, "abs(V_{z}) < 15");

    c->SaveAs(Form("Efficiency_vs_LeadingTrkPtCut_%s_HFEMax%i.png", inFileshort, HFEmaxcut));
}

double CalculatePurity(vector<double> xsec, 
                    vector<double> Eff,
                    vector<double> backgroundfraction) {
    double denominator =  0;
    for (int i = 0; i < xsec.size(); ++i) {
        denominator += xsec[i] * Eff[i] * backgroundfraction[i];
    }
    double Purity = (xsec[0] * Eff[0]) / denominator;
    return Purity;
}


// Draws two ROC curves with user-supplied legends as a vector<string>
void drawROCsame(
    const vector<vector<double>>& xVecs,
    const vector<vector<double>>& yVecs,
    EfficiencyResults results,
    HistVar1D hvar,
    const vector<std::string>& legends // legends[0], legends[1]
) {
    if (xVecs.size() != 2 || yVecs.size() != 2 || legends.size() != 2) {
        cerr << "drawROCsame: Need exactly 2 sets of points and 2 legends." << endl;
        return;
    }

    string outputDir = hvar.outFolderName;
    PlotUtils::setgstyle();

    cout << "------- Drawing ROC Curves -------" << endl;
    TCanvas* c = new TCanvas("cBoth", "", 800, 800);

    TGraph *gr1 = new TGraph(xVecs[0].size(), xVecs[0].data(), yVecs[0].data());
    TGraph *gr2 = new TGraph(xVecs[1].size(), xVecs[1].data(), yVecs[1].data());

    c->cd();
    gr1->SetTitle(hvar.histTitle.c_str());
    gr1->GetXaxis()->SetTitle(hvar.xLabel);
    gr1->GetYaxis()->SetTitle(hvar.yLabel);
    gr1->Draw("AP");
    gr2->Draw("P SAME");

    gr1->GetYaxis()->SetRangeUser(hvar.ymin, hvar.ymax);
    gr1->GetXaxis()->SetRangeUser(hvar.xmin, hvar.xmax);

    c->SetGrid();

    PlotUtils::setgraphstyle(gr1, 0.0, 0.0);
    PlotUtils::setgraphstyle(gr2, 0.0, 0.0);
    gr1->SetMarkerColor(kBlue);
    gr2->SetMarkerColor(kRed);

    PlotUtils::GraphLegend({
        {gr1, legends[0]},
        {gr2, legends[1]}
    }, 0.2, 0.55, 0.25, 0.65);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.5, 0.85, Form("HIJING OO(%.1f b)", results.xsec_had));
    latex.DrawLatex(0.5, 0.825, Form("Single Diffractive (%.1f b)", results.xsec_SD));
    latex.DrawLatex(0.5, 0.80, Form("Double Diffractive (%.5f b)", results.xsec_DD));
  // latex.DrawLatex(0.5, 0.775, Form("AlphaO (%.5f b)", results.xsec_alphaO));
    latex.DrawLatex(0.23, 0.40, "PVFilter == 1");
    latex.DrawLatex(0.23, 0.36, "ClusterCompatibilityFilter == 1");
    latex.DrawLatex(0.23, 0.32, "abs(V_{z}) < 15");
    c->SaveAs((outputDir + Form("PuritycurveBoth%i",results.N) + "_" + hvar.outFileName + ".png").c_str());
}

void drawROC(const char* fileNameShort, vector<double> x1, vector<double> y1, 
    EfficiencyResults results,
    HistVar1D hvar) {

    string outputDir = "/home/xirong/OOAnalysis_2025/Plots/FilterEfficiencyNewest";

    PlotUtils::setgstyle();

    cout << "------- Drawing ROC Curves -------" << endl;
    TCanvas* c = new TCanvas("cBoth", "", 800, 800);

    TGraph *gr1 = new TGraph(results.N, x1.data(), y1.data());

    c->cd();
    gr1->SetTitle(hvar.histTitle.c_str());
    gr1->GetXaxis()->SetTitle(hvar.xLabel);
    gr1->GetYaxis()->SetTitle(hvar.yLabel);
    gr1->Draw("AP");
   // gr1->GetYaxis()->SetRangeUser(hvar.ymin, hvar.ymax);
    gr1->GetXaxis()->SetRangeUser(hvar.xmin, hvar.xmax);
    c->SetGrid();

    PlotUtils::setgraphstyle(gr1, 0.0, 0.0);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.7, 0.85, Form("HIJING OO (%.1f b)", results.xsec_had));
    latex.DrawLatex(0.7, 0.82, Form("Starlight SD (%.1f b)", results.xsec_SD));
    latex.DrawLatex(0.7, 0.79, Form("Starlight DD (%.5f b)", results.xsec_DD));
    latex.DrawLatex(0.5, 0.5, Form("HFE_{max}^{+} and HFE_{max}^{-} > X"));

    latex.DrawLatex(0.23, 0.40, "PVFilter == 1");
    latex.DrawLatex(0.23, 0.36, "ClusterCompatibilityFilter == 1");
    latex.DrawLatex(0.23, 0.32, "abs(V_{z}) < 15");
    c->SaveAs((outputDir + Form("PuritycurveONE_Coinc%i_%i",results.coincidence,results.N) + "_" + hvar.outFileName + ".png").c_str());
}

void drawROCmulti(
    const vector<vector<double>>& xVecs,
    const vector<vector<double>>& yVecs,
    const vector<string>& labels,
    const map<string,float>& xsec,
    const HistVar1D& hvar,
    bool npartcut = false) {
    
    vector<int> colors = getRainbow2(xVecs.size());
    vector<int> markerStyles = {20, 21, 22, 23, 29,33,34,39,41,43,45,47};
    PlotUtils::setgstyle();

    TCanvas* c = new TCanvas("cROCmulti", "", 800, 800);
    vector<TGraph*> graphs;

    for (int i = 0; i < xVecs.size(); ++i) {
        int N = xVecs[i].size();
        TGraph* gr = new TGraph(N, xVecs[i].data(), yVecs[i].data());
        PlotUtils::setgraphstyle(gr, 0.0, 0.0);
        gr->SetMarkerColor(colors[i]);
        gr->SetLineColor(colors[i]);
       // gr->SetMarkerStyle(markerStyles[i % markerStyles.size()]);
        gr->SetMarkerSize(1.0);
        graphs.push_back(gr);
    }

    // Draw the first graph with axes, others with "P SAME"
    graphs[0]->SetTitle(hvar.histTitle.c_str());
    graphs[0]->GetXaxis()->SetTitle(hvar.xLabel);
    graphs[0]->GetYaxis()->SetTitle(hvar.yLabel);
    graphs[0]->GetYaxis()->SetRangeUser(hvar.ymin, hvar.ymax);
    graphs[0]->GetXaxis()->SetLimits(hvar.xmin, hvar.xmax); // or SetRangeUser
    c->cd();
    graphs[0]->Draw("AP");
    for (int i = 1; i < xVecs.size(); ++i) {
        graphs[i]->Draw("P SAME");
    }
    c->SetGrid();

    // Legend
    vector<pair<TObject*, string>> legendEntries;
    for (int i = 0; i < xVecs.size(); ++i) {
        legendEntries.push_back({graphs[i], labels[i]});
    }
    PlotUtils::GraphLegend(legendEntries, 0.7, 0.65, 0.9, 0.90,hvar.legendTitle.c_str());

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    double blockFilexsecx = 0.6;
    double blockFilexsecy = 0.85;
    double blockDY = 0.03;
    double xsec_had = xsec.at("Hijing");
    double xsec_SD = xsec.at("Starlight_SD");
    double xsec_DD = xsec.at("Starlight_DD");
    double xsec_alphaO = xsec.at("AlphaO");
    if (xsec_had != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy, Form("HIJING OO(%.1f b)", xsec_had));
    if (xsec_SD != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy- blockDY, Form("Starlight SD (%.1f b)", xsec_SD));
    if (xsec_DD != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy- 2*blockDY, Form("Starlight DD (%.5f b)", xsec_DD));
    if (xsec_alphaO != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy- 3*blockDY, Form("AlphaO (%.1f b)", xsec_alphaO));
    // Block for PVFilter etc. with adjustable y start
    double blockX = 0.23;
    double blockY = 0.5;
    latex.DrawLatex(blockX, blockY, "PVFilter == 1");
    latex.DrawLatex(blockX, blockY - blockDY, "ClusterCompatibilityFilter == 1");
    latex.DrawLatex(blockX, blockY - 2*blockDY, "abs(V_{z}) < 15");
    if (npartcut == true){
        latex.DrawLatex(blockX, blockY - 3*blockDY, "npart > 1 for HIJING");
    }
    // Save
    c->SaveAs((hvar.outFolderName + "Puritycurve_AllTrackptcuts_multi_" + hvar.outFileName + ".png").c_str());
}
void drawROCmultiShapes(
    const vector<EfficiencyResults>& resultsVec,
    const vector<vector<double>>& effVecs,
    const vector<vector<double>>& purityVecs,
    const vector<string>& labels,
    const map<string,float>& xsec,
    const HistVar1D& hvar,
    const char* varName,
    string andor = "and",
    bool npartcut = false) 
{
    vector<int> markerStyles = {20, 21, 22, 23, 29, 33, 34, 39, 41, 43, 45, 47};

    // 1. Collect all unique X values
    set<int> uniqueX;
    for (const auto& res : resultsVec) {
        for (double xval : res.x) {
            uniqueX.insert(static_cast<int>(xval));
        }
    }
    cout<< "Unique X values: ";
    for (int xval : uniqueX) {
        cout << xval << " ";
    }
    cout << endl;

    // 2. Generate rainbow palette
    int colorOffset = 4000 + (rand() % 1000);
    vector<int> colors = getRainbow(uniqueX.size(), colorOffset);

    // 3. Map each unique X value to a color
    map<int, int> x2color;
    size_t idx = 0;
    for (int xval : uniqueX) x2color[xval] = colors[idx++];
    
    PlotUtils::setgstyle();
    TCanvas* c = new TCanvas("cROCmulti", "", 800, 800);
    TMultiGraph* mg = new TMultiGraph();

    // 4. Plot, assigning color by X value
    for (int i = 0; i < resultsVec.size(); ++i) {
        const auto& res = resultsVec[i];
        for (int j = 0; j < res.x.size(); ++j) {
            int xint = static_cast<int>(res.x[j]);
            if (uniqueX.count(xint)) { // Only plot if X is in uniqueX (should always be true)
                TGraph* ptGraph = new TGraph(1, &effVecs[i][j], &purityVecs[i][j]);
                cout << "Adding point for X value: " << xint << " with efficiency: " << effVecs[i][j] << " and purity: " << purityVecs[i][j] << endl;
                int color = x2color[xint];
                ptGraph->SetMarkerColor(color);
                ptGraph->SetMarkerStyle(markerStyles[i % markerStyles.size()]);
                ptGraph->SetMarkerSize(1.5);
                mg->Add(ptGraph, "P");
            }
        }
    }

    mg->SetTitle(hvar.histTitle.c_str());
    mg->Draw("A");
    mg->GetXaxis()->SetTitle(hvar.xLabel);
    mg->GetYaxis()->SetTitle(hvar.yLabel);
    mg->GetXaxis()->SetLimits(hvar.xmin, hvar.xmax);
    mg->SetMinimum(hvar.ymin);
    mg->SetMaximum(hvar.ymax);
    c->SetGrid();

    // Color legend for X values
    mg->SetTitle(hvar.histTitle.c_str());
    mg->Draw("A");
    mg->GetXaxis()->SetTitle(hvar.xLabel);
    mg->GetYaxis()->SetTitle(hvar.yLabel);
    mg->GetXaxis()->SetLimits(hvar.xmin, hvar.xmax);
    mg->SetMinimum(hvar.ymin);
    mg->SetMaximum(hvar.ymax);
    c->SetGrid();

    vector<pair<TObject*, string>> colorLegendEntries;
    for (const auto& kv : x2color) {
        int xval = kv.first;
        cout << "Adding color legend entry for X value: " << xval << endl;
        int color = kv.second;
        TMarker* m = new TMarker(0, 0, 20);
        m->SetMarkerColor(color);
        m->SetMarkerSize(1.3);
        if (xval %1 == 0) {
            colorLegendEntries.push_back({m, Form("%s + %s - > %d", varName, andor.c_str(), xval)});
        }
    }
    PlotUtils::GraphLegend(colorLegendEntries, 0.5, 0.15, 0.8, 0.40, hvar.legendTitle.c_str());

    vector<pair<TObject*, string>> markerLegendEntries;
    for (size_t i = 0; i < labels.size(); ++i) {
        TMarker* m = new TMarker(0, 0, markerStyles[i % markerStyles.size()]);
        m->SetMarkerColor(kBlack);
        m->SetMarkerSize(1.5);
        markerLegendEntries.push_back({m, labels[i]});
    }
    PlotUtils::GraphLegend(markerLegendEntries, 0.25, 0.15, 0.45, 0.35);

    const auto& res1 = resultsVec[0];
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    double blockFilexsecx = 0.6;
    double blockFilexsecy = 0.85;
    double blockDY = 0.03;
    double xsec_had = xsec.at("Hijing");
    double xsec_SD = xsec.at("Starlight_SD");
    double xsec_DD = xsec.at("Starlight_DD");
    double xsec_alphaO = xsec.at("AlphaO");
    if (xsec_had != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy, Form("HIJING OO (%.1f b)", xsec_had));
    if (xsec_SD != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy- blockDY, Form("Starlight SD (%.1f b)", xsec_SD));
    if (xsec_DD != 0)
        latex.DrawLatex(blockFilexsecx, blockFilexsecy- 2*blockDY, Form("Starlight DD (%.5f b)", xsec_DD));
    //if (xsec_alphaO != 0)
       // latex.DrawLatex(blockFilexsecx, blockFilexsecy- 3*blockDY, Form("AlphaO (%.1f b)", xsec_alphaO));
    // Block for PVFilter etc. with adjustable y start
    double blockX = 0.20;
    double blockY = 0.85;
    latex.DrawLatex(blockX, blockY, "PVFilter == 1");
    latex.DrawLatex(blockX, blockY - blockDY, "ClusterCompatibility == 1");
    latex.DrawLatex(blockX, blockY - 2*blockDY, "abs(V_{z}) < 15");
    if (npartcut == true){
        latex.DrawLatex(blockX, blockY - 3*blockDY, "npart > 1 for HIJING");
    }

    c->SaveAs((hvar.outFolderName + "Puritycurve_AllTrackptcuts_shapes_" + andor.c_str() + hvar.outFileName +  ".png").c_str());
}

void drawROCasymm(const EfficiencyResults& res,
    HistVar1D hvar,
    string andor = "and") {
    
    const auto& xplus  = res.xplus;
    const auto& xminus  = res.xminus;
    const std::vector<double>* eff = nullptr;
    const std::vector<double>* purity = nullptr;

    if (andor == "and" || andor == "AND") {
        eff     = &res.EfficiencyHijingAND;
        purity  = &res.PurityAND;
    } else if (andor == "or" || andor == "OR") {
        eff     = &res.EfficiencyHijingOR;
        purity  = &res.PurityOR;
    } else {
        cerr << "Error: Invalid 'andor' parameter. Use 'and' or 'or'." << endl;
        return;
    }

    // ---- build colour map keyed to min(xplus,xminus) ------------------------
    set<int> uniqueMin;
    for (size_t i=0;i<xplus.size();++i)
        uniqueMin.insert(static_cast<int>(min(xplus[i], xminus[i])));

    vector<int> palette = getRainbow(uniqueMin.size());
    map<int,int> min2col;
    size_t idx = 0;
    for (int v : uniqueMin)  min2col[v] = palette[idx++];
    PlotUtils::setgstyle();
    // ---- draw points --------------------------------------------------------
    TCanvas* c  = new TCanvas("cEffPur", "Eff_Hijing(OR) vs Purity(OR)", 800, 700);
    TMultiGraph* mg = new TMultiGraph();

    for (size_t i=0;i<xplus.size();++i)
    {
        double x = (*eff)[i];
        double y = (*purity)[i];        
        int style;
        
        if (xplus[i] == xminus[i]) {
            style = 20;
        } else if (xminus[i] == xplus[i] - 1 || xminus[i] == xplus[i] + 1) {
            style = 21;
        } else if (xminus[i] == xplus[i] + 2 || xminus[i] == xplus[i] - 2) {
            style = 22;
        }
        
        int color = min2col[static_cast<int>(min(xplus[i],xminus[i])) ];

        TGraph* g = new TGraph(1,&x,&y);
        g->SetMarkerStyle(style);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetMarkerSize(1.3);

        mg->Add(g,"P");
    }

    mg->Draw("A");
    mg->SetTitle(hvar.histTitle.c_str());
    mg->GetXaxis()->SetTitle(hvar.xLabel);
    mg->GetYaxis()->SetTitle(hvar.yLabel);
    mg->GetXaxis()->SetLimits(hvar.xmin, hvar.xmax);
    mg->SetMinimum(hvar.ymin);
    mg->SetMaximum(hvar.ymax);

    c->SetGrid();

    // ---- legend & label ------------------------------------------------------
    vector<pair<TObject*, string>> legendEntries;

    // Marker for xplus > xminus
    TMarker* m1 = new TMarker(0, 0, 20);
    m1->SetMarkerColor(kBlack);
    m1->SetMarkerSize(1.3);
    legendEntries.push_back({m1, "x_{1} = x_{2}"});

    // Marker for xminus > xplus
    TMarker* m2 = new TMarker(0, 0, 21);
    m2->SetMarkerColor(kBlack);
    m2->SetMarkerSize(1.3);
    legendEntries.push_back({m2, "x_{1} = x_{2} #pm 1"});

    TMarker* m3 = new TMarker(0, 0, 22);
    m3->SetMarkerColor(kBlack);
    m3->SetMarkerSize(1.3);
    legendEntries.push_back({m3, "x_{1} = x_{2} #pm 2"});

    PlotUtils::GraphLegend(legendEntries, 0.25, 0.25, 0.45, 0.35);



    vector<pair<TObject*, string>> colorLegendEntries;
    for (const auto& kv : min2col) {
        TMarker* m = new TMarker(0, 0, 20);
        m->SetMarkerColor(kv.second);
        m->SetMarkerSize(1.3);
        colorLegendEntries.push_back({m, Form("x = %d", kv.first)});
    }
    PlotUtils::GraphLegend(colorLegendEntries, 0.65, 0.25, 0.95, 0.45);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.6, 0.85, Form("HIJING (%.1f b)", res.xsec_had));
    latex.DrawLatex(0.6, 0.82, Form("Starlight SD (%.1f b)", 0.6));
    latex.DrawLatex(0.6, 0.79, Form("Starlight DD (%.5f b)", res.xsec_DD));
    latex.DrawLatex(0.23, 0.55, "PVFilter == 1");
    latex.DrawLatex(0.23, 0.52, "ClusterCompatibilityFilter == 1");
    latex.DrawLatex(0.23, 0.49, "abs(V_{z}) < 15");
    latex.DrawLatex(0.23, 0.46, "npart > 1 for HIJING");
    // Use manual outFolderName if set, else current directory
    std::string outFolder = hvar.outFolderName.empty() ? "./" : hvar.outFolderName;
    c->SaveAs((outFolder + "Puritycurve_Asymm_" + hvar.outFileName + ".png").c_str());
}
