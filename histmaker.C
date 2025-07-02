/// begin aslkdfj;aslkjf;lkj;lkjelja;shgasufandskm
using namespace std;

void Centrality(TTree* T, TCut t, TH1D* h){
     T->Project(h->GetName(), "hiBin", t);
}

void HFEMaxPlus(TTree * T, TCut t, TH1D* h){
     T->Project(h->GetName(), "HFEMaxPlus", t);
}

void HFEMaxMinus(TTree * T, TCut t, TH1D* h){
     T->Project(h->GetName(), "HFEMaxMinus", t);
}

void HFEMax_Maximum(TTree* T, TCut t, TH1D* h) {
    TString expr = Form("TMath::Max(%s,%s)", "HFEMaxPlus", "HFEMaxMinus");
    T->Project(h->GetName(), expr, t);
}

void HFEMax_Minimum(TTree* T, TCut t, TH1D* h) {
    TString expr = Form("TMath::Min(%s,%s)", "HFEMaxPlus", "HFEMaxMinus");
    T->Project(h->GetName(), expr, t);
}

void ADCPlus(TTree* T, TCut t, TH1D* h) {
    T->Project(h->GetName(), "mMaxL1HFAdcPlus", t);
}

void ADCMinus(TTree* T, TCut t, TH1D* h) {
    T->Project(h->GetName(), "mMaxL1HFAdcMinus", t);
}  

void ADC_Maximum(TTree* T, TCut t, TH1D* h) {
    TString expr = Form("TMath::Max(%s,%s)", "mMaxL1HFAdcPlus", "mMaxL1HFAdcMinus");
    T->Project(h->GetName(), expr, t);
}   

void ADC_Minimum(TTree* T, TCut t, TH1D* h) {
    TString expr = Form("TMath::Min(%s,%s)", "mMaxL1HFAdcPlus", "mMaxL1HFAdcMinus");
    T->Project(h->GetName(), expr, t);
}

void HFEMaxMinScatter(TTree* T, TCut t, TH2D* h) {
    TString exprX = "HFEMaxPlus";
    TString exprY = "HFEMaxMinus";
    T->Project(h->GetName(), exprX + ":" + exprY, t);
}

void HFEPlusMinusScatter(TTree* T, TCut t, TH2D* h) {
    TString exprX = "HFEMaxPlus";
    TString exprY = "HFEMaxMinus";
    T->Project(h->GetName(), exprX + ":" + exprY, t);
}

Long64_t HowManyPass(TTree* T, TCut t){
    Long64_t nPass = T->GetEntries(t);
    return nPass;
}

TCut HFONLINE(int threshold, string andor){
    
    string expr = 
        "mMaxL1HFAdcPlus > " +
        std::to_string(threshold) + 
        " " +
        andor +  
        " " +
        "mMaxL1HFAdcMinus > " + 
        std::to_string(threshold);

    const char* cutExpr = expr.c_str();
    TCut ezcut(cutExpr);

    return ezcut;
}

TCut HFOFFLINE(int threshold, string andor){

    string expr = 
        "HFEMaxPlus > " + 
        std::to_string(threshold) + 
        " " + 
        andor + 
        " " +  
        "HFEMaxMinus > " + 
        std::to_string(threshold);

    const char* cutExpr = expr.c_str();
    TCut ezcut(cutExpr);

    return ezcut;
}

TCut TRKPTCUT(float threshold){
    string expr = "leadingPtEta1p0_sel > " + std::to_string(threshold);
    const char* cutExpr = expr.c_str();
    TCut ezcut(cutExpr);
    return ezcut;
}

TH1D* Divide(TH1D*h1, TH1D*h2){
    TH1D* h3 = (TH1D*)h1->Clone();
    h3->Divide(h2);
    return h3;
}

void histmaker(){

    /////////////////////////////
    /////////////////////////////
    ////                     ////
    ////  INPUT FILES HERE   ////
    ////                     ////
    /////////////////////////////
    /////////////////////////////

    TFile* f = new TFile("INSERT HERE");
    TTree* T = (TTree*)f->Get("Tree"); 

    /////////////////////////////
    /////////////////////////////
    ////                     ////
    ////  Event Selection    ////
    ////                     ////
    /////////////////////////////
    /////////////////////////////

    TCut VZcut("VZ > -15 && VZ < 15");
    TCut PVcut("PVFilter == 1");
    TCut CCcut("ClusterCompatibilityFilter == 1");
    TCut Nvtxcut("nVtx > 0");
    // use HFONLINE and HFOFFLINE functions for a tcut
    // use TRKPT for track pt cuts

    // EXAMPLE EVENT SELECTIONs
    TCut ESel_Example1 = VZcut && PVcut && CCcut && Nvtxcut && HFONLINE(14, "&&");
    TCut ESel_Example2 = VZcut && PVcut && CCcut && Nvtxcut && HFOFFLINE(10, "||") && HFONLINE(14, "||");

    /////////////////////////////
    /////////////////////////////
    ////                     ////
    ////  Generate Hists     ////
    ////                     ////
    /////////////////////////////
    /////////////////////////////

    TH1D* example1 = new TH1D("example", "Example Histogram", 201, -0.5, 200.5);
    TH1D* example2 = new TH1D("example2", "Example Histogram 2", 201, -0.5, 200.5);

    Centrality(T, ESel_Example1, example1);
    Centrality(T, ESel_Example2, example2);

    TH1D* cratio = Divide(example1, example2);
    

    /////////////////////////////
    /////////////////////////////
    ////                     ////
    ////  Save to output     ////
    ////                     ////
    /////////////////////////////
    /////////////////////////////

    TFile* outfile = new TFile("output.root", "RECREATE");
    example1->Write();
    example2->Write();
    cratio->Write();
    outfile->Close();
    

}

