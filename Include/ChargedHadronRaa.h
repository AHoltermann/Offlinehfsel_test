

class ChargedHadronRAATreeMessenger
{
public:
   TTree *Tree;
   int Run;
   long long Event;
   int Lumi;
   int hiBin;
   float VX, VY, VZ, VXError, VYError, VZError; //best vertex from track tree
   bool isFakeVtx;                              //best vertex from track tree
   int nTracksVtx, bestVtxIndx;                 //best vertex from track tree
   float chi2Vtx, ndofVtx;                      //best vertex from track tree
   float ptSumVtx;
   int nVtx;
   float HFEMaxPlus;
   float HFEMaxPlus2;
   float HFEMaxPlus3;
   float HFEMaxMinus;
   float HFEMaxMinus2;
   float HFEMaxMinus3;
   float ZDCsumPlus;
   float ZDCsumMinus;
   int ClusterCompatibilityFilter;
   int PVFilter;
   int mMaxL1HFAdcPlus, mMaxL1HFAdcMinus;
   float hiHF_pf;
   float Npart;
   float Ncoll;
   float leadingPtEta1p0_sel;
   int sampleType;

   std::vector<float> *trkPt;
   std::vector<float> *trkPhi;
   std::vector<float> *trkPtError;
   std::vector<float> *trkEta;
   std::vector<bool> *highPurity;
   std::vector<float> *trkDxyAssociatedVtx;
   std::vector<float> *trkDzAssociatedVtx;
   std::vector<float> *trkDxyErrAssociatedVtx;
   std::vector<float> *trkDzErrAssociatedVtx;
   std::vector<int> *trkAssociatedVtxIndx;
   std::vector<char> *trkCharge;
   std::vector<char> *trkNHits;
   std::vector<char> *trkNPixHits;
   std::vector<char> *trkNLayers;
   std::vector<float> *trkNormChi2;
   std::vector<float> *pfEnergy;

   // weighting properties
   std::vector<float> *trackWeight;

   // Debug mode quantities
   std::vector<float> *AllxVtx;
   std::vector<float> *AllyVtx;
   std::vector<float> *AllzVtx;
   std::vector<float> *AllxVtxError;
   std::vector<float> *AllyVtxError;
   std::vector<float> *AllzVtxError;
   std::vector<bool> *AllisFakeVtx;
   std::vector<int> *AllnTracksVtx;
   std::vector<float> *Allchi2Vtx;
   std::vector<float> *AllndofVtx;
   std::vector<float> *AllptSumVtx;

public:   // Derived quantities
   //bool GoodPhotonuclear; //FIXME: currently not implemented

private:
   bool WriteMode;
   bool Initialized;
   bool DebugMode;

public:
   ChargedHadronRAATreeMessenger(TFile &File, std::string TreeName = "tree", bool Debug = false);
   ChargedHadronRAATreeMessenger(TFile *File, std::string TreeName = "tree", bool Debug = false);
   ChargedHadronRAATreeMessenger(TTree *ChargedHadRAATree = nullptr, bool Debug = false);
   ~ChargedHadronRAATreeMessenger();
   bool Initialize(TTree *ChargedHadRAATree, bool Debug = false);
   bool Initialize(bool Debug = false);
   int GetEntries();
   bool GetEntry(int iEntry);
   bool SetBranch(TTree *T, bool Debug = false);
   void Clear();
   //void CopyNonTrack(ChargedHadronRAATreeMessenger &M);
   bool FillEntry();

};
