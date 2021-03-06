// Class: ReadBDT_50_50_03_20
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : BDT::BDT_50_50_03_20
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.06/01       [394753]
Creator        : matthies
Date           : Tue Sep 18 22:10:49 2018
Host           : Linux vocms0301.cern.ch 2.6.32-642.4.2.el6.x86_64 #1 SMP Wed Aug 24 09:19:54 CEST 2016 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb
Training events: 274422
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
H: "False" [Print method-specific help message]
NTrees: "50" [Number of trees in the forest]
MaxDepth: "3" [Max depth of the decision tree allowed]
MinNodeSize: "5%" [Minimum percentage of training events required in a leaf node (default: Classification: 5%, Regression: 0.2%)]
nCuts: "20" [Number of grid points in variable range used in finding optimal cut in node splitting]
BoostType: "AdaBoost" [Boosting type for the trees in the forest (note: AdaCost is still experimental)]
UseBaggedBoost: "True" [Use only a random subsample of all events for growing the trees in each boost iteration.]
AdaBoostBeta: "5.000000e-01" [Learning rate  for AdaBoost algorithm]
BaggedSampleFraction: "3.000000e-01" [Relative size of bagged event sample to original size of the data sample (used whenever bagging is used (i.e. UseBaggedBoost, Bagging,)]
SeparationType: "giniindex" [Separation criterion for node splitting]
# Default:
VerbosityLevel: "Default" [Verbosity level]
VarTransform: "None" [List of variable transformations performed before training, e.g., "D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)"]
CreateMVAPdfs: "False" [Create PDFs for classifier outputs (signal and background)]
IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are included for testing and performance evaluation)]
AdaBoostR2Loss: "quadratic" [Type of Loss function in AdaBoostR2]
Shrinkage: "1.000000e+00" [Learning rate for GradBoost algorithm]
UseRandomisedTrees: "False" [Determine at each node splitting the cut variable only as the best out of a random subset of variables (like in RandomForests)]
UseNvars: "3" [Size of the subset of variables used with RandomisedTree option]
UsePoissonNvars: "True" [Interpret "UseNvars" not as fixed number but as mean of a Possion distribution in each split with RandomisedTree option]
UseYesNoLeaf: "True" [Use Sig or Bkg categories, or the purity=S/(S+B) as classification of the leaf node -> Real-AdaBoost]
NegWeightTreatment: "inverseboostnegweights" [How to treat events with negative weights in the BDT training (particular the boosting) : IgnoreInTraining;  Boost With inverse boostweight; Pair events with negative and positive weights in traning sample and *annihilate* them (experimental!)]
Css: "1.000000e+00" [AdaCost: cost of true signal selected signal]
Cts_sb: "1.000000e+00" [AdaCost: cost of true signal selected bkg]
Ctb_ss: "1.000000e+00" [AdaCost: cost of true bkg    selected signal]
Cbb: "1.000000e+00" [AdaCost: cost of true bkg    selected bkg ]
NodePurityLimit: "5.000000e-01" [In boosting/pruning, nodes with purity > NodePurityLimit are signal; background otherwise.]
DoBoostMonitor: "False" [Create control plot with ROC integral vs tree number]
UseFisherCuts: "False" [Use multivariate splits using the Fisher criterion]
MinLinCorrForFisher: "8.000000e-01" [The minimum linear correlation between two variables demanded for use in Fisher criterion in node splitting]
UseExclusiveVars: "False" [Variables already used in fisher criterion are not anymore analysed individually for node splitting]
DoPreselection: "False" [and and apply automatic pre-selection for 100% efficient signal (bkg) cuts prior to training]
SigToBkgFraction: "1.000000e+00" [Sig to Bkg ratio used in Training (similar to NodePurityLimit, which cannot be used in real adaboost]
PruneMethod: "nopruning" [Note: for BDTs use small trees (e.g.MaxDepth=3) and NoPruning:  Pruning: Method used for pruning (removal) of statistically insignificant branches ]
PruneStrength: "0.000000e+00" [Pruning strength]
PruningValFraction: "5.000000e-01" [Fraction of events to use for optimizing automatic pruning.]
nEventsMin: "0" [deprecated: Use MinNodeSize (in % of training events) instead]
UseBaggedGrad: "False" [deprecated: Use *UseBaggedBoost* instead:  Use only a random subsample of all events for growing the trees in each iteration.]
GradBaggingFraction: "3.000000e-01" [deprecated: Use *BaggedSampleFraction* instead: Defines the fraction of events to be used in each iteration, e.g. when UseBaggedGrad=kTRUE. ]
UseNTrainEvents: "0" [deprecated: Use *BaggedSampleFraction* instead: Number of randomly picked training events used in randomised (and bagged) trees]
NNodesMax: "0" [deprecated: Use MaxDepth instead to limit the tree size]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 8
tmva_deltaphi_bottomlepton    tmva_deltaphi_bottomlepton    tmva_deltaphi_bottomlepton    #Delta#phi(l,b)               rad                               'F'    [1.04904174805e-05,3.14159083366]
tmva_deltaphi_bottomtoptag    tmva_deltaphi_bottomtoptag    tmva_deltaphi_bottomtoptag    #Delta#phi(t,b)               rad                               'F'    [1.51991844177e-06,3.14156723022]
tmva_pseudotop_mass           tmva_pseudotop_mass           tmva_pseudotop_mass           l#nub mass                    GeV                               'F'    [86.8835525513,6687.15380859]
tmva_ptbalance                tmva_ptbalance                tmva_ptbalance                pT balance                                                      'F'    [-0.99961078167,4.64822101593]
tmva_wass_pt                  tmva_wass_pt                  tmva_wass_pt                  asso. W pT                    GeV                               'F'    [0.213053271174,1770.44580078]
tmva_lepton_pt                tmva_lepton_pt                tmva_lepton_pt                lepton pT                     GeV                               'F'    [50.0001945496,1194.01318359]
tmva_lepton_eta               tmva_lepton_eta               tmva_lepton_eta               lepton #eta                                                     'F'    [-2.39999055862,2.39996170998]
tmva_deltaphi_leptonnextjet   tmva_deltaphi_leptonnextjet   tmva_deltaphi_leptonnextjet   #Delta#phi(l,jet)             rad                               'F'    [6.31809234619e-06,3.14115262032]
NSpec 0


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#define NN new BDT_50_50_03_20Node
   
#ifndef BDT_50_50_03_20Node__def
#define BDT_50_50_03_20Node__def
   
class BDT_50_50_03_20Node {
   
public:
   
   // constructor of an essentially "empty" node floating in space
   BDT_50_50_03_20Node ( BDT_50_50_03_20Node* left,BDT_50_50_03_20Node* right,
                          int selector, double cutValue, bool cutType, 
                          int nodeType, double purity, double response ) :
   fLeft         ( left         ),
   fRight        ( right        ),
   fSelector     ( selector     ),
   fCutValue     ( cutValue     ),
   fCutType      ( cutType      ),
   fNodeType     ( nodeType     ),
   fPurity       ( purity       ),
   fResponse     ( response     ){
   }

   virtual ~BDT_50_50_03_20Node();

   // test event if it decends the tree at this node to the right
   virtual bool GoesRight( const std::vector<double>& inputValues ) const;
   BDT_50_50_03_20Node* GetRight( void )  {return fRight; };

   // test event if it decends the tree at this node to the left 
   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;
   BDT_50_50_03_20Node* GetLeft( void ) { return fLeft; };   

   // return  S/(S+B) (purity) at this node (from  training)

   double GetPurity( void ) const { return fPurity; } 
   // return the node type
   int    GetNodeType( void ) const { return fNodeType; }
   double GetResponse(void) const {return fResponse;}

private:

   BDT_50_50_03_20Node*   fLeft;     // pointer to the left daughter node
   BDT_50_50_03_20Node*   fRight;    // pointer to the right daughter node
   int                     fSelector; // index of variable used in node selection (decision tree)   
   double                  fCutValue; // cut value appplied on this node to discriminate bkg against sig
   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise
   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal 
   double                  fPurity;   // Purity of node from training
   double                  fResponse; // Regression response value of node
}; 
   
//_______________________________________________________________________
   BDT_50_50_03_20Node::~BDT_50_50_03_20Node()
{
   if (fLeft  != NULL) delete fLeft;
   if (fRight != NULL) delete fRight;
}; 
   
//_______________________________________________________________________
bool BDT_50_50_03_20Node::GoesRight( const std::vector<double>& inputValues ) const
{
   // test event if it decends the tree at this node to the right
   bool result;
     result = (inputValues[fSelector] > fCutValue );
   if (fCutType == true) return result; //the cuts are selecting Signal ;
   else return !result;
}
   
//_______________________________________________________________________
bool BDT_50_50_03_20Node::GoesLeft( const std::vector<double>& inputValues ) const
{
   // test event if it decends the tree at this node to the left
   if (!this->GoesRight(inputValues)) return true;
   else return false;
}
   
#endif
   
#ifndef IClassifierReader__def
#define IClassifierReader__def

class IClassifierReader {

 public:

   // constructor
   IClassifierReader() : fStatusIsClean( true ) {}
   virtual ~IClassifierReader() {}

   // return classifier response
   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;

   // returns classifier status
   bool IsStatusClean() const { return fStatusIsClean; }

 protected:

   bool fStatusIsClean;
};

#endif

class ReadBDT_50_50_03_20 : public IClassifierReader {

 public:

   // constructor
   ReadBDT_50_50_03_20( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadBDT_50_50_03_20" ),
        fNvars( 8 ),
        fIsNormalised( false )
   {      
      // the training input variables
      const char* inputVars[] = { "tmva_deltaphi_bottomlepton", "tmva_deltaphi_bottomtoptag", "tmva_pseudotop_mass", "tmva_ptbalance", "tmva_wass_pt", "tmva_lepton_pt", "tmva_lepton_eta", "tmva_deltaphi_leptonnextjet" };

      // sanity checks
      if (theInputVars.size() <= 0) {
         std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
         fStatusIsClean = false;
      }

      if (theInputVars.size() != fNvars) {
         std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
                   << theInputVars.size() << " != " << fNvars << std::endl;
         fStatusIsClean = false;
      }

      // validate input variables
      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
         if (theInputVars[ivar] != inputVars[ivar]) {
            std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                      << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
            fStatusIsClean = false;
         }
      }

      // initialize min and max vectors (for normalisation)
      fVmin[0] = 0;
      fVmax[0] = 0;
      fVmin[1] = 0;
      fVmax[1] = 0;
      fVmin[2] = 0;
      fVmax[2] = 0;
      fVmin[3] = 0;
      fVmax[3] = 0;
      fVmin[4] = 0;
      fVmax[4] = 0;
      fVmin[5] = 0;
      fVmax[5] = 0;
      fVmin[6] = 0;
      fVmax[6] = 0;
      fVmin[7] = 0;
      fVmax[7] = 0;

      // initialize input variable types
      fType[0] = 'F';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';
      fType[5] = 'F';
      fType[6] = 'F';
      fType[7] = 'F';

      // initialize constants
      Initialize();

   }

   // destructor
   virtual ~ReadBDT_50_50_03_20() {
      Clear(); // method-specific
   }

   // the classifier response
   // "inputValues" is a vector of input values in the same order as the 
   // variables given to the constructor
   double GetMvaValue( const std::vector<double>& inputValues ) const;

 private:

   // method-specific destructor
   void Clear();

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   const bool fIsNormalised;
   bool IsNormalised() const { return fIsNormalised; }
   double fVmin[8];
   double fVmax[8];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[8];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)
   std::vector<BDT_50_50_03_20Node*> fForest;       // i.e. root nodes of decision trees
   std::vector<double>                fBoostWeights; // the weights applied in the individual boosts
};

double ReadBDT_50_50_03_20::GetMvaValue__( const std::vector<double>& inputValues ) const
{
   double myMVA = 0;
   double norm  = 0;
   for (unsigned int itree=0; itree<fForest.size(); itree++){
      BDT_50_50_03_20Node *current = fForest[itree];
      while (current->GetNodeType() == 0) { //intermediate node
         if (current->GoesRight(inputValues)) current=(BDT_50_50_03_20Node*)current->GetRight();
         else current=(BDT_50_50_03_20Node*)current->GetLeft();
      }
      myMVA += fBoostWeights[itree] *  current->GetNodeType();
      norm  += fBoostWeights[itree];
   }
   return myMVA /= norm;
};

void ReadBDT_50_50_03_20::Initialize()
{
  // itree = 0
  fBoostWeights.push_back(0.371248088332866);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 149.522, 0, 1, 0.694436,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.56062,-99) , 
NN(
0, 
0, 
-1, -0.149781, 0, -1, 0.32166,-99) , 
5, 189.953, 0, 0, 0.347627,-99) , 
7, 1.49549, 0, 0, 0.491612,-99)    );
  // itree = 1
  fBoostWeights.push_back(0.321704);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 167.667, 0, 1, 0.674201,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.547707,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.371616,-99) , 
4, 253.103, 0, 0, 0.428063,-99) , 
NN(
0, 
0, 
-1, 1.0468, 0, -1, 0.301235,-99) , 
2, 234.563, 0, 0, 0.375717,-99) , 
7, 1.49572, 0, 0, 0.499743,-99)    );
  // itree = 2
  fBoostWeights.push_back(0.292667);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 158.956, 0, 1, 0.672495,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.610677,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477693,-99) , 
2, 382.626, 0, 0, 0.526085,-99) , 
NN(
0, 
0, 
-1, 228.709, 0, -1, 0.35215,-99) , 
3, -0.12776, 0, 0, 0.396474,-99) , 
7, 1.64537, 0, 0, 0.500411,-99)    );
  // itree = 3
  fBoostWeights.push_back(0.248728);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.772481,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.612573,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.472656,-99) , 
3, -0.139165, 0, 0, 0.560622,-99) , 
5, 158.958, 0, 0, 0.625519,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509076,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.410046,-99) , 
3, -0.240316, 0, 0, 0.452862,-99) , 
NN(
0, 
0, 
-1, 0.89744, 0, -1, 0.312727,-99) , 
2, 231.1, 0, 0, 0.398083,-99) , 
7, 1.34589, 0, 0, 0.502889,-99)    );
  // itree = 4
  fBoostWeights.push_back(0.241256);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 124.939, 0, 1, 0.6418,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.599692,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.432309,-99) , 
5, 169.626, 0, 0, 0.464465,-99) , 
NN(
0, 
0, 
-1, 2.39396, 1, -1, 0.346426,-99) , 
2, 228.095, 0, 0, 0.419821,-99) , 
7, 1.64483, 0, 0, 0.504053,-99)    );
  // itree = 5
  fBoostWeights.push_back(0.198481);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 158.962, 0, 1, 0.60013,-99) , 
NN(
0, 
0, 
-1, 232.059, 0, -1, 0.404291,-99) , 
7, 1.19648, 0, 0, 0.503898,-99)    );
  // itree = 6
  fBoostWeights.push_back(0.219432);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 158.956, 0, 1, 0.644779,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.604573,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.445632,-99) , 
5, 108.1, 0, 0, 0.517562,-99) , 
NN(
0, 
0, 
-1, 1.0472, 0, -1, 0.394682,-99) , 
7, 1.04678, 0, 0, 0.450672,-99) , 
3, -0.0454891, 0, 0, 0.508296,-99)    );
  // itree = 7
  fBoostWeights.push_back(0.191383);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.72134,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.571605,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481625,-99) , 
5, 91.5166, 0, 0, 0.532604,-99) , 
5, 158.965, 0, 0, 0.589508,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.578745,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.450909,-99) , 
1, 0.747953, 0, 0, 0.480937,-99) , 
NN(
0, 
0, 
-1, 1.49612, 1, -1, 0.373417,-99) , 
2, 259.75, 0, 0, 0.430217,-99) , 
7, 1.34599, 0, 0, 0.503665,-99)    );
  // itree = 8
  fBoostWeights.push_back(0.169426);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 208.257, 0, 1, 0.686559,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.561622,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.478457,-99) , 
7, 1.34588, 0, 0, 0.539266,-99) , 
5, 158.96, 0, 0, 0.596441,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.541911,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.451494,-99) , 
5, 110.414, 0, 0, 0.49017,-99) , 
NN(
0, 
0, 
-1, 106.736, 1, -1, 0.372293,-99) , 
2, 219.82, 0, 0, 0.453022,-99) , 
3, -0.102194, 0, 0, 0.508509,-99)    );
  // itree = 9
  fBoostWeights.push_back(0.17272);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.672776,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.550658,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.411741,-99) , 
2, 209.373, 0, 0, 0.526887,-99) , 
5, 158.962, 0, 0, 0.566415,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.514616,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.424423,-99) , 
3, -0.235585, 0, 0, 0.462061,-99) , 
NN(
0, 
0, 
-1, 0.647934, 0, -1, 0.361022,-99) , 
0, 1.0472, 0, 0, 0.426193,-99) , 
7, 1.04689, 0, 0, 0.505758,-99)    );
  // itree = 10
  fBoostWeights.push_back(0.101632);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 158.965, 0, 1, 0.601264,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.578459,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.488786,-99) , 
5, 179.522, 0, 0, 0.498592,-99) , 
NN(
0, 
0, 
-1, -0.352079, 0, -1, 0.397741,-99) , 
0, 0.897605, 0, 0, 0.481312,-99) , 
3, -0.0308692, 0, 0, 0.512433,-99)    );
  // itree = 11
  fBoostWeights.push_back(0.135422);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.698044,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.545465,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.478896,-99) , 
5, 78.3063, 0, 0, 0.525077,-99) , 
5, 168.875, 0, 0, 0.56544,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.543003,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.459028,-99) , 
1, 0.897497, 0, 0, 0.47696,-99) , 
NN(
0, 
0, 
-1, 2.24423, 1, -1, 0.42153,-99) , 
2, 237.404, 0, 0, 0.453879,-99) , 
7, 1.04681, 0, 0, 0.516923,-99)    );
  // itree = 12
  fBoostWeights.push_back(0.102006);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.662004,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.539562,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.482728,-99) , 
3, -0.132611, 0, 0, 0.507895,-99) , 
NN(
0, 
0, 
-1, 0.116349, 0, -1, 0.390675,-99) , 
0, 0.897621, 0, 0, 0.491304,-99) , 
5, 213.432, 0, 0, 0.507999,-99)    );
  // itree = 13
  fBoostWeights.push_back(0.109671);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.68062,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477718,-99) , 
7, 1.19664, 0, 0, 0.604532,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.612438,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.495551,-99) , 
1, 1.19671, 0, 0, 0.509654,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505332,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.39204,-99) , 
1, 1.3471, 1, 0, 0.435222,-99) , 
2, 235.97, 0, 0, 0.485161,-99) , 
5, 158.955, 0, 0, 0.512971,-99)    );
  // itree = 14
  fBoostWeights.push_back(0.0933501);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.65683,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.550593,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.485931,-99) , 
6, 0.799941, 1, 0, 0.534738,-99) , 
NN(
0, 
0, 
-1, 207.714, 0, -1, 0.469144,-99) , 
7, 1.79368, 0, 0, 0.490382,-99) , 
5, 212.027, 0, 0, 0.506952,-99)    );
  // itree = 15
  fBoostWeights.push_back(0.102458);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 257.533, 0, 1, 0.600651,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.518695,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.421981,-99) , 
6, -1.25676, 0, 0, 0.507478,-99) , 
NN(
0, 
0, 
-1, 197.731, 0, -1, 0.416674,-99) , 
0, 1.04721, 0, 0, 0.491834,-99) , 
5, 158.955, 0, 0, 0.517774,-99)    );
  // itree = 16
  fBoostWeights.push_back(0.102611);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.675036,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.545895,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.485415,-99) , 
7, 0.897361, 0, 0, 0.525488,-99) , 
NN(
0, 
0, 
-1, 111.793, 1, -1, 0.430717,-99) , 
2, 203.817, 0, 0, 0.505485,-99) , 
5, 235.502, 0, 0, 0.516525,-99)    );
  // itree = 17
  fBoostWeights.push_back(0.0845463);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.672865,-99) , 
NN(
NN(
0, 
0, 
-1, 136.89, 0, 1, 0.530097,-99) , 
NN(
0, 
0, 
-1, -0.406529, 0, -1, 0.464085,-99) , 
7, 0.897477, 0, 0, 0.505179,-99) , 
5, 232.483, 0, 0, 0.516964,-99)    );
  // itree = 18
  fBoostWeights.push_back(0.0796812);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.638858,-99) , 
NN(
NN(
0, 
0, 
-1, 1.04717, 0, 1, 0.51998,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.512486,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.429196,-99) , 
1, 1.34651, 1, 0, 0.465605,-99) , 
2, 256.948, 0, 0, 0.499566,-99) , 
5, 213.437, 0, 0, 0.512642,-99)    );
  // itree = 19
  fBoostWeights.push_back(0.0711433);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.687565,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.605804,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.499426,-99) , 
1, 1.19672, 0, 0, 0.511688,-99) , 
NN(
0, 
0, 
-1, 1.49635, 1, -1, 0.438964,-99) , 
2, 232.622, 0, 0, 0.490661,-99) , 
5, 249.043, 0, 0, 0.501472,-99)    );
  // itree = 20
  fBoostWeights.push_back(0.0856462);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0.0688372, 0, 1, 0.539171,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.504172,-99) , 
NN(
0, 
0, 
-1, 153.092, 1, -1, 0.420961,-99) , 
1, 1.34652, 1, 0, 0.45209,-99) , 
2, 240.173, 0, 0, 0.513123,-99)    );
  // itree = 21
  fBoostWeights.push_back(0.0717815);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.65706,-99) , 
NN(
NN(
0, 
0, 
-1, 1.04713, 0, 1, 0.517258,-99) , 
NN(
0, 
0, 
-1, 1.19694, 1, -1, 0.455826,-99) , 
2, 258.289, 0, 0, 0.495093,-99) , 
5, 235.501, 0, 0, 0.506123,-99)    );
  // itree = 22
  fBoostWeights.push_back(0.0710166);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.652227,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.561606,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481806,-99) , 
7, 1.04705, 0, 0, 0.53316,-99) , 
NN(
0, 
0, 
-1, -0.800018, 0, -1, 0.478318,-99) , 
5, 135.301, 0, 0, 0.494689,-99) , 
5, 249.041, 0, 0, 0.502957,-99)    );
  // itree = 23
  fBoostWeights.push_back(0.0810087);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.651705,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.528388,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.471391,-99) , 
6, 0.767419, 1, 0, 0.512385,-99) , 
NN(
0, 
0, 
-1, 311.678, 0, -1, 0.446643,-99) , 
6, -1.02835, 0, 0, 0.501807,-99) , 
5, 235.502, 0, 0, 0.511811,-99)    );
  // itree = 24
  fBoostWeights.push_back(0.0594599);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.610503,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.58128,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.493807,-99) , 
1, 1.49597, 0, 0, 0.505403,-99) , 
NN(
0, 
0, 
-1, -0.573412, 0, -1, 0.459323,-99) , 
2, 206.67, 0, 0, 0.495201,-99) , 
5, 213.432, 0, 0, 0.505725,-99)    );
  // itree = 25
  fBoostWeights.push_back(0.0625166);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.619002,-99) , 
NN(
NN(
0, 
0, 
-1, 1.19671, 0, 1, 0.515204,-99) , 
NN(
0, 
0, 
-1, 1.79548, 1, -1, 0.459323,-99) , 
2, 224.077, 0, 0, 0.499973,-99) , 
5, 213.433, 0, 0, 0.51112,-99)    );
  // itree = 26
  fBoostWeights.push_back(0.0663018);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.620437,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.515278,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.458841,-99) , 
2, 198.405, 0, 0, 0.504117,-99) , 
NN(
0, 
0, 
-1, 0.747979, 0, -1, 0.436257,-99) , 
6, -1.0279, 0, 0, 0.493256,-99) , 
5, 236.603, 0, 0, 0.501135,-99)    );
  // itree = 27
  fBoostWeights.push_back(0.0729553);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.647534,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.554351,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.493313,-99) , 
7, 1.94418, 0, 0, 0.526489,-99) , 
4, 361.265, 0, 0, 0.552683,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.57976,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.491868,-99) , 
0, 2.54317, 1, 0, 0.520596,-99) , 
NN(
0, 
0, 
-1, -0.351875, 0, -1, 0.472318,-99) , 
2, 366.625, 0, 0, 0.487256,-99) , 
3, -0.0307666, 0, 0, 0.504502,-99)    );
  // itree = 28
  fBoostWeights.push_back(0.0588141);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.577231,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.51597,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.432358,-99) , 
6, -1.48584, 0, 0, 0.5098,-99) , 
NN(
0, 
0, 
-1, 0.427334, 1, -1, 0.448719,-99) , 
0, 0.748005, 0, 0, 0.502681,-99) , 
3, 0.0833908, 0, 0, 0.510887,-99)    );
  // itree = 29
  fBoostWeights.push_back(0.052304);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.626952,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.548545,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.498086,-99) , 
3, -0.107036, 0, 0, 0.526326,-99) , 
NN(
0, 
0, 
-1, 91.8902, 1, -1, 0.485006,-99) , 
5, 129.972, 0, 0, 0.498011,-99) , 
5, 236.604, 0, 0, 0.506423,-99)    );
  // itree = 30
  fBoostWeights.push_back(0.0511306);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.573147,-99) , 
NN(
NN(
0, 
0, 
-1, 0.447229, 0, 1, 0.524639,-99) , 
NN(
0, 
0, 
-1, 142.31, 0, -1, 0.482971,-99) , 
7, 1.94421, 0, 0, 0.493895,-99) , 
3, 0.0777092, 0, 0, 0.503441,-99)    );
  // itree = 31
  fBoostWeights.push_back(0.059816);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.597286,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.527518,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.450246,-99) , 
2, 200.003, 0, 0, 0.516432,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.504113,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.47398,-99) , 
1, 0.897545, 0, 0, 0.485802,-99) , 
5, 112.259, 0, 0, 0.498554,-99) , 
5, 213.432, 0, 0, 0.507741,-99)    );
  // itree = 32
  fBoostWeights.push_back(0.0457457);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.598937,-99) , 
NN(
NN(
0, 
0, 
-1, 147.364, 0, 1, 0.528601,-99) , 
NN(
0, 
0, 
-1, -1.02845, 0, -1, 0.488267,-99) , 
7, 1.79453, 0, 0, 0.501551,-99) , 
4, 385.733, 0, 0, 0.508192,-99)    );
  // itree = 33
  fBoostWeights.push_back(0.0471643);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.589934,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509175,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.473163,-99) , 
6, 0.919924, 1, 0, 0.501472,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.437879,-99) , 
6, -1.4852, 0, 0, 0.49664,-99) , 
5, 213.432, 0, 0, 0.505398,-99)    );
  // itree = 34
  fBoostWeights.push_back(0.064977);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.617824,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.568527,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.493999,-99) , 
5, 136.349, 0, 0, 0.538264,-99) , 
0, 2.5432, 1, 0, 0.560392,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.518126,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470659,-99) , 
2, 226.871, 0, 0, 0.504532,-99) , 
NN(
0, 
0, 
-1, 2.24396, 1, -1, 0.46788,-99) , 
4, 165.66, 0, 0, 0.495107,-99) , 
2, 453.596, 0, 0, 0.508776,-99)    );
  // itree = 35
  fBoostWeights.push_back(0.038601);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.606586,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.53559,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.498897,-99) , 
1, 1.19678, 0, 0, 0.508685,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.457834,-99) , 
4, 124.448, 0, 0, 0.504054,-99) , 
5, 235.502, 0, 0, 0.510738,-99)    );
  // itree = 36
  fBoostWeights.push_back(0.0449813);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.593703,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.581386,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.499414,-99) , 
1, 1.19679, 0, 0, 0.511514,-99) , 
NN(
0, 
0, 
-1, 169.295, 0, -1, 0.459904,-99) , 
2, 200.864, 0, 0, 0.501507,-99) , 
5, 236.604, 0, 0, 0.507117,-99)    );
  // itree = 37
  fBoostWeights.push_back(0.0562803);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.621323,-99) , 
NN(
NN(
0, 
0, 
-1, 1.04713, 0, 1, 0.51472,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.5688,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.463494,-99) , 
2, 151.417, 1, 0, 0.483713,-99) , 
2, 234.648, 0, 0, 0.50552,-99) , 
5, 249.042, 0, 0, 0.511512,-99)    );
  // itree = 38
  fBoostWeights.push_back(0.0450262);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.571407,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.511936,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476237,-99) , 
2, 205.017, 0, 0, 0.505427,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505067,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.458229,-99) , 
5, 55.6169, 1, 0, 0.470696,-99) , 
5, 73.5919, 0, 0, 0.496881,-99) , 
4, 374.603, 0, 0, 0.502514,-99)    );
  // itree = 39
  fBoostWeights.push_back(0.0525509);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.578033,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.549705,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.530197,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481311,-99) , 
5, 185.896, 0, 0, 0.48668,-99) , 
7, 0.149575, 1, 0, 0.490357,-99) , 
4, 372.632, 0, 0, 0.497457,-99)    );
  // itree = 40
  fBoostWeights.push_back(0.0351056);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.560003,-99) , 
NN(
NN(
0, 
0, 
-1, 0.149602, 1, 1, 0.509621,-99) , 
NN(
0, 
0, 
-1, 1.63787, 1, -1, 0.463551,-99) , 
6, 1.25694, 1, 0, 0.504221,-99) , 
3, 0.0935272, 0, 0, 0.509866,-99)    );
  // itree = 41
  fBoostWeights.push_back(0.0405649);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.566278,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509224,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.477669,-99) , 
6, -0.799874, 0, 0, 0.502203,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.437283,-99) , 
0, 0.448841, 0, 0, 0.498283,-99) , 
3, 0.0940732, 0, 0, 0.50504,-99)    );
  // itree = 42
  fBoostWeights.push_back(0.0559048);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.607323,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.542726,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.493862,-99) , 
7, 0.897477, 0, 0, 0.527336,-99) , 
0, 2.39399, 1, 0, 0.542558,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.51413,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.472604,-99) , 
2, 197.809, 0, 0, 0.503954,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.436267,-99) , 
4, 118.229, 0, 0, 0.497976,-99) , 
2, 413.545, 0, 0, 0.509849,-99)    );
  // itree = 43
  fBoostWeights.push_back(0.0591869);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.566553,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.558055,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.452036,-99) , 
7, 0.897058, 0, 0, 0.52243,-99) , 
NN(
0, 
0, 
-1, 112.372, 0, -1, 0.482345,-99) , 
2, 455.097, 0, 0, 0.489844,-99) , 
3, 0.0777092, 0, 0, 0.498946,-99)    );
  // itree = 44
  fBoostWeights.push_back(0.0411346);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.551483,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.50945,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.469468,-99) , 
6, -1.25715, 0, 0, 0.504119,-99) , 
NN(
0, 
0, 
-1, 2.69277, 1, -1, 0.462561,-99) , 
6, 1.02855, 1, 0, 0.497096,-99) , 
3, 0.0802795, 0, 0, 0.50338,-99)    );
  // itree = 45
  fBoostWeights.push_back(0.0390548);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.569475,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.519434,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.490308,-99) , 
5, 83.9753, 0, 0, 0.507203,-99) , 
NN(
0, 
0, 
-1, 0.598262, 0, -1, 0.481357,-99) , 
4, 247.379, 1, 0, 0.498473,-99) , 
4, 399.543, 0, 0, 0.502571,-99)    );
  // itree = 46
  fBoostWeights.push_back(0.0471613);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.580581,-99) , 
NN(
NN(
0, 
0, 
-1, 1.04707, 1, 1, 0.530968,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.507638,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.4718,-99) , 
1, 0.719505, 0, 0, 0.489948,-99) , 
1, 0.149603, 1, 0, 0.495898,-99) , 
5, 236.603, 0, 0, 0.501273,-99)    );
  // itree = 47
  fBoostWeights.push_back(0.0440597);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.541472,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.53692,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.488151,-99) , 
5, 129.751, 0, 0, 0.507974,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506323,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.475336,-99) , 
1, 2.19411, 0, 0, 0.484744,-99) , 
1, 0.299198, 1, 0, 0.491027,-99) , 
3, 0.0940732, 0, 0, 0.495992,-99)    );
  // itree = 48
  fBoostWeights.push_back(0.0613439);
  fForest.push_back( 
NN(
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.542109,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.487078,-99) , 
6, 1.22414, 1, 0, 0.529989,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.507342,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.4646,-99) , 
5, 133.442, 1, 0, 0.495764,-99) , 
7, 1.04682, 0, 0, 0.514654,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.527358,-99) , 
NN(
0, 
0, 
-1, 1.04716, 0, -1, 0.459099,-99) , 
3, -0.00184957, 0, 0, 0.47393,-99) , 
6, -0.342796, 0, 0, 0.499713,-99)    );
  // itree = 49
  fBoostWeights.push_back(0.0388974);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.557441,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.508518,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.468196,-99) , 
6, -1.02796, 0, 0, 0.502062,-99) , 
4, 374.204, 0, 0, 0.506956,-99) , 
NN(
0, 
0, 
-1, 1.79522, 0, -1, 0.454055,-99) , 
4, 133.323, 0, 0, 0.501261,-99)    );
   return;
};
 
// Clean up
inline void ReadBDT_50_50_03_20::Clear() 
{
   for (unsigned int itree=0; itree<fForest.size(); itree++) { 
      delete fForest[itree]; 
   }
}
   inline double ReadBDT_50_50_03_20::GetMvaValue( const std::vector<double>& inputValues ) const
   {
      // classifier response value
      double retval = 0;

      // classifier response, sanity check first
      if (!IsStatusClean()) {
         std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
                   << " because status is dirty" << std::endl;
         retval = 0;
      }
      else {
         if (IsNormalised()) {
            // normalise variables
            std::vector<double> iV;
            iV.reserve(inputValues.size());
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));
            }
            retval = GetMvaValue__( iV );
         }
         else {
            retval = GetMvaValue__( inputValues );
         }
      }

      return retval;
   }
