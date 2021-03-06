// Class: ReadBDT_50_50_03_10
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : BDT::BDT_50_50_03_10
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.06/01       [394753]
Creator        : matthies
Date           : Tue Sep 18 21:15:16 2018
Host           : Linux vocms0301.cern.ch 2.6.32-642.4.2.el6.x86_64 #1 SMP Wed Aug 24 09:19:54 CEST 2016 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb
Training events: 189370
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
H: "False" [Print method-specific help message]
NTrees: "50" [Number of trees in the forest]
MaxDepth: "3" [Max depth of the decision tree allowed]
MinNodeSize: "5%" [Minimum percentage of training events required in a leaf node (default: Classification: 5%, Regression: 0.2%)]
nCuts: "10" [Number of grid points in variable range used in finding optimal cut in node splitting]
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

NVar 9
tmva_n_btags                  tmva_n_btags                  tmva_n_btags                  N of b-tags                                                     'I'    [2,6]
tmva_deltaphi_bottomlepton    tmva_deltaphi_bottomlepton    tmva_deltaphi_bottomlepton    #Delta#phi(l,b)               rad                               'F'    [1.07288360596e-06,3.14150905609]
tmva_deltaphi_bottomtoptag    tmva_deltaphi_bottomtoptag    tmva_deltaphi_bottomtoptag    #Delta#phi(t,b)               rad                               'F'    [0.000285387039185,3.14159035683]
tmva_pseudotop_mass           tmva_pseudotop_mass           tmva_pseudotop_mass           l#nub mass                    GeV                               'F'    [86.5497436523,3728.6262207]
tmva_ptbalance                tmva_ptbalance                tmva_ptbalance                pT balance                                                      'F'    [-0.998895287514,3.84762001038]
tmva_wass_pt                  tmva_wass_pt                  tmva_wass_pt                  asso. W pT                    GeV                               'F'    [0.287396162748,1448.25354004]
tmva_lepton_pt                tmva_lepton_pt                tmva_lepton_pt                lepton pT                     GeV                               'F'    [50.0000190735,1128.80273438]
tmva_lepton_eta               tmva_lepton_eta               tmva_lepton_eta               lepton #eta                                                     'F'    [-2.39984273911,2.39958810806]
tmva_deltaphi_leptonnextjet   tmva_deltaphi_leptonnextjet   tmva_deltaphi_leptonnextjet   #Delta#phi(l,jet)             rad                               'F'    [1.07288360596e-06,3.13939738274]
NSpec 0


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#define NN new BDT_50_50_03_10Node
   
#ifndef BDT_50_50_03_10Node__def
#define BDT_50_50_03_10Node__def
   
class BDT_50_50_03_10Node {
   
public:
   
   // constructor of an essentially "empty" node floating in space
   BDT_50_50_03_10Node ( BDT_50_50_03_10Node* left,BDT_50_50_03_10Node* right,
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

   virtual ~BDT_50_50_03_10Node();

   // test event if it decends the tree at this node to the right
   virtual bool GoesRight( const std::vector<double>& inputValues ) const;
   BDT_50_50_03_10Node* GetRight( void )  {return fRight; };

   // test event if it decends the tree at this node to the left 
   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;
   BDT_50_50_03_10Node* GetLeft( void ) { return fLeft; };   

   // return  S/(S+B) (purity) at this node (from  training)

   double GetPurity( void ) const { return fPurity; } 
   // return the node type
   int    GetNodeType( void ) const { return fNodeType; }
   double GetResponse(void) const {return fResponse;}

private:

   BDT_50_50_03_10Node*   fLeft;     // pointer to the left daughter node
   BDT_50_50_03_10Node*   fRight;    // pointer to the right daughter node
   int                     fSelector; // index of variable used in node selection (decision tree)   
   double                  fCutValue; // cut value appplied on this node to discriminate bkg against sig
   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise
   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal 
   double                  fPurity;   // Purity of node from training
   double                  fResponse; // Regression response value of node
}; 
   
//_______________________________________________________________________
   BDT_50_50_03_10Node::~BDT_50_50_03_10Node()
{
   if (fLeft  != NULL) delete fLeft;
   if (fRight != NULL) delete fRight;
}; 
   
//_______________________________________________________________________
bool BDT_50_50_03_10Node::GoesRight( const std::vector<double>& inputValues ) const
{
   // test event if it decends the tree at this node to the right
   bool result;
     result = (inputValues[fSelector] > fCutValue );
   if (fCutType == true) return result; //the cuts are selecting Signal ;
   else return !result;
}
   
//_______________________________________________________________________
bool BDT_50_50_03_10Node::GoesLeft( const std::vector<double>& inputValues ) const
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

class ReadBDT_50_50_03_10 : public IClassifierReader {

 public:

   // constructor
   ReadBDT_50_50_03_10( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadBDT_50_50_03_10" ),
        fNvars( 9 ),
        fIsNormalised( false )
   {      
      // the training input variables
      const char* inputVars[] = { "tmva_n_btags", "tmva_deltaphi_bottomlepton", "tmva_deltaphi_bottomtoptag", "tmva_pseudotop_mass", "tmva_ptbalance", "tmva_wass_pt", "tmva_lepton_pt", "tmva_lepton_eta", "tmva_deltaphi_leptonnextjet" };

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
      fVmin[8] = 0;
      fVmax[8] = 0;

      // initialize input variable types
      fType[0] = 'I';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';
      fType[5] = 'F';
      fType[6] = 'F';
      fType[7] = 'F';
      fType[8] = 'F';

      // initialize constants
      Initialize();

   }

   // destructor
   virtual ~ReadBDT_50_50_03_10() {
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
   double fVmin[9];
   double fVmax[9];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[9];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)
   std::vector<BDT_50_50_03_10Node*> fForest;       // i.e. root nodes of decision trees
   std::vector<double>                fBoostWeights; // the weights applied in the individual boosts
};

double ReadBDT_50_50_03_10::GetMvaValue__( const std::vector<double>& inputValues ) const
{
   double myMVA = 0;
   double norm  = 0;
   for (unsigned int itree=0; itree<fForest.size(); itree++){
      BDT_50_50_03_10Node *current = fForest[itree];
      while (current->GetNodeType() == 0) { //intermediate node
         if (current->GoesRight(inputValues)) current=(BDT_50_50_03_10Node*)current->GetRight();
         else current=(BDT_50_50_03_10Node*)current->GetLeft();
      }
      myMVA += fBoostWeights[itree] *  current->GetNodeType();
      norm  += fBoostWeights[itree];
   }
   return myMVA /= norm;
};

void ReadBDT_50_50_03_10::Initialize()
{
  // itree = 0
  fBoostWeights.push_back(0.389074909527247);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 261.754, 0, 1, 0.795071,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.493351,-99) , 
5, 186.302, 0, 0, 0.713064,-99) , 
NN(
0, 
0, 
-1, 1.11561, 0, -1, 0.336839,-99) , 
1, 1.42793, 0, 0, 0.492341,-99)    );
  // itree = 1
  fBoostWeights.push_back(0.369097);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 203.68, 0, 1, 0.719703,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.537052,-99) , 
NN(
0, 
0, 
-1, -0.271853, 0, -1, 0.328424,-99) , 
8, 1.10382, 0, 0, 0.367949,-99) , 
1, 1.42796, 0, 0, 0.513744,-99)    );
  // itree = 2
  fBoostWeights.push_back(0.443513);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.850721,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.740119,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.415357,-99) , 
5, 175.282, 0, 0, 0.636198,-99) , 
NN(
0, 
0, 
-1, 197.166, 0, -1, 0.329981,-99) , 
3, 221.208, 0, 0, 0.441615,-99) , 
3, 379.918, 0, 0, 0.514422,-99)    );
  // itree = 3
  fBoostWeights.push_back(0.365804);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 232.496, 0, 1, 0.703137,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.727609,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.526791,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.321751,-99) , 
3, 229.549, 0, 0, 0.35825,-99) , 
3, 311.271, 0, 0, 0.409212,-99) , 
1, 1.71352, 0, 0, 0.504104,-99)    );
  // itree = 4
  fBoostWeights.push_back(0.293615);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.723823,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.699638,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.52321,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.357702,-99) , 
3, 218.83, 0, 0, 0.409701,-99) , 
3, 377.7, 0, 0, 0.441728,-99) , 
8, 1.4262, 0, 0, 0.51293,-99)    );
  // itree = 5
  fBoostWeights.push_back(0.349519);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, -0.350541, 0, 1, 0.772195,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.669278,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.423738,-99) , 
4, -0.281452, 0, 0, 0.55503,-99) , 
NN(
0, 
0, 
-1, 132.343, 1, -1, 0.360328,-99) , 
3, 212.486, 0, 0, 0.428946,-99) , 
3, 317.495, 0, 0, 0.511724,-99)    );
  // itree = 6
  fBoostWeights.push_back(0.261666);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.759941,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.65103,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.478194,-99) , 
6, 139.637, 0, 0, 0.562985,-99) , 
2, 1.9992, 1, 0, 0.688292,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.644653,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.500768,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.34817,-99) , 
3, 204.999, 0, 0, 0.408955,-99) , 
3, 347.166, 0, 0, 0.441052,-99) , 
4, -0.116278, 0, 0, 0.514722,-99)    );
  // itree = 7
  fBoostWeights.push_back(0.23448);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.715136,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.639951,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.574793,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.423683,-99) , 
3, 273.454, 0, 0, 0.446711,-99) , 
3, 379.278, 0, 0, 0.467132,-99) , 
8, 1.427, 0, 0, 0.527686,-99)    );
  // itree = 8
  fBoostWeights.push_back(0.246629);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, -0.344044, 0, 1, 0.728466,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.669741,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.44069,-99) , 
2, 2.57033, 1, 0, 0.60536,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.577825,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.41985,-99) , 
4, 0.0135359, 0, 0, 0.434239,-99) , 
3, 247, 0, 0, 0.47327,-99) , 
3, 337.98, 0, 0, 0.525268,-99)    );
  // itree = 9
  fBoostWeights.push_back(0.198727);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 293.559, 0, 1, 0.687649,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.561366,-99) , 
NN(
0, 
0, 
-1, 2.28502, 1, -1, 0.441741,-99) , 
3, 304.847, 0, 0, 0.464818,-99) , 
4, -0.113641, 0, 0, 0.532672,-99)    );
  // itree = 10
  fBoostWeights.push_back(0.206628);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 1.42794, 0, 1, 0.672876,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.586763,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.516903,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.423978,-99) , 
3, 252.722, 0, 0, 0.438845,-99) , 
3, 314.532, 0, 0, 0.466334,-99) , 
4, -0.116278, 0, 0, 0.529874,-99)    );
  // itree = 11
  fBoostWeights.push_back(0.17296);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.680384,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.610778,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.505401,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.427183,-99) , 
4, -0.26537, 0, 0, 0.455865,-99) , 
4, 0.00835653, 0, 0, 0.469597,-99) , 
8, 1.427, 0, 0, 0.522781,-99)    );
  // itree = 12
  fBoostWeights.push_back(0.203639);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.762842,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.45529,-99) , 
4, -0.349655, 0, 0, 0.673191,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.57431,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.593253,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.43351,-99) , 
3, 129.874, 1, 0, 0.44847,-99) , 
2, 1.71371, 1, 0, 0.479348,-99) , 
3, 318.731, 0, 0, 0.524646,-99)    );
  // itree = 13
  fBoostWeights.push_back(0.155673);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.718543,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.616724,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501263,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.402865,-99) , 
8, 0.620125, 0, 0, 0.464265,-99) , 
3, 347.166, 0, 0, 0.485748,-99) , 
8, 1.70625, 0, 0, 0.525057,-99)    );
  // itree = 14
  fBoostWeights.push_back(0.166644);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.691594,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.615774,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.572178,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.454756,-99) , 
3, 251.471, 0, 0, 0.485176,-99) , 
3, 135.256, 1, 0, 0.498188,-99) , 
3, 348.403, 0, 0, 0.533086,-99)    );
  // itree = 15
  fBoostWeights.push_back(0.129262);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 291.497, 0, 1, 0.635082,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.562765,-99) , 
NN(
0, 
0, 
-1, 183.599, 0, -1, 0.469612,-99) , 
3, 349.31, 0, 0, 0.482359,-99) , 
4, -0.112738, 0, 0, 0.527225,-99)    );
  // itree = 16
  fBoostWeights.push_back(0.145456);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.668534,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.571276,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.573914,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.459991,-99) , 
3, 133.934, 1, 0, 0.472535,-99) , 
2, 1.42854, 1, 0, 0.494181,-99) , 
3, 347.166, 0, 0, 0.527128,-99)    );
  // itree = 17
  fBoostWeights.push_back(0.12863);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.774579,-99) , 
NN(
NN(
0, 
0, 
-1, 271.471, 0, 1, 0.555295,-99) , 
NN(
0, 
0, 
-1, 0.155704, 0, -1, 0.457686,-99) , 
8, 0.570832, 0, 0, 0.526159,-99) , 
6, 241.823, 0, 0, 0.540121,-99)    );
  // itree = 18
  fBoostWeights.push_back(0.124137);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 272.224, 0, 1, 0.618297,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.527234,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.483411,-99) , 
2, 1.99919, 0, 0, 0.511223,-99) , 
NN(
0, 
0, 
-1, -0.598366, 0, -1, 0.451876,-99) , 
8, 0.85271, 0, 0, 0.47666,-99) , 
4, -0.117711, 0, 0, 0.518834,-99)    );
  // itree = 19
  fBoostWeights.push_back(0.144498);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, -0.350866, 0, 1, 0.655483,-99) , 
NN(
NN(
0, 
0, 
-1, 149.089, 0, 1, 0.539344,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.587429,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.445824,-99) , 
3, 127.872, 1, 0, 0.462346,-99) , 
3, 197.806, 0, 0, 0.497009,-99) , 
3, 329.702, 0, 0, 0.530123,-99)    );
  // itree = 20
  fBoostWeights.push_back(0.121402);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.644701,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.554336,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.490026,-99) , 
2, 2.28487, 1, 0, 0.515084,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.410837,-99) , 
5, 116.182, 0, 0, 0.501805,-99) , 
8, 1.42621, 0, 0, 0.534751,-99)    );
  // itree = 21
  fBoostWeights.push_back(0.0915562);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.686849,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.625121,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481147,-99) , 
8, 0.856099, 0, 0, 0.559406,-99) , 
NN(
0, 
0, 
-1, 1.71394, 0, -1, 0.483865,-99) , 
6, 132.527, 0, 0, 0.50325,-99) , 
5, 352.326, 0, 0, 0.519806,-99)    );
  // itree = 22
  fBoostWeights.push_back(0.0955709);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.657656,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.623289,-99) , 
NN(
0, 
0, 
-1, 127.933, 0, -1, 0.482419,-99) , 
5, 349.195, 0, 0, 0.494429,-99) , 
8, 1.69968, 0, 0, 0.520642,-99)    );
  // itree = 23
  fBoostWeights.push_back(0.101356);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.725466,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.592648,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.515382,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.471736,-99) , 
8, 1.1416, 0, 0, 0.484574,-99) , 
4, -0.00501118, 0, 0, 0.49867,-99) , 
5, 366.653, 0, 0, 0.518904,-99)    );
  // itree = 24
  fBoostWeights.push_back(0.0897668);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.704958,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.596286,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.498759,-99) , 
3, 311.271, 0, 0, 0.51762,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.34764,-99) , 
5, 101.891, 0, 0, 0.50458,-99) , 
5, 367.778, 0, 0, 0.521692,-99)    );
  // itree = 25
  fBoostWeights.push_back(0.117892);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.701965,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.632384,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.447182,-99) , 
8, 0.84561, 0, 0, 0.553861,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.522489,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470576,-99) , 
2, 2.57073, 0, 0, 0.494148,-99) , 
6, 130.389, 0, 0, 0.510086,-99) , 
5, 367.278, 0, 0, 0.526139,-99)    );
  // itree = 26
  fBoostWeights.push_back(0.096625);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.682227,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.605463,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.514674,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.401289,-99) , 
5, 100.134, 0, 0, 0.504561,-99) , 
4, -0.04103, 0, 0, 0.521614,-99) , 
5, 366.392, 0, 0, 0.534724,-99)    );
  // itree = 27
  fBoostWeights.push_back(0.0880896);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.656209,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.530807,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.488758,-99) , 
8, 0.795169, 0, 0, 0.511014,-99) , 
1, 1.42798, 0, 0, 0.599092,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.519291,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.491153,-99) , 
1, 1.06453, 1, 0, 0.505766,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.439959,-99) , 
1, 0.285599, 0, 0, 0.497901,-99) , 
4, -0.119975, 0, 0, 0.527211,-99)    );
  // itree = 28
  fBoostWeights.push_back(0.0893994);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.668714,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.587818,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481243,-99) , 
2, 2.71718, 0, 0, 0.53428,-99) , 
1, 1.1424, 0, 0, 0.605713,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.568217,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.497635,-99) , 
3, 293.463, 0, 0, 0.51163,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.40835,-99) , 
5, 97.7611, 0, 0, 0.502483,-99) , 
6, 148.074, 0, 0, 0.527525,-99)    );
  // itree = 29
  fBoostWeights.push_back(0.0961851);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.63714,-99) , 
NN(
NN(
0, 
0, 
-1, 0.216093, 1, 1, 0.54179,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.525424,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470714,-99) , 
8, 0.495014, 1, 0, 0.49438,-99) , 
8, 1.08908, 0, 0, 0.506894,-99) , 
8, 1.71145, 0, 0, 0.527318,-99)    );
  // itree = 30
  fBoostWeights.push_back(0.0855951);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 1.14222, 0, 1, 0.61121,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.562895,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.492926,-99) , 
5, 300.407, 0, 0, 0.49953,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.418657,-99) , 
5, 107.012, 0, 0, 0.490474,-99) , 
6, 145.914, 0, 0, 0.520445,-99)    );
  // itree = 31
  fBoostWeights.push_back(0.0869081);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.626881,-99) , 
NN(
NN(
0, 
0, 
-1, 317.495, 0, 1, 0.524751,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.449045,-99) , 
5, 118.345, 0, 0, 0.515305,-99) , 
8, 1.71145, 0, 0, 0.5326,-99)    );
  // itree = 32
  fBoostWeights.push_back(0.0990056);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.637052,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.525497,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.449171,-99) , 
1, 0.285458, 0, 0, 0.517806,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.41205,-99) , 
4, -0.637409, 0, 0, 0.50856,-99) , 
8, 1.71206, 0, 0, 0.528442,-99)    );
  // itree = 33
  fBoostWeights.push_back(0.072828);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.645678,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.582666,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.535101,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.486389,-99) , 
2, 2.85602, 0, 0, 0.497824,-99) , 
6, 185.169, 0, 0, 0.504524,-99) , 
5, 343.992, 0, 0, 0.519151,-99)    );
  // itree = 34
  fBoostWeights.push_back(0.0866038);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.635031,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.615681,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.498694,-99) , 
4, -0.574263, 1, 0, 0.513,-99) , 
3, 304.392, 0, 0, 0.544026,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501958,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.436207,-99) , 
3, 283.622, 1, 0, 0.486763,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.374759,-99) , 
8, 0.467031, 1, 0, 0.465384,-99) , 
8, 0.570818, 0, 0, 0.521007,-99)    );
  // itree = 35
  fBoostWeights.push_back(0.0942679);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.621885,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.591708,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.534045,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.484041,-99) , 
2, 2.28505, 1, 0, 0.507057,-99) , 
3, 136.238, 1, 0, 0.515538,-99) , 
3, 349.31, 0, 0, 0.534217,-99)    );
  // itree = 36
  fBoostWeights.push_back(0.0699579);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.643802,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.56721,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.499874,-99) , 
3, 301.74, 0, 0, 0.515695,-99) , 
NN(
0, 
0, 
-1, -0.32636, 0, -1, 0.445642,-99) , 
2, 2.856, 1, 0, 0.501773,-99) , 
6, 203.4, 0, 0, 0.515222,-99)    );
  // itree = 37
  fBoostWeights.push_back(0.0856006);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.611958,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.529557,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476858,-99) , 
7, 0.811951, 1, 0, 0.516258,-99) , 
NN(
0, 
0, 
-1, -0.354282, 1, -1, 0.449956,-99) , 
7, -1.09123, 0, 0, 0.507295,-99) , 
4, 0.0056194, 0, 0, 0.522497,-99)    );
  // itree = 38
  fBoostWeights.push_back(0.0862056);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.600179,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.597599,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.52163,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.425826,-99) , 
7, 0.813741, 1, 0, 0.507774,-99) , 
7, 1.52759, 0, 0, 0.515169,-99) , 
4, 0.0197146, 0, 0, 0.525761,-99)    );
  // itree = 39
  fBoostWeights.push_back(0.0729605);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.682793,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.55891,-99) , 
NN(
0, 
0, 
-1, 319.001, 1, -1, 0.482274,-99) , 
8, 1.4195, 0, 0, 0.50004,-99) , 
5, 396.637, 0, 0, 0.510784,-99)    );
  // itree = 40
  fBoostWeights.push_back(0.0617725);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.632914,-99) , 
NN(
NN(
0, 
0, 
-1, 379.608, 0, 1, 0.523243,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.444149,-99) , 
4, -0.707798, 0, 0, 0.519027,-99) , 
6, 229.281, 0, 0, 0.525216,-99)    );
  // itree = 41
  fBoostWeights.push_back(0.0609264);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.588205,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.57143,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.494801,-99) , 
1, 0.285529, 1, 0, 0.504008,-99) , 
3, 338.263, 0, 0, 0.518579,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.421064,-99) , 
4, -0.631269, 0, 0, 0.510671,-99)    );
  // itree = 42
  fBoostWeights.push_back(0.0703737);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.606592,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.510626,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.436634,-99) , 
0, 3, 1, 0, 0.503193,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.440127,-99) , 
2, 1.14258, 0, 0, 0.493847,-99) , 
8, 1.70641, 0, 0, 0.510852,-99)    );
  // itree = 43
  fBoostWeights.push_back(0.0735051);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.608695,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.543944,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.460318,-99) , 
5, 141.747, 0, 0, 0.52025,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.502671,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.412788,-99) , 
1, 1.14231, 0, 0, 0.480073,-99) , 
2, 2.57057, 0, 0, 0.497412,-99) , 
5, 385.681, 0, 0, 0.504152,-99)    );
  // itree = 44
  fBoostWeights.push_back(0.060709);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.599477,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.572595,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.487698,-99) , 
6, 195.811, 0, 0, 0.49318,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.391927,-99) , 
5, 95.7777, 0, 0, 0.486031,-99) , 
5, 346.845, 0, 0, 0.49731,-99)    );
  // itree = 45
  fBoostWeights.push_back(0.107532);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.650315,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.532071,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.433447,-99) , 
8, 0.466368, 1, 0, 0.48558,-99) , 
8, 0.855039, 0, 0, 0.574007,-99) , 
NN(
NN(
0, 
0, 
-1, 1.14194, 1, 1, 0.561585,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.523405,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.469309,-99) , 
6, 112.412, 0, 0, 0.482155,-99) , 
7, 0.653938, 0, 0, 0.505899,-99) , 
6, 148.075, 0, 0, 0.52161,-99)    );
  // itree = 46
  fBoostWeights.push_back(0.0651557);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.592455,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.584415,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.504438,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.467711,-99) , 
8, 0.853601, 1, 0, 0.488306,-99) , 
4, 0.0104284, 0, 0, 0.500624,-99) , 
0, 3, 0, 0, 0.510447,-99)    );
  // itree = 47
  fBoostWeights.push_back(0.077861);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.606584,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.515812,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.463716,-99) , 
1, 0.571138, 0, 0, 0.502282,-99) , 
NN(
0, 
0, 
-1, -1.44786, 0, -1, 0.43718,-99) , 
7, -1.09092, 0, 0, 0.494152,-99) , 
8, 1.70646, 0, 0, 0.511689,-99)    );
  // itree = 48
  fBoostWeights.push_back(0.0818478);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.588811,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.526801,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470366,-99) , 
7, -0.655879, 0, 0, 0.51323,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509594,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.334231,-99) , 
8, 1.00906, 1, 0, 0.465113,-99) , 
3, 178.436, 0, 0, 0.49615,-99) , 
3, 419.153, 0, 0, 0.505894,-99)    );
  // itree = 49
  fBoostWeights.push_back(0.0721538);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.610338,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.534147,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.492522,-99) , 
4, -0.221341, 1, 0, 0.516013,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.526358,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.431119,-99) , 
5, 233.548, 0, 0, 0.459999,-99) , 
1, 0.571172, 0, 0, 0.503743,-99) , 
5, 366.653, 0, 0, 0.511897,-99)    );
   return;
};
 
// Clean up
inline void ReadBDT_50_50_03_10::Clear() 
{
   for (unsigned int itree=0; itree<fForest.size(); itree++) { 
      delete fForest[itree]; 
   }
}
   inline double ReadBDT_50_50_03_10::GetMvaValue( const std::vector<double>& inputValues ) const
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
