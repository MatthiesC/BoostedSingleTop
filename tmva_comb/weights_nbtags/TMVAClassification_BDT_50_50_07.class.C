// Class: ReadBDT_50_50_07
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : BDT::BDT_50_50_07
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.06/01       [394753]
Creator        : matthies
Date           : Wed Sep 12 15:12:02 2018
Host           : Linux vocms0301.cern.ch 2.6.32-642.4.2.el6.x86_64 #1 SMP Wed Aug 24 09:19:54 CEST 2016 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb
Training events: 463792
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
AdaBoostBeta: "7.000000e-01" [Learning rate  for AdaBoost algorithm]
BaggedSampleFraction: "5.000000e-01" [Relative size of bagged event sample to original size of the data sample (used whenever bagging is used (i.e. UseBaggedBoost, Bagging,)]
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
GradBaggingFraction: "5.000000e-01" [deprecated: Use *BaggedSampleFraction* instead: Defines the fraction of events to be used in each iteration, e.g. when UseBaggedGrad=kTRUE. ]
UseNTrainEvents: "0" [deprecated: Use *BaggedSampleFraction* instead: Number of randomly picked training events used in randomised (and bagged) trees]
NNodesMax: "0" [deprecated: Use MaxDepth instead to limit the tree size]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 9
tmva_n_btags                  tmva_n_btags                  tmva_n_btags                  N of b-tags                                                     'I'    [1,6]
tmva_deltaphi_bottomlepton    tmva_deltaphi_bottomlepton    tmva_deltaphi_bottomlepton    #Delta#phi(l,b)               rad                               'F'    [2.32458114624e-06,3.14158797264]
tmva_deltaphi_bottomtoptag    tmva_deltaphi_bottomtoptag    tmva_deltaphi_bottomtoptag    #Delta#phi(t,b)               rad                               'F'    [4.29153442383e-06,3.14158964157]
tmva_pseudotop_mass           tmva_pseudotop_mass           tmva_pseudotop_mass           l#nub mass                    GeV                               'F'    [86.2544937134,6195.96142578]
tmva_ptbalance                tmva_ptbalance                tmva_ptbalance                pT balance                                                      'F'    [-0.997143745422,3.99755883217]
tmva_wass_pt                  tmva_wass_pt                  tmva_wass_pt                  asso. W pT                    GeV                               'F'    [0.646613299847,2188.19360352]
tmva_lepton_pt                tmva_lepton_pt                tmva_lepton_pt                lepton pT                     GeV                               'F'    [50.0001106262,1194.46325684]
tmva_lepton_eta               tmva_lepton_eta               tmva_lepton_eta               lepton #eta                                                     'F'    [-2.39999055862,2.39999628067]
tmva_deltaphi_leptonnextjet   tmva_deltaphi_leptonnextjet   tmva_deltaphi_leptonnextjet   #Delta#phi(l,jet)             rad                               'F'    [2.32458114624e-06,3.14143514633]
NSpec 0


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#define NN new BDT_50_50_07Node
   
#ifndef BDT_50_50_07Node__def
#define BDT_50_50_07Node__def
   
class BDT_50_50_07Node {
   
public:
   
   // constructor of an essentially "empty" node floating in space
   BDT_50_50_07Node ( BDT_50_50_07Node* left,BDT_50_50_07Node* right,
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

   virtual ~BDT_50_50_07Node();

   // test event if it decends the tree at this node to the right
   virtual bool GoesRight( const std::vector<double>& inputValues ) const;
   BDT_50_50_07Node* GetRight( void )  {return fRight; };

   // test event if it decends the tree at this node to the left 
   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;
   BDT_50_50_07Node* GetLeft( void ) { return fLeft; };   

   // return  S/(S+B) (purity) at this node (from  training)

   double GetPurity( void ) const { return fPurity; } 
   // return the node type
   int    GetNodeType( void ) const { return fNodeType; }
   double GetResponse(void) const {return fResponse;}

private:

   BDT_50_50_07Node*   fLeft;     // pointer to the left daughter node
   BDT_50_50_07Node*   fRight;    // pointer to the right daughter node
   int                     fSelector; // index of variable used in node selection (decision tree)   
   double                  fCutValue; // cut value appplied on this node to discriminate bkg against sig
   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise
   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal 
   double                  fPurity;   // Purity of node from training
   double                  fResponse; // Regression response value of node
}; 
   
//_______________________________________________________________________
   BDT_50_50_07Node::~BDT_50_50_07Node()
{
   if (fLeft  != NULL) delete fLeft;
   if (fRight != NULL) delete fRight;
}; 
   
//_______________________________________________________________________
bool BDT_50_50_07Node::GoesRight( const std::vector<double>& inputValues ) const
{
   // test event if it decends the tree at this node to the right
   bool result;
     result = (inputValues[fSelector] > fCutValue );
   if (fCutType == true) return result; //the cuts are selecting Signal ;
   else return !result;
}
   
//_______________________________________________________________________
bool BDT_50_50_07Node::GoesLeft( const std::vector<double>& inputValues ) const
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

class ReadBDT_50_50_07 : public IClassifierReader {

 public:

   // constructor
   ReadBDT_50_50_07( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadBDT_50_50_07" ),
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
   virtual ~ReadBDT_50_50_07() {
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
   std::vector<BDT_50_50_07Node*> fForest;       // i.e. root nodes of decision trees
   std::vector<double>                fBoostWeights; // the weights applied in the individual boosts
};

double ReadBDT_50_50_07::GetMvaValue__( const std::vector<double>& inputValues ) const
{
   double myMVA = 0;
   double norm  = 0;
   for (unsigned int itree=0; itree<fForest.size(); itree++){
      BDT_50_50_07Node *current = fForest[itree];
      while (current->GetNodeType() == 0) { //intermediate node
         if (current->GoesRight(inputValues)) current=(BDT_50_50_07Node*)current->GetRight();
         else current=(BDT_50_50_07Node*)current->GetLeft();
      }
      myMVA += fBoostWeights[itree] *  current->GetNodeType();
      norm  += fBoostWeights[itree];
   }
   return myMVA /= norm;
};

void ReadBDT_50_50_07::Initialize()
{
  // itree = 0
  fBoostWeights.push_back(0.594933950548316);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 158.956, 0, 1, 0.739957,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.614899,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.41257,-99) , 
5, 225.541, 0, 0, 0.506132,-99) , 
NN(
0, 
0, 
-1, 219.894, 0, -1, 0.267748,-99) , 
3, 261.181, 0, 0, 0.365169,-99) , 
8, 1.4957, 0, 0, 0.500564,-99)    );
  // itree = 1
  fBoostWeights.push_back(0.431685);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 131.372, 0, 1, 0.699956,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.648077,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476953,-99) , 
8, 1.01859, 0, 0, 0.535825,-99) , 
NN(
0, 
0, 
-1, 1.94486, 1, -1, 0.329896,-99) , 
3, 264.704, 0, 0, 0.417032,-99) , 
8, 1.64552, 0, 0, 0.506594,-99)    );
  // itree = 2
  fBoostWeights.push_back(0.32292);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.736337,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.625955,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.480226,-99) , 
1, 1.94485, 0, 0, 0.586003,-99) , 
6, 158.955, 0, 0, 0.646221,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.602972,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476632,-99) , 
6, 133.443, 0, 0, 0.511111,-99) , 
NN(
0, 
0, 
-1, 2, 1, -1, 0.340265,-99) , 
3, 215.959, 0, 0, 0.440706,-99) , 
4, -0.100289, 0, 0, 0.512143,-99)    );
  // itree = 3
  fBoostWeights.push_back(0.253006);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 179.524, 0, 1, 0.580684,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.546546,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.441264,-99) , 
1, 0.833407, 0, 0, 0.483139,-99) , 
NN(
0, 
0, 
-1, 192.006, 0, -1, 0.350002,-99) , 
0, 2, 1, 0, 0.408496,-99) , 
1, 1.3464, 0, 0, 0.525356,-99)    );
  // itree = 4
  fBoostWeights.push_back(0.251926);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 213.432, 0, 1, 0.574612,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.588657,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.417467,-99) , 
4, -0.342948, 0, 0, 0.52709,-99) , 
NN(
0, 
0, 
-1, 2, 1, -1, 0.387388,-99) , 
3, 218.594, 0, 0, 0.429989,-99) , 
3, 270.909, 0, 0, 0.505234,-99)    );
  // itree = 5
  fBoostWeights.push_back(0.206739);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 158.964, 0, 1, 0.595563,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.550189,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.450855,-99) , 
6, 98.5768, 0, 0, 0.498111,-99) , 
NN(
0, 
0, 
-1, 2, 1, -1, 0.427516,-99) , 
8, 1.04715, 0, 0, 0.455036,-99) , 
4, -0.0734563, 0, 0, 0.49776,-99)    );
  // itree = 6
  fBoostWeights.push_back(0.146708);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 159.013, 0, 1, 0.583853,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.535583,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.471297,-99) , 
8, 1.19664, 0, 0, 0.505467,-99) , 
NN(
0, 
0, 
-1, 243.3, 0, -1, 0.453478,-99) , 
4, -0.224795, 0, 0, 0.473152,-99) , 
4, -0.0436745, 0, 0, 0.501786,-99)    );
  // itree = 7
  fBoostWeights.push_back(0.164375);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.67269,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.655861,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.478661,-99) , 
2, 1.0472, 0, 0, 0.51382,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.51021,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.411208,-99) , 
0, 2, 1, 0, 0.464845,-99) , 
3, 261.181, 0, 0, 0.490501,-99) , 
6, 213.495, 0, 0, 0.506859,-99)    );
  // itree = 8
  fBoostWeights.push_back(0.103612);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.636712,-99) , 
NN(
NN(
0, 
0, 
-1, 120.067, 0, 1, 0.548306,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503117,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.466759,-99) , 
3, 261.559, 0, 0, 0.482322,-99) , 
8, 1.79485, 0, 0, 0.499077,-99) , 
6, 213.497, 0, 0, 0.511032,-99)    );
  // itree = 9
  fBoostWeights.push_back(0.0912709);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.662185,-99) , 
NN(
NN(
0, 
0, 
-1, 119.387, 0, 1, 0.543489,-99) , 
NN(
0, 
0, 
-1, -0.800023, 0, -1, 0.482449,-99) , 
8, 1.94452, 0, 0, 0.496021,-99) , 
6, 232.102, 0, 0, 0.507068,-99)    );
  // itree = 10
  fBoostWeights.push_back(0.0971905);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.62793,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.554797,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.510722,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.446833,-99) , 
3, 201.744, 0, 0, 0.49135,-99) , 
5, 337.751, 0, 0, 0.495845,-99) , 
6, 213.495, 0, 0, 0.507295,-99)    );
  // itree = 11
  fBoostWeights.push_back(0.0928903);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 243.783, 0, 1, 0.612363,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.540424,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.486321,-99) , 
4, 0.0122683, 0, 0, 0.49463,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.433678,-99) , 
7, -1.25716, 0, 0, 0.487691,-99) , 
6, 200.717, 0, 0, 0.501432,-99)    );
  // itree = 12
  fBoostWeights.push_back(0.0627357);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 452.604, 0, 1, 0.587907,-99) , 
NN(
NN(
0, 
0, 
-1, 189.203, 0, 1, 0.509438,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.451395,-99) , 
7, -1.25716, 0, 0, 0.502914,-99) , 
5, 337.945, 0, 0, 0.512725,-99)    );
  // itree = 13
  fBoostWeights.push_back(0.0699817);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.593965,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.55217,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.496107,-99) , 
4, -0.0089794, 0, 0, 0.512039,-99) , 
NN(
0, 
0, 
-1, 0.598406, 0, -1, 0.480988,-99) , 
6, 120.066, 0, 0, 0.491646,-99) , 
6, 213.496, 0, 0, 0.500383,-99)    );
  // itree = 14
  fBoostWeights.push_back(0.0995442);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.602104,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.616864,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.499191,-99) , 
2, 1.19679, 0, 0, 0.521878,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.522891,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.436285,-99) , 
0, 2, 1, 0, 0.482134,-99) , 
3, 240.641, 0, 0, 0.505656,-99) , 
5, 385.62, 0, 0, 0.511869,-99)    );
  // itree = 15
  fBoostWeights.push_back(0.0787081);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.594024,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.521265,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481414,-99) , 
7, 1.25719, 1, 0, 0.515506,-99) , 
NN(
0, 
0, 
-1, 0.897769, 0, -1, 0.472826,-99) , 
7, -1.02838, 0, 0, 0.508606,-99) , 
6, 213.432, 0, 0, 0.51593,-99)    );
  // itree = 16
  fBoostWeights.push_back(0.14098);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.62818,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.559553,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.47836,-99) , 
6, 158.998, 0, 0, 0.501795,-99) , 
2, 1.04719, 0, 0, 0.527107,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.554827,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.521897,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.441502,-99) , 
2, 1.34645, 1, 0, 0.466608,-99) , 
3, 141.461, 1, 0, 0.477608,-99) , 
3, 248.768, 0, 0, 0.506357,-99)    );
  // itree = 17
  fBoostWeights.push_back(0.0611299);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 2.69356, 0, 1, 0.548157,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.578888,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.498386,-99) , 
5, 285.825, 0, 0, 0.512617,-99) , 
NN(
0, 
0, 
-1, 152.745, 1, -1, 0.484463,-99) , 
8, 0.897527, 0, 0, 0.497845,-99) , 
8, 2.09424, 0, 0, 0.507687,-99)    );
  // itree = 18
  fBoostWeights.push_back(0.0910922);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.601317,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.53481,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.498807,-99) , 
3, 283.543, 0, 0, 0.519215,-99) , 
6, 167.395, 0, 0, 0.537691,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.525668,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.468784,-99) , 
5, 109.412, 0, 0, 0.517885,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.503809,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.433948,-99) , 
2, 0.897627, 0, 0, 0.47293,-99) , 
5, 253.852, 1, 0, 0.506119,-99) , 
8, 1.19671, 0, 0, 0.520111,-99)    );
  // itree = 19
  fBoostWeights.push_back(0.0865626);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.577907,-99) , 
NN(
NN(
0, 
0, 
-1, 2, 0, 1, 0.51965,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.536846,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.454204,-99) , 
3, 149.671, 1, 0, 0.474976,-99) , 
3, 207.326, 0, 0, 0.505135,-99) , 
3, 594.767, 0, 0, 0.510124,-99)    );
  // itree = 20
  fBoostWeights.push_back(0.10212);
  fForest.push_back( 
NN(
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.605393,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.484014,-99) , 
5, 169.199, 0, 0, 0.56804,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.527131,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.468978,-99) , 
0, 2, 1, 0, 0.497546,-99) , 
3, 238.776, 0, 0, 0.521891,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.526,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.488102,-99) , 
5, 209.201, 1, 0, 0.502243,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.4353,-99) , 
7, 1.02862, 1, 0, 0.490995,-99) , 
2, 0.748003, 0, 0, 0.507715,-99)    );
  // itree = 21
  fBoostWeights.push_back(0.0966667);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.58905,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.567399,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.488629,-99) , 
8, 1.47417, 1, 0, 0.507627,-99) , 
6, 167.395, 0, 0, 0.525622,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.549242,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.497624,-99) , 
4, -0.3143, 0, 0, 0.513727,-99) , 
NN(
0, 
0, 
-1, -0.109549, 0, -1, 0.444109,-99) , 
5, 224.746, 1, 0, 0.489111,-99) , 
8, 1.19663, 0, 0, 0.505443,-99)    );
  // itree = 22
  fBoostWeights.push_back(0.080036);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.580094,-99) , 
NN(
NN(
0, 
0, 
-1, 2, 0, 1, 0.519829,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.513966,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.446117,-99) , 
2, 2.84244, 0, 0, 0.465243,-99) , 
3, 187.914, 0, 0, 0.506869,-99) , 
3, 619.973, 0, 0, 0.511288,-99)    );
  // itree = 23
  fBoostWeights.push_back(0.0586273);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.587783,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.554193,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.49693,-99) , 
0, 2, 0, 0, 0.509139,-99) , 
NN(
0, 
0, 
-1, 145.957, 0, -1, 0.477943,-99) , 
1, 1.04721, 0, 0, 0.501165,-99) , 
5, 382.086, 0, 0, 0.507069,-99)    );
  // itree = 24
  fBoostWeights.push_back(0.072712);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, -0.0354401, 0, 1, 0.534003,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.606483,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.491155,-99) , 
3, 309.983, 0, 0, 0.508313,-99) , 
NN(
0, 
0, 
-1, 226.348, 1, -1, 0.476546,-99) , 
2, 0.748, 0, 0, 0.49655,-99) , 
8, 1.9447, 0, 0, 0.505535,-99)    );
  // itree = 25
  fBoostWeights.push_back(0.0800646);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.615669,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.539027,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.47931,-99) , 
6, 158.998, 0, 0, 0.50848,-99) , 
1, 2.5434, 1, 0, 0.544002,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.550948,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.499362,-99) , 
0, 2, 0, 0, 0.512584,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.530012,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.439383,-99) , 
6, 89.7464, 1, 0, 0.486047,-99) , 
3, 192.824, 0, 0, 0.50476,-99) , 
3, 454.369, 0, 0, 0.511187,-99)    );
  // itree = 26
  fBoostWeights.push_back(0.0788256);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 2.54318, 1, 1, 0.535697,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.511769,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.479103,-99) , 
1, 2.54321, 1, 0, 0.500543,-99) , 
7, 0.723816, 1, 0, 0.524376,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.528691,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.492381,-99) , 
2, 2.45768, 1, 0, 0.512695,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.46885,-99) , 
2, 0.747994, 0, 0, 0.491869,-99) , 
7, -0.799954, 0, 0, 0.51726,-99)    );
  // itree = 27
  fBoostWeights.push_back(0.0501537);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 213.495, 0, 1, 0.515986,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.529096,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.483618,-99) , 
1, 1.6456, 1, 0, 0.50493,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.472806,-99) , 
5, 225.936, 1, 0, 0.491021,-99) , 
7, 0.800004, 1, 0, 0.510321,-99)    );
  // itree = 28
  fBoostWeights.push_back(0.0453173);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.559888,-99) , 
NN(
NN(
0, 
0, 
-1, 0.748169, 1, 1, 0.510727,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.500272,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.427478,-99) , 
5, 243.73, 1, 0, 0.481321,-99) , 
8, 0.598324, 0, 0, 0.502858,-99) , 
3, 593.129, 0, 0, 0.506768,-99)    );
  // itree = 29
  fBoostWeights.push_back(0.0907816);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.627005,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.561272,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.478788,-99) , 
2, 0.106858, 1, 0, 0.497316,-99) , 
2, 0.747986, 0, 0, 0.524091,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.537912,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.484756,-99) , 
1, 0.29919, 0, 0, 0.51899,-99) , 
NN(
0, 
0, 
-1, 196.433, 0, -1, 0.47896,-99) , 
8, 0.448656, 1, 0, 0.487868,-99) , 
3, 340.721, 0, 0, 0.500359,-99)    );
  // itree = 30
  fBoostWeights.push_back(0.0512568);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 2.24354, 0, 1, 0.517606,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501715,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.455409,-99) , 
8, 1.19665, 1, 0, 0.479916,-99) , 
7, -1.02869, 0, 0, 0.511537,-99)    );
  // itree = 31
  fBoostWeights.push_back(0.0745351);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.556253,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.55099,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.48553,-99) , 
2, 1.1968, 0, 0, 0.500335,-99) , 
NN(
0, 
0, 
-1, 0.114117, 1, -1, 0.471776,-99) , 
3, 217.61, 0, 0, 0.490734,-99) , 
5, 385.438, 0, 0, 0.494908,-99)    );
  // itree = 32
  fBoostWeights.push_back(0.0518473);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.55197,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.541605,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.510888,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.479626,-99) , 
6, 135.001, 1, 0, 0.502914,-99) , 
4, 0.0800765, 0, 0, 0.506659,-99) , 
3, 593.129, 0, 0, 0.509764,-99)    );
  // itree = 33
  fBoostWeights.push_back(0.0562272);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 315.058, 0, 1, 0.51596,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.513131,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.481028,-99) , 
3, 205.762, 1, 0, 0.494705,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.425808,-99) , 
6, 152.745, 1, 0, 0.481234,-99) , 
8, 0.747848, 0, 0, 0.503871,-99)    );
  // itree = 34
  fBoostWeights.push_back(0.0779248);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.585241,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.525394,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.470507,-99) , 
8, 1.34635, 0, 0, 0.501423,-99) , 
2, 0.598392, 0, 0, 0.522782,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.518683,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.491332,-99) , 
6, 79.476, 0, 0, 0.509904,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.506891,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.450088,-99) , 
0, 2, 1, 0, 0.475319,-99) , 
2, 1.7952, 1, 0, 0.494067,-99) , 
3, 377.193, 0, 0, 0.501893,-99)    );
  // itree = 35
  fBoostWeights.push_back(0.0662842);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.538787,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.536128,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.486716,-99) , 
0, 2, 1, 0, 0.511489,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.525521,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.479728,-99) , 
8, 0.448692, 1, 0, 0.486526,-99) , 
1, 1.94479, 1, 0, 0.498068,-99) , 
6, 235.501, 0, 0, 0.500354,-99)    );
  // itree = 36
  fBoostWeights.push_back(0.0526761);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 2.69509, 0, 1, 0.519163,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.52814,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.488611,-99) , 
0, 2, 0, 0, 0.498526,-99) , 
NN(
0, 
0, 
-1, -0.42711, 0, -1, 0.475628,-99) , 
1, 1.0472, 0, 0, 0.490502,-99) , 
8, 1.79506, 0, 0, 0.49839,-99)    );
  // itree = 37
  fBoostWeights.push_back(0.0435721);
  fForest.push_back( 
NN(
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.58191,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.492834,-99) , 
3, 261.574, 0, 0, 0.523148,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.543934,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.497524,-99) , 
4, -0.5478, 1, 0, 0.501871,-99) , 
0, 2, 0, 0, 0.507667,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.463929,-99) , 
1, 0.4488, 0, 0, 0.503696,-99)    );
  // itree = 38
  fBoostWeights.push_back(0.0572592);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.548098,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.512445,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.479916,-99) , 
8, 1.49554, 1, 0, 0.503538,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509216,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.44699,-99) , 
8, 0.747992, 0, 0, 0.488353,-99) , 
6, 136.416, 1, 0, 0.499462,-99) , 
4, 0.079509, 0, 0, 0.50477,-99)    );
  // itree = 39
  fBoostWeights.push_back(0.0631832);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.539884,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509126,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.454153,-99) , 
1, 0.897648, 0, 0, 0.503853,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.547326,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.467803,-99) , 
3, 141.239, 1, 0, 0.479195,-99) , 
3, 230.589, 0, 0, 0.494068,-99) , 
3, 591.47, 0, 0, 0.49723,-99)    );
  // itree = 40
  fBoostWeights.push_back(0.064058);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 2.09446, 1, 1, 0.533241,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.530983,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.486507,-99) , 
1, 2.09441, 1, 0, 0.499044,-99) , 
0, 2, 0, 0, 0.507729,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.514659,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476042,-99) , 
8, 0.59783, 0, 0, 0.498432,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.462093,-99) , 
5, 124.052, 0, 0, 0.488011,-99) , 
4, -0.304486, 0, 0, 0.500813,-99)    );
  // itree = 41
  fBoostWeights.push_back(0.0348207);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.524712,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.508119,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.490113,-99) , 
8, 1.19657, 0, 0, 0.498472,-99) , 
NN(
0, 
0, 
-1, 0.712429, 1, -1, 0.468961,-99) , 
7, 1.25715, 1, 0, 0.49507,-99) , 
2, 2.99199, 0, 0, 0.496956,-99)    );
  // itree = 42
  fBoostWeights.push_back(0.0460635);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.523582,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.509569,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.478701,-99) , 
3, 201.391, 0, 0, 0.500754,-99) , 
NN(
0, 
0, 
-1, 0.594162, 0, -1, 0.47575,-99) , 
4, -0.41413, 0, 0, 0.495003,-99) , 
2, 0.149606, 1, 0, 0.498267,-99)    );
  // itree = 43
  fBoostWeights.push_back(0.0669337);
  fForest.push_back( 
NN(
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.522505,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.495584,-99) , 
3, 261.635, 1, 0, 0.511914,-99) , 
NN(
0, 
0, 
-1, 1.3464, 0, -1, 0.468353,-99) , 
6, 152.744, 1, 0, 0.503852,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.55148,-99) , 
NN(
0, 
0, 
-1, 2.04405, 0, -1, 0.470096,-99) , 
5, 314.284, 0, 0, 0.487433,-99) , 
8, 1.49569, 1, 0, 0.498255,-99)    );
  // itree = 44
  fBoostWeights.push_back(0.0709841);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.564211,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.521604,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476538,-99) , 
7, -0.114705, 1, 0, 0.497913,-99) , 
3, 261.447, 0, 0, 0.51911,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.520418,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.489502,-99) , 
5, 139.762, 0, 0, 0.511234,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.523024,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.476786,-99) , 
2, 1.19679, 0, 0, 0.486095,-99) , 
5, 209.201, 1, 0, 0.497032,-99) , 
0, 2, 0, 0, 0.503767,-99)    );
  // itree = 45
  fBoostWeights.push_back(0.0346484);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.532439,-99) , 
NN(
NN(
0, 
0, 
-1, 257.759, 1, 1, 0.52636,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.501364,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.466724,-99) , 
3, 184.286, 0, 0, 0.493524,-99) , 
4, -0.538111, 1, 0, 0.497833,-99) , 
3, 593.129, 0, 0, 0.50017,-99)    );
  // itree = 46
  fBoostWeights.push_back(0.0433703);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0.190361, 1, 1, 0.521746,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.517441,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.489802,-99) , 
3, 154.438, 1, 0, 0.492853,-99) , 
3, 435.091, 0, 0, 0.498019,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.451974,-99) , 
7, 1.48571, 1, 0, 0.494655,-99)    );
  // itree = 47
  fBoostWeights.push_back(0.0618055);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.546792,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.523679,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.459588,-99) , 
1, 0.897511, 1, 0, 0.493837,-99) , 
3, 240.649, 0, 0, 0.513627,-99) , 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.507309,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.461778,-99) , 
7, 1.4857, 1, 0, 0.50277,-99) , 
NN(
0, 
0, 
-1, 1.19664, 0, -1, 0.466174,-99) , 
7, -0.800006, 0, 0, 0.494673,-99) , 
0, 2, 0, 0, 0.500357,-99)    );
  // itree = 48
  fBoostWeights.push_back(0.0486704);
  fForest.push_back( 
NN(
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.532197,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.518634,-99) , 
NN(
0, 
0, 
-1, 0, 1, -1, 0.491526,-99) , 
0, 2, 0, 0, 0.49813,-99) , 
4, -0.533998, 1, 0, 0.501313,-99) , 
NN(
NN(
0, 
0, 
-1, 0, 1, 1, 0.526064,-99) , 
NN(
0, 
0, 
-1, 2.94924, 1, -1, 0.454061,-99) , 
4, -0.295912, 0, 0, 0.47756,-99) , 
2, 2.69279, 1, 0, 0.497335,-99)    );
  // itree = 49
  fBoostWeights.push_back(0.0393245);
  fForest.push_back( 
NN(
NN(
0, 
0, 
-1, 1.04716, 1, 1, 0.528335,-99) , 
NN(
NN(
0, 
0, 
-1, 2.87089, 0, 1, 0.50689,-99) , 
NN(
0, 
0, 
-1, 205.949, 0, -1, 0.485696,-99) , 
7, 0.114308, 1, 0, 0.497272,-99) , 
2, 0.299203, 1, 0, 0.503778,-99)    );
   return;
};
 
// Clean up
inline void ReadBDT_50_50_07::Clear() 
{
   for (unsigned int itree=0; itree<fForest.size(); itree++) { 
      delete fForest[itree]; 
   }
}
   inline double ReadBDT_50_50_07::GetMvaValue( const std::vector<double>& inputValues ) const
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
