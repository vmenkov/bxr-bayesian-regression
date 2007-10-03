/*
3.03 -R option fix; (-R is only a switcharg)
3.04 sparse topicFeats output;
3.05 -R changed to valueArg
3.06 bbrformat fix.
3.07 modelname. 
3.08 inf as a value for -V
3.09 allow comment line in ind prior file
3.10 warning if ind prior file has no content
3.11 -R classid bug (if classid not existing, or not the final one)
3.12 0 in ind prior file
3.13 -R class info; add -R classid into DataFactory; overwrite the ver3.11 changes.
3.14 --bbrtrainformat error check
3.15 add try catch block
3.16 skip reference class when reading ind prior file setting. 
3.17 modify train data reading;
3.18 linux 64 bit change
3.19 -bmrtrainformat 
3.20 -I class 
3.21 remove standardization code: commandline option -s; avgsqunorm calc in tuneModel; standardization in train; WriteModel. 
3.22 nonzero mode features;
4.0  read ind prior file with ABS; put ind prior file reading into train() instead of in Commandline.cpp
4.1  symbolic feature reading;

todo list:
1. consider +1 and 1 different; (should be solved after the symbolic modification)
2. -R 1 and -b 1 generates different results for penalty.data; (waiting for Dave's notes)


 */

#ifdef _MSC_VER
#define VERSION "3.21 Windows Executable"
#endif 
#ifdef USE_GCC 
#define VERSION "3.21 Linux Executable"
#endif



#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdexcept>
#define  _USE_MATH_DEFINES
#include <math.h>

#include <sys/resource.h>

using namespace std;

#include "Compiler.h"
#include "logging.h"
#include "FeatSelectParamManager.h"
#include "BayesParamManager.h"
#include "ModelTypeParam.h"
#include "Data.h"
#include "PolyZO.h"
#include "IndividualPriors.h"
#include "Design.h"
#include "DataFactory.h"
#include "ModelFile.h"


#include <tclap/CmdLine.h>
using namespace TCLAP;

// TODO: These should be "somewhere" but maybe not here.
logging Log(15);

bool readFromStdin(const char* fname) {
    return 0==stricmp( "-", fname ) || 0==stricmp( "=", fname );    
}
// 


struct Config {
    HyperParamPlan* m_hyperParamPlan;

//    IndividualPriorsHolder m_individualPriorsHolder;  // ver 4.0
    DesignParameter m_designParameter;
    
//    bool m_standardize;  // v3.21
    
    double m_convergenceLimit;
    int m_iterationLimit;

    ModelType::optimizer m_optimizerType;
    
    ModelType::zoStoppingRule m_stopType;
    
    bool m_highAccuracy;
    double m_probabilityThreshold;
    string m_stringArg;
    string m_trainPlainFile;
    string m_testPlainFile;
    string m_modelWriteFileName;
    string m_modelReadFileName;
    string m_resultFileName; 
    string m_indpriorFile;
    string m_initFile;

    int m_logLevel;
    
    string m_referenceClassId; // ver 3.05
    bool m_allZeroes;
    double m_stopThreshold;
    double m_stopRatio;
    bool m_computeLogLikelihood;
    
    bool m_legacyWriteModel;
    
    bool m_binary;
    
    int m_format;  // ver 3.06

    string m_modelname; // ver 3.07

    bool m_symbolic; // ver 4.1

public:
    Config() : m_hyperParamPlan(0), // m_standardize(false), // v3.21
	       m_convergenceLimit(-1), m_iterationLimit(-1), m_optimizerType(ModelType::ZO),
	       m_stopType(ModelType::zoStoppingRule_linearScores),
	       m_highAccuracy(false),
	       m_logLevel(0),
	       m_legacyWriteModel(false),
	       m_binary(false),
	       m_symbolic(false)
	{}


};


// Obsolete args:
static ValueArg<int> nclassesArg("y","yvalues","OBSOLETE: Number of classes",false,2,"[2..]");

// Args for stopping criteria:
static ValueArg<double> stopRatioArg("", "stopCountRatio", "Ratio above stopping criterion threshold", false, 0.01, "");
static ValueArg<double> stopThresholdArg("", "stopCountThreshold", "Prob of change stopping criterion threshold", false, 0.1, "");
static ValueArg<int> stopTypeArg("", "stopcriterion", "Stopping criterion, 1=linear scores, 2=prob of change", false, 1, "[1,2]");
static ValueArg<unsigned> iterLimitArg("","iter","Max number of iterations",false,iterDefault,"integer");
static ValueArg<double>  convergeArg("e","eps","Convergence threshold",false,convergeDefault,"float");

// Args for logging:
static ValueArg<int> computeLogLikelihoodArg("", "loglikelihood", "Display log likelihood, 1=yes, 0=no", false, 0, "");
static ValueArg <int>  logArg("l","log","Log verbosity level",false,0, "[0..2]");

// Args for data standardization and normalization:
static SwitchArg cosnormArg("c","cosnorm","Cosine normalization",false); 
// static SwitchArg standardizeArg("s","standardize","Standardize variables",false);  // v3.21

// Args for parameter searching 
//back-compatibility only -->
static ValueArg <string>  searchArg("S","search",
    "DEPRECATED Search for hyperparameter value",false,"","list of floats, comma-separated, no spaces"); 
static ValueArg <double>  hypArg("H","hyperparameter",
    "DEPRECATED Hyperparameter, depends on the type of prior", false,0,"float"); 
//<--back-compatibility only


// Args configuring output file:
static SwitchArg resScoreArg("","rscore","Scores, not probabilities, in the result file", false);
static ValueArg<string> resfileArg("r","resultfile","Results file",false,"","resultfile");
static ValueArg<int> formatArg("","outputformat",
			       "format: \n\t 1-bbr \n\t 2-bmr \n\t 3-sparse", false, 2, "[1,3]"); 
static ValueArg<string> modelnameArg("","modelname","name of the model",false,"","model name");    // ver 3.07
static SwitchArg symbolicArg("","symbolic","string feature or class labels",false);  // ver 4.1

// Args configuring optimizer:
static ValueArg <int> optArg("o","opt",
    "Optimizer: \n\t 1-ZO \n\t 2-quasi-Newton, smoothed penalty \n\t 3-quasi-Newton, double coordinate",
    false,1,"[1,3]");
static SwitchArg highAccuracyArg("","accurate","High accuracy mode",false);


// Binary versus non:
static ValueArg<int> binaryArg("b", "binary","1=Binary, 0=Poly", false, 0, "[0,1]");


// Bayes args:
static ValueArg <int>  priorArg("p","prior","Type of prior, 1-Laplace 2-Gaussian",false,2,"[1,2]"); 
static ValueArg<string> indPriorsArg("I", "individualPriorsFile", "Individual Priors file", false, "", "indPriorsFile"); 
static ValueArg<string> initArg("", "init", "init file", false, "", "initfile");
static SwitchArg allZeroArg("z","zerovars","Exclude all-zero per class", false);
//static SwitchArg refClassArg("R","refClass","Use Reference Class", false);
static ValueArg <string>  refClassArg("R","refClass","Use Reference Class", false, "", "reference class");  // ver 3.05
static SwitchArg pvarSearchArg("","autosearch","Auto search for prior variance, no grid required", false);
static SwitchArg negOnlyArg("","neg","Negative only model coefficients", false);
static SwitchArg posOnlyArg("","pos","Positive only model coefficients", false);


// Args for cross validation:
static SwitchArg errBarArg("","errbar","Error bar rule for cross-validation", false); 
static ValueArg <string>  cvArg("C","cv","Cross-validation",false,"10,10","#folds[,#runs]"); 
static ValueArg <string>  priorVarArg("V","variance",
    "Prior variance values; if more than one, cross-validation will be used",false,"",
    "number[,number]*"); 


// Args specifying test data file, train data file, and model data file:
static ValueArg<string> testfileArg("", "testfile", "Test file", false, "", "");
static UnlabeledValueArg<string>  datafileArg("x","Data file; '-' signifies standard input","","data file");
static UnlabeledValueArg<string>  modelfileArg("modelfile","Model file","","model file");

//
static ValueArg<double> probThrArg("","probthres","Probability threshold",false,0.5,"0<=p<=1");

static ValueArg<int> legacyWriteModelArg("", "legacywm", "Legacy write model", false, 0, "[0,1]");


static Arg* args[] = {&nclassesArg,
		      &legacyWriteModelArg,
		      &negOnlyArg,
		      &posOnlyArg,
		      &probThrArg,
		      &highAccuracyArg,
		      &pvarSearchArg,
		      &stopRatioArg,
		      &stopThresholdArg,
		      &stopTypeArg,
		      &iterLimitArg,
		      &convergeArg,
		      &computeLogLikelihoodArg,
		      &logArg,
		      &cosnormArg,
		      // &standardizeArg,  // v3.21
		      &searchArg,
		      &hypArg,
		      &resScoreArg,
		      &resfileArg,
		      &formatArg, // ver 3.06
		      &modelnameArg, // ver 3.07
		      &symbolicArg, // ver 4.1
		      &optArg,
		      &priorArg,
		      &indPriorsArg,
		      &initArg,
		      &allZeroArg,
		      &refClassArg,
		      &errBarArg,
		      &cvArg,
		      &priorVarArg,
		      &testfileArg,
		      &datafileArg,
		      &modelfileArg,
		      &binaryArg,
		      0
};

static bool parseArgs(int argc, char** argv) {
//    try {  
	CmdLine cmd( "Bayesian Multinomial Logistic Regression - Training", ' ', VERSION );
	for (Arg** a = args; *a; a++) {
	    cmd.add(*a);
	}
	
	// Parse the args.
        cmd.parse( argc, argv );
	/*   } catch( std::exception& e){
	Log(0)<<std::endl<<"Exception parsing parameters: "<<e.what();
        cerr<<std::endl<<"Exception parsing parameters: "<<e.what();
        return false;
    } catch(...){
        Log(0)<<std::endl<<"***Unrecognized exception";
        cerr<<std::endl<<"***Unrecognized exception";
        return false;
    }
	*/
    return true;
    
}

static ModelType::optimizer getOptimizerType() {
    if (optArg.getValue() == 2) {
	return ModelType::QuasiNewtonSmooth;
    } else if (optArg.getValue() == 3) {
	return ModelType::QuasiNewtonDoubleCoord;
    } else {
	return ModelType::ZO;
    }
}

static HyperParamPlan getHyperParamPlan() {
    
    HyperParamPlan plan;
    
    enum PriorType prior = PriorType(priorArg.getValue());
    if( prior!=1 && prior!=2 )
        throw runtime_error("Illegal prior type; should be 1-Laplace or 2-Gaussian");
    int skew = posOnlyArg.getValue() ? 1 : negOnlyArg.getValue() ? -1 : 0;
    
    
    if( priorVarArg.isSet() ) { // new mode
        if( pvarSearchArg.isSet() ) throw runtime_error("Incompatible arguments: auto search and grid");
        plan = HyperParamPlan( prior, skew, priorVarArg.getValue(), cvArg.getValue(),
			       HyperParamPlan::AsVar, errBarArg.getValue() );
    } else if( searchArg.isSet() ) { //back-compatibility
        plan = HyperParamPlan( prior, skew, searchArg.getValue(), cvArg.getValue(),
			       HyperParamPlan::Native, errBarArg.getValue() );
    } else if( hypArg.isSet() ) { //fixed hyperpar - back-compatibility
        plan = HyperParamPlan( prior, hypArg.getValue(), skew );
    } else if( pvarSearchArg.isSet() ) { // auto search, no grid
        plan = HyperParamPlan( prior, skew, cvArg.getValue() );
    }
    else //auto-select hyperpar
        plan = HyperParamPlan( prior, skew );
    
    return plan;
}


static ModelType* createModelType(Config* config) {
    
    ModelType* modelType = new ModelType(ModelType::logistic,
					 config->m_optimizerType,
					 config->m_probabilityThreshold,
					 config->m_convergenceLimit,
					 config->m_iterationLimit, 
					 config->m_highAccuracy,
					 config->m_referenceClassId, // ver 3.05
					 config->m_allZeroes,
					 config->m_modelname, // ver 3.07
					 config->m_stringArg, 
					 config->m_stopType,
					 config->m_stopThreshold,
					 config->m_stopRatio,
					 config->m_computeLogLikelihood,
					 config->m_binary
        );
    Log(2)<<endl<<*modelType;

    return modelType;
}


static DataFactory* createDataFactory(Config* config) {

    DataFactory* dataFactory = new DataFactory(config->m_symbolic);
    
    dataFactory->setTestAndTrainFileSpec(
	config->m_trainPlainFile,
	config->m_testPlainFile,
	config->m_indpriorFile,
	config->m_referenceClassId
	);

    dataFactory->readFiles();

    // ver 3.14; --bbrtrainformat error check
    if(config->m_format == 1) {
	// --bbrtrainformat only allow with two classes
	const RowSetMetadata& metadata = dataFactory->getRowSetMetadata();
	if (metadata.getClassCount()!=2)
	    throw logic_error("--bbrtrainformat is allowed only with two classes.\n");
	// the two classes could only be -1 and +1
        #ifdef SYMBOLIC
	if (metadata.getClassId(0)!="-1" || (metadata.getClassId(1)!="1"&&metadata.getClassId(1)!="+1") )
        #else
	if (metadata.getClassId(0)!=-1 || metadata.getClassId(1)!=1)    
        #endif
	    throw logic_error("--bbrtrainformat is allowed only when the two classes are -1 and +1.\n");
    }

    return dataFactory;
}

static WriteModel* createWriteModel(Config* config, const RowSetMetadata& rowSetMeta) {
    ofstream* resultFileStream = new ofstream(config->m_resultFileName.c_str());

    WriteModel* writeModel = new WriteModel(config->m_modelWriteFileName, config->m_format,  // ver 3.06
					    rowSetMeta, 
					    resultFileStream, config->m_legacyWriteModel);

    return writeModel;
}

static void run(DataFactory& dataFactory, ModelType& modelType, WriteModel& writeModel, Config *config) {

    LRModel model(*config->m_hyperParamPlan, 
		  dataFactory.getIndPriorHolder(), // ver 4.0
		  config->m_designParameter, 
		  modelType,
		  dataFactory.getRowSetMetadata(),
		  config->m_initFile );

    model.train(*dataFactory.getTrainData(), writeModel, dataFactory.getTestData());

    bool test = (dataFactory.getTestData() != 0);
    if (test)  {
	model.test(*dataFactory.getTestData(), writeModel);
    }

    writeModel.closeResultsFile();
}

static void logHeader(int argc, char** argv) {
    Log(0)<<endl<<"Bayesian Multinomial Logistic Regression - Training \tVer. "<<VERSION;
    Log(2)<<"\nCommand line: ";
    for( int i=0; i<argc; i++ ) {
        Log(2)<<" "<<argv[i];
    }

    Log(2)<<"\nLog Level: "<<Log.level()-5;
}


static Config* configureLogger() {
    Config* config = new Config();
    config->m_logLevel = logArg.getValue()+5;
    Log.setLevel(config->m_logLevel);
    return config;
}

static void readConfig(Config* config) {
    
    config->m_hyperParamPlan = new HyperParamPlan(getHyperParamPlan());
    // config->m_standardize = standardizeArg.getValue();   // v3.21
    config->m_convergenceLimit = convergeArg.isSet() ? convergeArg.getValue() : convergeDefault;
    config->m_iterationLimit = iterLimitArg.isSet() ? iterLimitArg.getValue() : iterDefault;
    config->m_optimizerType = getOptimizerType();
    config->m_stopType = stopTypeArg.getValue() == 1 ? ModelType::zoStoppingRule_linearScores : ModelType::zoStoppingRule_changeProb;
    config->m_stopRatio = stopRatioArg.getValue();
    config->m_stopThreshold = stopThresholdArg.getValue();
    config->m_referenceClassId = refClassArg.isSet() ? refClassArg.getValue() : "" ;  // ver 3.05

    config->m_modelname = modelnameArg.isSet() ? modelnameArg.getValue() : datafileArg.getValue()+".model"; // ver 3.07

    config->m_symbolic = symbolicArg.isSet() ? true : false;  // ver 4.1

    config->m_allZeroes = allZeroArg.getValue();
    
    config->m_trainPlainFile = datafileArg.getValue();
    config->m_testPlainFile = testfileArg.getValue();
    config->m_modelWriteFileName = modelfileArg.getValue();  
    
    config->m_modelReadFileName = ""; // TODO This should be "" when training, set when testing to: modelfileArg.getValue();  
    config->m_resultFileName = resfileArg.getValue();
    
    config->m_designParameter = DesignParameter((enum DesignType)designPlain);
    
    config->m_indpriorFile = indPriorsArg.isSet() ? indPriorsArg.getValue() : "";

    config->m_initFile = initArg.isSet() ? initArg.getValue() : "";

    /* // ver 4.0
    if( indPriorsArg.isSet() ) {
	config->m_individualPriorsHolder = IndividualPriorsHolder( indPriorsArg.getValue(), IND_PRIORS_MODE_RELATIVE );
    }
    */
    
    config->m_binary = binaryArg.getValue() == 1 ? true : false;
    config->m_highAccuracy = highAccuracyArg.getValue();
    config->m_probabilityThreshold = probThrArg.getValue();
    
    config->m_legacyWriteModel = legacyWriteModelArg.getValue() == 1 ? true : false;

    // Does this do anything?
    config->m_stringArg = string("");

    config->m_format = formatArg.getValue();  // ver 3.06

}

void logConfig(Config* config) {
    Log(2)<<"\nData file for training: "
	  <<( readFromStdin(config->m_trainPlainFile.c_str()) ? "stdin" : config->m_trainPlainFile);
    Log(2)<<"\nData file for testing: "<<config->m_testPlainFile.c_str();
    Log(2)<<"\nWrite Model file: "<<config->m_modelWriteFileName;
    Log(2)<<endl; // <<config->m_individualPriorsHolder; // ver 4.0
}

void postSanityCheck(Config* config){
    // ver 3.14; if --bbrtrainformat specified, and -R is only allowed with value +1;
    if(config->m_format == 1){
	if(config->m_referenceClassId!="1")
	    throw logic_error("--bbrtrainformat is allowed only when '-R -1' is specified, i.e. class -1 is forced to be the reference class.\n");
    }
}

int main(int argc, char** argv) {
    try{
	bool ok = parseArgs(argc, argv);
	
	if (!ok) {
	    exit(1);
	}
	

	Config* config = configureLogger();
	logHeader(argc, argv);

	readConfig(config);
	logConfig(config);
	postSanityCheck(config);  // ver 3.14
	
	// ver 3.13; add referenceclassid into Y vector;
	DataFactory* dataFactory = createDataFactory(config);

	ModelType* modelType = createModelType(config);
	
	WriteModel* writeModel = createWriteModel(config, dataFactory->getRowSetMetadata());
	run(*dataFactory, *modelType, *writeModel, config);
	
	struct rusage rsg;
	if(getrusage(RUSAGE_SELF, &rsg)==0)
	    cout<<"max memory used:"<<rsg.ru_maxrss<<endl;

	delete modelType;
	delete writeModel;
	delete dataFactory;
    }
    catch (ArgException e)  {
	cerr << "***Command line error: " << e.error() 
	     << " for arg " << e.argId() << endl;
	return 1;
    }
    catch(std::exception& e){
	cerr<<std::endl<<"***Exception: "<<e.what();
	return 1;
    }
    catch(...){
	cerr<<std::endl<<"***Unrecognized exception";
	return 1;
    }

}
