// 3.21 remove -s related code in ModelType

#ifndef _MODEL_TYPE_PARAMETER_HPP
#define _MODEL_TYPE_PARAMETER_HPP

#include <ostream>
#include <string>
#include <stdexcept>

#define convergeDefault 0.001
#define iterDefault 10000

class ModelType {
    friend class ReadModel;
 public:
    //types
    enum link { probit=0, logistic=1, SVM=2, linkUndef };
    enum optimizer { EM=0, ZO=1, SVMlight=2, QuasiNewtonSmooth=3, QuasiNewtonDoubleCoord=4, optUndef };
    enum thrtune { thrNo=0, thrSumErr=1, thrT11U=2, thrF1=3, thrBER=4, thrT13U=5,  thrUndef };
    enum zoStoppingRule { zoStoppingRule_linearScores=0, zoStoppingRule_changeProb=1 };

    ModelType() : m_link(logistic), m_opt(ZO), m_thr(thrSumErr), m_tuneEst(false), 
	m_thrConverge(convergeDefault), m_iterLimit(iterDefault), m_highAccuracy(false),
	m_zoStoppingRule(zoStoppingRule_linearScores),
	m_stoppingRuleRScoreThreshold(0),
	m_stoppingRuleRowLimit(0)
	{}

	ModelType( enum link link_, 
		   enum optimizer opt_, 
		   double probThr_, 
		   double thrConverge, 
		   unsigned iterLimit, 
		   bool highAccuracy, 
		   string referenceClassId_, // ver 3.05
		   bool allZeroClassMode_,
		   string modelname_ ="", // ver 3.07
		   string strParam_="",
		   enum zoStoppingRule zoStoppingRule_ = zoStoppingRule_linearScores,
		   double stoppingRuleRScoreThreshold_ = 0,
		   double stoppingRuleRowLimit_ = 0,
		   bool computeLogLikelihood_ = false,
		   bool binary_ = false
	    )
	    : m_link(link_), m_opt(opt_), 
	    m_thr(thrNo), m_probThr(probThr_),
	    m_tuneEst(false), 
	    m_thrConverge(thrConverge), m_iterLimit(iterLimit), m_highAccuracy(highAccuracy),
	    m_referenceClassId(referenceClassId_), // ver 3.05
	    m_modelname(modelname_), // ver 3.07
	    m_allZeroClassMode(allZeroClassMode_),
	    strParam(strParam_),
	    m_zoStoppingRule(zoStoppingRule_),
	    m_stoppingRuleRScoreThreshold(stoppingRuleRScoreThreshold_),
	    m_stoppingRuleRowLimit(stoppingRuleRowLimit_),
	    m_computeLogLikelihood(computeLogLikelihood_),
	    m_binary(binary_)
	    {}

	enum link getLinkType() const { return m_link; }
	enum optimizer getOptimizerType() const { return m_opt; }
	enum thrtune getTuneThreshold() const { return m_thr; }
	double getProbThreshold() const { return  m_probThr; }
	bool isBinary() const { return m_binary; }
	bool isTuneEst() const { return m_tuneEst; }					// ?
	// bool isStandardize() const { return m_standardize; }			// v3.21
	string getStringParam() const { return strParam; }
	double getConvergenceThreshold() const { return m_thrConverge; }			// ?
	bool isHighAccuracy() const { return m_highAccuracy; } //bbr v2.04
	unsigned getIterLimit() const { return m_iterLimit; } //bbr v2.04
	string getReferenceClassId() const { return m_referenceClassId; }  // ver 3.05 
	string getModelName() const { return m_modelname; } // ver 3.07
	bool isAllZeroClassMode() const { return m_allZeroClassMode; }
	bool isComputeLogLikelihood() const { return m_computeLogLikelihood; }
	
	enum zoStoppingRule getZoStoppingRule() const {return m_zoStoppingRule; }
	double getZoStoppingRuleRScoreThreshold() const { return m_stoppingRuleRScoreThreshold; }
	double getZoStoppingRuleRowLimit() const { return m_stoppingRuleRowLimit; }
	
	//input
private:
	enum link m_link;
	enum optimizer m_opt;
	enum thrtune m_thr;
	double m_probThr;
	/** If m_tuneEst==true, tuning is done by using estimated
	    probabilities on the (large, unlabeled) Training Population,
	    rather than by using actual labels on the (small, labeled)
	    training set */
	bool m_tuneEst; 
	// bool m_standardize;  // v3.21
	string strParam;
	double m_thrConverge;
	bool m_highAccuracy;  //bbr v2.04
	unsigned m_iterLimit; //bbr v2.04
	string m_referenceClassId; // ver 3.05
	string m_modelname; //ver 3.07
	bool m_allZeroClassMode;
	
	// Whether to compute and display log-likelihood each time
	bool m_computeLogLikelihood;


	// Valid only when using the ZO optimizer:
	enum zoStoppingRule m_zoStoppingRule;

	// Valid only when using the ZO optimizer and the changeProb stopping rule:
	double m_stoppingRuleRScoreThreshold;
	double m_stoppingRuleRowLimit;

	//
	bool m_binary;
};

inline std::ostream& operator<<( std::ostream& o, const ModelType& mt ) {
    o << "Model - link function: "<<( 
        mt.getLinkType()==ModelType::probit ? "Probit "
        : mt.getLinkType()==ModelType::logistic ? "Logistic " 
        : mt.getLinkType()==ModelType::SVM ? "SVMlight " 
        : "<undefined> " )
      << "\tOptimizer: "<<( 
        mt.getOptimizerType()==ModelType::EM ? "EM"
        : mt.getOptimizerType()==ModelType::ZO ? "ZO" 
        : mt.getOptimizerType()==ModelType::SVMlight ? "SVMlight" 
        : mt.getOptimizerType()==ModelType::QuasiNewtonSmooth ? "quasi-Newton, smoothed penalty" 
        : mt.getOptimizerType()==ModelType::QuasiNewtonDoubleCoord ? "quasi-Newton, double coordinate" 
        : "<undefined> " );
    
	o << (mt.isTuneEst() ? ", on TP estimates" : "" )
	    // << "\nStandardize: "<<( mt.isStandardize() ? "Yes" : "No" )  // v3.21
	  << "\nConvergence threshold: "<<mt.getConvergenceThreshold()<<"\tIterations limit: "<<mt.getIterLimit()
	  <<( mt.isHighAccuracy()? "\tHigh accuracy mode" : "" );
	
	if( mt.getReferenceClassId()!="" ) {
	    o<<"Using Reference Class"<<endl;
	}
	
	if( mt.isAllZeroClassMode() ) {
        o<<"\nAllZeroClassMode: Variable/class const zero gets zero beta";
	}
	
	if( mt.getStringParam().size()>0 ) {
        o<<"\nString parameter: '"<<mt.getStringParam()<<"'";
	}
	
	if (mt.getOptimizerType()==ModelType::ZO) {
		o <<"\nStopping rule: ";
		if (mt.getZoStoppingRule()== ModelType::zoStoppingRule_linearScores) {
			o << "linear scores";
		} else if (mt.getZoStoppingRule()== ModelType::zoStoppingRule_changeProb) {
			o << "change probabilities (" << mt.getZoStoppingRuleRScoreThreshold() << "/" << mt.getZoStoppingRuleRowLimit() << ")";
		}
	}
    return o;
}

#endif //_MODEL_TYPE_PARAMETER_HPP
