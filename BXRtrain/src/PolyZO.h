#ifndef POLY_H
#define POLY_H

#include "Matrix.h"
#include "Data.h"
#include "BayesParamManager.h"
#include "ModelTypeParam.h"
#include "Design.h"
#include "ParamMatrix.h"
#include "IndividualPriors.h"

// Forward Def
class RowSetStats;

class LRModel {
    string m_initFName;

    bool m_bTrained;
    const class ModelType& m_modelType;

    BayesParameter m_bayesParam;
    ParamMatrix m_beta;
    
    ParamMatrix m_priorMode; // ver 4.0 ; to move ind prior reading to ZOLR
    ParamMatrix m_priorSkew; // ver 4.0
    ParamMatrix m_penalty; // ver 4.0

    double m_threshold;
    
    vector<double> m_lastDiff;
    unsigned m_dScoreAboveThreshold;
	
    double m_convergenceScore;
    unsigned m_trainingSetAccuracy;
    
    const class DesignParameter& m_designParameter;
    const class HyperParamPlan &m_hyperParamPlan;

    // const class IndividualPriorsHolder& m_individualPriorsHolder;  // ver 4.0
    class IndividualPriorsHolder* m_individualPriorsHolder; // ver 4.0
    const class RowSetMetadata& m_metadata;

 public:
    LRModel(const class HyperParamPlan& hyperParamPlan_,
	    class IndividualPriorsHolder* const individualPriorsHolder_,  // ver 4.0
	    const class DesignParameter& designParameter_,
	    const class ModelType& modelType_,
	    const class RowSetMetadata& metadata_,
	    const string initfilename_ )
	: m_bTrained(false), m_hyperParamPlan(hyperParamPlan_),
	m_modelType(modelType_), m_designParameter(designParameter_), 
        m_individualPriorsHolder(individualPriorsHolder_), 
	m_metadata(metadata_), m_initFName(initfilename_)
	/*m_pDesign(0), */
	// m_holdoutRowSet(0), m_holdoutRowSetIterator(0) 
	{	
	    /*
	    if (m_holdoutRowSet != 0) {
		m_holdoutRowSetIterator = m_holdoutRowSet->iterator();
	    }
	    */

	}
    
    ~LRModel() { /*delete m_holdoutRowSetIterator;*/ }

    void train(MemRowSet& trainData, class WriteModel& modelFile, IRowSet* heldout);
    void restore( class ReadModel& modelFile, const RowSetMetadata& metadata);
    void test( PlainYRowSet & TestRowSet, WriteModel& modelFile);

private:
    void tuneModel(MemRowSetIterator & drs,
		   const vector<int>& wordIndex,
		   const HyperParamPlan& hyperParamPlan,
		   const class FixedParams& fixedParams,
		   RowSetStats& stats );
    
    void invokeZOLR(
	MemRowSetIterator& rowSet,
	const class InvData& invData, //input - generates design data
	const class BayesParameter& bayesParam,
	const class FixedParams& fixedParams,  //input
	ParamMatrix & beta,  // i/0: init/result values of model coeffs, 
	bool binary);
    

    void ZOLR(
	MemRowSetIterator& rowSet,
        const class InvData& invData, //input - generates design data
	const class BayesParameter& bayesParam,
        const class FixedParams& fixedParams,  //input
        ParamMatrix& beta  // i/0: init/result values of model coeffs, 
	);
    
    void LRModel::ZOLRBinary(
	MemRowSetIterator& rowSet,
	const InvData& invData, //input - generates design data
	const BayesParameter& bayesParam, //input
	//const vector<double>& priorMode,  //input
	//const vector<double>& priorScale,  //input
	//const vector<double>& priorSkew,  //input
	vector<bool>& fixedParams,
	vector<double>& beta
	);
    

    bool stopTest(MemRowSetIterator& rowset, 
		int nRows, 
		  int nClasses, 
		  const ModelType& modelType, 
		  const ParamMatrix& w, 
		  double sum_abs_dr, 
		  const vector< vector<double> >& wTX, 
		  vector<double>& lastDiff);
    
    vector<double>* hyperParamLoop(
	MemRowSetIterator& drs,
	unsigned dim, unsigned classes,
	RowSetIterator& drsTest,
	ModelType modelType,
	const HyperParamPlan& hpPlan,
	const class FixedParams& fixedParams
	);
    
    
    vector<int> buildWordIndexForSample(RowSetIterator& drs, unsigned dim) const;
    pair<int,int> countCorrect(RowSetIterator& rowset, const ParamMatrix& beta) const;
    void iterationLogger(MemRowSetIterator& rowSet, 
			 const BayesParameter& bayesParam, 
			 const ParamMatrix& beta, 
			 const ParamMatrix & penalty, unsigned itr) const;
    double AvgSquNorm( RowSetStats& stats ) const;
};

class InvData {

    vector< pair<unsigned,double> > data;

    unsigned nrows;
    unsigned m_c;
    vector<unsigned> m_y;
    const vector<int> m_wordIndex;
public:

    void print() const {
	for(int i=0; i<data.size(); i++) {
	    cout<<"["<<i<<"]:"<<data.at(i).first<<":"<<data.at(i).second<<endl;
	}
    }

    unsigned n() const { return nrows; }
    unsigned c() const { return m_c; } // #classes in 'y'
    unsigned y( unsigned row ) const { return m_y[row]; }

    typedef pair<unsigned, double> pair_type;
    typedef vector< pair<unsigned, double> > pair_vector_type;
    typedef pair<pair_vector_type::const_iterator, pair_vector_type::const_iterator > pair_range_type;

    InvData(RowSetIterator& it, unsigned c, unsigned dim, const vector<int>& wordIndex ) : m_wordIndex(wordIndex), m_c(c)
    {

	Log(3)<<std::endl<<"Starting InvData initialization 1, - Time "<<Log.time();

	try{	
	    data = vector< pair<unsigned, double> >(wordIndex.at(wordIndex.size()-1));
	    //data.reserve(wordIndex.at(wordIndex.size()-1));
	}
	catch(std::exception & e){
	    throw logic_error("Insufficient memory to build inverted index to training data.  Aborting execution.");
	}

	vector<unsigned> wordCounter = vector<unsigned>(dim, 0);
	Log(20)<<endl; cout<<data.size()<<","<<data.capacity()<<","<<data.max_size()<<endl;

	Log(3)<<std::endl<<"Starting InvData initialization 2, - Time "<<Log.time();
	it.rewind();
        for( nrows=0; it.next(); nrows++ ) {
            const SparseVector& x=it.xsparse();
	    for( SparseVector::const_iterator xitr=x.begin(); xitr!=x.end(); xitr++ ) {
		int index = xitr->first;
		data[wordIndex.at(index) + wordCounter.at(index)] = pair<unsigned, double>(nrows, xitr->second);
		wordCounter.at(index)++;
	    }
            m_y.push_back( it.y() );
        }
	
	Log(3)<<"Starting InvData initialization 3, - Time "<<Log.time()<<endl;
    }
    
    const pair_range_type getRange(unsigned var) const {
	pair_vector_type::const_iterator i1 = data.begin();
	return pair<pair_vector_type::const_iterator, pair_vector_type::const_iterator >
	    ( data.begin() + m_wordIndex.at(var), 
	      data.begin() + m_wordIndex.at(var+1)
		);
    }
    
};

class WakeUpPlan {
    vector<bool> plan;        
 public:
    bool operator()(unsigned int step ) {
	if( step<plan.size() ) return plan.at(step);
	else return 0==step%100;
    }
    WakeUpPlan(unsigned size=1000) : plan( size+1, false ) { //ctor
	for(unsigned i=0;i*i*i<size;i++) plan[i*i*i] = true;     }
};



#endif // POLY_H

/*
    Copyright 2005, Rutgers University, New Brunswick, NJ.

    All Rights Reserved

    Permission to use, copy, and modify this software and its documentation for any purpose 
    other than its incorporation into a commercial product is hereby granted without fee, 
    provided that the above copyright notice appears in all copies and that both that 
    copyright notice and this permission notice appear in supporting documentation, and that 
    the names of Rutgers University, DIMACS, and the authors not be used in advertising or 
    publicity pertaining to distribution of the software without specific, written prior 
    permission.

    RUTGERS UNIVERSITY, DIMACS, AND THE AUTHORS DISCLAIM ALL WARRANTIES WITH REGARD TO 
    THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
    ANY PARTICULAR PURPOSE. IN NO EVENT SHALL RUTGERS UNIVERSITY, DIMACS, OR THE AUTHORS 
    BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER 
    RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, 
    NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR 
    PERFORMANCE OF THIS SOFTWARE.
*/
