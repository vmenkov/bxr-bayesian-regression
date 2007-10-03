#define  _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <limits>
#include <iomanip>

#include "logging.h"
#include "Compiler.h"
#include "Matrix.h"
#include "BayesParamManager.h"
#include "IndividualPriors.h"
#include "Design.h"
#include "PolyZO.h"
#include "ModelTypeParam.h"
#include "RowSetStats.h"
#include "ModelFile.h"
#include "Likelihood.h"
#include "QuasiNewton.h"
#include "OutputFormat.h"
#include "Data.h"

using namespace std;



/*
 * Zhang & Oles regularized logistic regression
 * applied to polytomous here
 *
 * 1.01   Feb 23 05     error bar rule for hyperparameter cv (option --errbar)
 *                      'wakeUpPlan' does not depend on iter.limit
 * 1.05   Jun 20, 05    one std error rule: fixed bug with stderr denom calc //PolyZO.cpp
 * 1.60   Oct 24, 05    sparse beta in model file
 *                      Individual priors
 *                      fixed log.prior with infinite prior var
 * 3.01   Jan 31, 06    fixed bug: Gaussian penalty should be 1/(2*var), not 1/var
 */


inline double OneCoordStep(
    double dw_numerator,
    double dw_denominator,
    const InvData::pair_range_type range,
    const InvData& invData, //for 'y' only
    unsigned k, //class
    const vector< vector<double> >& linearScores, //added for exp overflow
    const vector<double>& classLinearScores,
    const vector<double>& classExponentialScores,
    const vector<double>& classSumExponentialScores,
    double trustregion,
    bool highAccuracy
    )
{

    vector<double> dnums; vector<double> ddenoms;
    if(highAccuracy) {
        dnums.push_back(dw_numerator);
        ddenoms.push_back(dw_denominator);
    }

    for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ ) {
        unsigned i = ix->first; //case #
        double x = ix->second; //value
        if( 0.0==x ) continue; //just in case - shouldn't ever hit this
        unsigned target = invData.y(i)==k ? 1 : 0; // y_ik in MGLF05

        double D;
        double delta = trustregion * fabs(x);
        bool numbranch = false;
	double dwNumeratorChange = 0.0; // DELTAv_kj in MGLF05
        //Log(10)<<"\ndelta = trustregion * fabs(x) "<<delta<<" "<<trustregion<<" "<<x;
        if( finite(classExponentialScores[i]) && finite(classSumExponentialScores[i]) 
	    && classSumExponentialScores[i]>0 ) 
	{ 
            const double exp_wkTX_i = classExponentialScores[i]; // exp(r_ik) in MGLF05
            const double K = classSumExponentialScores[i] - exp_wkTX_i; // E_ik in MGLF05

            const double phat_ik = exp_wkTX_i / classSumExponentialScores[i];

            // contribution to DELTAv_kj from example i 
            dwNumeratorChange += target ? (x - x * phat_ik) : (-x * phat_ik); 

            const double e_delta = exp(delta); /// 1+delta; 
            const double e_rlo = exp_wkTX_i / e_delta;
            const double e_rhi = exp_wkTX_i * e_delta;

	    // if exp(r_ik - delta) <= K <= exp(r_ik+delta), F(B, x_i; delta) = 2; MGLF05
	    // D = F(B, x_i, DELTAv_kj) + 2; in MGLF05
	    if( e_rlo<=K && K<=e_rhi )
                D = 4;
            else
                D = min( e_rlo/K+K/e_rlo, e_rhi/K+K/e_rhi ) + 2;
        }
        else
        {
	    //std::cout << "here";
            numbranch = true;
            double invPminus1 = 0.0;
            double qlo=0, qhi=0; // K/e_rlo, K/e_rhi
            for( unsigned kk=0; kk<invData.c(); kk++ ) {
                if( kk!=k ) {
                    invPminus1 += exp( linearScores[kk][i] - linearScores[k][i] );
                    qlo += exp( linearScores[kk][i] - (linearScores[k][i]-delta) ); //wkTX[i]
                    qhi += exp( linearScores[kk][i] - (linearScores[k][i]+delta) ); //wkTX[i]
                }
            }
            if( 1<=qlo && qhi<=1 ) //e_rlo<=K && K<=e_rhi )
                D = 4;
            else
                D = min( qlo+1/qlo, qhi+1/qhi ) + 2;
            dwNumeratorChange = (target ? x : 0) - x/(1+invPminus1);  //ver 0.24
        }

        //Log(10)<<setw(23)<<( x * ( target ? invPminus1 : -1 ) / (1+invPminus1) )<<setw(23)<< x * x / D<<setw(23)<<dw_numerator<<setw(23)<<dw_denominator;
        //Log(10)<<"\ninvPhat "<<setw(23)<< invPminus1;
        //Log(7)<<setw(14)<< x *( target - p_hat )<<setw(14)<< x * x / D<<setw(14)<<dw_numerator*(ii+1)<<setw(14)<<dw_denominator*(ii+1);

        if( highAccuracy ){
            dnums.push_back(dwNumeratorChange);
            ddenoms.push_back(x * x / D);
        }else{
            dw_numerator += dwNumeratorChange;
            dw_denominator += x * x / D;
        }

    }
    
    // Log(3) << "\nOneCoordStep: " << k << ' ' << Log.time()<<endl;

    //update w
    double dw = 0;
    if( highAccuracy ){  //ver 2.04
        sort( dnums.begin(), dnums.end() );
        double dnumerator=0;
        for( vector<double>::const_iterator itr=dnums.begin(); itr!=dnums.end(); itr++ )
            dnumerator += *itr;

        sort( ddenoms.begin(), ddenoms.end() );
        double ddenominator=0;
        for( vector<double>::const_iterator itr=ddenoms.begin(); itr!=ddenoms.end(); itr++ )
            ddenominator += *itr;
	
        dw = 0.0==dnumerator ? 0.0 : dnumerator / ddenominator;  // ver 3.11; compare with 0.0 first
    }
    else {
        dw = 0.0==dw_numerator ? 0.0 : dw_numerator / dw_denominator;
    }

    dw = min( max(dw,-trustregion), trustregion );
    if( !finite(dw) )
        Log(1)<<"\nErrDW k dw_numerator dw_denominator dw "<<k<<" "<<dw_numerator<<" "<<dw_denominator<<" "<<dw;
    return dw;
}

inline double OneCoordStepBinary(
    double dw_numerator,
    double dw_denominator,
    const InvData::pair_range_type range,
    const class InvData& invData, //for 'y' only
    const Vector& wTX,
    double trustregion,
    bool highAccuracy
    )
{ 
    vector<double> dnums; vector<double> ddenoms;
    if(highAccuracy) {
        dnums.push_back(dw_numerator);
        ddenoms.push_back(dw_denominator);
    }

    for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ )
    {
        unsigned i = ix->first; //case #
        double x = ix->second; //value
        if( 0.0==x ) continue; //just in case - shouldn't ever hit this
        int target = invData.y(i) ? 1 : -1;
        double r = 1==target ? wTX[i] : -wTX[i]; //wTX[i] * target;
        //double q = 1 / (1 + exp(r));  //p = 1 / (1 + exp(- r))
        //try this approximation
        double a = fabs(r);
        double b = fabs(trustregion * x);
        double F; // = a<=b ? 0.25 : 1.0/( 2.0 + exp(a-b) + exp(b-a) );
        if( a<=b ) F = 0.25;
        else {
            double e = exp(a-b);
            F = 1.0/( 2.0 + e + 1.0/e );
        }
        /*/original Z&O
        double F = min( 0.25, exp(fabs(trustregion * x)) * q * (1 - q)); //* p * (1 - p))
            */
        if( highAccuracy ){
            dnums.push_back( x*target/(1 + exp(r)) );  //ver 2.04  accurate summation
            ddenoms.push_back( F*x*x );  //ver 2.04  accurate summation
        }else{
            dw_numerator += x * target / (1 + exp(r)) ;  //q * x * target; //(1 - p) * x * target
            dw_denominator += F * x * x;
        }
	}

    //update w
    double dw;
    if( highAccuracy ){  //ver 2.04
        sort( dnums.begin(), dnums.end() );
        double dnumerator=0;
        for( vector<double>::const_iterator itr=dnums.begin(); itr!=dnums.end(); itr++ )
            dnumerator += *itr;

        sort( ddenoms.begin(), ddenoms.end() );
        double ddenominator=0;
        for( vector<double>::const_iterator itr=ddenoms.begin(); itr!=ddenoms.end(); itr++ )
            ddenominator += *itr;

        dw = dnumerator / ddenominator;
    }
    else
        dw = dw_numerator / dw_denominator;

    dw = min( max(dw,-trustregion), trustregion );
        //if(trace) cout<<"\nOUT dnumerator ddenominator trustrgn dw "<<dnumerator<<" "<<ddenominator<<" "<<trustregion<<" "<< dw ;
    return dw;
}

pair<double,double>  //numer, denom
OneCoordStepBinary2(
	const InvData::pair_range_type range,
    const InvData& invData, //for 'y' only
    const Vector& wTX,
    double trustregion,
    bool highAccuracy,
    int group=-1 //disabled if <0
    )
{ 
    double dw_numerator = 0;
    double dw_denominator = 0;
    vector<double> dnums; vector<double> ddenoms;  // highAccuracy mode

    for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ )
    {
        unsigned i = ix->first;
//        if( group>=0 && !invData.group( i, group ) ) continue; //--->>---
        double x = ix->second;
        if( 0.0==x ) continue; //just in case - shouldn't ever hit this
        int target = invData.y(i) ? 1 : -1;
        double r = 1==target ? wTX[i] : -wTX[i]; //wTX[i] * target;
        double a = fabs(r);
        double b = fabs(trustregion * x);
        double F; // = a<=b ? 0.25 : 1.0/( 2.0 + exp(a-b) + exp(b-a) );
        if( a<=b ) F = 0.25;
        else {
            double e = exp(a-b);
            F = 1.0/( 2.0 + e + 1.0/e );
        }
        /*/original Z&O
        double F = min( 0.25, exp(fabs(trustregion * x)) * p * (1 - p))
            */
        double dd_numer = x * target / (1 + exp(r)) ;  //q * x * target; //(1 - p) * x * target
        double dd_denom = F * x * x;
        if( highAccuracy ){  //ver 2.04  accurate summation
            dnums.push_back( dd_numer );
            ddenoms.push_back( dd_denom );
        }else{
            dw_numerator += dd_numer;
            dw_denominator += dd_denom;
        }
        //if(trace) cout<<"\n x y xy r expr ddnum "<<x<<" "<<target<<" "<<x*target<<" "<<r<<" "<<exp(r)<<" "<<x * target / (1 + exp(r));
        //if(trace) cout<<"\n i x y r dw_numerator F dw_denominator "<<i<<" "<<x<<" "<<target<<" "<<r
          //  <<"\t  "<<dw_numerator<<" "<<F<<" "<<dw_denominator;
   }

    if( highAccuracy ){  //ver 2.04
        sort( dnums.begin(), dnums.end() );
        for( vector<double>::const_iterator itr=dnums.begin(); itr!=dnums.end(); itr++ )
            dw_numerator += *itr;

        sort( ddenoms.begin(), ddenoms.end() );
        for( vector<double>::const_iterator itr=ddenoms.begin(); itr!=ddenoms.end(); itr++ )
            dw_denominator += *itr;
    }

    return pair<double,double>( dw_numerator, dw_denominator );
}




pair<int,int> LRModel::countCorrect(RowSetIterator& rowset, const ParamMatrix& beta) const {

	rowset.rewind();

	int correct = 0;
	int count = 0;
	for(unsigned i=0; rowset.next(); i++ ) {
	    double bestScore = -1E10;
	    unsigned bestCategory = numeric_limits<unsigned>::infinity();
	    for( unsigned k=0; k<beta.numClasses(); k++ ) {
		const SparseVector& v = rowset.xsparse();
		double score = v.dot(beta.classparam(k));
		if ((score > bestScore)) {
		    bestScore = score;
		    bestCategory = k;
		}
	    }
	    unsigned y = rowset.y();
	    if (bestCategory == y) {
		correct++;
	    }
	    count++;
	}

	return pair<int,int>(correct, count);
}

bool LRModel::stopTest(MemRowSetIterator& rowset, int nRows, int nClasses, const ModelType& modelType, const ParamMatrix& beta, 
						  double sum_abs_dr, const vector< vector<double> >& linearScores, vector<double>& lastDiff) 
{
	if (modelType.getZoStoppingRule() == ModelType::zoStoppingRule_linearScores) {
		double sum_abs_r = 0;
		for( int i=0; i<nRows; i++ ) { 
			for( int k=0; k<nClasses; k++ ) {
				if( !finite(linearScores[k][i]) )
					Log(1)<<"\nNAN wTX[k][i] k i "<<linearScores[k][i]<<" "<<k<<" "<<i;
				sum_abs_r += fabs(linearScores[k][i]);
			}
		}

		/* commented for speed test, Shenzhi 07May07
		pair<int, int> trainingSetAccuracy =  countCorrect(rowset, beta);
		m_trainingSetAccuracy = trainingSetAccuracy.first;
		*/


		double rel_sum_abs_dr = sum_abs_dr/(1+sum_abs_r);

		m_convergenceScore = rel_sum_abs_dr;

		bool converge = (rel_sum_abs_dr < modelType.getConvergenceThreshold());
		return converge;

	} else if (modelType.getZoStoppingRule() == ModelType::zoStoppingRule_changeProb) {
		m_dScoreAboveThreshold = 0;
		m_trainingSetAccuracy = 0;

	    for(int i=0; i<nRows; i++ ) {
		    double bestScore = -1E10;
		    double secondBestScore = -1E10;
		    double dScore;

			unsigned bestCategory = numeric_limits<unsigned>::infinity();
			for( int k=0; k<nClasses; k++ ) {
				double score = linearScores[k][i];
				if ((score > bestScore) && (score > secondBestScore)) {
					secondBestScore = bestScore;
					bestScore = score;
					bestCategory = k;
				}
				else if (score > secondBestScore) {
					secondBestScore = score;
				}
			}
			if (bestCategory == rowset.getY(i)) {
				m_trainingSetAccuracy++;
			}

			dScore = (bestScore-secondBestScore-lastDiff[i])/lastDiff[i];
			lastDiff[i] = (bestScore-secondBestScore);

			if (fabs(dScore) > modelType.getZoStoppingRuleRScoreThreshold()) {
				m_dScoreAboveThreshold++;
			}
		}

		double ratioAboveThreshold = (double)m_dScoreAboveThreshold / (double)nRows;
		m_convergenceScore = ratioAboveThreshold;
		return ratioAboveThreshold <= modelType.getZoStoppingRuleRowLimit();
	}

	return false;

}

void LRModel::invokeZOLR(
    MemRowSetIterator& drsIterator,
    const InvData& invData, //input - generates design data
    const BayesParameter& bayesParam,
    const class FixedParams& fixedParams,  //input
    ParamMatrix & beta,  // i/0: init/result values of model coeffs, 
    bool binary)
{
    if (binary) {
	// TODO: Fix up const cast
	// vector<double> skew = vector<double>(drsIterator.dim(), 0);
	ZOLRBinary( drsIterator, invData,
		    bayesParam, 
// priorMean.classparam(0), priorScale.classparam(0), priorSkew.classparam(0), 
 		    const_cast<vector<bool>& >(fixedParams.classparam(0) ),
	 	    const_cast<vector<double>& >(beta.classparam(0)) );

    } else {
	ZOLR( drsIterator, invData, bayesParam, fixedParams, beta );
    }
}


void LRModel::ZOLRBinary(
    MemRowSetIterator& rowSet,
    const InvData& invData, //input - generates design data
    const BayesParameter& bayesParam, //input
//    const vector<double>& priorMode,  //input
//    const vector<double>& priorScale,  //input
//    const vector<double>& priorSkew,  //input
    vector<bool>& fixedParams,
    vector<double>& beta
    )
{
    Log(3)<<std::endl<<"Starting ZO model, Time "<<Log.time();
    //Log(3)<<std::endl<<bayesParam;

    unsigned nDesignVars = rowSet.dim();
    unsigned nRows = invData.n();
    if( nDesignVars==m_priorMode.numFeatures() ) ;//ok
    else    throw DimensionConflict(__FILE__,__LINE__);

    double common_variance = bayesParam.getPriorVar(); 
    Log(11) << " common_variance = "         << common_variance << endl;
    double p;
    if( Normal==bayesParam.getPriorType() ) {
	// analogous to 1/(2*tau_j) from Equation 8 of GLM2007
	p = 1 / (2 *common_variance);
    }
    else if( Laplace==bayesParam.getPriorType() ) {
	// analogous to lambda_j from Equation 8 of GLM2007
	double lambda = sqrt(2.0) / sqrt(common_variance); 
	p = lambda;
    }
    else
	throw runtime_error("ZO only allows Normal or Laplace prior");

    //convert prior scale into penalty parameter
    m_penalty.reset(nDesignVars, 1, p);

    if(m_individualPriorsHolder!=0) {
	m_individualPriorsHolder->reset();

	for(unsigned j=0; j<nDesignVars; j++) {
	    double final_variance;
	    
	    // type: 0 - not found; 1 - feature-level; 2 - coefficient-level; 
	    int type;
	    ModeVarSkew p = m_individualPriorsHolder->hasIndPrior(j,0,type);

	    // if not found, continue
	    if(type==0) continue;
	    
	    // if the prior is specified in indprior file, overwrite the current value
	    m_priorMode(j,0) = p.mode;
	    m_priorSkew(j,0) = p.skew;
	    
	    final_variance = type==0 ? common_variance : p.abs ? p.var : common_variance*p.var;
	    
	    // penalty(j,k) is multiplier of the beta portion of the actual penalty
	    if( Normal==bayesParam.getPriorType() ) {
		// analogous to 1/(2*tau_j) from Equation 8 of GLM2007
		m_penalty(j,0) = 1 / (2 * final_variance); //'2*' fixes bug, ver2.01
	    }
	    else if( Laplace==bayesParam.getPriorType() ) {
		// analogous to lambda_j from Equation 8 of GLM2007
		double lambda = sqrt(2.0) / sqrt(final_variance); 
		m_penalty(j,0) = lambda; 
	    }
	    else
		throw runtime_error("ZO only allows Normal or Laplace prior");
	}
    }

    Log(11)<<m_penalty<<endl;

    /*
    for( unsigned j=0; j<nDesignVars; j++ )
        if( Normal==bayesParam.getPriorType() )
            penalty[j] = 1.0 /( 2*bayesParam.getVar()*priorScale[j]*priorScale[j] );  //'2*' fixes bug, ver3.01
        else if( Laplace==bayesParam.getPriorType() )
            penalty[j] = bayesParam.getGamma() / priorScale[j];
        else
            throw runtime_error("ZO only allows Normal or Laplace prior");
    */

    // initialize the search
    Vector wTX( 0.0, nRows );
    for( unsigned j=0; j<nDesignVars; j++ ) {
        if( beta[j] != 0.0 ) {
	    //pair<InvDataItr,InvDataItr> range = invData.getRange(j);
	    //for( InvDataItr ix=range.first; ix!=range.second; ix++ ) // i such that xi,j != 0.0
	    const InvData::pair_range_type range = invData.getRange(j);
	    for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ )
	    {
		wTX[ ix->first ] += beta[j] * ix->second;
	    }
	}
    }
    
    
    Vector wTXprev = wTX;
    
    Vector trustregion( 1.0, nDesignVars );
    unsigned ndropped;
    
    vector<unsigned> sleep( nDesignVars, 0 );
    class WakeUpPlan {
        vector<bool> plan;        
    public:
        bool operator()(unsigned int step ) {
            if( step<plan.size() ) return plan.at(step);
            else return 0==step%100;
        }
        WakeUpPlan(unsigned size=1000) : plan( size+1, false ) { //ctor
            for(unsigned i=0;i*i<size;i++) plan[i*i] = true;     }
    };
    WakeUpPlan wakeUpPlan;
    
    bool converge = false;
    unsigned k;
    bool highAccuracy = m_modelType.isHighAccuracy();
    unsigned iterLimit = m_modelType.getIterLimit();
    double thrConverge = m_modelType.getTuneThreshold();
    
    int n1coordsteps=0, ndoubltries=0, nslept=0;
    for( k=0; !converge; k++ ) //---iterations loop--
    {
        vector<double> wPrev = beta;
        double dot_product_total_change = 0.0;
        ndropped = 0;
        for( unsigned j=0; j<nDesignVars; j++ ) //--design vars--
        {
	    double dw, dw_numerator, dw_denominator;
	    
	    const InvData::pair_range_type range = invData.getRange(j);
            //pair<InvDataItr,InvDataItr> range = invData.getRange(j);
            //Log(10)<<"\nVar "<<j<<" range: "<<range.first.i()<<"/"<<range.first.x()<<" - "<<range.second.i()<<"/"<<range.second.x();
	    

            if( Normal==bayesParam.getPriorType() )
            {
		
                dw_numerator = (beta[j]==m_priorMode(j,0)) ? 0.0  //in case of infinite penalty // ver 2.57
                    : - 2 * ( beta[j] - m_priorMode(j,0) ) * m_penalty(j,0);
		
                dw_denominator = 2 * m_penalty(j,0);
                dw = OneCoordStepBinary( dw_numerator, dw_denominator, range, invData, wTX, trustregion[j], highAccuracy );
	    }
            else { //Laplace
                if( beta[j]!=m_priorMode(j,0) || wakeUpPlan( sleep[j] ) ){
                    pair<double,double> numer_denom = OneCoordStepBinary2( range, invData, wTX, trustregion[j], highAccuracy );
                    n1coordsteps++;
                    if( beta[j]-m_priorMode(j,0)>0 ) {
                        dw = (numer_denom.first - m_penalty(j,0)) / numer_denom.second;
                    }
                    else if( beta[j]-m_priorMode(j,0)<0 ) {
                        dw = (numer_denom.first + m_penalty(j,0)) / numer_denom.second;
                        //dw = min( max(dw,-trustregion), trustregion );
                    }
                    else // at mean right now
                    { // try both directions
                        //if(k>0) Log(9)<<"\nTry Wake-up j/sleep[j] "<<j<<"/"<<sleep[j];
                        dw = 0;
                        if( m_priorSkew(j,0) > -1 ) { //positive allowed
                            double dwPlus = (numer_denom.first - m_penalty(j,0)) / numer_denom.second;
                            //dwPlus = min( max(dwPlus,-trustregion), trustregion );
                            if( dwPlus>0 )
                                dw = dwPlus;
                        }
                        if( dw==0 ) //positive step not allowed or unsuccessful
                        {
                            if( m_priorSkew(j,0) < 1 ) { //negative allowed
                                double dwMinus = (numer_denom.first + m_penalty(j,0)) / numer_denom.second;
                                //dwMinus = min( max(dwMinus,-trustregion), trustregion );
                                if( dwMinus<0 )
                                    dw = dwMinus; //othw dw stays at 0
                            }
                        }
                    }
                }
		else {
		    std::cout << "SLEEP" << wakeUpPlan(sleep[j]) << nslept << "\n";
                    nslept++;//Log(9)<<"\nNo Wake-up j/sleep[j] "<<j<<"/"<<sleep[j]
		}
		
                dw = min( max(dw,-trustregion[j]), trustregion[j] ); //for Normal prior this is done in OneCoordStep
		
                //check if we crossed the break point
                if(( beta[j] < m_priorMode(j,0) && m_priorMode(j,0) < beta[j]+dw )
		   || ( beta[j]+dw < m_priorMode(j,0) && m_priorMode(j,0) < beta[j] )
                    )
                    dw = m_priorMode(j,0) - beta[j]; //stay at mean
            } //Laplace
	    
            beta[j] += dw;
            if( beta[j]==m_priorMode(j,0) )   ndropped++;
            //update local data
	    for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ )
	        //for( InvTripletItr ix=range.first; ix!=range.second; ix++ ) // i such that xi,j != 0.0
	        //for( InvDataItr ix=range.first; ix!=range.second; ix++ ) // i such that xi,j != 0.0
            {
                unsigned i = ix->first;
                double x = ix->second;
                if( 0.0==x ) continue; //just in case - shouldn't ever hit this
		//double old = wTX[i];
                wTX[i] += dw * x;
		dot_product_total_change += fabs(dw * x); 
            }
            //trust region update
	    trustregion[j] = max( 2*fabs(dw), trustregion[j]/2 );
	    
            if( 0==beta[j]-m_priorMode(j,0) ) sleep[j]++;
            else {
                //if( sleep[j]>0 && k>0 ) Log(9)<<"\nWake-up j/sleep[j] "<<j<<"/"<<sleep[j];
                sleep[j]=0;
            }
	    
        }//--design vars-- j --
	
        double sum_abs_r = 0;
        for( unsigned i=0; i<nRows; i++ )
            sum_abs_r += fabs(wTX[i]);
        //ZO stopping
        double sum_abs_dr = 0.0;
        for( unsigned i=0; i<nRows; i++ )
            sum_abs_dr += fabs(wTX[i]-wTXprev[i]);
	
        double dot_product_rel_change = dot_product_total_change/(1+sum_abs_r);//DL stopping
        double rel_sum_abs_dr = sum_abs_dr/(1+sum_abs_r);//ZO stopping
	
        converge =( k>=iterLimit  ||
		    rel_sum_abs_dr<thrConverge );//ZO stopping
	//DL stopping || dot_product_rel_change<thrConverge );
	
        Log(7)<<"\nZO iteration "<<k+1
	      <<"  Dot product abs sum "<<sum_abs_r
	      <<"  Rel change "<<rel_sum_abs_dr //ZO stopping //DL stopping dot_product_rel_change
	      <<"  Beta relative increment = "<<normdiff(beta,wPrev)/norm(wPrev)
	      <<"  Beta components dropped: "<<ndropped
	      <<"  Time "<<Log.time()<<endl;
        /*dbg
	  Log(4)<<"\nBeta sparse:";
	  for( unsigned i=0; i<w.size(); i++ )  if(0!=w[i]) Log(4)<<" "<<i<<":"<<w[i];  */
	
        wTXprev = wTX;
    }//---iter loop---
    
    Log(7)<<std::endl<<"1-coord steps: "<<n1coordsteps<<"  Double tries: "<<ndoubltries<<" Slept: "<<nslept;
    Log(5)<<std::endl<<( k>=iterLimit ? "Stopped by iterations limit"  : "Stopped by original ZO rule" );
    Log(3)<<std::endl<<"Built ZO model "<<k<<" iterations, Time "<<Log.time() << endl; 
    Log(0).flush();
//    return w;
}



void LRModel::ZOLR(
    MemRowSetIterator& rowSet,
    const InvData& invData, //input - generates design data
    const BayesParameter& bayesParam,
    const class FixedParams& fixedParams,  //input
    ParamMatrix & beta  // i/0: init/result values of model coeffs, 
    )
{
    
    unsigned nDesignVars = rowSet.dim();
    unsigned nClasses = rowSet.c();
    unsigned nRows = invData.n();
    
    // common_variance is the variance of prior distribution used
    // for coefficients not mentioned in the individual prior file

    double common_variance = bayesParam.getPriorVar(); 
    Log(11) << " common_variance = "         << common_variance << endl;
    double p;
    if( Normal==bayesParam.getPriorType() ) {
	// analogous to 1/(2*tau_j) from Equation 8 of GLM2007
	p = 1 / (2 *common_variance);
    }
    else if( Laplace==bayesParam.getPriorType() ) {
	// analogous to lambda_j from Equation 8 of GLM2007
	double lambda = sqrt(2.0) / sqrt(common_variance); 
	p = lambda;
    }
    else
	throw runtime_error("ZO only allows Normal or Laplace prior");
    
    if(m_individualPriorsHolder!=0) 
	m_individualPriorsHolder->reset();

    for(unsigned k=0; k<nClasses; k++) {
	for(unsigned j=0; j<nDesignVars; j++) {
	    double final_variance;

	    if(m_individualPriorsHolder==0) 
		m_penalty(j,k) = p;
	    else {
		// type: 0 - not found; 1 - feature-level; 2 - coefficient-level; 
		int type;
		ModeVarSkew p = m_individualPriorsHolder->hasIndPrior(j,k,type);

		// if not found, continue
		if(type==0) continue;

		// if the prior is specified in indprior file, overwrite the current value
		m_priorMode(j,k) = p.mode;
		m_priorSkew(j,k) = p.skew;
		
		final_variance = type==0 ? common_variance : p.abs ? p.var : common_variance*p.var;
	    
		// penalty(j,k) is multiplier of the beta portion of the actual penalty
		if( Normal==bayesParam.getPriorType() ) {
		    // analogous to 1/(2*tau_j) from Equation 8 of GLM2007
		    m_penalty(j,k) = 1 / (2 * final_variance); //'2*' fixes bug, ver2.01
		}
		else if( Laplace==bayesParam.getPriorType() ) {
		    // analogous to lambda_j from Equation 8 of GLM2007
		    double lambda = sqrt(2.0) / sqrt(final_variance); 
		    m_penalty(j,k) = lambda; 
		}
		else
		    throw runtime_error("ZO only allows Normal or Laplace prior");
	    }
	}
    }
    Log(11)<<m_penalty<<endl;

    if( nDesignVars==m_priorMode.numFeatures() && nDesignVars==m_penalty.numFeatures() ) ;//ok
    else    throw DimensionConflict(__FILE__,__LINE__);

    Log(3)<<std::endl<<"Starting PolyZO model 2, - Time "<<Log.time()<<endl;
    
    // initialize temp data
    vector< vector<double> > linearScores( nClasses, vector<double>(nRows,0.0) ); //dot products, by classes by cases
    bool highAccuracy = m_modelType.isHighAccuracy();
    for( unsigned j=0; j<nDesignVars; j++ ) //--design vars--
    {
	const InvData::pair_range_type range = invData.getRange(j);
	for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ ) 
	    for( unsigned k=0; k<nClasses; k++ ) 
		if( beta(j,k) != 0.0 )
		    linearScores[k][ ix->first ] += beta(j,k) * ix->second;
    }
    Log(3)<<std::endl<<"Starting PolyZO model 3, - Time "<<Log.time() << endl;
    
    ParamMatrix trustregion(nDesignVars, nClasses, 1.0);
    unsigned ndropped;
    
    // TODO: Check that I didn't change this.
    const unsigned IterLimit = m_modelType.getIterLimit();
    
    ParamMatrix sleep( nDesignVars, nClasses,  0.0);
    WakeUpPlan wakeUpPlan;
    
    bool converge = false;
    unsigned itr;
    unsigned steps=0, stepsskipped=0, stepsslept=0;
    ParamMatrix lastDiff = ParamMatrix(beta.numFeatures(), beta.numClasses());
    Log(3)<<std::endl<<"Starting PolyZO model loop, - Time "<<Log.time()<<endl;
    for( itr=0; !converge; itr++ ) //---iterations loop--
    {
        double sum_abs_dr = 0.0;
        ndropped = 0;
        // unsigned classParams = m_modelType.getReferenceClass() ? nClasses-1 : nClasses;

        // for( unsigned k=0; k<classParams; k++ ) //--classes--
	for( unsigned k=0; k<nClasses; k++)  // ver 3.11
        {

	    // skip if this is the reference class
	    int classid = m_metadata.getRefClassId();
	    if (k==classid) continue;
	    
            //this could be outside all loops, but then error cumulates. ver 0.23
            vector<double> S( nRows, 0.0 );
            vector<double> exponentialScores( nRows );
	    for( unsigned i=0; i<nRows; i++ ) {
                exponentialScores[i] = exp( linearScores[k][i] );
		for( unsigned k=0; k<nClasses; k++ ) {
                    S[i] += exp( linearScores[k][i] );
		}
	    }
	    
	    
	    double t = Log.time();
	    double lastT = Log.time();
	    unsigned stepssk = 0;
	    
            vector<double> dwkTX( nRows, 0.0 );
            for( unsigned j=0; j<nDesignVars; j++ ) //--design vars--
            {
                if( fixedParams(j,k) ){  
		    Log(11)<<"skipped:" << j << " " << k << endl; 
		    stepsskipped++;
		    continue;        //----->>--
                }
                else{
		    Log(11)<<"steps: " << j << " " << k << endl;  
		    steps++;
		}
                if( beta(j,k)==m_priorMode(j,k) )
                    if( ! wakeUpPlan( (int)sleep(j,k) ) ) {
                        sleep(j,k) ++;
                        stepsslept++;
                        ndropped++;
                        continue;       //----->>--
                    }
		
                double dw;
		
		const InvData::pair_range_type range = invData.getRange(j);
		Log(11) << "invdata range: "<<range.second-range.first << endl; 
		
                if( Normal==bayesParam.getPriorType() )
                {
		    double dw_numerator = (beta(j,k)==m_priorMode(j,k)) ? 0.0 //in case of inf penalty; ver 3.11
			: - 2 * ( beta(j,k) - m_priorMode(j,k) ) * m_penalty(j,k);
                    double dw_denominator = 2 * m_penalty(j,k);
                    dw = OneCoordStep(
                        dw_numerator, dw_denominator, range, invData, k, linearScores, 
			linearScores[k], exponentialScores, S, trustregion(j,k), highAccuracy );
		    
		    Log(11) <<"dw_numerator="<<dw_numerator<<"; dw_denominator=" << dw_denominator 
		    <<"; dw="<<dw<<endl; 
                }
                else { //Laplace
                    if( beta(j,k)-m_priorMode(j,k)>0 ) { //step for positive
                        dw = OneCoordStep( -m_penalty(j,k), 0, range, invData, k, linearScores, linearScores[k],
					   exponentialScores, S, trustregion(j,k), highAccuracy );
                    }
                    else if( beta(j,k)-m_priorMode(j,k)<0 ) { //step for negative
                        dw = OneCoordStep( m_penalty(j,k), 0, range, invData, k, linearScores, linearScores[k],
					   exponentialScores, S, trustregion(j,k), highAccuracy );
                    }
                    else // at mean right now
                    {
			bool skewAllowsPositiveDirection = (m_priorSkew(j,k) >= -1);
			bool skewAllowsNegativeDirection = (m_priorSkew(j,k) <= 1);
			
			double dwPlus = 0;
			double dwMinus = 0;
			
			if (skewAllowsPositiveDirection) { 
                            dwPlus = OneCoordStep( -m_penalty(j,k), 0, range, invData, k, linearScores, linearScores[k], exponentialScores, S, trustregion(j,k), highAccuracy );
                            if(!finite(dwPlus))   Log(1)<<"\nInfDW zero neg "<<j<<" "<<k;
			    else if (dwPlus < 0) dwPlus = 0;
			}

			bool positiveDirectionSucceeded = (dwPlus > 0);
			if (!positiveDirectionSucceeded && skewAllowsNegativeDirection) {
			    dwMinus = OneCoordStep( m_penalty(j,k), 0, range, invData, k, linearScores, linearScores[k], exponentialScores, S, trustregion(j,k), highAccuracy );
                            if(!finite(dwMinus))   Log(1)<<"\nInfDW zero neg "<<j<<" "<<k;
			    else if (dwMinus > 0) dwMinus = 0;
			}
			
			dw = dwPlus + dwMinus;
                    }
		    
                    //check if we crossed the break point
                    if(( beta(j,k) < m_priorMode(j,k) && m_priorMode(j,k) < beta(j,k)+dw )
		       || ( beta(j,k)+dw < m_priorMode(j,k) && m_priorMode(j,k) < beta(j,k) )
                        )
                        dw = m_priorMode(j,k) - beta(j,k); //stay at mean
		}
		
		beta(j,k) += dw;
                if( beta(j,k)==m_priorMode(j,k) )   ndropped++;
		
                //update local data
                for( InvData::pair_vector_type::const_iterator ix=range.first; ix!=range.second; ix++ )
		    //for( vector< pair<unsigned, double> >::const_iterator ix=range.begin(); ix!=range.end(); ix++ )
                {
                    unsigned i = ix->first;
                    double x = ix->second;
                    if( 0.0==x ) continue; //just in case - shouldn't ever hit this
                    linearScores[k][i] += dw * x;
                    dwkTX[i] += dw * x;
                    double e_r = exp( linearScores[k][i] );
                    S[i] -= exponentialScores[i]; 
		    S[i] += e_r;
                    exponentialScores[i] = e_r;                
		    Log(11)<<"i="<<i<<";x="<<x<<"; linearScores[k][i]="<<linearScores[k][i]
			   <<";dwkTX[i]="<<dwkTX[i]<<";e_r="<<e_r<<";S[i]="<<S[i]
			   <<";expScore[i]="<<exponentialScores[i]<<endl;  
                }
                //trust region update
		trustregion(j,k) = max( 2*fabs(dw), trustregion(j,k)/2 );   //2*fabs(dw)*1.1
		
                if( beta(j,k)==m_priorMode(j,k) ) sleep(j,k)++;
                else                     sleep(j,k)=0;
                
		
            }//--design vars-- j --
	    
	    for( unsigned i=0; i<nRows; i++ ) {
                sum_abs_dr += fabs( dwkTX[i] );
	    }


	    Log(7)<<"iteration loop " << itr << ", Class "<<k<<" - Time "<<Log.time()<<endl;

	}//--classes-- k --
	
	converge = stopTest(rowSet, nRows, nClasses, m_modelType, beta, sum_abs_dr, linearScores, m_lastDiff);
	converge |= (itr >= IterLimit);
	
	iterationLogger(rowSet, bayesParam, beta, m_penalty, itr);
	
	Log(8)<<"iteration loop " << itr << ", - Time "<<Log.time()<<endl;

    }//---iter loop---
    
    Log(3)<<std::endl<<"Built ZO model "<<itr<<" iterations, Time "<<Log.time()<<endl;
    Log(5)<<"\nTotal steps "<<steps<<", skipped "<<stepsskipped<<", slept "<<stepsslept<<endl;
}


void LRModel::iterationLogger(MemRowSetIterator& rowSet, const BayesParameter& bayesParam, const ParamMatrix& beta, 
							  const ParamMatrix& penalty, unsigned itr) const {

    double likelihood = m_modelType.isComputeLogLikelihood() ? 
	LogLikelihoodPenalty(rowSet, beta.numClasses(), m_modelType, beta, bayesParam, penalty) : 0;
    
    int holdoutAccuracy = 0;
    int holdoutRows = 0;

    /*
    if (m_holdoutRowSetIterator != 0) {
	pair<int, int> holdoutCounts = countCorrect(*m_holdoutRowSetIterator, beta);
	holdoutAccuracy =  holdoutCounts.first;
	holdoutRows = holdoutCounts.second;
    }
    */

    Log(8) << "Iteration " << (itr+1) << '\t' <<
	m_convergenceScore << '\t' <<
	likelihood << '\t' <<
	Log.time() << '\t' <<
	m_trainingSetAccuracy << '/' << rowSet.n() << '\t' <<
	holdoutAccuracy << '/' << holdoutRows<<endl;
}

vector<int> LRModel::buildWordIndexForSample(RowSetIterator &drs, unsigned dim) const {
	drs.rewind();
	
	vector<int> wordIndex = vector<int>(dim+1, 0);
	for (const SparseVector& v = drs.xsparse(); drs.next(); ) {
		for (SparseVector::const_iterator i = v.begin(); i < v.end(); i++) {
			wordIndex[1 + i->first]++;
		}
	}

	int runningsum = 0;
	for (vector<int>::iterator i = wordIndex.begin(); i < wordIndex.end(); i++) {
		runningsum += *i;
		*i = runningsum;
	}

	drs.rewind();

	return wordIndex;
}

/* Perform loop of hyperparameter values for a given train/validation split */
vector<double>* LRModel::hyperParamLoop( MemRowSetIterator& drs,
					 unsigned dim, unsigned classes,
					 RowSetIterator& drsTest,
					 ModelType modelType,
					 const HyperParamPlan& hpPlan,
					 const class FixedParams& fixedParams
    )

{  
    vector<double>* lglkli = new vector<double>();
	
    unsigned planSize = hpPlan.getPlan().size();
    ParamMatrix localbeta( dim, classes );
    
    InvData* p_invData = 0;
    if (modelType.getOptimizerType() == ModelType::ZO) {
	vector<int> wordIndex = buildWordIndexForSample(drs, dim);
	p_invData = new InvData(drs, classes, dim, wordIndex);
    }
    
    for( unsigned iparam=0; iparam<planSize; iparam++ )
    {
        Log(5)<<"\nHyperparameter plan #"<<iparam+1<<" value="<<(hpPlan.getPlan()[iparam]);
        const BayesParameter& localBayesParam = hpPlan.getInstBayesParam( iparam );
	
        // build the model
	if( ModelType::ZO==modelType.getOptimizerType() ) {
	    MemRowSetIterator drsIterator = MemRowSetIterator(drs);
	    invokeZOLR(drsIterator, *p_invData, localBayesParam, fixedParams, localbeta, m_modelType.isBinary());
	    
	} else if( ModelType::QuasiNewtonDoubleCoord==modelType.getOptimizerType() || ModelType::QuasiNewtonSmooth==modelType.getOptimizerType() ) {
	    QN( drs, localBayesParam, modelType, fixedParams, localbeta, m_individualPriorsHolder );
	} else 
            throw logic_error("Unsupported optimizer");
	
        double eval = LogLikelihood( drsTest, classes, modelType, localbeta );
        lglkli->push_back(eval);
    }
    
    delete p_invData;
    return lglkli;
}


double LRModel::AvgSquNorm(RowSetStats& stats) const { //separated from class Stats' to exclude constant term
    double s = 0.0;
    for(unsigned j=0; j<stats.Means().size(); j++) {
        if (j==0) //HACK: drop constant term
            continue;
        s += stats.Means()[j]*stats.Means()[j] + stats.Stddevs()[j]*stats.Stddevs()[j];
	Log(16)<<"AvgSquNorm: "<<stats.Means()[j]<<";"<<stats.Stddevs()[j]<<endl;
    }
    return s;
}

void LRModel::tuneModel(   MemRowSetIterator & drs,
			   const vector<int>& rowIndex,
			   const HyperParamPlan& hyperParamPlan,
			   const class FixedParams& fixedParams,
			   RowSetStats& stats )

{
    // build the model
    if( ! hyperParamPlan.hasGrid() ) {    //==paramPlan.size() ) {
        BayesParameter localBayesParam;
        if( hyperParamPlan.getAutoHPar() ) {
	    double avgsqunorm = AvgSquNorm(stats)/(double)drs.dim();    // v3.21
            double hpar = Normal==hyperParamPlan.getPriorType() ? 
		1.0/avgsqunorm  : sqrt(avgsqunorm*2.0);
            localBayesParam = BayesParameter( hyperParamPlan.getPriorType(), hpar, hyperParamPlan.getSkew() );

            Log(5)<<"\nAverage square norm (no const term) "<<avgsqunorm<<" Hyperparameter "<<hpar<<endl;
        }
        else if( hyperParamPlan.isFixed() )
            localBayesParam = hyperParamPlan.getBp();
        else throw logic_error("Inconsistent 'HyperParamPlan' setting");


	MemRowSetIterator& drsIterator = drs;
        if( ModelType::ZO==m_modelType.getOptimizerType() ) {
	    Log(0) << "before invData initiation " << Log.time() << endl;  // for time track
            InvData invData(drsIterator, drs.c(), drs.dim(), rowIndex );
	    // invData.print();
	    Log(0) << "after invData initiation  " << Log.time() << endl; // for time track

            invokeZOLR(drsIterator, invData, localBayesParam, fixedParams, m_beta, m_modelType.isBinary());
        }
        else if( ModelType::QuasiNewtonDoubleCoord==m_modelType.getOptimizerType() || ModelType::QuasiNewtonSmooth==m_modelType.getOptimizerType() )
	    QN( drsIterator, localBayesParam, m_modelType, fixedParams, m_beta, m_individualPriorsHolder);
        else 
            throw logic_error("Unsupported optimizer");

        m_bayesParam = localBayesParam;
    }
    else {    //model averaging/selection

        unsigned planSize = hyperParamPlan.getPlan().size();
        vector<double> lglkliMean( planSize, 0.0 );
        vector<double> lglkliMeanSqu( planSize, 0.0 );
        const double randmax = RAND_MAX;

        const std::vector<double>& hpplan = hyperParamPlan.getPlan();

        //prepare CV
        unsigned nfolds = hyperParamPlan.getNfolds(); //10
        if( nfolds>drs.n() ) {
            nfolds = drs.n();
            Log(1)<<"\nWARNING: more folds requested than there are data. Reduced to "<<nfolds;
        }
        unsigned nruns = hyperParamPlan.getNruns(); //2
        if( nruns>nfolds )  nruns = nfolds;

	RowSetSampler sampler = RowSetSampler(drs, nfolds);
        vector<bool> foldsAlreadyRun( nfolds, false );

        //cv loop
        for( unsigned irun=0; irun<nruns; irun++ )
        {
            //next fold to run
            unsigned x = unsigned( rand()*(nfolds-irun)/randmax );
            unsigned ifold, y=0;
            for( ifold=0; ifold<nfolds; ifold++ )
                if( !foldsAlreadyRun[ifold] ) {
                    if( y==x ) break;
                    else y++;
                }
            foldsAlreadyRun[ifold] = true;
            Log(5)<<"\nCross-validation "<<nfolds<<"-fold; Run "<<irun+1<<" out of "<<nruns<<", fold "<<ifold;
	    
            //training
	    MemRowSetIterator* cvtrain = sampler.getSample(ifold, false);
            //Log(5)<<"\ncv training sample "<<cvtrain.n()<<" rows";
            MemRowSetIterator* cvtest = sampler.getSample(ifold, true); // ( drs, rndind, ifold, ifold+1, true ); //int(randmax*ifold/nfolds), int(randmax*(ifold+1)/nfolds)
            //Log(5)<<"\ncv test sample "<<cvtest.n()<<" rows";
	    
            vector<double>* loglikelihoods = hyperParamLoop( *cvtrain, drs.dim(), drs.c(), *cvtest, //invData,
							     m_modelType, hyperParamPlan, //bayesParameter, hpplan, 
							     fixedParams);
	    
            Log(5)<<"\nCross-validation test log-likelihood: ";
            for( unsigned iparam=0; iparam<planSize; iparam++ ) {
                double eval = (*loglikelihoods)[iparam];
                Log(5)<<eval<<" ";
                //arith ave of loglikeli, i.e. geometric mean of likelihoods
                lglkliMean[iparam] = lglkliMean[iparam]*irun/(irun+1) + eval/(irun+1);
                lglkliMeanSqu[iparam] = lglkliMeanSqu[iparam]*irun/(irun+1) + eval*eval/(irun+1);
            }
	    
	    delete loglikelihoods;
	    delete cvtrain;
	    delete cvtest;
	    
        }//cv loop
	
        vector<double> lglkliStdErr; //( hpplan.size(), 0.0 );
        for(unsigned i=0; i<planSize; i++ )
            lglkliStdErr.push_back( sqrt( lglkliMeanSqu[i] - lglkliMean[i]*lglkliMean[i] )//stddev
                    / sqrt(double(nruns)) ); 
	
        // best by cv
        double bestEval = - numeric_limits<double>::max();
        unsigned bestParam = unsigned(-1);
        Log(6)<<"\nCross-validation results - hyperparameter values, "
	      <<(Laplace==hyperParamPlan.getPriorType()?"prior var, ":"")
	      <<"cv mean loglikelihood, std.error:";
        for( unsigned i=0; i<planSize; i++ ) {
            if( lglkliMean[i]>bestEval ) {
                bestEval = lglkliMean[i];
                bestParam = i;
            }
            Log(6)<<"\n\t"<<hpplan[i];
            if(Laplace==hyperParamPlan.getPriorType()) 
                Log(6)<<"\t"<<hyperpar2var(hyperParamPlan.getPriorType(),hyperParamPlan.getSkew(),hpplan[i]);
            Log(6)<<"\t"<<lglkliMean[i]<<"\t"<<lglkliStdErr[i];
        }
        if( bestParam==unsigned(-1) )
            throw runtime_error("No good hyperparameter value found");
        if( hyperParamPlan.isErrBarRule() ) {
            double valueneeded = lglkliMean[bestParam] - lglkliStdErr[bestParam];
            unsigned oldBestParam = bestParam;
            for( unsigned i=0; i<planSize; i++ )
                if( lglkliMean[i]>=valueneeded ) {
                    bestParam = i;
                    break;
                }
            Log(6)<<"\nError bar rule selected hyperparameter value "<<hpplan[bestParam]<<" instead of "<<hpplan[oldBestParam];
        }
        Log(3)<<"\nBest hyperparameter value "<<hpplan[bestParam]<<" cv-average loglikely "<<bestEval;

        //build final model
        m_bayesParam = hyperParamPlan.getInstBayesParam( bestParam );

        Log(3)<<std::endl<<"Starting final model after cv, Time "<<Log.time();
	MemRowSetIterator drsIterator = MemRowSetIterator(drs);
        if( ModelType::ZO==m_modelType.getOptimizerType() ) {
            InvData invData(drsIterator, drs.c(), drs.dim(), rowIndex );
	    invokeZOLR(drsIterator, invData, m_bayesParam, fixedParams, m_beta, m_modelType.isBinary());
        }
        else if( ModelType::QuasiNewtonDoubleCoord==m_modelType.getOptimizerType() || ModelType::QuasiNewtonSmooth==m_modelType.getOptimizerType() )
	    QN( drsIterator, m_bayesParam, m_modelType, fixedParams, m_beta, m_individualPriorsHolder );
        else 
            throw logic_error("Unsupported optimizer");

    }//model averaging/selection

}


void LRModel::train(MemRowSet& trainData, WriteModel& modelFile, IRowSet* heldout)
{
    m_lastDiff = vector<double>(trainData.n(), 0);

    const BayesParameter& bayesParameter = m_hyperParamPlan.getBp();
    m_bTrained = false;
    
    MemRowSet& drs = trainData;
    // drs.print();

    int numClasses = drs.c();
    m_beta = ParamMatrix( drs.dim(), numClasses );
    m_priorMode = ParamMatrix( drs.dim(), numClasses);
    m_priorSkew = ParamMatrix( drs.dim(), numClasses, m_hyperParamPlan.getSkew() );
    m_penalty = ParamMatrix( drs.dim(), numClasses);
    //FixedParams
    const vector<vector<bool> > dummy; //(drs.dim(),vector<bool>(drs.c(),false));
    int classid = m_metadata.getRefClassId();
    FixedParams fixedParams( m_modelType.isAllZeroClassMode() ? trainData.allzeroes() : dummy,
			     m_individualPriorsHolder,
			     classid, drs.dim(), numClasses );
    Log(11)<<"FIXEDPARAMS "<<fixedParams<<endl;

    Log(3)<<"\nTotal parameters: "<<drs.dim()*numClasses<<"\t Of them fixed: "
	  <<fixedParams.count()<<endl;

    // m_beta initialization, read in the init file
    InitModel initfilemodel(m_initFName, m_metadata, bayesParameter, 
			fixedParams, *m_individualPriorsHolder, m_beta);
    Log(11)<<"coef after reading init file:" <<endl;
    Log(11)<< m_beta << endl;

    MemRowSetIterator* p_workingRs = new MemRowSetIterator(drs);

    RowSetStats stats(drs);
    RowSetIterator* it = drs.iterator();

    Log(0) << "\n Check 1 " << Log.time() << endl;

    //this sets m_beta and m_bayesParam
    tuneModel( *p_workingRs, //bayesParameter, 
	       trainData.getWordIndex(),
	       m_hyperParamPlan,
	       fixedParams, 
	       stats );
    
    cout<<m_beta<<endl;
    Log(0) << "\n Check 2 " << Log.time() << endl;
    

    //Beta components dropped finally
    unsigned bdropped = m_priorMode.countequal( m_beta );
    Log(3)<<"\nBeta components dropped finally: "<<bdropped
	  <<" Left: "<< m_beta.numClasses()*m_beta.numFeatures()-bdropped << endl;

    vector< vector<double> > resubstScore( drs.n(), vector<double>(numClasses,0.0) );
    vector<unsigned> y( drs.n() );
    for( unsigned j=0; it->next(); j++ ) {
        for( unsigned k=0; k<numClasses; k++ )
            resubstScore[j][k] = it->xsparse().dot(m_beta.classparam(k));
        y[j] = it->y();
    }
    
    // report model
    Log(10)<<std::endl<<"Beta:";
    m_beta.displayModel( Log(10), drs );

    // resubstitution - evaluate
    vector<unsigned> resubst( drs.n() );
    for( unsigned i=0; i<drs.n(); i++ )
        resubst[i] = argmax( resubstScore[i] );
    Log(3)<<"\n\n---Resubstitution results---";
    Log(3)<<"\nConfusion Table:";   //TODO? make it Log(6)
	makeConfusionTable( Log(3), drs.getRowSetMetadata(), y, resubst );

    double trainLogLikeli = LogLikelihood( drs, m_modelType, m_beta );
    Log(3)<<"Training set loglikelihood "<<setprecision(12)<<trainLogLikeli
	  <<" Average "<<trainLogLikeli/drs.n()<<setprecision(6)<<endl;
    double logPrior = LogPrior( m_bayesParam, m_priorMode, m_penalty, m_beta );
    Log(3)<<"Log prior (penalty) "<<logPrior<<" Log posterior "<<trainLogLikeli+logPrior<<endl;

    makeCT2by2( Log(3), drs.getRowSetMetadata(), y, resubstScore, resubst );

    if( Log.level()>=12 ) {
        Log()<<"Resubstitution Scores";
        for( unsigned i=0; i<drs.n(); i++ )
            Log()<<endl<<resubstScore[i]<<":"<<resubst[i];
    }

    m_bTrained = true;
    modelFile.WriteParams( m_modelType, m_designParameter, trainData.getRowSetMetadata(),  m_beta );

    // write results file

    it->rewind();
    for( unsigned r=0; it->next(); r++ ) {
		modelFile.resultsFile() << drs.getRowSetMetadata().getClassId(y[r])
            <<" "<< estprob(resubstScore[r]) //p_hat
            <<" "<< drs.getRowSetMetadata().getClassId(argmax(resubstScore[r]))
	        <<endl;
    }

    delete it;

}


void LRModel::test(PlainYRowSet& testData, WriteModel& modelFile)
{
    if( !m_bTrained )
        throw logic_error("Model not trained, unable to test");

    //DesignRowSet testDesignData( m_pDesign, testData );

    vector<unsigned> TP(testData.c(),0), FP(testData.c(),0), FN(testData.c(),0);
    vector< vector<unsigned> > CT( testData.c(), vector<unsigned>( testData.c(), 0 ) );
    vector< vector<double> > allScores;
    vector<unsigned> allYs;
    double logLhood = 0;
    unsigned nAlien = 0; //cases with untrained class

    unsigned n = 0;
    PlainYRowSetIterator* it = testData.plainYiterator();
    while( it->next() ) 
    {
        vector<double> predictScore( testData.c() );

	for( unsigned k=0; k<testData.c(); k++ ) {
            predictScore[k] = it->xsparse().dot(m_beta.classparam(k));
	}

	allScores.push_back( predictScore );
        allYs.push_back( it->y() );

        // get the classes on the test data
        unsigned prediction = argmax( predictScore );
        if( it->ygood() ) {
            CT [it->y()] [prediction] ++;
            if( prediction==it->y() )  
                TP[it->y()]++;
            else {
                FN[it->y()]++;
                FP[prediction]++;
            }
            logLhood += PointLogLikelihood( predictScore, it->y(), m_modelType );
            n++;
        }
        else 
            nAlien++;

		

	modelFile.resultsFile() << it->yid() //testData.getRowSetMetadata().getClassId(it->y())
				<<" "<< estprob( predictScore) //p_hat predictScore
				<<" "<< testData.getRowSetMetadata().getClassId(argmax(predictScore))
				<<endl;

    }
	delete it;

    Log(1)<<"\n\n---Validation results---";
    Log(1)<<"\nCases of trained classes "<<n<<"  Other "<<nAlien<<"  Total "<<n+nAlien;
    Log(1)<<"\nConfusion Table:";


	displayConfusionTable(Log(1), testData.getRowSetMetadata(), CT);
    Log(1)<<"\nTest set loglikelihood "<<logLhood<<" Average "<<logLhood/n;

	for( unsigned k=0; k<testData.c(); k++ ) {
		Log(1)<<"\n\nOne-vs-All view: Class "<<testData.getRowSetMetadata().getClassId(k);
        double roc = calcROC( allScores, allYs, k );
        displayCT2by2( Log(1), TP[k], FP[k], FN[k], n-TP[k]-FP[k]-FN[k] , roc );
    }

    Log(3)<<endl<<"Time "<<Log.time();

}





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
