#define  _USE_MATH_DEFINES
#include <math.h>

#include "Likelihood.h"

double PointLogLikelihood(
                const vector<double> & linscores,
                unsigned y,
                const class ModelType& modelType
            )
{
    double pointLogLhood;
//    if( ModelType::logistic==modelType.getLinkType() ) {
        double invProb = 1.0;
        for( unsigned c=0; c<linscores.size(); c++ )
            if( c!=y )
                invProb += exp( linscores[c] - linscores[y] );
        pointLogLhood = -log( invProb );
        //Log(11)<<"\nscore/y/invProb/lhood "<<linscores<<" | "<<y<<" "<<invProb<<" "<<pointLogLhood;
//    }
//    else
//        throw logic_error("only logistic link supported for likelihood calculations");

    return pointLogLhood;
}

double PointLogLikelihoodAndGradient(
                const vector<double> & linscores,
                unsigned y,
				ParamMatrix& grad,
				const SparseVector& data
            )
{

    double pointLogLhood;
    double invProb = 1.0;
    vector<double> invProbs( linscores.size() , 1);
	for( unsigned c=0; c<linscores.size(); c++ ) {
		if( c!=y ) {
			double expScore = exp( linscores[c] - linscores[y] );
			invProbs[c] += expScore;
            invProb += invProbs[c];
		}
	}
    pointLogLhood = -log( invProb );

	for( SparseVector::const_iterator ix=data.begin(); ix!=data.end(); ix++ ) {
		for( unsigned k=0; k<linscores.size(); k++ ) {
//            if( ! fixedParams(ix->first,k) ) {
                grad(ix->first,k) -= ( k==y ? ix->second : 0 ) - ix->second/invProbs[k];
//            }
		}
	}


    return pointLogLhood;
}

double LogLikelihood(
                IRowSet & rowset,
                const class ModelType& modelType,
                const ParamMatrix& beta
            )
{
    RowSetIterator* it = rowset.iterator();
    double logLhood = 0;
    double i;

    // fit term
    for( i=0; it->next(); i++ ) {
        vector<double> scores;
        for( unsigned k=0; k<rowset.c(); k++ )
            scores.push_back( it->xsparse().dot(beta.classparam(k) ));
        logLhood = logLhood*(i/(i+1)) + PointLogLikelihood( scores, it->y(), modelType )/(i+1);
    }

    delete it;

    return logLhood*i;
}

double LogLikelihood(
                RowSetIterator& it,
				int classes,
                const class ModelType& modelType,
                const ParamMatrix& beta
            )
{
    double logLhood = 0;
    double i;

    // fit term
    for( i=0; it.next(); i++ ) {
        vector<double> scores;
        for( unsigned k=0; k<classes; k++ )
            scores.push_back( it.xsparse().dot(beta.classparam(k) ));
        logLhood = logLhood*(i/(i+1)) + PointLogLikelihood( scores, it.y(), modelType )/(i+1);
    }

    return logLhood*i;
}

double LogLikelihoodAndGradient(IRowSet & rowset,
                const class ModelType& modelType,
                const ParamMatrix& beta,
				ParamMatrix& grad)
 {

    RowSetIterator* it = rowset.iterator();
    double logLhood = 0;
    double i;

    // fit term
    for( i=0; it->next(); i++ ) {
        vector<double> scores;
		const SparseVector& row(it->xsparse());
        for( unsigned k=0; k<rowset.c(); k++ )
            scores.push_back( row.dot(beta.classparam(k) ));
        logLhood = logLhood*(i/(i+1)) + PointLogLikelihoodAndGradient( scores, it->y(), grad, row )/(i+1);
    }
	delete it;

    return logLhood*i;


}


double LogLikelihoodPenalty(
                RowSetIterator& rowset,
		int classes,
                const class ModelType& modelType,
                const ParamMatrix& beta,
		const BayesParameter& bayesParam,
		const ParamMatrix & penalty
    )
{
    double logLhood = 0;
    double i;
    
    // fit term
    for( i=0; rowset.next(); i++ ) {
        vector<double> scores;
        for( unsigned k=0; k<classes; k++ )
            scores.push_back( rowset.xsparse().dot(beta.classparam(k) ));
        logLhood = logLhood*(i/(i+1)) + PointLogLikelihood( scores, rowset.y(), modelType )/(i+1);
    }
    
    
    // penalty term
    double penaltyTerm = 0;
    double p = 0;
    double s = 0;
    
    {//---penalty term---
	for( unsigned j=0; j<beta.numFeatures(); j++ ) {
	    for( unsigned k=0; k<beta.numClasses(); k++ ) {
		p = s = 0;
		if( Normal==bayesParam.getPriorType() || Laplace==bayesParam.getPriorType() )
		    p = penalty(j,k);
		else
		    throw runtime_error("Only Normal or Laplace priors allowed");
		if( j==beta.numFeatures()-1 ) p = 0.0; //HACK - no penalty for intercept
		
		
		double s = beta(j,k); //  - priorMean(j,k);
		if( Normal==bayesParam.getPriorType() ) {
		    penaltyTerm += s * s * p;
		}
		else if( Laplace==bayesParam.getPriorType() ) {
		    if( ModelType::QuasiNewtonSmooth==modelType.getOptimizerType() ) {
			/*
			  double p = sqrt( s*s + esmooth );
			  penaltyTerm += p * penalty[i];
			*/
		    }else{ //2-coord
			penaltyTerm += fabs(s) * p;
		    }
		}
	    }
	}
    }
    
    
    //printf("\nLL/penalty: %f %f\n", -logLhood*i, penaltyTerm);
    
    return logLhood*i - penaltyTerm;
}


double LogPrior(
                const BayesParameter& bayesParameter,
                const ParamMatrix& priorMode,
                const ParamMatrix& penalty,
                const ParamMatrix& beta
            )
{
    double penaltyLogLhood = 0;
    for( unsigned j=0; j<beta.numFeatures(); j++ ) 
	for( unsigned k=0; k<beta.numClasses(); k++ ) 
	{
	    double betaLogLhood;
	    if( penalty(j,k) ==0 )
		betaLogLhood = 0;
	    else if( Normal==bayesParameter.getPriorType() ) {
		double u = beta(j,k)-priorMode(j,k);
		betaLogLhood = - log(sqrt(0.5/penalty(j,k))) -log(2*M_PI)/2 - u*u/(2*penalty(j,k));
	    }
	    else if( Laplace==bayesParameter.getPriorType() ){
		betaLogLhood = - log(2.0) + log(penalty(j,k)) - penalty(j,k)*fabs( beta(j,k)-priorMode(j,k) ); 
	    }
	    penaltyLogLhood += betaLogLhood;
	}
    
    return penaltyLogLhood;
}

// Some kind of area under the curve thing.
double calcROC( const vector< vector<double> >& allScores,
    const vector<unsigned>& allYs, 
    unsigned k )
{
    unsigned n=allScores.size();
    vector< std::pair<double,bool> > forROC;
    for( unsigned i=0; i<n; i++ )
        forROC.push_back( pair<double,bool>( allScores.at(i).at(k), k==allYs.at(i)) );

    std::sort( forROC.begin(), forROC.end() );
    double area = 0;
    double x=0, xbreak=0;
    double y=0, ybreak=0;
    double prevscore = - numeric_limits<double>::infinity();
    for( vector< pair<double,bool> >::reverse_iterator ritr=forROC.rbegin(); ritr!=forROC.rend(); ritr++ )
    {
        double score = ritr->first;
        bool label = ritr->second;
        if( score != prevscore ) {
            area += (x-xbreak)*(y+ybreak)/2.0;
            xbreak = x;
            ybreak = y;
            prevscore = score;
        }
        if( label )  y ++;
        else     x ++;
    }
    area += (x-xbreak)*(y+ybreak)/2.0; //the last bin
    if( 0==y || x==0 )   area = 0.0;   // degenerate case
    else        area = 100.0 * area /( x*y );
    return area;
}



vector<double> estprob( const vector<double>& score ) {
    vector<double> ret;
    double denom = 0;
    for( vector<double>::const_iterator itr=score.begin(); itr!=score.end(); itr++ )  {
        double e = exp(*itr);
        denom += e;
        ret.push_back( e );
    }
    for( vector<double>::iterator itr=ret.begin(); itr!=ret.end(); itr++ )
        *itr /= denom;
    return ret;
}

