#ifndef _LIKELIHOOD_H
#define _LIKELIHOOD_H

#include "Compiler.h"
#include "Data.h"
#include "PolyZO.h"

double PointLogLikelihood(
                const vector<double> & linscores,
                unsigned y,
                const class ModelType& modelType
            );


Vector Stddevs( IRowSet & drs );

double LogLikelihood(
                IRowSet & trainDesignData, // HACK!!! DesignRowSet needed: colId's are 0, 1, 2... 
                const class ModelType& modelType,
                const ParamMatrix& beta
            );

double LogLikelihood(
                RowSetIterator& it,
				int classes,
                const class ModelType& modelType,
                const ParamMatrix& beta
            );

double LogLikelihoodPenalty(
                RowSetIterator& it,
				int classes,
                const class ModelType& modelType,
                const ParamMatrix& beta,
                const BayesParameter& bayesParameter,
				const ParamMatrix & priorScale
            );

double LogPrior(
                const BayesParameter& bayesParameter,
                const ParamMatrix& priorMean,
                const ParamMatrix& priorScale,
                const ParamMatrix& beta
            );

double calcROC( const vector< vector<double> >& allScores,
    const vector<unsigned>& allYs, 
    unsigned k );

vector<double> estprob( const vector<double>& score );


#endif // _LIKELIHOOD_H

