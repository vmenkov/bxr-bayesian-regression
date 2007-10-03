#ifndef QUASI_NEWTON_H
#define QUASI_NEWTON_H

void QN(
                MemRowSetIterator& xdata,
                const BayesParameter& bayesParam, //input
                //const ParamMatrix & priorMean,  //input
                //const ParamMatrix & priorScale,  //input
                const ModelType& modelType,  //input
                const class FixedParams& fixedParams,  //input
                ParamMatrix & w,  // i/0: init/result values of model coeffs, 
		IndividualPriorsHolder* indprior
            );




#endif // QUASI_NEWTON_H

