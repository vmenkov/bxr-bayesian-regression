/*
 *  quasi-Newton optimization wrapper
 */
#define  _USE_MATH_DEFINES
#include <math.h>
//#include <algorithm>
//#include <limits>
//#include <iomanip>

#include "logging.h"
//#include "Matrix.h"
#include "BayesParamManager.h"
//#include "PriorTerms.h"
#include "Design.h"
#include "PolyZO.h"
#include "ModelTypeParam.h"
//#include "knitro/knitro.h"
#include "IndividualPriors.h"

using namespace std;

#ifdef _MSC_VER
#define finite(x) (_finite(x))
#endif


//prototyped from blmvm-1.0\example-c\combustionv1.c
extern "C" {
#include "Blmvmabc.h"
}

double LogLikelihoodAndGradient(IRowSet & rowset,
                const class ModelType& modelType,
                const ParamMatrix& beta,
				ParamMatrix& grad);

    //int BLMVMFunctionGradient(void*,double[],double[],int,double*);
    //int BLMVMConverge(void*,double[],int,double,int*);
    //extern int BLMVMSolve(void*,double[],double[],double[],int,int);

typedef struct {
    double f,fmin;             /* Current function value, stop is less than     			     */
    double pgnorm,pgtol;       /* Current residual norm, stop if less than      			     */
    int iter, maxiter;         /* Current number of iterations, stop if greater than                     */
    int fgeval, maxfgeval;     /* Current number of objective/gradient evaluations, stop ig greater than */
    int whystop;               /* Reason why we stopped BLMVM                                            */
	double ll;
} CConverge;

class Context {     /* This user-defined context contains application specific information */
    MemRowSetIterator& data;
    const BayesParameter& bayesParam;
    ParamMatrix priorMode;
    //const ParamMatrix & priorScale;
    const ModelType& modelType;
    const class FixedParams& fixedParams;
    class IndividualPriorsHolder* individualPriorsHolder;
    unsigned c, d, nparams, nPseudoParams;
    bool dupcoord;
    vector<double> penalty;
    ParamMatrix & w; //fixed terms never change
    ParamMatrix *wCopy;
    ParamMatrix *lastDiff;
    ParamMatrix grad; //only to avoid reallocation
    ParamMatrix hessInvDiagonal;
    double esmooth; //constant for smoothing L1
    double *px, *pxl, *pxu;
    int history;
public:
    CConverge cc;
    Context(
        MemRowSetIterator& data_,
        const BayesParameter& bayesParam_,
        //const ParamMatrix & priorMode_,
        //const ParamMatrix & priorScale_,
        const ModelType& modelType_,
        const FixedParams& fixedParams_,
        ParamMatrix & w,
	IndividualPriorsHolder* indprior_);
    ~Context() {
        free(px);
        free(pxl);
        free(pxu);
    }

    int FuncGrad( const double *x, double *gFlat, int n, double *f, double *hessianInverseDiagonal );
    // int FuncGrad( const double *x, double *g, int nparams, double *f );
    int Converge( double *residual, int nparams, double stepsize, int *whystop );

    unsigned NParams() const { return  nPseudoParams; }
    double* Low() {return pxl;} 
    double* X() {return px;} 
    double* Up() {return pxu;} 
    int History() {return history;} 
    double Esmooth() { return esmooth; }
    double set_esmooth( double e ) { return esmooth=e; }

    void square2flat( const ParamMatrix & w, double* x ) {
        unsigned i=0;
        for( unsigned j=0; j<w.numFeatures(); j++ )
            for( unsigned k=0; k<w.numClasses(); k++ )
                if( ! fixedParams(j,k) )
                    x[i++] = w(j,k);
    }
    void flat2square( double* x, ParamMatrix & w ) {
        unsigned i=0;
        for( unsigned j=0; j<w.numFeatures(); j++ )
            for( unsigned k=0; k<w.numClasses(); k++ )
                if( ! fixedParams(j,k) )
                    w(j,k) = x[i++];
    }
};
Context:: Context(
        MemRowSetIterator& data_,
        const BayesParameter& bayesParam_,
        //const ParamMatrix & priorMode_,
        //const ParamMatrix & priorScale_,
        const ModelType& modelType_,
        const FixedParams& fixedParams_,
        ParamMatrix & w_ , //CConverge& cc_ 
	IndividualPriorsHolder* indprior_)
        : data(data_), bayesParam(bayesParam_), 
	  /*priorMode(priorMode_), priorScale(priorScale_),*/ modelType(modelType_), fixedParams(fixedParams_), individualPriorsHolder(indprior_),
        w(w_), 
        d( data.dim() ), c( data.c() ),
        grad( d, c, 0.0 ),
        hessInvDiagonal( d, c, 0.0 ),
        px(0), pxl(0), pxu(0)
{
    Log(0) << "\n Context 1 " << Log.time();
    
    wCopy = new ParamMatrix(w_.numFeatures(), w_.numClasses());
    lastDiff = new ParamMatrix(w_.numFeatures(), w_.numClasses());
    priorMode = ParamMatrix(w_.numFeatures(),w_.numClasses());

    Log(0) << "\n Context 2 " << Log.time();
    if( data.dim()==priorMode.numFeatures() );
    else    throw DimensionConflict(__FILE__,__LINE__);
    
    dupcoord = ( Laplace==bayesParam.getPriorType() && ModelType::QuasiNewtonDoubleCoord==modelType.getOptimizerType() );

    nparams = c*d - fixedParams.count();
    nPseudoParams = dupcoord ? 2*nparams : nparams;

    Log(0) << "\n Context 3 " << Log.time();
    /* initial vector, lower bound, upperbound */
    px=(double*)malloc(nPseudoParams*sizeof(double)); 
    pxl=(double*)malloc(nPseudoParams*sizeof(double));
    pxu=(double*)malloc(nPseudoParams*sizeof(double));

    Log(0) << "\n Context 4 " << Log.time();
    /* set bounds: no bounds */
    for( unsigned i=0; i<nPseudoParams; i++ ) {
        pxl[i] = dupcoord ? 0 : -numeric_limits<double>::max();
        pxu[i] = numeric_limits<double>::max();
    }

    Log(0) << "\n Context 5 " << Log.time();
    /* initial point */
    square2flat( w, px );
    if(dupcoord)
        for( unsigned i=0; i<nparams; i++ )
            if( px[i]>=0 ) px[nparams+i]=0;
            else{ px[nparams+i]=-px[i]; px[i]=0; }
    
    Log(0) << "\n Context 6 " << Log.time();
    // init convergence data
    cc.fmin = 0; //-numeric_limits<double>::max(); -1e30;
    cc.pgtol = 1E-8; // modelType.ThrConverge(); //.001;
    cc.iter = 0;
    cc.maxiter = 55; // numeric_limits<int>::max(); //100;
    cc.fgeval = 0;
    cc.maxfgeval = numeric_limits<int>::max(); //10000; 

    history = 10;

    //convert prior scale into penalty parameter
    // Log(3)<<std::endl<<"No penalty for the intercept!";
    penalty.reserve(nparams);
    if(individualPriorsHolder!=0) 
	individualPriorsHolder->reset();

    // calculate the value if no ind prior file is specified
    double common_variance = bayesParam.getPriorVar(); 
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
	throw runtime_error("only allows Normal or Laplace prior");
    
    // then push the penalty data into the vector, over write if -I specified
    for(unsigned k=0; k<c; k++) {
	for(unsigned j=0; j<d; j++) {
	    // skip if it is locked
	    if( fixedParams(j,k) )  continue;

	    double final_variance;

	    if(individualPriorsHolder!=0) {
		// type: 0 - not found; 1 - feature-level; 2 - coefficient-level; 
		int type;
		ModeVarSkew mvs = individualPriorsHolder->hasIndPrior(j,k,type);

		// if not found, continue
		if(type==0) continue;
		
		// if the prior is specified in indprior file, overwrite the current value
		priorMode(j,k) = mvs.mode;
		
		final_variance = type==0 ? common_variance 
		    : mvs.abs ? mvs.var : common_variance*mvs.var;
		
		// penalty(j,k) is multiplier of the beta portion of the actual penalty
		if( Normal==bayesParam.getPriorType() ) {
		    // analogous to 1/(2*tau_j) from Equation 8 of GLM2007
		    p = 1.0 / (2.0 * final_variance); //'2*' fixes bug, ver2.01
		}
		else if( Laplace==bayesParam.getPriorType() ) {
		    // analogous to lambda_j from Equation 8 of GLM2007
		    double lambda = sqrt(2.0) / sqrt(final_variance); 
		    p = lambda; 
		}
		else
		    throw runtime_error("ZO only allows Normal or Laplace prior");
		
	    }

	    // keep the penalty value
	    if( j==0 ) p = 0.0; //HACK - no penalty for intercept
	    penalty.push_back( p );
	}
    }

    /*
    for( unsigned j=0; j<d; j++ )
	for( unsigned k=0; k<c; k++ )
	{
	    if( fixedParams(j,k) )  continue;
	    double p;
	    if( Normal==bayesParam.getPriorType() )
		p = 1.0 /( 2.0*bayesParam.getVar()*priorScale(j,k)*priorScale(j,k) );
	    else if( Laplace==bayesParam.getPriorType() )
		p = bayesParam.getGamma() / priorScale(j,k);
	    else
		throw runtime_error("Only Normal or Laplace priors allowed");
	    if( j==d-1 ) p = 0.0; //HACK - no penalty for intercept
	    penalty.push_back( p );
	    // Log(11)<<"\n k/j/priorScale(j,k)/priorMode(j,k)/penalty(j,k)/"<<k<<" "<<j<<" "<<priorScale(j,k)<<" "<<priorMean(j,k)<<" "<<p;
	}
    */

    Log(3)<<std::endl<<"No penalty for the intercept!";
    if( penalty.size() != nparams )
        throw logic_error("Penalty array size exception");
    
    Log(0) << "\n Context 7 " << Log.time();
    
}

void /*returns model coeffs thru the last arg*/
QN_BLMVM(
    MemRowSetIterator& data,
    const BayesParameter& bayesParam, //input
    //const ParamMatrix & priorMode,  //input
    //const ParamMatrix & priorScale,  //input
    const ModelType& modelType,  //input
    const class FixedParams& fixedParams,  //input
    ParamMatrix & w,  // i/0: init/result values of model coeffs, 
    IndividualPriorsHolder* indprior
    )
{
    //throw logic_error("Quasi-Newton NOT IMPLEMENTED YET");
    if( Laplace!=bayesParam.getPriorType()
        && Normal!=bayesParam.getPriorType() )
            throw runtime_error("Unsupported prior type; only Normal or Laplace currently supported");

    if( ModelType::QuasiNewtonSmooth==modelType.getOptimizerType() );
    else if( ModelType::QuasiNewtonDoubleCoord==modelType.getOptimizerType() );
    else 
        throw logic_error("Unsupported optimizer for quasi-Newton");

    /*w(1,1)=.005;dbg*/

//	const IParamMatrix& priorMean2 = IParamMatrix();
//	priorMean2.
    Context context(data, bayesParam, /*priorMode, priorScale,*/
		    modelType, fixedParams, w, indprior );

    Log(3)<<endl<<"Starting BLMVM optimization, "
	  <<( Laplace==bayesParam.getPriorType() ? 
	      ( ModelType::QuasiNewtonSmooth==modelType.getOptimizerType() ? "Smoothed, " : "Double coordinates, ") //Crude L1
	      : "" )
	  <<"Time "<<Log.time();
    
    if( Laplace==bayesParam.getPriorType() && ModelType::QuasiNewtonSmooth==modelType.getOptimizerType() )
        Log(3)<<"\nesmooth "<<context.set_esmooth( 1e-25 ); //1e+200*numeric_limits<double>::min() );
    
    /*for(unsigned i=0;i<nparams;i++) Log(8)<<"\nx "<<x[i];dbg*/
    int info=BLMVMSolveIt( (void*)&context,
			   context.Low(), context.X(), context.Up(), context.NParams(), context.History() );

    // opt point (context has modified 'w' thru its reference)
    Log(9)<<"\nFinal Beta "<<w;
    string diagnos;
    switch(context.cc.whystop) {
        case 1: diagnos = "Solution found"; break;
        case 2: diagnos = "Max iterations reached"; break;
        case 3: diagnos = "Max  reached"; break;
        case 4: diagnos = "Min Objective reached"; break;
        default: diagnos = "Stopped due to numerical difficulties";
    }
    Log(0)<<endl<<"Whystop " << context.cc.whystop<<endl;
    Log(3)<<endl<<"Optimization finished "<<context.cc.iter<<" iterations, "<<diagnos<<", Time "<<Log.time();
    Log(4)<<endl<<"Function/Gradient evaluations "<<context.cc.fgeval;
    
}

//comments from blmvm-1.0\example-c\combustionv1.c:
/* --------------------------------------------------------------------------
 *    BLMVMFunctionGradient -- User must implement this routine that evaluates 
 *    the objective function and its gradient at current point x 
 *
 *    Input:   ctx - pointer to user-defined context.
 *             x - current variables (array of length n).
 *             g - empty (array of length n)
 *             n - number of variables.
 *             f - address of current objective function value;
 *
 *    Output:  ctx - pointer to user-defined context.
 *             x - unchanged.
 *             g - the gradient of objective function f(x) evaluated at x.
 *             n - unchanged.
 *             f - objective function evaluated at x.
 *
 *             return:  0 for Normal return and nonzero to abort mission.
 *                                                                         */

//extern int BLMVMFunctionGradient(void*,double[],double[],int,double*,double*);

int BLMVMFunctionGradient(void *ctx,double x[],double g[],int n, double *f, double *hessInvDiag)
{
    Context* pCtx = (Context*)ctx;
    return pCtx->FuncGrad( x, g, n, f, hessInvDiag );
}

extern double LogLikelihoodPenalty(
                IRowSet & trainDesignData, // HACK!!! DesignRowSet needed: colId's are 0, 1, 2... 
                const class ModelType& modelType,
                const ParamMatrix& beta,
                const BayesParameter& bayesParameter,
				const ParamMatrix & priorScale
            );


/*
int Context::ApplyInvHessianDiagonal(const double *x, double *g, int nparams) {
	if( n != nPseudoParams ) {
        throw logic_error("Wrong parameter vector size returned from BLMVMFunctionGradient");
	}

    unsigned i=0;
	for( unsigned j=0; j<w.D(); j++ ) {
		for( unsigned k=0; k<w.C(); k++ ) {



			i++;
		}
	}


}
*/


int Context::FuncGrad( const double *x, double *gFlat, int n, double *f, double *hessianInverseDiagonal )
{
    if( n != nPseudoParams )
        throw logic_error("Wrong parameter vector size returned from BLMVMFunctionGradient");

    Log(0) << "\nCheck 1 "<< Log.time();

    //Log(0)<<"\nComputing gradient at X "; for(unsigned i=0;i<(dupcoord?2*nparams:nparams);i++) Log(0)<<" "<<x[i];
    {//flat2square( wFlat, w ); //fixed terms already there
	unsigned i=0;
	for( unsigned j=0; j<w.numFeatures(); j++ )
	    for( unsigned k=0; k<w.numClasses(); k++ )
		if( ! fixedParams(j,k) ) {
		    if(dupcoord)  w(j,k) = x[i] - x[nparams+i];
		    else  w(j,k) = x[i];
		    i++;
		}
    }
    //Log(0)<<"\nBeta "<<w;
    Log(0) << "\nCheck 2 "<< Log.time();
    
    for( unsigned j=0; j<d; j++ ) {
	for( unsigned k=0; k<c; k++ ) {
            grad( j, k ) = 0.0;
	    hessInvDiagonal( j, k ) = 0.0;
	}
    }
    
//    double** m2 = new double*[grad.numFeatures()];
//    for (int i = 0; i < grad.numFeatures(); i++) m2[i] = new double[grad.numClasses()];
    
    Log(0) << "\nCheck 3 "<< Log.time();
    
    double lossTerm = 0;
    
    Log(0) << "\nCheck 3.7 "<< Log.time() << " " << lossTerm;
    
    //---loss term---
    MemRowSetIterator i = MemRowSetIterator(data);
    i.rewind();
    int index = 0;
    while( i.next() ) {
	
        vector<double> linscores( c, 0.0 );
	
	const SparseVector& x = i.xsparse();
        for( unsigned k=0; k<c; k++ )
	    linscores[k] += x.dot(w.classparam(k));
	  /*
	  try{ //pass one
	  for( SparseVector::const_iterator ix=data.xsparse().begin(); ix!=data.xsparse().end(); ix++ )
	  for( unsigned k=0; k<c; k++ )
	  linscores[k] += w(ix->first,k) * ix->second;
	  }catch(...){
	  throw logic_error("sparse vector index beyond dense vector size"); }
	  //Log(9)<<"\nlinscore "<<linscores;
	  */
	
	vector<double> expscores(c);
	double sumExpscores = 0;
	for (unsigned k = 0; k < c; k++) {
	    expscores[k] = exp(linscores[k]);
	    sumExpscores += expscores[k];
	}
	
        vector<double> invPhat( c );
        for( unsigned k=0; k<c; k++ )  {
            invPhat[k] = sumExpscores/expscores[k];
        }
	
	/*
	  for( unsigned k=0; k<c; k++ )  {
	  invPhat[k] = invPhat[k]/expscores[k];
	  }
	  
	  vector<double> invPhat2( c );
	  for( unsigned k=0; k<c; k++ )  {
	  invPhat2[k] = 1.0;
	  for( unsigned kk=0; kk<c; kk++ )
	  if( kk!=k )
	  invPhat2[k] += exp( linscores[kk] - linscores[k] );
	  }*/
	
	
	lossTerm += log( invPhat[i.y()] );
	
        //pass two - gradient
	YType y = data.getY(index);
	
	
	for( SparseVector::const_iterator ix=x.begin(); ix!=x.end(); ix++ ) {
	    double value = ix->second;
	    for( unsigned k=0; k<c; k++ ) {
//                if( ! fixedParams(ix->first,k) ) {
		double& r = grad(ix->first,k);
                r -= ( k==y ? value : 0 ) - value/invPhat[k];
//                    hessInvDiagonal(ix->first,k) -=  (ix->second*ix->second)* (1.0/invPhat[k]) * (1.0-1.0/invPhat[k]);
		//Log(8)<<"\n j x y k p_hat[k] d_g "<<ix->first<<" "<<ix->second<<" "<<data.y()<<" "<<p_hat[data.y()]<<" "<<( k==data.y() ? ix->second : 0 ) - ix->second*p_hat[k];
//                }
	    }
	    // grad(ix->first, k) -= ix->second where k==y
	}
	
	index++;
	
    }
    /*
      lossTerm = LogLikelihoodAndGradient(data,
      modelType,
      w,
      grad
      );
    */
    //Log(9)<<"\nLoss-only Gradient "<<grad;
    Log(0) << "\nCheck 4 "<< Log.time();
    flush(std::cout);
    
    square2flat( grad, gFlat );
    square2flat( hessInvDiagonal, hessianInverseDiagonal);
    if( dupcoord )
        for( unsigned i=0; i<nparams; i++ )
            gFlat[i+nparams] = -gFlat[i];
    Log(0) << "\nCheck 5 "<< Log.time();
    
    
    double penaltyTerm = 0;
    {//---penalty term---
	unsigned i=0;
	for( unsigned j=0; j<w.numFeatures(); j++ )
	    for( unsigned k=0; k<w.numClasses(); k++ )
		if( ! fixedParams(j,k) ) {
		    double s = x[i] - priorMode(j,k);
		    if( Normal==bayesParam.getPriorType() ) {
			penaltyTerm += s * s * penalty[i];
			gFlat[i] += 2 * s * penalty[i];
			hessianInverseDiagonal[i] -= 2 * penalty[i];
		    }
		    else if( Laplace==bayesParam.getPriorType() ) {
			if( ModelType::QuasiNewtonSmooth==modelType.getOptimizerType() ) {
			    double p = sqrt( s*s + esmooth );
			    penaltyTerm += p * penalty[i];
			    gFlat[i] += s * penalty[i] / p;
			    
			}else{ //2-coord
			    if( ! dupcoord ) throw logic_error("dupcoord inconsistent");
			    s = fabs(x[i] + x[nparams+i]);
			    penaltyTerm += s * penalty[i];
			    //Log(0) << "Penalty term " << fabs(x[i] - x[nparams+i]) << " " << s;
			    gFlat[i] += penalty[i];
			    gFlat[nparams+i] += penalty[i];
			}
		    }
		    i++ ;
		}
    }
    Log(0) << "\nCheck 6 "<< Log.time();
    
    
    /*
      FILE* fg = fopen("c:\\out3.csv", "a");
      fprintf(fg, "---\n");
      for ( unsigned i = 0; i < nparams; i++) {
      hessianInverseDiagonal[i] = 1.0/hessianInverseDiagonal[i];
      if (hessianInverseDiagonal[i] < -50.0) {
      hessianInverseDiagonal[i] = -50.0;
      }
      //printf("%f\n", hessianInverseDiagonal[i]); 
      //fprintf(fg, "%d %f %f \n", i, hessianInverseDiagonal[i], gFlat[i]); 
      }
      fflush(fg);
      fclose(fg);
    */
    
    std::cout << "Loss term penalty term " << lossTerm << "/" << penaltyTerm << '\n';
    *f = lossTerm + penaltyTerm;
    
    // update convergence data
    cc.fgeval++;
    cc.f = *f;
    cc.ll = lossTerm;
    
    
//    Log(0)<<"\nFGeval: loss penalty f "<<lossTerm<<" "<<penaltyTerm<<" "<<*f;
//    Log(0)<<"\nGradient"; for(unsigned i=0;i<nparams;i++) Log(0)<<" "<<gFlat[i];
//	Log(0)<<"\n";
    //Log(0)<<"\nGradient at X is"; for(unsigned i=0;i<(dupcoord?2*nparams:nparams);i++) Log(0)<<" "<<gFlat[i];
    /*/suggested by S.Benson:*/
    if( ! finite(*f) ) {
	Log(0)<<"Infinite F\n";
        *f = numeric_limits<double>::max();
        return 1;
    }
    Log(0)<<"Finite F " << f << "\n";
    return 0;
}

//comments from blmvm-1.0\example-c\combustionv1.c:
/* --------------------------------------------------------------------------
 *    BLMVMConverge -- User must implement this routine that tells the solver
 *    the current solution is sufficiently accurate and the solver should 
 *    terminate.
 *
 *    Input:   ctx - Pointer to user-defined structure 
 *             residual - residuals of each variable (array of length n) (=0 at solution).
 *             n - number of variables.
 *             stepsize - change is the solution (2 norm).
 *             whystop - .
 *
 *    Output:  ctx - Pointer to user-defined structure 
 *    residual - unchanged.
 *    whystop - flag set to zero if BLMVM should continue, and nonzero if the solver should stop.
 *    return: 0 for Normal return and nonzero to abort mission.
 *                                                                            */

// New function signature:
//int BLMVMConverge(void*ctx, int iter, double stepsize, double residual[], int n){


int BLMVMConverge( void* ctx, int iter, double stepsize, double* residual, int n )
{
    Context* pCtx = (Context*)ctx;

    int whystop;
    return pCtx->Converge( residual, n, stepsize, &whystop );
}

int Context::Converge( double *residual, int n, double stepsize, int *whystop )
{ 
    if( n != nPseudoParams )
        throw logic_error("Wrong parameter vector size returned from BLMVMConverge");
    
    {
//		Log(0) << "\n";
        double pgnorm=0.0;
	for( int i=0; i<n; i++ ) {
//			Log(0) << residual[i] << " ";
            pgnorm += residual[i]*residual[i];
	}
        cc.pgnorm = sqrt(pgnorm);
    }
    Log(0)<<"\nIter: "<<cc.iter<<", F: "<<cc.f<<",  pgnorm: "<<cc.pgnorm<<", Time "<<Log.time();
    
    double diffsum = 0.0;
    double diffmax = 0.0;
    
    double dotprod = 0.0;
    double normNew = 0.0;
    double normOld = 0.0;
    
    double diff = 0.0;
    double cos = 0.0;
    
    for (unsigned i = 0; i < w.numFeatures(); i++) {
	for (unsigned j = 0; j < w.numClasses(); j++) {
	    diff = (*wCopy)(i, j) - w(i, j);
	    diffsum += fabs(diff);
	    diffmax = max(diffmax, fabs(diff));
	    
	    dotprod += diff * (*lastDiff)(i, j);
	    normNew += diff * diff;
	    normOld += (*lastDiff)(i, j) * (*lastDiff)(i, j);
	}
    }
    cos = dotprod / sqrt(normNew * normOld);
    //printf("%f\n", diffsum);
    
    for (unsigned i = 0; i < w.numFeatures(); i++) {
	for (unsigned j = 0; j < w.numClasses(); j++) {
	    (*lastDiff)(i, j) = (*wCopy)(i, j) - w(i, j);
	    (*wCopy)(i, j) = w(i, j);
	}
    }
    
    
    //FILE *f = fopen("c:\\out.csv", "a");
    printf("%d\t%0.20f\t%f\t%0.20f\t%0.20f\t%0.20f\t%0.20f\n", cc.iter, cc.f, cc.pgnorm, cc.ll, diffsum, diffmax, cos);
    flush(std::cout);
    //fclose(f);
    
    int finished;
    if     ( cc.pgnorm <= cc.pgtol        ) finished=1;
    else if( cc.iter   >= cc.maxiter      ) finished=2;
    else if( cc.fgeval >= cc.maxfgeval    ) finished=3;
    else if( cc.f      <= cc.fmin         ) finished=4;
    else if( cc.iter > 0 && stepsize <= 0 ) finished=5;
    else                                    finished=0;
 
    cc.iter++;
    cc.whystop = *whystop = finished;
    
    return finished;
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

/*------------------------------------------------------------------*/
/*     FUNCTION getProblemSizes                                     */
/*------------------------------------------------------------------*/
/** Define sizes for problem HS15.
 */

static Context *ctx = 0;
//static KTR_context  *kc;


void  getProblemSizes (int *  const  n,
                       int *  const  m,
                       int *  const  nnzJ,
                       int *  const  nnzH)
{
    *n = ctx->NParams();
    *m = 0;
    *nnzJ = 0;
    *nnzH = 0;

    return;
}

/*------------------------------------------------------------------*/
/*     FUNCTION getProblemData                                   */

void /*returns model coeffs thru the last arg*/
QN(
                MemRowSetIterator& data,
                const BayesParameter& bayesParam, //input
                //const ParamMatrix & priorMean,  //input
                //const ParamMatrix & priorScale,  //input
                const ModelType& modelType,  //input
                const class FixedParams& fixedParams,  //input
                ParamMatrix & w,  // i/0: init/result values of model coeffs, 
		IndividualPriorsHolder* indprior
            )
{
    QN_BLMVM(data, bayesParam, /*priorMean, priorScale,*/ modelType, fixedParams, w, indprior);
	//QN_KNITRO(data, bayesParam, priorMean, priorScale, modelType, fixedParams, w);
}
