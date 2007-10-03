// updated:
// 2.60     Oct 05, 05  hyperparameter search without a grid

#ifndef _BAYES_PARAMETER_HPP
#define _BAYES_PARAMETER_HPP

#include <ostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
using namespace std;

enum PriorType { Jeffreys=0, Laplace=1, Normal=2, hyperbolic=3, priorUndef=4 };

inline double var2hyperpar( PriorType priorType, double skew, double var ) {
    return Laplace==priorType ? 0==skew ? sqrt(2/var) : sqrt(1/var)
        : var; }
inline double hyperpar2var( PriorType priorType, double skew, double hyperpar ) {
    return Laplace==priorType ? 0==skew ? 2/hyperpar/hyperpar : 1/hyperpar/hyperpar
        : hyperpar; }

class BayesParameter {
public:
    //access members
    enum PriorType getPriorType() const { return m_prior; }
    double getPriorVar() const { 
        return hyperpar2var( m_prior, m_skew, m_hyperpar ); }
    double getGamma() const {
        if( m_prior!=Laplace ) throw logic_error("Irrelevant Laplace parameter requested");
        return m_hyperpar; }
    double getVar() const {
        if( m_prior!=Normal ) throw logic_error("Irrelevant Normal parameter requested");
        return m_hyperpar; 
	}
    double getQ() const {
        if( m_prior!=hyperbolic ) throw logic_error("Irrelevant hyperbolic parameter requested");
        return m_hyperpar; }
    double getSkew() const { return m_skew; }
    bool isValid() const { return m_validPar; } // hyperparameter value is set
    //ctors
    BayesParameter() : m_validPar(false), m_prior(priorUndef) {} //default ctor
    BayesParameter(enum PriorType priorType_, double param_, double skew_ )
        : m_prior(priorType_), m_validPar(true)/*autoHpar(false)*/, m_skew( skew_>0? 1 : skew_<0? -1 : 0 )
    {
        m_hyperpar = param_;
        is_consistent();
    }
private:
    bool m_validPar; // hyperparameter value is set
    enum PriorType m_prior;
    double m_hyperpar;
    //double m_hyperpar2; hierarchical
    double m_skew;
    void is_consistent()
    {
        if( Laplace==m_prior ) {
            if( m_validPar && m_hyperpar<0 )
                throw std::runtime_error("Laplacian gamma wrong or undefined");
        }else if( Normal==m_prior ) {
            if( m_validPar && m_hyperpar<0 )
                throw std::runtime_error("Normal variance wrong or undefined");
        }else if( hyperbolic==m_prior ) {
            if( m_validPar && m_hyperpar<=0 )
                throw std::runtime_error("Hyperbolic \'q\' wrong or undefined");
        }else if( Jeffreys==m_prior ) {
            ;
        }else
            throw std::runtime_error("Undefined prior type");
    }

    friend class HyperParamPlan;
    BayesParameter(enum PriorType priorType, double skew_ ) // auto hpar
        : m_prior(priorType), m_validPar(false)/*autoHpar(true)*/, m_skew( skew_>0? 1 : skew_<0? -1 : 0 ) 
    { is_consistent(); }
};

inline std::ostream& operator<<( std::ostream& o, const BayesParameter& bp ) {
    o << "PriorType: ";
    if( bp.getPriorType()==Laplace ) {
        o<<"Laplace\tLambda="<<bp.getGamma()<<" Prior var="<<hyperpar2var(bp.getPriorType(),bp.getSkew(),bp.getGamma());
    }else if( bp.getPriorType()==Normal ){
        o<<"Normal diagonal\tVar="<<bp.getVar();
    }else if( bp.getPriorType()==hyperbolic ){
        o<<"Hyperbolic\tq="<<bp.getQ();
    }else if( bp.getPriorType()==Jeffreys )
        o<<"Jeffreys";
    else
        o<<"<undefined>";
    if( bp.getSkew() != 0.0 )
        o<<"\t Skew="<<bp.getSkew();
    return o;
}

class HyperParamPlan {
    BayesParameter m_bp;
    std::vector<double> m_plan;
    unsigned m_nfolds, m_nruns;
    bool errBarRule;
    bool m_search; //search for hyperpar with cv without a grid
    void sortPlan() {   //sort hpar plan by penalty decreasing
        if( ! hasGrid() )   return;
        if(Laplace==PriorType()) 
            sort( m_plan.rbegin(), m_plan.rend() );
        else
            sort( m_plan.begin(), m_plan.end() );
    }
public:
    enum PriorType getPriorType() const { return m_bp.getPriorType(); }
    double getSkew() const { return m_bp.getSkew(); }
    enum HyperParMode { Native, AsVar };
    const BayesParameter& getBp() const { return m_bp; }

    unsigned getPlansize() const { return m_plan.size(); }
    const std::vector<double>& getPlan() const { return m_plan; }
    double getPlanAsVar(unsigned i) const { return hyperpar2var( getPriorType(), getSkew(), m_plan.at(i) ); }

    unsigned getNfolds() const { return m_nfolds; }
    unsigned getNruns() const { return m_nruns; }

    bool isValid() const { return PriorType()<priorUndef; }
    bool getNormDefault() const { return isValid() && (m_plan.size()==0) && !m_search; } //norm-based default
    bool getAutoHPar() const { return getNormDefault(); } //this alias should fade out
    bool hasGrid() const { return isValid() && m_plan.size()>=2; }
    bool isFixed() const { return isValid() && m_plan.size()==1; }
    bool isSearch() const { return isValid() && m_search; }

    bool isErrBarRule() const { return hasGrid() && errBarRule; }

    BayesParameter getInstBayesParam(unsigned iplan) const { 
        return BayesParameter( m_bp.m_prior, m_plan.at(iplan), m_bp.getSkew() ); }
    BayesParameter getInstByVar(double var) const { 
        return BayesParameter( m_bp.m_prior, 
            var2hyperpar( m_bp.m_prior, m_bp.getSkew(), var ), 
            m_bp.getSkew() ); 
    }

    //ctor's
    HyperParamPlan() {} //default - disabled, since prior type is undef
    HyperParamPlan( enum PriorType priorType, double skew, std::string strHP, std::string strCV, 
        HyperParMode hyperParMode = Native,
        bool errBarRule_=false)
        : errBarRule(errBarRule_), m_search(false)
    {
        m_nruns = m_nfolds = 0;
        char comma;
        {
	    /*
            std::istringstream istr(strHP);
            double d;
            while( istr>>d ) {
                m_plan.push_back( Native==hyperParMode ? d : var2hyperpar(priorType,skew,d) );
                istr>>comma;
            }
	    */

            char* var = strtok((char*)strHP.c_str(),", ");
            while( var ) {
		double d;
                if(var=="inf")
		    d = std::numeric_limits<double>::infinity();
                else 
                    d = atof(var);
		m_plan.push_back( Native==hyperParMode ? d : var2hyperpar (priorType,\
									       skew,d));
                var = strtok(NULL,", ");
            }
	
        }
        {
            std::istringstream istr2(strCV);
	    double d;
            istr2>>d>>comma;
	    m_nfolds = (unsigned)d;
	    if( !(istr2>>d) ) {
                m_nruns = m_nfolds;
	    } else {
		m_nruns = (unsigned)d;
	    }
	    if( m_nfolds<2 || m_nruns>m_nfolds )
                throw std::runtime_error(std::string("Wrong cross-validation parameters: ")+strCV);
        }
        if( m_plan.size()==1 )
            m_bp = BayesParameter( priorType, m_plan[0], skew );
        else
            m_bp = BayesParameter( priorType, skew ); // auto hpar
	
        if( m_plan.size()>1 ) {
            sortPlan();
            Log(6)<<"\nHyperparameter plan sorted by penalty decreasing:\n";
            for( unsigned i=0; i<m_plan.size(); i++ ) Log(6)<<" "<<m_plan[i];
        }
    }
    HyperParamPlan( enum PriorType priorType, double hypval, double skew_  ) : m_search(false) { //fixed hyperpar
        m_plan.push_back( hypval );
        m_bp = BayesParameter( priorType, m_plan[0], skew_ );
    }
    HyperParamPlan( enum PriorType priorType, double skew_  ) : m_search(false) { //auto-select norm-based hyperpar
        m_bp = BayesParameter( priorType, skew_ );
    }
    HyperParamPlan( enum PriorType priorType, double skew, std::string strCV ) //search without a grid
        : m_search(true)
    {
        m_nruns = m_nfolds = 0;
        char comma;
        {
            std::istringstream istr(strCV);
            istr>>m_nfolds>>comma;
            if( !(istr>>m_nruns) )
                m_nruns = m_nfolds;
            if( m_nfolds<2 || m_nruns>m_nfolds )
                throw std::runtime_error(std::string("Wrong cross-validation parameters: ")+strCV);
        }
        m_bp = BayesParameter( priorType, skew );
    }
};

inline std::ostream& operator<<( std::ostream& o, const HyperParamPlan& p ) {
    //o<<p.bp();
    o << "PriorType: ";
    if( p.getPriorType()==Laplace ) {
        o<<"Laplace\tLambda=";
        if(p.isFixed()) o<<p.getBp().getGamma();
        else o<<" auto-select";
    }else if( p.getBp().getPriorType()==Normal ){
        o<<"Normal diagonal\tVar=";
        if(p.isFixed()) o<<p.getBp().getVar();
        else o<<" auto-select";
    }else if( p.getBp().getPriorType()==hyperbolic ){
        o<<"Hyperbolic\tq=";
        if(p.isFixed()) o<<p.getBp().getQ();
        else o<<" auto-select";
    }else if( p.getBp().getPriorType()==Jeffreys )
        o<<"Jeffreys";
    else
        o<<"<undefined>";
    if( p.getBp().getSkew() != 0.0 )
        o<<"\t Skew="<<p.getBp().getSkew();
    o << "\nHyperParameter";
    if( p.hasGrid() ) {
        o<<" plan:";
        for( unsigned i=0; i<p.getPlan().size(); i++ )
            o<<" "<<p.getPlan()[i];
        o<<"\tCross-validation: "<<p.getNfolds()<<" folds, "<<p.getNruns()<<" runs";
    }
    else if( p.getNormDefault() )
        o<<": Norm-based default";
    else if( p.isSearch() )
        o<<": Search"<<"\tCross-validation: "<<p.getNfolds()<<" folds, "<<p.getNruns()<<" runs";
    else
        o<<": Fixed";
    if( p.hasGrid() && p.isErrBarRule() )
        o<<"\nOne-standard error rule in effect";
    return o;
}

#endif //_BAYES_PARAMETER_HPP

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
