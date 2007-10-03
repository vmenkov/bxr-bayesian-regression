// This header includes code to handle compiler dependencies.

#ifndef _PARAM_MATRIX_H
#define _PARAM_MATRIX_H

#include "IndividualPriors.h"

class ParamMatrix {
    vector< vector<double> > m_matrix; //replace by sparse row storage when needed
    unsigned m_features;
    unsigned m_classes;
public:

    ParamMatrix() {};
    ParamMatrix( unsigned features_, unsigned classes_, double v=0.0)
        : m_features(features_), m_classes(classes_)//, m( c, vector<double>(d,v) )  
    {
        m_matrix = vector< vector<double> >( m_classes, vector<double>(m_features,v) );
    }

    unsigned numFeatures() const { return m_features; }
    unsigned numClasses() const { return m_classes; }

    inline double& operator()(unsigned j, unsigned k) {
//		if( j>=m_features || k>=m_classes )  throw DimensionConflict(__FILE__,__LINE__);
        return m_matrix[k][j];
    }

    inline double operator()(unsigned j, unsigned k) const {
//        if( j>=m_features || k>=m_classes )  throw DimensionConflict(__FILE__,__LINE__);
        return m_matrix[k][j];
    }

    const vector<double>& classparam( unsigned k ) const { return m_matrix[k]; }
    
    bool active(unsigned j, unsigned k) const { return true; } //temp

    void reset( unsigned features_, unsigned classes_, double v=0.0)
	{
	    m_features = features_; m_classes = classes_;
	    m_matrix = vector< vector<double> >( m_classes, vector<double>(m_features,v) );
	}
    
    unsigned countequal( const ParamMatrix& pm2 ) const {
	unsigned nz=0;
	for( unsigned j=0; j<this->numFeatures(); j++ ) {
	    for( unsigned k=0; k<pm2.numClasses(); k++ ) {
		if( (*this)(j,k) == pm2(j,k) ) nz++;
	    }
	}
	return nz;
    }
    
    void displayModel( ostream& o, const IRowSet& datasource ) const
	{
	    o<<std::endl<<std::setw(20)<<"Feature\\Class";
	    for( unsigned k=0; k<this->numClasses(); k++ )
		o<<std::setw(20)<<k;
	    for( unsigned j=0; j<this->numFeatures(); j++ ) {
		o<<std::endl<<std::setw(20)<<j; // datasource.colName(j);
		for( unsigned k=0; k<this->numClasses(); k++ )
		    o<<std::setw(20)<<(*this)(j,k);
	    }
	}
    
    
    
    
};

inline std::ostream& operator<<( std::ostream& s, const ParamMatrix& m ) {
    s<<endl<<m.numFeatures()<<" "<<m.numClasses()<<endl;
    for( size_t j=0; j<m.numFeatures(); j++ ) {
        s<<"["<<j<<"]: ";
        for( size_t k=0; k<m.numClasses(); k++ ) 
            s<<" "<<m(j,k);
	s<<endl;
    }
    return s;
}

class FixedParams {
    mutable vector< vector<bool> > m_locked;  //jk;  ===> change this
    mutable int nfixed;

public:

    int dim() const {
	if(active())
	    return m_locked.at(0).size();
	else 
	    return 0;
    }

    int c() const {
	return m_locked.size();
    }

    bool active() const { return m_locked.size()>0; }

    bool operator() (unsigned j, unsigned k) const {  // j-->feature; k-->class
        if( k>=m_locked.size() ) return false; // btw takes care if intercept is unaccounted for
        if( j>=m_locked[k].size() ) return false;
        return m_locked[k][j];
    }

    const vector<bool>& classparam( unsigned k ) const { return m_locked[k]; }

    FixedParams( const vector< vector<bool> >& allzeroes, 
		 IndividualPriorsHolder* indpriorholder,
		 int referenceClassId,
		 unsigned featnum, unsigned clsnum) 
        : nfixed(-1), m_locked(allzeroes)
    {
	if(indpriorholder!=0) {
	    if(!active()) 
		m_locked = vector<vector<bool> > (clsnum,vector<bool>(featnum,false));
	    indpriorholder->reset();
	
	    for(int j=0; j<clsnum ; j++) {
		for(int i=0; i<featnum ; i++) {	
		    // level 6 && 4
		    int type;
		    ModeVarSkew level6and4 = indpriorholder->hasIndPrior(i,j,type);
		    if(type == 2 )
			m_locked[j][i] = level6and4.var == 0 ? true : false;
		    // level 4 
		    else if (type==1 && level6and4.var==0)
			m_locked[j][i] = true;
		}
	    }
	}
    }

    unsigned count() const {
        if( nfixed>=0 ) return nfixed; //already set
        nfixed=0;
        if( !active() ) ;
        else
            for( unsigned j=0; j<m_locked.size(); j++ )
                for( unsigned k=0; k<m_locked[j].size(); k++ )
                    if( (*this)(k,j) ) nfixed++;
        return nfixed;
    }
};


inline std::ostream& operator<<( std::ostream& s, const FixedParams& m ) {
    s<<endl<<m.dim()<<" "<<m.c()<<endl;
    for( size_t j=0; j<m.dim(); j++ ) {
        s<<"["<<j<<"]: ";
        for( size_t k=0; k<m.c(); k++ ) 
            s<<" "<<m(j,k);
	s<<endl;
    }
    return s;
}

#endif // _PARAM_MATRIX_H
