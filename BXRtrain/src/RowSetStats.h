//  v.2.04a     May 30, 05  stddev is zero when it seems that meansqu < mean^2

#ifndef ROW_SET_STATS_H
#define ROW_SET_STATS_H

#include "Data.h"

class RowSetStats {
    IRowSet &m_source;
    
    bool m_ready; // lazy initialization

	Vector m_stddevs;
    Vector m_means;
    double m_avgSquNorm;
    unsigned m_nrows;

public:
    RowSetStats(  IRowSet & drs_ ) : m_source(drs_), m_ready(false) {}

    const unsigned NRows() { 
        if( !m_ready ) getready();
        return m_nrows; 
    }

    const Vector& Means() { 
        if( !m_ready ) getready();
        return m_means; 
    }

    const Vector& Stddevs() { 
        if( !m_ready ) getready();
        return m_stddevs; 
    }

    double AvgSquNorm() {  //for default bayes parameter: avg x*x
        if( !m_ready ) getready();
        return m_avgSquNorm; 
    }

private:
    void getready() { //laziness
        m_means.resize( m_source.dim(), 0.0 );
        m_avgSquNorm = 0.0;
        vector<double> meansqu( m_source.dim(), 0.0 );
        vector<unsigned> xn( m_source.dim(), 0 ); //non-zeroes

        Log(8)<<"\nCompute stddevs start Time "<<Log.time();

        // compute m_means and stddevs for non-zeroes only first
        RowSetIterator* it = m_source.iterator();
        unsigned r;
        for( r=0; it->next(); r++ )
        {
            const SparseVector& x=it->xsparse();
            double squNorm = 0;
            for( SparseVector::const_iterator xitr=x.begin(); xitr!=x.end(); xitr++ )
            {
                unsigned iFeat = xitr->first; // m_source.colNumById( xitr->first );
                double val = xitr->second;
                double fNew = 1.0 /( xn[iFeat] + 1 );
                double fPrev = xn[iFeat] * fNew;
                m_means[iFeat] = m_means[iFeat] * fPrev + val * fNew;
                meansqu[iFeat] = meansqu[iFeat] * fPrev + val*val * fNew;
                xn[iFeat] ++;
                squNorm += val*val;
                Log(16)<<"- "<<iFeat<<" "<<val<<" "<<m_means[iFeat]<<" "<<meansqu[iFeat]<<" "<<squNorm<<endl;
            }
            m_avgSquNorm = m_avgSquNorm*r/(r+1) + squNorm/(r+1);
            //Log(6)<<"\nm_avgSquNorm "<<m_avgSquNorm;
        }
		delete it;
        //Log(5)<<"\nr "<<r;

        //now adjust for zeroes
        m_stddevs.resize( m_source.dim(), 0.0 );
    
	for( unsigned f=0; f<m_source.dim(); f++ ) {
            double adjust = double(xn[f]) / double(r);
            m_means[f] *= adjust;
            meansqu[f] *= adjust;
            double stddev_squ = meansqu[f] - m_means[f]*m_means[f];
            m_stddevs[f] = stddev_squ<=0 ? 0 : sqrt( stddev_squ );
        }
	
        m_nrows = r;
        m_ready = true;
    }
};

#endif // ROW_SET_STATS_H
