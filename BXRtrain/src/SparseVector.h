
#ifndef SPARSE_VECTOR_H
#define SPARSE_VECTOR_H

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <limits>
#include <stdexcept>

#include "Matrix.h"
#include "logging.h"
#include "Compiler.h"
using namespace std;

typedef pair<int,double> SparseItem;

class SparseVector : public vector< SparseItem > {
public:
    typedef vector< SparseItem >::const_iterator const_iterator;
    typedef vector< SparseItem >::iterator iterator;

    //
    SparseVector() {}
    
    SparseVector( const vector<SparseItem>& v ) : vector<SparseItem>( v ) {
        std::sort( begin(), end() );
    }
	
    //
    const_iterator begin() const { return vector<SparseItem>::begin(); }
    const_iterator end() const { return vector<SparseItem>::end(); }

    iterator begin() { return vector<SparseItem>::begin(); }
    iterator end() { return vector<SparseItem>::end(); }
    
    //
    void sort() {
        std::sort( begin(), end() );  
    }
    
    void insert( int i, double d ) {
	push_back( SparseItem( i, d ) ); 
    }

    inline const_iterator find( int i ) const { 
	SparseItem pattern( i, - numeric_limits<double>::max() );
        vector< SparseItem >::const_iterator found = lower_bound( begin(), end(), pattern );
	if( found != end() && found->first==i ) {
            return found;
	} else {
            return end();
	}
    }
    
    const double dot(const vector<double>& b) const {
	double score = 0.0;
	for( SparseVector::const_iterator ix=begin(); ix!=end(); ix++ )
	    score += b[ix->first] * ix->second;
	return score;
    }
  
    // SL ver3.02
    void resize() {
	reserve(size()+1);
    }  
    
};
/*class SparseVector : private map<int,double> {};*/
inline std::ostream& operator<<( std::ostream& o, const SparseVector& sv ) {
    for( SparseVector::const_iterator isv=sv.begin(); isv!=sv.end(); isv++ )
        o << isv->first <<":"<<isv->second<<"  ";
    return o;
}

/*
inline double TFWeight(double rawTF, TFMethod tf ) {
    if (tf == RAWTF) {
        return (rawTF);
    } else if (tf == LOGTF) {
        return (rawTF > 0) ?  log(rawTF) +1 : 0;
    } else {  // default to raw TF
        throw logic_error("unknown TF method");
    }
}
*/

#endif //SPARSE_VECTOR_H

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
