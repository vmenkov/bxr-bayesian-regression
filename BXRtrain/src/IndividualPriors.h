 // TODO: Rename IndPriors -> IndividualPriors
// TODO: Fix capitalization outside file
// TODO: what does iVar mean?
// TODO: Rename IndividualPriorsMap m variable
// TODO: Rename this file to something consistent with IndividualPriors

#ifndef INDIVIDUALPRIORS_H
#define INDIVIDUALPRIORS_H

#include <map>
#include <vector>
#include <set>
#include <string>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "Compiler.h"
#include "Data.h"

using namespace std;

struct ModeVarSkew
{
    double mode, var, skew;
    bool abs;

    ModeVarSkew(): mode(0.0), var(1.0), skew(0.0), abs(false) {}
    ModeVarSkew(double mode_, double var_, double skew_, bool abs_): mode(mode_), var(var_), skew(skew_), abs(abs_) {}

    static ModeVarSkew* INVALID() {
	return new ModeVarSkew(0, -1, 0, false);
    };
};

typedef map <unsigned, ModeVarSkew> IndividualPriorsMap;


class IndividualPriorsHolder
{
    // the values are already mapped
    const map<int, IndividualPriorsMap>& m_byClass;  // coefficent-level
    const IndividualPriorsMap& m;   // feature level

    // iterators for the above 2 maps; used when calculating the final variance in ZOLR()
    mutable IndividualPriorsMap::const_iterator fiter, subfiter;
    mutable map<int, IndividualPriorsMap>::const_iterator citer;

    string m_file;

 public:

    bool isPriorvar0CoefLevel(int fid, int cid) const {
	map<int,IndividualPriorsMap>::const_iterator firstlevel = m_byClass.find(cid);
	if (firstlevel == m_byClass.end()) return false;
	IndividualPriorsMap::const_iterator seclevel = firstlevel->second.find(fid);
	if (seclevel == firstlevel->second.end()) return false;
	if (seclevel -> second.var == 0) return true;
	return false; 
    }

    // set the iterators to the beginning of the maps; SL Jul 07
    void reset() const {
	// iterator for the feature-level data 
	fiter=m.begin();

	// iterator for the coefficient-level map
	
	citer = m_byClass.begin();
	if(m_byClass.size()>0)    subfiter = citer->second.begin();
    }

    // SL Jul 07
    // check for the <feature,class> pair, whether ind prior file has specifications;
    // type: 0 - not found; 1 - feature-level;
    //                      2 - coefficient-level; 
    // NOTE: when iterating through the maps, we could also get the non-active features with non-zero modes
    ModeVarSkew hasIndPrior(int f, int c, int& type) {
	ModeVarSkew p(0,1,0,false);
	bool found;

	// first check the coefficient-level
	p=_checkCoefLevel(m_byClass, citer, subfiter, c, f, found);
	if(found) { type = 2; return p; }

	// then, check the feature-level
	p=_checkFeatLevel(m, fiter, f, found);
	if(found) { type = 1; return p; }

	// if nothing found, 
	type = 0;   return p;
    }
                       
    
    bool valid(const ModeVarSkew &p) const
    {
        return p.var >= 0;
    } 

	
    bool active() const {
	return m.size() > 0 || m_byClass.size()>0;
    }
    
    IndividualPriorsHolder(
	const map<int, IndividualPriorsMap>& byClass_, 
	const IndividualPriorsMap& m_, 
	const string filename_)
	: m_byClass(byClass_), m(m_), m_file(filename_)
	{}

    string individualPriorsFile() { return m_file; }

 private:


    // SL Jul 07
    // if current_class < c, go to next entry
    // if current_class > c, do nothing
    // if current_class = c, check the features
    // if the class or the (class, feature) not found, go check scale value coefficient-level
    ModeVarSkew _checkCoefLevel(const map<int,IndividualPriorsMap>& cmap,
			 map<int,IndividualPriorsMap>::const_iterator& citer, 
			 IndividualPriorsMap::const_iterator& fiter,
			 int c, int f, bool& found) 
	{
	    found = false; bool reset = false;
	    ModeVarSkew p(0,1,0,false);
	    while(citer!=cmap.end() && citer->first < c) {
		citer++;  reset = true;
	    }

	    if(citer==cmap.end()) return p;

	    if (reset)	    fiter = citer->second.begin(); 	

	    if(citer->first == c ) {
		while (fiter->first < f && fiter!=citer->second.end() ) fiter++;
		if(fiter->first==f && fiter!=citer->second.end() ) {
		    p = fiter->second;  found = true;  return p;
		}
	    }
	    return p;
	}


    // SL Jul 07
    // if current_feat < f, go to next entry; NOTE: this feature is non-active feature in training examples; 
    // if current_feat > f, do nothing
    // if current_feat = f, return the ModeVarSkew struct
    ModeVarSkew _checkFeatLevel(const IndividualPriorsMap& fmap,
				IndividualPriorsMap::const_iterator& fiter, 
				int f, bool& found) 
	{
	    found = false;
	    ModeVarSkew p(0,1,0,false);

	    while(fiter!=fmap.end() && fiter->first < f) 
		fiter++;

	    if(fiter!=fmap.end() && fiter->first==f) {
		p = fiter->second;  found = true;  return p;
	    }
	    
	    return p;
	}
    
};


inline std::ostream &operator << (std::ostream &o, const IndividualPriorsHolder &pt)
{
    return o;
}



#endif // INDIVIDUALPRIORS_H

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
