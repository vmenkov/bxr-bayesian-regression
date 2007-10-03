#ifndef DATA_FACTORY_
#define DATA_FACTORY_

#ifdef _MSC_VER
  #include <hash_map>
  #include <hash_set>
  #include <sys/stat.h>
  #if _MSC_VER >= 1300
    using namespace stdext;
  #endif
#endif //_MSC_VER
#ifdef USE_GCC 
  #include <ext/hash_map>
  #include <ext/hash_set>
  using namespace __gnu_cxx;
#endif

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <map>

#include <iostream>
#include <fstream>
#include <set>
#include <stdlib.h>
#include "Data.h"
#include "IndividualPriors.h"

using namespace std;


class DataFactory  {

    string trainFName;
    string testFName;
    string indpriorFName;
    string m_refClass;
    int m_refClassId;

    vector<SparseVector> m_Mtrain;

    SIHM featmap;  // to keep the <feature str, feature id> hash_map if SYMBOLIC is defined
    Vec m_wordRestrict;
    vector<int> m_wordIndex;  // to keep the number of appearances of each feature in m_featSelect
    Vec m_featSelect, m_addFeats;  // to keep the actual unique features (string or integer, depends on SYMBOLIC)

    vector<YType> m_yTrain;  // keep the index of the class 
    Vec m_classes, m_addClasses;  // to keep the actual unique class labels
    SIHM clsmap;  // to keep the <class str, class id> map if SYMBOLIC is defined

    map<int, IndividualPriorsMap> m_byClass;  // coefficent-level
    IndividualPriorsMap m;   // feature level
    
    // keep the modes for those non active class/feats
    map<MyType,double> nonactiveFeats;
    map<MyType, map<MyType,double> > nonactiveCoeflevelFeats;

    MemRowSet* trainData;
    PlainYRowSet* testData;
    IndividualPriorsHolder* indpriorHolder;
    RowSetMetadata* m_metadata;

 public:

    DataFactory(bool symbolic_) {
    }

    void setTestAndTrainFileSpec( 
        string trainFile_,
        string testFile_,
	string indpriorFile_,
	string strRefClassId_ // ver3.13
        );

    void readFiles();

    MemRowSet* const getTrainData() {
	return trainData;
    }

    PlainYRowSet* const getTestData() {
	return testData;
    }

    IndividualPriorsHolder* const getIndPriorHolder() {
	return indpriorHolder;
    }

    // should be called after MakeSparseData and createIndPriors
    const RowSetMetadata& getRowSetMetadata() const {
	return *m_metadata;
    }

    
 private:

    MemRowSet* createTrainRowSet();
    PlainYRowSet* createTestRowSet();
    IndividualPriorsHolder* createIndPriorHolder();
    
    void MakeSparseMatrix( string fName);
    void relabelFeatures();
    void collectRecodeClasses( const vector<int>& clids );
    void addIntercept();

    int getIntID(const string& stringID, SIHM& keys, Vec& reversekeys, int& maxid);    

    ModeVarSkew _parsePriorsLine(istringstream& s);
};

#endif //DATA_FACTORY_

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
