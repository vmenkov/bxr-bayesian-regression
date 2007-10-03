// author: Bing Bai
// email: bbai@cs.rutgers.edu
// 2006/06/29
#ifndef MULTI_MODEL_H
#define MULTI_MODEL_H


#ifdef _MSC_VER
#define VERSION "0.10, Windows Executable"
  #include <hash_map>
  #include <hash_set>
#include <sys/stat.h>
  #if _MSC_VER >= 1300
    using namespace stdext;
  #endif
#endif //_MSC_VER
#ifdef USE_GCC 
#define VERSION "0.10, Linux Executable"
  #include <ext/hash_map>
  #include <ext/hash_set>
  using namespace __gnu_cxx;
#endif

#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <list>
#include <string.h>
#include <tclap/CmdLine.h>

using namespace TCLAP;
#ifdef _MSC_VER
  #include <direct.h>
  #include <windows.h>
  #include <io.h>
  #ifndef MAXPATHLEN
  #define MAXPATHLEN 256
  #endif
#endif
#ifdef USE_GCC
  #include <dirent.h>
#endif
#include <errno.h>
#include <math.h>

#define KW_betaSparse "betaClassSparse"
#define BBRheader "Bayesian Binary Regression"
#define BPRheader "Bayesian Polytomous Regression"
#define BMRheader "Bayesian Multinomial Regression"
#define MLRheader "Multinomial logistic regression model format"
#define SPSheader "BXRclassify Sparse Model File"

using namespace std;

struct eqstr
{
  bool operator()(const string s1, const string s2) const{
    return s1==s2;
  }
}; 

#ifdef USE_GCC
namespace __gnu_cxx{
  // DJB hashing. 
  // "An algorithm produced by Professor Daniel J. Bernstein and shown first 
  // to the world on the usenet newsgroup comp.lang.c. It is one of the most 
  // efficient hash functions ever published."
  // -- http://www.partow.net/programming/hashfunctions/
  template<> struct hash<string> {
    size_t operator()(const std::string& str) const
    {
      unsigned int hash = 5381;
      for(unsigned int i = 0; i < str.length(); i++){
	hash = ((hash << 5) + hash) + str[i];
      }
      return (hash & 0x7FFFFFFF);
    }  
  };
}
#endif // USE_GCC
#ifdef _MSC_VER  // visual c++ hash_map has different form
class stringhasher : public hash_compare <std::string>
{
public:
  size_t operator()(const string& str) const
  {
      unsigned int hash = 5381;
      for(unsigned int i = 0; i < str.length(); i++){
	hash = ((hash << 5) + hash) + str[i];
      }
      return (hash & 0x7FFFFFFF);
  }  

  bool operator() (const string& s1, const string& s2) const
  {
    return s1 < s2;
  }
};
#endif

typedef struct feature_element{
  int modelID;
  int cateID;
  double weight;
  feature_element(){}
  feature_element(int mid, int cid, double w)
    :modelID(mid), cateID(cid), weight(w){}
} TFE;

typedef struct score_tuple{
  int mid;
  int cid;
  double innprod;
  double odds;
  double prob;
  double logodds; // log of odd
  double logprob; // log probability, called log likelihood in Dave lewis's email
  //  score_tuple(int m, int c, int i=0):mid(m),cid(c), innprod(i){}
  bool used;
} TST;

// artument to pass into model DB constructor
typedef struct arguments{
  bool classic; // -k
  string criterion; // -c
  double threshold; // -t
  string scale; // -s
  bool printTrue; // -r
  bool printPred; // -p
    bool integerclasses;
  string omit; // -o
/*
  string modelType; // -m
  bool outputActiveClassesOnly; // --activeonly
  bool docID; // whether doc ID is included in data lines
  bool modelID; // whether model ID is included in data lines
  bool cateID; // whether category ID is included in data lines
*/
  string labelFile; // -l
  string modelFile;
  string dataFile;
  string outFile;
  string errorFile;
} Args;
  
// true ID, for each document, there could be at most 1 true category
// for each model.
typedef struct struct_tid{
  int modelID;
  int cateID;
  struct_tid(int mid, int cid):modelID(mid), cateID(cid){}
} TID;

typedef map<int, TST> TCM; 

// customized exception for invalid options
// THM is type of hash map
#ifdef USE_GCC
  //typedef hash_map<string, map<int, map<int, double> >, hash<string>, eqstr> THM;

// (commented by Shenzhi)
//  typedef hash_map<string, vector<TFE>, hash<string>, eqstr> THM;

// added by Shenzhi, in testing
typedef hash_map<string, map<pair<int,int>,double>, hash<string>, eqstr> THM;

  // Why this name: type, hash, key, lookup
  typedef hash_map<string, int, hash<string>, eqstr> THKL; 
  // Why this name: type, hash, label, lookup
  typedef hash_map<string, vector<TID>, hash<string>, eqstr> THLL;

  typedef hash_map<int, double> TSM; // the sum of values in each model
  typedef hash_map<int, TCM > TACM; // the accumulator type, for each model 
#endif // USE_GCC
#ifdef _MSC_VER 
#if _MSC_VER >= 1300
  // typedef stdext::hash_map<string, vector<TFE>, stringhasher> THM;
  typedef stdext::hash_map<string, map<pair<int,int>,double>, stringhasher> THM;
  typedef stdext::hash_map<string, int, stringhasher> THKL; 
  typedef stdext::hash_map<string, vector<TID>, stringhasher> THLL;
  typedef stdext::hash_map<int, double> TSM; // the sum of values in each model
  typedef stdext::hash_map<int, TCM > TACM; // the accumulator type, for each model 
#else
  //typedef std::hash_map<string, vector<TFE>, stringhasher> THM;
  typedef std::hash_map<string, map<pair<int,int>,double>, stringhasher> THM;
  typedef std::hash_map<string, int, stringhasher> THKL; 
  typedef std::hash_map<string, vector<TID>, stringhasher> THLL;
  typedef std::hash_map<int, double> TSM; // the sum of values in each model
  typedef std::hash_map<int, TCM > TACM; // the accumulator type, for each model 
#endif
#endif //_MSC_VER
typedef list<TCM::iterator> TMMM;

namespace std{
  template<> struct less<TID>{
    bool operator()(const TID & a, const TID & b) const{
      if(a.modelID < b.modelID) return true;
      return false;
    }
  };
}

// type, vector, key, reverse lookup
typedef vector<string> TVKRL;
// type, vector, feature element
typedef vector<TFE> TVFE; 

class MyException{
  string msg;
public:
  MyException(string & s){
    msg = s;
  }
  MyException(const char * s){
    msg = s;
  }
  string & error(){
    return msg;
  }
};


typedef struct struct_output{
  //  int mid;
  double maxP;
    int maxPCid;
    double maxLinear;
  double maxPPred; // in case printPred is specified.
  int maxPPredCid; 
  vector<TCM::iterator> v;
    bool used;  // added by shenzhi
//    stuct_output() {used = false; } // added by shenzhi
} TOUTPUT;

class ModelDB{

  // Inverted index of features over multiple models and multiple classes.
  // hash table of models. First layer is marked with feature ID, which
  // is a string. 
  THM hmModel;

  // hash map, key: string ID; value: integerID.
  // Model ID and Category ID share the same dictionary
  THKL modelKeys; 
  THKL cateKeys; 

  // vector to keep the integerID ~ stringID map;  
  TVKRL reverseModelKeys;
  TVKRL reverseCateKeys;

  // the table to keep all the labels <mid, cid> pairs for docID
  THLL labelTable;
  int maxlabelnumber;  // keep track the number of <mid,cid>s in labelTable when trigering a bad_alloc exception  
  int totalclasses;

  // when only the info for one docID is needed, keep the <mid, cid> pairs in this map;
  map<int, int> m_trueLabelForOneVector;

  int maxModelID;
  int maxCateID;
  Args args;
  TACM acm;
  TSM sumIP; 

  enum ENUM_Criteria{
      CT_all, CT_prob, CT_logProb, CT_odds, CT_logOdds, CT_linear, CT_maxlinear, CT_maxlinearall,
      CT_maxProbAll, CT_predicted
  };
  
  enum ENUM_Scales{
    SC_prob, SC_logProb, SC_odds, SC_logOdds, SC_linear, SC_all
  };

  // map for <BBRModelID, BBRmodelthresh>, will be filled when reading BBRtrain model files
    // the Thresh will be used when choosing classes: 
    // if the linear score of class "+1" >= Thresh, choose +1; otherwise choose -1; 
  map<int,double> m_BBRModelThresh;

  // error msg output
  ostream *errout;

public:  
  static string criteria[];
  static string scales[];
  static string modeltypes[];

  int getIntID(const string & stringID, THKL & keys, 
	       TVKRL & reverseKeys, int & max_ID);
  string & getModelStringID(int intID);
  string & getCateStringID(int intID);

  // parse model files
  void parseModelFileDir(string dirName); 
  void parseSingleModelFile(string fileName, string & modelID);

  #ifdef _MSC_VER 
  void SearchDirectory(char *pathname);
  #endif

  void parseSparseModel(const string & fileName);

  void addEntry(const string & featureID, const string & modelID,
		const string & cateID, const double weight);
  void addExtraClass(const string & extraClassName);

  ModelDB(const Args & args);
  ~ModelDB();

  // read the whole label file into memory
  void readLabelFile(const string & labelFile);
  void readLabelFile(ifstream &fin, int &lineindex, int &labelnumber, string &upperDocID);

  // read the label file until finding the upperDocID, keep all (mid,cid)pairs in the multimap.
  void cntReadLabelFile(ifstream &fin, int &lineindex, string &upperDocID);

  bool readLine(ifstream &fin, int &lineindex, int &labelnumber, string &upperDocID, string &currentDocID, bool &found);

  void getProbabilities(istringstream & sstr, TMMM & v);
  void processData();
  static int getIntCriterion(string & crit);
  static int getIntScale(string & sc);

  void sparseOutput(TMMM & v, string &docID, string &trueLabel, int crit, int sc, ostream * pout);

  void updateOutputVector(TCM::iterator & tit, int crit, map<int, TOUTPUT> & vd);
  void clearEntryMarks(TMMM & v);
  static void parseArgs(int argc, char ** argv, Args & args);

  void labelFileErrorCheck2(TACM::iterator itACM, int cid, string &upperModelID, string &upperCateID, string &upperDocID);

  // static check whether arguments are good
  // no return value, but exception throwing
  static void sanityCheck(Args & args);
  void postSanityCheck();

  string toUpperCase(const string &s);  

  bool intvalue(const string& s);

  char *emergencypool;

  bool lowmemory;

  bool bDir;
};
#endif //MULTI_MODEL_H
