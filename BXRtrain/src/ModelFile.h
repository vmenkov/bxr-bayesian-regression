// v3.21 remove standardize flag in model file

#ifndef STORED_MODEL_
#define STORED_MODEL_

#include <string>
#include <fstream>
#include <vector>
#include <set> //ver 3.04

#include "Compiler.h"
#include "DataFactory.h"
#include "BayesParamManager.h"
#include "ModelTypeParam.h"
#include "Data.h"
#include "Design.h"
#include "ParamMatrix.h"

using namespace std;


#define LEGACYTITLE "Bayesian Polytomous Regression model file ver "
static const double legacyVersionId = 2.5;

static const double versionId = 3.1;
#define TITLE1 "Model file produced by BMR train version "
#define TITLE2 "Multinomial logistic regression model format: sparse 3.0 "
#define TITLE3 "produced by BXRtrain version "

static const double bbrversionId = 2.0;
#define BBRTITLE "Bayesian Binary Regression ver "
#define SPARSETITLE "Sparse Model File ver "
static const double sparseVersionId = 1 ;

// key words
#define KW_endofheader "endofheader"
#define KW_endofbeta "endofbeta"
#define KW_beta "beta"
#define KW_classes "classes"
#define KW_nDesignParams "nDesignParams"
#define KW_featRestrict "featRestrict"
#define KW_featRestrictRange "featRestrictRange"
#define KW_idf "idf"
#define KW_tfMethod "tfMethod"
#define KW_idfMethod "idfMethod"
#define KW_cosNorm "cosineNormalize"

#define KW_topic "topic"
#define KW_modelType "modelType"
#define KW_design "design"
#define KW_topicFeats "topicFeats"
#define KW_betaDense "betaDense"
#define KW_betaClassSparse "betaClassSparse"
#define KW_intercept "intercept"
#define KW_threshold "threshold"


// error msg when reading init file
#define A8 "The optimization starting point file,"
#define A81 "contains coefficients for features which have no nonzero values in the training data and no nonzero modes in the individual priors file.  This is allowed and will not cause any difficulties for training, but be aware that these features will not occur in the trained model."


class InitModel {
    bool m_active;
    string m_filename;

    vector<int> m_feats;

    // the following two maps are used for original -> index exchange when NON-SYMBOLIC
    map<int,unsigned> m_featmap;
    map<int,unsigned> m_clsmap;

 public:
    //access
    bool isActive() const { return m_active; }

    const vector<int>& getFeatures() const { return m_feats; }

    InitModel(string filename_,const RowSetMetadata& metadata, 
	      const BayesParameter& bayesparam, 
	      const FixedParams& fixedparam, 
	      const IndividualPriorsHolder& indprior,
	      ParamMatrix& beta )   : m_filename(filename_)
    {
        m_active = 0<m_filename.size();
        if(!m_active) return;         //--->>--

        ifstream m_file(m_filename.c_str());
        if( !m_file.good() )
            throw runtime_error(string("Cannot open model file: ")+m_filename);

	// build the original class/feature -> index map 
	// for SYMBOLIC, it is already built in hash_map in metadata;
	#ifndef SYMBOLIC
	map<int, unsigned> featmap, clsmap;
	_buildmap(metadata.getFeatures(), m_featmap);
	_buildmap(metadata.getClasses(), m_clsmap);
	#endif

        string buf;
        //only read the beta lines; if bbrmode, read the topicFeats also;
        while( getline( m_file, buf ) ) { //---line-wise---
            istringstream rowbuf( buf );
            string kw;
            rowbuf>>kw;
            
	    if(kw==KW_topicFeats ) {
		string str;
                while( rowbuf>>str ) 
                    m_feats.push_back( _getFeatIndex(str,metadata) );
		
                if( !rowbuf.eof() ) //the only good reason to end the loop
                    throw runtime_error("Corrupt model file: '" KW_topicFeats "'");
	    }
            else if( kw==KW_beta ) {

		// check class id, active and -R
		int classid = _activeClsId("-1", metadata);
		if(classid==-1) return;

		// if the class -1 is active, read in the beta
		double d;
		vector<double> localbeta;
                while( rowbuf>>d ) 
		    localbeta.push_back( d );

		// the intercept
		beta(0,classid) = localbeta.at(localbeta.size()-1);

		// overwrite the beta values
		bool hasunknownfeats = false;
		for(int i=0; i<m_feats.size(); i++) {
		    if(m_feats.at(i)!=-1)
			_checkbetavalue(m_feats.at(i),classid, localbeta.at(i), bayesparam, 
					metadata, fixedparam, indprior, beta);
		    else
			hasunknownfeats = true;
		}

		// if there are unactive features, issue warning
		if(hasunknownfeats)
		    cout<<"WARNING: "<<A8<<m_filename<<" , "<<A81<<endl;
 
                if( !rowbuf.eof() ) //the only good reason to end the loop
                    throw runtime_error("Corrupt model file: 'beta'");
            }
            else if( kw==KW_betaClassSparse ) {// ver 2.5
		// first read the class id
		string str;
                rowbuf>>str;
		int classid = _activeClsId(str, metadata);
		if (classid==-1) continue;

                string featpair;
		bool hasunknownfeats = false;
	
                while( rowbuf>>featpair) {
		    if( rowbuf.fail() ) 
			throw runtime_error("Corrupt model file: " + buf);

		    string feat; double d;
		    string::size_type loc = featpair.find( ":", 0 );
		    if( loc != string::npos ) {
			feat = featpair.substr(0, loc);
			d = atof(featpair.substr(loc+1).c_str());
		    } else 
                        throw runtime_error(string("missing colon in line ")+ buf);

		    int featid = _getFeatIndex(feat,metadata);
		    if(featid!=-1)
			_checkbetavalue(featid,classid, d, bayesparam, metadata, fixedparam, indprior, beta);
		    else if (!hasunknownfeats)
			hasunknownfeats = true;
                }

		// if there are unactive features, issue warning
		if(hasunknownfeats)
		    cout<<"WARNING: "<<A8<<m_filename<<" , "<<A81<<endl;
	    }

            if( !rowbuf.good() && !rowbuf.eof() )
                throw runtime_error(string("Corrupt model file, line: ")+buf);
	
        }//---line-wise---

        if( !m_file.eof() && !m_file.good() )
            throw runtime_error(string("Corrupt model file after following line: ")+buf);

    }


 private:

    #ifndef SYMBOLIC
    void _buildmap(const Vec& featvec, map<int,unsigned>& map) {
	for(int i=0; i<featvec.size(); i++)
	    map[featvec.at(i)] = i;
    }
    #endif

    int _activeClsId(string id, const RowSetMetadata& metadata) {
	// check whether class -1 is active in current model
	int classid = _getClsIndex(id, metadata);

	// if class -1 is not active, return
	if (classid == -1) {
	    string A82 = "contains coefficients for classes which are not present in the training data or the individual priors file.  This is allowed and will not cause any difficulties for training, but be aware that these classes will not occur in the trained model.";
	    cout<<"WARNING: "<<A8<<m_filename<<" ,"<<A82<<endl;
	    return -1;
	}

	// check -R
	if(classid == metadata.getRefClassId()) {
	    string A831 = "contains one or more nonzero value for coefficients of class ";
	    string A832 = ", which was specified to be the reference class by -R.   Since all coefficients for the reference class are required to be 0, the suggested starting point coefficients for class CLASS will be ignored.  Training will otherwise proceed normally.";
	    cout<<"WARNING: "<<m_filename<<" , "<<A831<<metadata.getClassId(classid)<<A832<<endl;
	    return -1;
	}

	return classid;
    }

    int _getClsIndex(string str, const RowSetMetadata& metadata) {
	#ifdef SYMBOLIC
        return metadata.getClassIndex(str);
	#else
	map<int,unsigned>::iterator miter = m_clsmap.find(atoi(str.c_str()));
	if(miter!=m_clsmap.end()) return miter->second;
	else return -1;
	#endif
    }

    int _getFeatIndex(string str, const RowSetMetadata& metadata) {
	#ifdef SYMBOLIC
        return metadata.getFeatIndex(str);
	#else
	map<int,unsigned>::iterator miter = m_featmap.find(atoi(str.c_str()));
	if(miter!=m_featmap.end()) return miter->second;
	else return -1;
	#endif
    }

    void _checkbetavalue(int fid, int cid, double coef, const BayesParameter& bayesparam, 
			 const RowSetMetadata& metadata, const FixedParams& fixedparam, 
			 const IndividualPriorsHolder& indprior, ParamMatrix& beta) {

	// check skew
	double skew = bayesparam.getSkew();
	if( (skew==-1 && coef>0 ) || (skew==1 && coef<0) ) {
	    string A84 = "contains one or more values which fall outside the range of the prior probability distributions for the corresponding coefficients.  These coefficients will use the usual default starting point values instead. Training will otherwise proceed normally.";
	    cout<<"WARNING: "<<m_filename<<" , "<<A84<<endl;
	    return;
	}

	// check locked 
	if (fixedparam(fid,cid)) {
	    // check prior var 0 first
	    if(indprior.isPriorvar0CoefLevel(fid,cid)) {
		string A86 = "contains one or more values for coefficients whose priors have a variance of 0. The suggested starting point values for these coefficients will be ignored. Training will otherwise proceed normally.";
		cout<<"WARNING: "<<m_filename<<" , "<<A86<<endl;
		return;
	    }
	    else {
		string A85 = "contains one or more nonzero values for coefficients that were locked to 0 by -z. These suggested starting point values will be ignored, but training will otherwise proceed normally.";
		cout<<"WARNING: "<<m_filename<<" , "<<A85<<endl;
		return;
	    }
	}
	
	if(coef==0.0) {
	    string A87 = "contains explicit coefficients of 0.  These will be used as the starting point for optimization for the corresponding coefficients.  We remind the user that implicit coefficients of 0 in the OSPF, i.e. coefficients whose identifier is not mentioned, have no impact on the starting point of optimization.";   // for log file
	    Log(5)<<"WARNING: "<<m_filename<<" , "<<A87<<endl;
	}

	beta(fid,cid) = coef;
    }
};


class WriteModel {
    string m_filename;

    bool m_legacy;
    int m_format;  
    //const vector<MyType>& m_classes;
    //coanst vector<MyType>& m_features;
    ofstream m_file;
    unsigned m_ntopics; //state
    bool m_active;

    // Results file access
    ofstream* m_resultsFile;

    const RowSetMetadata& m_metadata;
public:
    // ctor
    WriteModel(string filename_, int format_, 
	       const RowSetMetadata& metadata_,
	       std::ofstream* resultsFile_, bool legacy_)
        : m_filename(filename_), m_format(format_), m_metadata(metadata_),
	m_resultsFile(resultsFile_), m_legacy(legacy_)  
	{
	    m_active = 0<m_filename.size();
	    if(!m_active) return;         //--->>--
	    m_file.open(m_filename.c_str());
	    
	    // BBR train format ; ver 3.06.  
	    if(m_format == 1){
		m_file<<BBRTITLE<<bbrversionId<<endl;
		m_file<<"tfMethod 0"<<endl;
		m_file<<"idfMethod 0"<<endl;
		m_file<<"cosineNormalize 0"<<endl;
		writeFeats(m_metadata.getFeatures(), m_metadata.getAddFeatures());
		return;
	    }

	    if(m_format == 3) {
		m_file<<SPARSETITLE<<sparseVersionId<<endl;
		return;
	    }

	    // BMRtrain format; ver 3.19
	    if (m_legacy) 
		m_file<<LEGACYTITLE<<legacyVersionId<<endl;
	    else 
		m_file<<TITLE2<<TITLE3<<versionId<<endl;  // ver 3.19
	    
	    if(m_legacy) { // ver 3.19 
		writeClasses(m_metadata.getClasses(), m_metadata.getAddClasses());
		m_file<<"tfMethod 0"<<endl;   // ver 3.19 
		m_file<<"idfMethod 0"<<endl;  // ver 3.19
		m_file<<"cosineNormalize 0"<<endl;  // ver 3.19
		writeFeats(m_metadata.getFeatures(), m_metadata.getAddFeatures());
	    }
	   
	    m_file.flush();
	    if( !m_file.good() )
		throw runtime_error(string("Error creating model file '")+m_filename+"'");
	}
    
    void WriteParams(
        const ModelType& modelType,
        const DesignParameter& design,
        const RowSetMetadata& metadata,
        const ParamMatrix& beta )         //double threshold
   {

       if(!m_active) return;         //--->>--

       if(m_format==3) {
	   writeCoefs(beta, metadata, modelType.getModelName());
	   return;
       }
	    
       if (m_legacy || m_format==1) {
	   m_file<<KW_nDesignParams<<" "<<beta.numFeatures()<<endl;
	   
	   m_file<<KW_modelType 
		 <<" "<<modelType.getLinkType()
		 <<" "<<modelType.getOptimizerType()
		 <<" "<<modelType.getTuneThreshold()
		 <<" "<<modelType.getReferenceClassId()
		 <<" "<<modelType.getStringParam()
		 <<endl;
	   
	   m_file<<KW_design<<" "<<design.DesignType()<<endl;
       }

       m_file<<KW_endofheader<<endl;       
       m_file<<"modelname "<<modelType.getModelName()<<endl; // ver 3.07
       
       set<unsigned> sparseTopicFeats;	    
       ostringstream bbrfeats, bbrvals;
       ostringstream betastr;
       for( unsigned k=0; k<beta.numClasses(); k++ ) { // ver 2.5
	   if(m_format==1)
	       betastr<<KW_beta<<" "<<metadata.getClassId(k)
		      <<setiosflags(ios_base::scientific)<<setprecision(12);
	   else
	       betastr<<KW_betaClassSparse<<" "<<metadata.getClassId(k )
		      <<setiosflags(ios_base::scientific)<<setprecision(12);
	   
	   for( unsigned j=0; j<beta.numFeatures(); j++ )
	       if( 0!=beta(j,k) ) {
		   sparseTopicFeats.insert(j);  // ver 3.04
		   if(m_format==1)
		       betastr<<" "<<beta(j,k);
		   else
		       betastr<<" "<<metadata.getFeatureId(j)<<":"<<beta(j,k);
	       }
	   
	   if(m_format==2)
	       betastr<<" "<<metadata.getNonactiveFeatsVarStr4bmr(metadata.getClassId(k)).str()<<endl; 
	   else {
	       metadata.getNonactiveFeatsVarStr4bbr(bbrfeats,bbrvals);
	       betastr<<" "<<bbrvals.str()<<endl;
	   }

	   if(m_format==1) break;
       }
       
       if(m_legacy || m_format==1){  // ver 3.19
	   // ver 3.04
	   m_file<<KW_topicFeats;
	   // first outout the active feats
	   for(set<unsigned>::iterator siter=sparseTopicFeats.begin(); siter!=sparseTopicFeats.end(); siter++)
	       m_file<<" "<< metadata.getFeatureId(*siter);
	   // then output the others
	   if(m_format==1)
	       m_file<<" "<<bbrfeats.str()<<endl;
	   else
	       m_file<<" "<<metadata.getNonactiveFeatsStr()<<endl;
       }
       
       m_file<<betastr.str();
       
       if(m_format==1) {
	   m_file<<"threshold "<<modelType.getProbThreshold()<<endl;
	   m_file<<"endoftopic"<<endl;
       }
       
   }
    
    std::ofstream& resultsFile() { return *m_resultsFile; }
    
    void closeResultsFile() { 
	if (m_resultsFile == 0)  return;
	
	m_resultsFile->close();
	delete m_resultsFile;
	m_resultsFile = 0;
    }
    
 private:
    void writeClasses( const Vec& classes1, const Vec& classes2 )
	{
	    if(!m_active) return;         //--->>--
	    m_file<<KW_classes;
	    for( unsigned i=0; i<classes1.size(); i++ )
		m_file<<" "<<classes1.at(i);
	    for( unsigned i=0; i<classes2.size(); i++ )
		m_file<<" "<<classes2.at(i);
	    m_file<<endl;
	    m_file.flush();
	}

    void writeFeats( const Vec& feats1, const Vec& feats2 )
	{
	    if(!m_active) return;         //--->>--
	    m_file<<KW_featRestrict;
	    for( unsigned i=0; i<feats1.size(); i++ )
		m_file<<" "<<feats1.at(i);
	    for( unsigned i=0; i<feats2.size(); i++ )
		m_file<<" "<<feats2.at(i);
	    m_file<<endl;
	    m_file.flush();
	}
    void writeCoefs(const ParamMatrix& beta, const RowSetMetadata& metadata, string modelname)
	{
	    for(unsigned k=0; k<beta.numClasses(); k++)
		for(unsigned j=0; j<beta.numFeatures(); j++) 
		    if(beta(j,k)!=0)
			m_file<<modelname<<" "<<metadata.getClassId(k)<<" "
			      <<metadata.getFeatureId(j)<<" "<<beta(j,k)<<endl; 
	    m_file<<endl;

	    map<MyType, map<MyType,double> >::const_iterator iter1;
	    map<MyType,double>::const_iterator iter2;
 
	    const map<MyType, double>& temp = metadata.getNonActiveFeats();
	    for(iter2=temp.begin();iter2!=temp.end();iter2++) {
		for(int i=0; i<metadata.getClasses().size();i++)
		    m_file<<modelname<<" "<<metadata.getClassId(i)<<" "<<iter2->first<<" "<<iter2->second<<endl;
		for(int i=0; i<metadata.getAddClasses().size();i++)
		    m_file<<modelname<<" "<<metadata.getAddClassId(i)<<" "<<iter2->first<<" "<<iter2->second<<endl;
	    }
	    m_file<<endl;    

	    const map<MyType, map<MyType,double> >& temp1 = metadata.getNonactiveCoeflevelFeats();
	    for(iter1=temp1.begin();iter1!=temp1.end();iter1++) 
		for(iter2=iter1->second.begin();iter2!=iter1->second.end();iter2++)
		    m_file<<modelname<<" "<<iter1->first<<" "<<iter2->first<<" "<<iter2->second<<endl;
	    m_file<<endl;
	    
	}

};

#endif //STORED_MODEL_

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
