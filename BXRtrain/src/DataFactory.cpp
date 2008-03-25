/*
bugs:

SL Mar25,08:  convert 0 in indprior file to @constant if defined SYMBOLIC
 */
#include "DataFactory.h"
#include "PolyZO.h"


void DataFactory::setTestAndTrainFileSpec( 
    string trainFile_,
    string testFile_,
    string indpriorFile_,
    string strRefClassId_ 
    ) 
{

    trainFName = trainFile_;
    testFName = testFile_;
    indpriorFName = indpriorFile_;
    m_refClass = strRefClassId_;

    Log(3) <<"\nTime: " << Log.time();

    MakeSparseMatrix( trainFName);  // ver 4.0

    Log(3) <<"\nTime: " << Log.time();
    Log(3)<<"\nCases in training data: "<<m_Mtrain.size();
    
    Log(3)<<"\nFeatures in training data: "<<m_featSelect.size();
    Log(3)<<"\nClasses in training data: "<<m_classes.size();
    
}


void DataFactory::readFiles() {

    if(!indpriorFName.empty())
	indpriorHolder = createIndPriorHolder();
    else
	indpriorHolder = 0;

    #ifdef SYMBOLIC
    m_metadata = new RowSetMetadata(m_classes, m_featSelect, m_addClasses, m_addFeats,nonactiveFeats,
			      nonactiveCoeflevelFeats, clsmap, featmap, m_refClass );
    #else
    m_metadata = new RowSetMetadata(m_classes, m_featSelect, m_addClasses, m_addFeats,nonactiveFeats,
			      nonactiveCoeflevelFeats, m_refClass );
    #endif
    if (!trainFName.empty()) 
	trainData = createTrainRowSet();
    else 
	trainData = 0;
   
    if (!testFName.empty()) 
	testData = createTestRowSet();
    else 
	testData = 0;
    
}

#ifndef SYMBOLIC
void DataFactory::collectRecodeClasses( const vector<int>& clids ) {
    m_classes.clear();
    set<int> cset;

    // collect all the unique classes appearred in the training examples
    for( vector<int>::const_iterator yitr=clids.begin(); yitr!=clids.end(); yitr++ )
	cset.insert( *yitr );
    // store the unique classes into m_classes
    for(set<int>::const_iterator sitr=cset.begin(); sitr!=cset.end(); sitr++ )
        m_classes.push_back( *sitr );

    m_yTrain.clear();
    // replace the actual class string with the mapped class id 
    for( vector<int>::const_iterator yitr=clids.begin(); yitr!=clids.end(); yitr++ ) //recode
        m_yTrain.push_back( lower_bound( m_classes.begin(), m_classes.end(), *yitr ) - m_classes.begin() );
}
#endif


/*
read the data file during training.

m_Mtrain: used to store the training data, declared as vector<sparsevector>;
          each vector in data file is stored as a SparseVector, which is a vector of SparseItem == pair(int, double);

words: keeps all the features, including repeated ones;

Y: keeps the class information for each train vector.
 */

void DataFactory::MakeSparseMatrix( string fName ) { 
    // clear the m_Mtrain, and reserve some space for training data
    m_Mtrain.clear();
    m_Mtrain.reserve(300000);
    vector<int> Y;//temp for collection of class label index

    istream *p_ifs;
    // could read from either stdin or the file
    if( readFromStdin( fName.c_str() ) ) {
        p_ifs = &cin;
    } else {
        p_ifs = new ifstream(fName.c_str() );
    }

    // if cannot open the file, throw error
    if( !*p_ifs )
        throw runtime_error(string("Cannot open training file ")+fName);


    // keep track of the number of rows read in
    int nrows = 0;
    vector<int> words;  // ver 4.1
    istreambuf_iterator<char> it(*p_ifs);
    const std::istreambuf_iterator<char, std::char_traits<char> > end;
    bool linebegin = true;
    int maxy=0, maxfeat=0;  // keep record of the number of unique feature or class labels

    // add intercept first, make sure the index is always 0
    #ifdef SYMBOLIC
    getIntID("@constant", featmap, m_featSelect, maxfeat);
    #else
    m_featSelect.push_back(0);
    #endif

    // if not the end of the stream
    string linebuf;

    while(p_ifs->good()) {
	vector<SparseItem> xprime; 
	SparseVector x(xprime);     // SL ver3.02
	x.insert(0,1.0);  // intercept

	// read in a line
	getline(*p_ifs, linebuf);

	// use istringstream to process the line, space as delim
	istringstream sstr;
	sstr.str(linebuf);

	// read in the first word of the line; if not empty, this is the class label
	string first;
	sstr >> first;

	// if line is empty
	if (sstr.fail()) continue;

	// if line is a comment line, skip it
	if(first[0]=='#')
	    continue;

	#ifdef SYMBOLIC
	int classs = getIntID(first, clsmap, m_classes, maxy);
	m_yTrain.push_back(classs);
	#else
	int classs = atoi(first.c_str());
	Y.push_back(classs);
	#endif

	while (!sstr.eof()){
	    string featurepair;
	    sstr >> featurepair;

	    if(sstr.fail()){
		continue; // end of the line
	    }

	    // find the location of the separator
	    size_t pcol = featurepair.find_first_of(":");
	    // get the featureID:weight pair
	    if(pcol==string::npos) {	
		string emsg = "Wrong data format: Missing colon";
		throw logic_error(emsg);
	    }

	    #ifdef SYMBOLIC
	    string featureStr = featurepair.substr(0,pcol);
	    int feature = getIntID(featureStr, featmap, m_featSelect, maxfeat);
	    #else
	    int feature = atoi(featurepair.substr(0,pcol).c_str());
	    #endif

	    double value = atof(featurepair.substr(pcol+1).c_str());
	    Log(25)<<" "<<feature<<":"<<value; // for test

	    if(value!=0.0) {
		x.insert(feature,value); 
		words.push_back(feature);
	    }
	} // end of the line

	m_Mtrain.push_back(x); //m_Mtrain.push_back( SparseVector(x)); // SL ver3.02
	nrows++;
	if(nrows % 1000 == 0) {
	    std::cout << ".";
	    flush(std::cout);
	}
    } // end of the file

    if( !readFromStdin( fName.c_str() ) ) {
        delete p_ifs;
    }

    /* 
       "words" keeps every appearances of all features or feature index
       m_wordIndex keeps the position (the start position and end position) of each feature.
       eg. if a feature occurs in "words" from postion 0 to 7, then the m_wordIndex will look like:
       0, 8, ....          
    */
    sort(words.begin(), words.end());
    m_wordIndex.clear();
    m_wordIndex.push_back(0);
    m_wordIndex.push_back(nrows);  // number of intercept
    for (unsigned i = 0; i < words.size()-1; i++) {
	if (words.at(i) != words.at(i+1)) {
	    m_wordIndex.push_back(i+1+nrows);
	}
    }
    m_wordIndex.push_back(words.size()+nrows);

    #ifndef SYMBOLIC
    words.erase(unique(words.begin(),words.end()),words.end());
    sort(words.begin(),words.end());
    m_featSelect.insert(m_featSelect.end(), words.begin(), words.end());
    collectRecodeClasses(Y);
    #endif

    // keep the index of reference class; 
    // if it does not exit in training data, keep it in m_addClasses (original form)
    if(m_refClass!="") {

        #ifdef SYMBOLIC
	int temp = maxy;
	m_refClassId = getIntID(m_refClass, clsmap, m_classes, maxy);
	if(maxy>temp) m_addClasses.push_back(m_refClass);
	#else
	m_refClassId = atoi(m_refClass.c_str());
	vector<int>::iterator viter = find(m_classes.begin(), m_classes.end(), m_refClassId);
	if(viter==m_classes.end()) m_addClasses.push_back(m_refClassId);
	else m_refClassId = viter - m_classes.begin();
	#endif
    }
}

IndividualPriorsHolder* DataFactory::createIndPriorHolder() {
    ifstream ifs(indpriorFName.c_str());
    if( !ifs.good() )
	throw runtime_error(string("Cannot open individual priors file ")
			    + indpriorFName);
    Log(10)<<"\nReading  individual priors file "<<indpriorFName;
    
    int nrows = 0;
    set<int> nonactivefeats;

    while(  ifs.good() )
    {
	string buf;
	getline( ifs, buf );

	//check for keyword
	std::istringstream rowbuf( buf);
	string keyword;
	rowbuf>>keyword;
	if( rowbuf.fail() ) //empty line 
	    continue;      //--->>-- ver 3.09
	
	if(keyword[0]=='#')  // ver 3.09
	    continue;
	

	if( keyword == "class" ) { //coefficient-specific prior
	    MyType iclass, feat;
	    int classid=-1, featid=-1;
	    rowbuf>>iclass>>feat;  

	    Vec::iterator viter = find(m_classes.begin(),m_classes.end(), iclass);
	    if(viter==m_classes.end()) 	m_addClasses.push_back(iclass);
	    else  		classid = viter-m_classes.begin();

          // if defined SYMBOLIC, then 0 represents @constant; SL Mar25,08              
          #ifdef SYMBOLIC                                                               
          if (feat=="0") feat = "@constant";                                            
          #endif 
	    viter = find(m_featSelect.begin(),m_featSelect.end(),feat);
	    if(viter == m_featSelect.end()) m_addFeats.push_back(feat); 
	    else 		featid = viter - m_featSelect.begin();

	    ModeVarSkew prior = _parsePriorsLine(rowbuf);
	    
	    // only store that (class,feature) combination if both of them are active in training examples
	    if(classid!=-1 && featid!=-1) {
		map<int,IndividualPriorsMap>::iterator citer = m_byClass.find(classid);
		if(citer!=m_byClass.end() && (citer->second).find(featid)!=citer->second.end() ) 
		    throw logic_error("repeated ind prior mode declaration\n");
		else
		    m_byClass[classid][featid] = prior;
	    }
	    else if (prior.mode !=0 )
		nonactiveCoeflevelFeats[iclass][feat] = prior.mode ;
	    
	    Log(10)<<"- "<<iclass<<" "<<feat<<" "
		   <<prior.mode<<" "<<prior.var<<" "<<prior.skew<<endl;
	}
	else {//feature-level prior
	    std::istringstream rowbuf2( buf); //read whole line from the begining
	    MyType feat;
	    int featid = -1;
	    rowbuf2>>feat;  

          // if defined SYMBOLIC, then 0 represents @constant; SL Mar25,08              
          #ifdef SYMBOLIC                                                               
          if (feat=="0") feat = "@constant";                                            
          #endif 
	    Vec::iterator viter = find(m_featSelect.begin(),m_featSelect.end(),feat);
	    if(viter==m_featSelect.end()) m_addFeats.push_back(feat);
	    else featid = viter-m_featSelect.begin();

	    ModeVarSkew prior = _parsePriorsLine(rowbuf2);
	    if(featid!=-1) {
		if(m.find(featid)!=m.end())
		    throw logic_error("repeated feature level declaration\n");
		else
		    m[featid] = prior;
	    }
	    else if (prior.mode!=0) 
		nonactiveFeats[feat] = prior.mode;
	    
	    Log(10)<<"\n- "<<feat<<" "<<prior.mode<<" "<<prior.var<<" "<<prior.skew;
	}
	nrows ++;
    } // end of while

    if(nrows==0) {  // ver 3.10
	cerr<<"WARNING: The individual prior file " << indpriorFName << " has no contents."<<endl;
    }

    m_addClasses.erase(unique(m_addClasses.begin(),m_addClasses.end()),m_addClasses.end());
    m_addFeats.erase(unique(m_addFeats.begin(),m_addFeats.end()),m_addFeats.end());
    
    if( ! ifs.eof() ) { //the only legal way to end the loop
	std::ostringstream buf; 
	buf<<"Corrupt individual priors file after line # " << nrows;
	throw runtime_error(buf.str());
    }

    return new IndividualPriorsHolder(m_byClass, m, indpriorFName);  
}


MemRowSet* DataFactory::createTrainRowSet() {
    #ifndef SYMBOLIC
    relabelFeatures();
    #endif

    #ifdef SYMBOLIC
    MemRowSet* trainRowSet = new MemRowSet( m_Mtrain, m_yTrain, featmap, m_featSelect, 
					    clsmap, m_classes, m_wordIndex, getRowSetMetadata());  
    #else
    MemRowSet* trainRowSet = new MemRowSet( m_Mtrain, m_yTrain, m_featSelect, m_classes, 
					    m_wordIndex, getRowSetMetadata());  
    #endif

    return trainRowSet;
}

#ifdef SYMBOLIC
int DataFactory::getIntID(const string& stringID, SIHM& keys, Vec& reversekeys, int& maxid) {
    int iIntID;
    SIHM::iterator keyPair;
    keyPair = keys.find(stringID);
    if(keyPair != keys.end()){
	iIntID = keyPair->second;
    }else{
	iIntID = maxid;
	reversekeys.push_back(stringID);
	keys[stringID] = maxid++;
    }
    return iIntID;
}
#endif

void DataFactory::addIntercept() {
     unsigned dim = m_featSelect.size(); 
     for (vector<SparseVector>::iterator i = m_Mtrain.begin(); i<m_Mtrain.end(); i++) {
	 i->insert(dim,1.0);
     }
}

PlainYRowSet* DataFactory::createTestRowSet() {

    return new PlainYRowSet(
        testFName.c_str(),
        m_featSelect,
	getRowSetMetadata()
        );
}

#ifndef SYMBOLIC
void DataFactory::relabelFeatures() {
	vector<int>::iterator featBegin = m_featSelect.begin();
	vector<int>::iterator featEnd = m_featSelect.end();
	for (vector<SparseVector>::iterator i = m_Mtrain.begin(); i < m_Mtrain.end(); i++) {
		SparseVector& v = *i;
		vector<int>::iterator searchPos = featBegin;
		for (SparseVector::iterator j = v.begin(); j < v.end(); j++) {
			vector<int>::iterator p = lower_bound(searchPos, featEnd, j->first);
			int newIndex = p-featBegin;
			j->first = newIndex;
			searchPos = p+1;
		}
	}
}
#endif

RowSetIterator* MemRowSet::iterator() {
	return new MemRowSetIterator(*this);
}

RowSetIterator* PlainYRowSet::iterator() {
	return new PlainYRowSetIterator(*this);
}

PlainYRowSetIterator* PlainYRowSet::plainYiterator() {
	return new PlainYRowSetIterator(*this);
}

ModeVarSkew DataFactory::_parsePriorsLine( istringstream& s)
{
    ModeVarSkew p(0,1,0,false); //skew is optional, so needs to be initialized here
    std::string varstr;

    s>>p.mode>>varstr;
    if( varstr == "inf" ) 
        p.var = std::numeric_limits<double>::infinity();
    else{
        std::istringstream varbuf(varstr);
        varbuf>>p.var;
    }

    if( s.fail() ) //feat, mode,var are mandatory
        throw runtime_error(string("Corrupt individual priors file line: ") + s.str());

    string absstr="";
    s>>absstr;

    if(absstr=="ABS") { 
	p.abs = true;
	s>>p.skew; //optional
    }
    else if(absstr!="") {
	p.skew = atoi(absstr.c_str());
    }

    if( p.skew!=0 && p.skew!=1 && p.skew!=-1 )
	throw runtime_error(string("Illegal skew value: should be 1, 0, or -1"));

    return p;
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
