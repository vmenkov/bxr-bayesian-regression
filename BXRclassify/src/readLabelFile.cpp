
#include "MultiModel.h"

/******************************************************************************************
readLabelFile.cpp     READ IN LABEL FILE WHEN --label IS SPECIFIED

FORMAT OF LABEL FILE:

           FILE := (COMMENTLINE | DATALINE | REFCLASSLINE)*
           COMMENTLINE := <whitespace>*(#<any characters>)?
           DATALINE := DOCID MODELNAME CLASSNAME
           REFCLASSLINE := DOCID @allothermodels REFCLASSNAME
           DOCID := NAME
           MODELNAME := NAME
           REFCLASSNAME := CLASSNAME
           CLASSNAME := NAME
 *****************************************************************************************/



/****************************************** 
readLabelFile(const string &labelFile)

Functionality:  This is the old version of label file reading. 
                It reads the label file into memory, keep the info in labelTable.
	    
Usage:          Call it in MultiDB constructor, after reading the model files. 
                The error check 2 for label file needs the info from model files parsing.
************************ */


void ModelDB::readLabelFile(const string & labelFile){

    // open the label file, 
    ifstream fin(labelFile.c_str());
    if(!fin.is_open()){
	// ERROR MESSAGE 48
	string emsg = "Can't open label file "+labelFile;
	if(args.errorFile!="") (*errout)<<emsg<<endl;
	throw MyException(emsg);
    }

    int lineindex = 0;
    // parse the label file line by line
    while(!fin.eof()){
	string dataLine;
	getline(fin, dataLine);
	lineindex++;

	istringstream sstr;
	sstr.str(dataLine);
	// parse each string in the line; use space as delim
	while(!sstr.eof()){
	    string docID;

	    // get the first string
	    sstr >> docID;

	    // if encountering errors during reading, quit with error message.
	    if(sstr.fail()) {
		// ERROR MESSAGE 49
		if(args.errorFile!="")
		    (*errout)<<"Line "<<lineindex<<" with invalid format encoutered when reading label file "<<labelFile<<endl;
		char buf[200];
		sprintf(buf, "Line %d with invalid format encountered when reading label file %s", lineindex, labelFile.c_str());
		throw MyException(buf);
	    }

	    // skip the COMMENTLINE
	    if(docID[0] == '#')
		continue;
	    
	    // otherwise, data line; get the modelname and classname
	    string modelID, cateID;
	    sstr >> modelID >> cateID;
	    // if encountering errors during reading, quit with error message.
	    if(sstr.fail()) {
		// ERROR MESSAGE 50
		if(args.errorFile!="")
		    (*errout)<<"Line "<<lineindex<<" with invalid format encountered when reading label file "<<labelFile<<endl;
		char buf[200];
		sprintf(buf, "Line %d with invalid format encounteered when reading label file %s", lineindex, labelFile.c_str());
		throw MyException(buf);
	    }

	    // the model name and classname should not be case sensitive
	    string trueModelID = toUpperCase(modelID);
	    string trueCateID = toUpperCase(cateID);
	    string trueDocID = toUpperCase(docID);

	    if(args.classic || args.integerclasses) {
		int tempcid = atoi(trueCateID.c_str());
		char buf[128];
		sprintf(buf,"%d",tempcid);
		trueCateID = buf;
	    }

	    // if this is REFCLASSLine
	    if(trueModelID == "@allothermodels") {

		int cid = getIntID(trueCateID, cateKeys, reverseCateKeys, maxCateID);

		THLL::iterator ltit = labelTable.find(trueDocID);

		// find the entry to the docID in the labelTable
		if(ltit == labelTable.end()) {
		    labelTable[trueDocID] = vector<TID>();
		    ltit = labelTable.find(trueDocID);
		}

		// get the vector<TID> reference
		vector<TID> & vTID = ltit->second;
		vector<TID>::iterator itTID = vTID.begin();
		
		// as sorting is performed right after each inserting, the vector is already sorted
		// compare with the mid list in TACM, which is a hash_map == sorted also
		// store the mid not in THLL, and add them later
		vector<int> tempMids;
		TACM::iterator itACM = acm.begin();
		while(itACM!=acm.end() && itTID!=vTID.end()) {
		    if(itACM->first == itTID->modelID) {
			itACM++; itTID++;
		    }
		    else if(itACM->first < itTID->modelID) {
			tempMids.push_back(itACM->first);
			labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
			itACM++;
		    }
		    // in this case, the model id in label file does not occur in model file; ignore it
		    else {
			itTID++;
		    }
		}
		// add the TID pairs into THLL
		for(int i=0; i<tempMids.size(); i++) {
		    vTID.push_back(TID(tempMids.at(i),cid));
		}
		// if more mids in acm, add them also
		while(itACM!=acm.end()) {
		    vTID.push_back(TID(itACM->first,cid));
		    labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
		    itACM++;
		}
		sort(vTID.begin(),vTID.end(),less<TID>());
		
	    } // end of processing REFCLASSLINE
	    
	    // otherwise, it is a DATALINE
	    else{
		int mid=getIntID(trueModelID, modelKeys, reverseModelKeys, 
				 maxModelID);
		int cid = getIntID(trueCateID, cateKeys, reverseCateKeys, 
				  maxCateID);
		THLL::iterator ltit = labelTable.find(trueDocID);
		if(ltit == labelTable.end()){
		    //	vector<TID> tvtid;
		    labelTable[trueDocID] = vector<TID>();
		    ltit = labelTable.find(trueDocID);
		}
/* 	  labelTable[trueDocID].push_back(TID(mid, cid)); */ // commented by shenzhi, use iterator instead
		// add the new TID pair, and sort
		vector<TID> & vtid = ltit->second;
		vtid.push_back(TID(mid,cid));
		TACM::iterator itACM = acm.find(mid);
		if(itACM!=acm.end())
		    labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
		sort(vtid.begin(),vtid.end(),less<TID>());
		
	    } // end of processing DATALINE
	    
	} // end of reading the line
    } // end of reading the label file
    
    
/* // commented by shenzhi; move the sorting right after inserting; so that when reading @allotherclasses, the TIDs are already sorted.
    // got all the IDs, sort them
    for(THLL::iterator ltit = labelTable.begin();
    ltit != labelTable.end();
    ltit++){
    vector<TID> & vtid = ltit->second;
    sort(vtid.begin(), vtid.end(), less<TID>());
    }
*/
}

void ModelDB::labelFileErrorCheck2(TACM::iterator itACM, int cid, string & upperModelID, string & upperCateID, string & upperDocID){
    TCM & mTCM = itACM->second;
    if(mTCM.find(cid)==mTCM.end()) {
	// ERROR MESSAGE 51
	(*errout) << "WARNING:Label file specifies one or more classes as being associated with model names where the loaded version of that model does not list the class. One example is " + upperDocID + " " + upperModelID + " " + upperCateID + " from label file " + args.labelFile + ". The " + upperModelID + " was loaded from file " + args.modelFile + ". If you specified --printlabel we will go ahead and include these classes in the output, but you should verify you really intended this." << endl;
    }
}


void ModelDB::cntReadLabelFile(ifstream &fin, int &lineindex, string &upperDocID) {

    // flag to indicate whether the upperDocID has been found or not in the label file
    bool found = false;

    // read the label file
    while(!fin.eof()) {
	// read the next line, keep track the line number 
	string dataline;
	getline(fin, dataline);
	lineindex++;

	// read the line into stringstream, space as delim
	istringstream sstr;
	sstr.str(dataline);
	// parse the line word by word, space is used as delim
	// the line should have "docid modelname classname" format
	while(!sstr.eof()) {
	    // read the first word: docID
	    string docID;
	    sstr >> docID;

	    // if empty line, skip it
	    if(sstr.fail()) {
		continue;
	    }

	    // skip the COMMENTLINE
	    if(docID[0] == '#')
		continue;

	    // otherwise, this is a dataline; check whether the docID matches the specified one 
	    string trueDocID = toUpperCase(docID);
	    // if the specified docID has not been found, continue reading;
	    // if already found, then return; 
	    if(trueDocID != upperDocID) {
		// if the specified docID has been found, this is the new docID which should not be read;
		// put back the dataline read 
		if (found) {
		    fin.unget();
		    long pos = fin.tellg();
		    fin.seekg(pos-dataline.length()-1);
		    return;
		}
		// if has not found the docID, read the next line
		else 
		    break; // break the loop for processing this line
	    }
	    // if it is a match, set the flag; and keep the (mid,cid)paris  
	    else{
		if (!found) found = true;

		string modelID, cateID;
		sstr >> modelID >> cateID;
		// if reading is not successful, quit with error
		if(sstr.fail()) {
		    // ERROR MESSAGE 51
		    if(args.errorFile!="")
			(*errout)<<"Line "<<lineindex<<" with invalid format encoutered when reading label file "<<args.labelFile<<endl;
		    char buf[200];
		    sprintf(buf, "Line %d with invalid format encountered when reading label file %s", lineindex, args.labelFile.c_str());
		    throw MyException(buf);
		}

		// the model name and classname should not be case sensitive
		string trueModelID = toUpperCase(modelID);
		string trueCateID = toUpperCase(cateID);

		// get the integer index for the classname
		int cid = getIntID(trueCateID, cateKeys, reverseCateKeys, maxCateID);

		// if this is REFCLASSLine
		if(trueModelID == "@allothermodels") {

		    map<int,int>::iterator mmit = m_trueLabelForOneVector.begin();

		    // check the rest modelnames to be added
		    vector<int> tempMids;
		    TACM::iterator itACM = acm.begin();
		    // compare the model names for this docID with the all the model names kept in acm 
		    while(itACM!=acm.end() && mmit!= m_trueLabelForOneVector.end()) {
			if(itACM->first == mmit->first) {
			    itACM++; mmit++;
			}
			else if(itACM->first < mmit->first) {
			    tempMids.push_back(itACM->first);
			    labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
			    itACM++;
			}
			// in this case, the model id in label file does not occur in model file; ignore it
			else {
			    mmit++;
			}
		    }
		    // add the (mid,cid) pairs into map
		    for(int i=0; i<tempMids.size(); i++) {
			m_trueLabelForOneVector[tempMids.at(i)] = cid;
		    }
		    // if more mids in acm, add them also
		    while(itACM!=acm.end()) {
			m_trueLabelForOneVector[itACM->first] = cid;
			labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
			itACM++;
		    }
		} // end of processing REFCLASSLINE
		
		// if not a reference line, it is a DATALINE
		else{
		    // get the integer index for the model name
		    int mid=getIntID(trueModelID, modelKeys, reverseModelKeys, maxModelID);
		    // add the (mid, cid) pair into the multimap
		    m_trueLabelForOneVector[mid] = cid;
		    // check whether the mid exists
		    // if not, fine, ignore this; if yes, check the classname status
		    TACM::iterator itACM = acm.find(mid);
		    if(itACM!=acm.end())
			labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
		} // end of dataline processing 
	
	    } // end of processing when the docID read in is a match
	
	} // end of parsing the line

	    
    } // end of reading the label file

	
}

void ModelDB::readLabelFile(ifstream &fin, int &lineindex, int &labelnumber, string &upperDocID){

    // if the docID already read in, do nothing
    if(labelTable.find(upperDocID)!=labelTable.end())
	return;

    // parse the label file line by line
    string currentDocID = "";
    bool found = false;

    while(!fin.eof()){

	try{

	    if(!lowmemory) {
		readLine(fin, lineindex, labelnumber, upperDocID, currentDocID, found );
	    }
	    else {
		bool stopreading = readLine(fin, lineindex, labelnumber, upperDocID, currentDocID, found);
		if(stopreading)
		    return;
	    }
	}
	catch(bad_alloc & bae) {
	    if(!lowmemory) {
		delete[] emergencypool;
		cout<<"Out of memory when reading line "<<lineindex<<endl;
		lowmemory = true;		
	    }
	    else {
		throw MyException("Out of memory");
	    }   
	} 
	
    } // end of reading the label file
	    
}


bool ModelDB::readLine(ifstream &fin, int &lineindex, int &labelnumber, string &upperDocID, string &currentDocID, bool& found){

    bool stopreading = false;

    string dataLine;
    getline(fin, dataLine);
    lineindex++;
	    
    istringstream sstr;
    sstr.str(dataLine);

    // parse each string in the line; use space as delim
    while(!sstr.eof()){
	string docID;
	// get the first string
	sstr >> docID;
	
	// if encountering errors during reading, empty line
	if(sstr.fail()) {
	    return false;
	}

	string trueDocID = toUpperCase(docID);

	if(lowmemory) {
	    if(trueDocID!=upperDocID) {
		if(found) {
		    fin.unget();
		    long pos = fin.tellg();
		    fin.seekg(pos-dataLine.length()-1);
		    return true;
		}
		else
		    return false;
	    }
	    else{
		if(!found) found = true;
	    }
	}
	else if(trueDocID==upperDocID && !found) {
	    found = true;
	}

	// skip the COMMENTLINE
	if(docID[0] == '#')
	    continue;
	
	// otherwise, data line; get the modelname and classname
	string modelID, cateID;
	sstr >> modelID >> cateID;
	// if encountering errors during reading, quit with error message.
	if(sstr.fail()) {
	    // ERROR MESSAGE 52
	    if(args.errorFile!="")
		(*errout)<<"Line "<<lineindex<<" with invalid format encountered when reading label file "<<args.labelFile<<endl;
	    char buf[200];
	    sprintf(buf, "Line %d with invalid format encounteered when reading label file %s", lineindex, args.labelFile.c_str());
	    throw MyException(buf);
	}
	
	// the model name and classname should not be case sensitive
	string trueModelID = toUpperCase(modelID);
	string trueCateID = toUpperCase(cateID);

	if(args.classic || args.integerclasses) {
	    int tempcid = atoi(trueCateID.c_str());
	    char buf[128];
	    sprintf(buf,"%d",tempcid);
	    trueCateID = buf;
	}

	currentDocID = trueDocID;
	
	// if this is REFCLASSLine
	if(trueModelID == "@allothermodels") {
	    
	    int cid = getIntID(trueCateID, cateKeys, reverseCateKeys, maxCateID);
	    
	    THLL::iterator ltit = labelTable.find(trueDocID);
	    
	    // find the entry to the docID in the labelTable
	    if(ltit == labelTable.end()) {
		labelTable[trueDocID] = vector<TID>();
		ltit = labelTable.find(trueDocID);
	    }
	    
	    // get the vector<TID> reference
	    vector<TID> & vTID = ltit->second;
	    vector<TID>::iterator itTID = vTID.begin();
	    
	    // as sorting is performed right after each inserting, the vector is already sorted
	    // compare with the mid list in TACM, which is a hash_map == sorted also
	    // store the mid not in THLL, and add them later
	    vector<int> tempMids;
	    TACM::iterator itACM = acm.begin();
	    while(itACM!=acm.end() && itTID!=vTID.end()) {
		if(itACM->first == itTID->modelID) {
		    itACM++; itTID++;
		}
		else if(itACM->first < itTID->modelID) {
		    tempMids.push_back(itACM->first);
		    labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
		    itACM++;
		}
		// in this case, the model id in label file does not occur in model file; ignore it
		else {
		    itTID++;
		}
	    }
	    // add the TID pairs into THLL
	    for(int i=0; i<tempMids.size(); i++) {
		vTID.push_back(TID(tempMids.at(i),cid));
	    }
	    // if more mids in acm, add them also
	    while(itACM!=acm.end()) {
		vTID.push_back(TID(itACM->first,cid));
		labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
		itACM++;
	    }
	    sort(vTID.begin(),vTID.end(),less<TID>());
	    
	} // end of processing REFCLASSLINE
	
	// otherwise, it is a DATALINE
	else{
	    int mid=getIntID(trueModelID, modelKeys, reverseModelKeys, 
			     maxModelID);
	    int cid = getIntID(trueCateID, cateKeys, reverseCateKeys, maxCateID);
	    THLL::iterator ltit = labelTable.find(trueDocID);
	    if(ltit == labelTable.end()){
		//	vector<TID> tvtid;
		labelTable[trueDocID] = vector<TID>();
		ltit = labelTable.find(trueDocID);
	    }
	    
	    // add the new TID pair, and sort
	    vector<TID> & vtid = ltit->second;
	    vtid.push_back(TID(mid,cid));
	    TACM::iterator itACM = acm.find(mid);
	    if(itACM!=acm.end())
		labelFileErrorCheck2(itACM,cid,trueModelID,trueCateID,trueDocID);
	    sort(vtid.begin(),vtid.end(),less<TID>());
	    
	} // end of processing DATALINE
	
    } // end of reading the line
    
    return false;    
}
