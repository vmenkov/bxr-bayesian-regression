
#include "MultiModel.h"

/*
  After the model is built. read in the test vector file.

  read in one test vector
      get the docID, trueclasslabel, (feature:weight) pairs
      if --label specified, 
           open labelfile
           get the (docid, mid, cid) tuples for this docid
	   for each (mid,cid)
	        calculate the score and print
      else  // only valid for one model in the system
           calculate the score and print 
      
  Each vector has the form as: 
      vector:=docidline? dataline
      docidline:=ID: docid
      docid:=NAME
      dataline:=classname (feature:weight)*
      feature:=NAME

  If --label is specified, docid must be available, otherwise, exit with error.    
  If no --label specified, and only one model read in, then no docid will be fine; treat as --classic situation.

  If --printlabel is specified, please refer genResultFile.cpp for the logic. In short,
     use the info from label file if --label is specified;
  otherwise, 
     use the classname in dataline as the true class label. (only valid when there is one model)
*/



void ModelDB::processData(){



  ostream * pout;
  // output to a file if --out is specified; otherwise, stdout.
  if(args.outFile != ""){
    pout = new ofstream(args.outFile.c_str());
    if(!((ofstream *)pout)->is_open()){
      cerr << "error";
      exit (1);
    }
  }
  else
    pout =  & cout;

  // open the test vector file; if fails, exit with error.
  istream* fin;
  if(args.dataFile=="-")
      fin = &cin;
  else
      fin = new ifstream(args.dataFile.c_str());

  if(!*fin) {
      // ERROR MESSAGE 41
      if(args.errorFile!="")
	  (*errout) <<"Can't open specified data file: " << args.dataFile << endl;
      throw MyException("Can't open specified data file.");
  }

			 
  // open the label file if --label is specified
  ifstream labelfin;
  int linenum = 0;
  int labelnumber = 0;
  if (args.labelFile != "") {
      labelfin.open(args.labelFile.c_str());
      if(!labelfin.is_open()) {
	  // ERROR MESSAGE 42
	  if(args.errorFile!="")
	      (*errout) << "Can't open specified label file: " << args.labelFile << endl;
	  throw MyException("Cannot open specified label file.");
      }
  }

  string linebuf;   // store each line from the test vector file;
  int lineindex = 0;   // the number of the line in the test vector file
  int crit = getIntCriterion(args.criterion); // the integer index for the parameter of --criterion
  int sc = getIntScale(args.scale); // the integer index for the parameter of --scale 
  bool idFlag = false;    // flip flop flag, make sure the id line and data line happens alternatingly.
  string docID;  // keep the docID of the test vector

  // process each line of the test vector file:
  // if the line is a comment line, skip it;
  // if the line is a docID line, keep the docID;
  // if the line is a data line, getProbabilities() will deal with the line; and sparseOutput() will write the results.
  while(fin->good()){

    // read in a line, keep track of the line number;  
    getline(*fin, linebuf);
    lineindex++;

    // use istringstream to process the line, using space as delim.
    istringstream sstr;
    sstr.str(linebuf);

    // read in the first word of the line
    string first;
    sstr >> first;

    // if the line is empty, then, read the next line.
    if (sstr.fail()) 	continue; 
    
    // if the line is comment line, skip it;
    // if the line is docID line, keep the docID.

    // Here is the logic used for docID reading:
    /*
      If --label specified,           docID is required;
      Else	  set docID to empty, but still outputs

      change to:
      if --classic specified, no docID required;
      else  docID is mandatory.
    */
    if(first[0] == '#'){

	bool idline = false;

        // if the line is a docID line with #ID:
	if(first[1] == 'I' && first[2] == 'D' && first[3] == ':')
	    idline = true;
	// if the line begins with # ID:
	else { // read the next word and it should be ID:
	    string docidflag;
	    sstr >> docidflag;
	    // if fails, then this is a comment line 
	    if (sstr.fail())
		continue;
	    // if the word is ID: then this is a idline
	    if(docidflag[0] == 'I' && docidflag[1] == 'D' && docidflag[2] == ':' )
		idline = true;
	    else // still a comment line
		continue;
	}
	
	{
	    // if there are multiple docIDs for one dataline, quit with error.
	    if(idFlag){
		// ERROR MESSAGE 43
		string emsg = "Wrong format for data file: "+args.dataFile+
		    ". Exactly one ID line should be in front of a dataline."; 
		if(args.errorFile!="")
		    (*errout) << emsg << endl;
		throw MyException(emsg);
	    }

	    // get the docID
	    string originDocID;
	    sstr >> originDocID; 

	    // if --label specified and no docID read in, quit with error; 
	    // change to if --classic not specified and no docID read in, quit with error.
	    // otherwise, keep reading.
	    if(sstr.fail()){
		if(!args.classic) {
		    // ERROR MESSAGE 44
		    if(args.errorFile!="")
			(*errout)<<"Invalid doc ID at line " << lineindex << " in data file: " << args.dataFile << endl;
		    char buf[300];
		    sprintf(buf,"Invalid document ID at line %d in %s.", lineindex, args.dataFile.c_str());
		    throw MyException(buf);
		}
		else 
		    continue;
	    }

	    // docID should not be case sensitive
	    docID = toUpperCase(originDocID);

	    // no colon is allowed in docID 
	    if(docID.find(":") != string::npos){
		// ERROR MESSAGE 45
		string emsg = "Wrong data format: doc ID contains ':' : "+ linebuf; 
		if(args.errorFile!="")
		    (*errout) << emsg << endl;
		throw MyException(emsg);
	    }

	    // pass all test, set idFlag to true
	    idFlag = true;
	} // end of line parsing starting with #

      continue;
    }

    // if not an empty line, not a line beginning with #, then it should be a data line.
    // before processing data line, make sure the docID is existing.
    if(!idFlag) {
	// when --label is specified, docID is required for each test vector; if missing, quit with error
	if(args.labelFile!=""){
	    // ERROR MESSAGE 46
	    if(args.errorFile!="")
		(*errout)<<"All vectors must have DocIDs. Vector ID missing on line "
			 << lineindex << " of test vector file " << args.dataFile << endl;
	    char buf[300];
	    sprintf(buf,"All vectors must have DOCIDs when label file is used. \
Vector ID missing for vector on line %d of test vector file %s.", lineindex, args.dataFile.c_str());
	    throw MyException(buf);
	}
    }


    TMMM modelMarkerMap; 

    // attach "@contant:1" to the dataline, read the line into istringstream
    linebuf += " @constant:1";
    sstr.clear();
    sstr.str(linebuf);
    // the data line is sent and parsed in getProbabilities() method
    getProbabilities(sstr, modelMarkerMap );

    // if --label if specified, and --printtrue specified
    if(args.labelFile!="" && args.printTrue) {
	readLabelFile(labelfin, linenum, labelnumber, docID); 
    }

    // generate results lines for this test vector; the "first" represents the true class id
    sparseOutput(modelMarkerMap, docID, first, crit, sc, pout);

    // set the flag to false after succefully parsing the data line; set docID to empty.
    idFlag = false;
    docID = "";    
    clearEntryMarks(modelMarkerMap);

  }

  // close the output file
  if(args.outFile != ""){
    ((ofstream *)pout) -> close();
    delete pout;
  }

  // close the label file
  if(args.labelFile != "") 
      labelfin.close();


  // data file handle
  if(args.dataFile!="-")
      delete fin;

} // END OF PROCESSDATA()

// accumulators. Two layers, the first layer contains model IDs; the second layer contains class IDs.
// parameters: 
//   line-- the data line in test vector file; "trueClassID <feature:value>*"
//   accumulator -- The probability table. Should be empty when input, 
//     contains probability value when return;
void ModelDB::getProbabilities(istringstream & sstr, TMMM & v)
{    
    // The first word in the data line is the true class label 
    // as the first word has been read in processData(), thus, no test needed here.
    string trueCate;
    sstr >> trueCate;

    // The rest of the line is featureID:weight pairs
    // calculate the sum of inner products of all classes for each model, the key is model ID.
  
    // processing the featureID:weight pairs, get the linear score 
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
	    // ERROR MESSAGE 47
	    string emsg = "Wrong data format for data line:"+sstr.rdbuf()->str()+". Missing colon";
	    if(args.errorFile!="")
		(*errout) << emsg << endl;
	    throw MyException(emsg);
	}

	string featureID = toUpperCase(featurepair.substr(0,pcol));
	double count = atof(featurepair.substr(pcol+1).c_str());

	// THM : the inverted index of featureID; 
	// if the featureID could not be found, the beta for this feature is 0.
	THM:: iterator hmit;
	hmit = hmModel.find(featureID);
	if(hmit == hmModel.end()){
	    continue;
	}
    
	// if the feature is found, get the beta, and multiply with the weight
	map<pair<int,int>, double> &betaMap = hmit->second;
	map<pair<int,int>, double>::iterator feit;

	// for each entry of (<mid,cid>,beta) 
	for(feit = betaMap.begin(); feit != betaMap.end(); feit++){

	    int mid = (feit->first).first;
	    int cid = (feit->first).second;
	    double weight = feit->second;
    
	    // get the accumulator entry for corresponding model
	    TACM::iterator acmit; // acumulator iterator
	    acmit = acm.find(mid);

	    // finally we got model_id, class(category)_id and weight for this 
	    // feature. now update the accumulator
	    TCM & acCates = acmit->second;
	    TCM::iterator accit;
	    accit = acCates.find(cid);

	    double singleProd = weight * count;
	    accit->second.innprod += singleProd;
	    if(!accit->second.used){
		v.push_back(accit);
		accit->second.used = true;
	    }
	}
    } // end of while loop

  
    // all lines are processed, the table has been updated. now we need to
    // divide each entry by the number in sumIP to calculate the probability
    TMMM::iterator vit;
  
    for(vit = v.begin(); vit != v.end(); vit++){
	TCM::iterator tit = *vit;
	TST & t = tit->second;
	double expscore = exp(t.innprod);
	t.prob = expscore;
	sumIP[t.mid] += expscore;
    }
    
    for(vit = v.begin(); vit != v.end(); vit++){
	TCM::iterator tit = *vit;
	TST & t = tit->second;
	t.prob /= sumIP[t.mid];
	t.logprob = log(t.prob);
	t.odds = t.prob / (1-t.prob);
	t.logodds = log(t.odds);  
    }
}


