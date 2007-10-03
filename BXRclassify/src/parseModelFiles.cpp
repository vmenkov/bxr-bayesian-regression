
#include "MultiModel.h"
#include <algorithm>

/* =====================================================================================================================
Shenzhi Li, Jan07 

parseModelFileDir(string dirName)

   Parameters: dirName == could be a directory or a file path

   Functionality: if the dirName is a file path, read in the model file;
                  if the dirName is a directory, read each model file under that directory. if not a model file, skip it; 

   How to use: the method will be called in ModelDB constructor, when the -m directory is specified on command line. 
               ./BXRclassify -m directory [OPTIONs] MODELPATH DATAPATH > RESULTFILEPATH

   Background information: 
         this method is designed to read in the model format #2, which has three (or four, or six, depending how you 
         count) formats, the choice of which depends on the pathname specified by MODELPATH on the command line: 

     2a. Format 2a applies if MODELPATH is the pathname of a regular file (not a directory), the first line of the 
         file has the form: 

            Bayesian Binary Regression TOKEN TOKEN+ 

	 We parse the file as a classic BBRtrain format model file containing a single binary logistic regression model.  
	 This will be the only model applied on this run.  

     2b. Format 2b applies if MODELPATH is the pathname of a regular file, and the first line of that file has one of 
         these forms:  

            Bayesian Polytomous Regression TOKEN TOKEN+
            Bayesian Multinomial Regression TOKEN TOKEN+
            Multinomial logistic regression model format sparse VERSION TOKEN+
            Multinomial logistic regression model format sparse-symbolic VERSION TOKEN+ 
 
	 Then we parse the file as a BMRtrain format model file containing a polytomous logistic regression model for 2 
         or more classes. This will be the only model applied on this run, and the MODELNAME for the model will be the 
         filename. 

     2c. Format 2c applies if MODELPATH is the pathname of a directory.  In that case each file in the directory should 
         be a BBRtrain format model file or a BMRtrain format model file. Any combination of BBRtrain and BMRtrain 
         models can be present.  
	 We read, store, and apply all of these models to the test vectors.  
	 The MODELNAME for each model is the corresponding file name, and the classes for the model are the classes 
         specified in the corresponding file.  

     - From "bmrclassify-specs-20070116-change-tracked-from-prev.doc".	 
===================================================================================================================== */ 

void ModelDB::parseModelFileDir(string dirName){
#ifdef USE_GCC
    DIR * dirp = opendir(dirName.c_str());
    struct dirent * dp;
    ifstream fin(dirName.c_str());

    // if the path is a directory, parse the each model in this directory
    // ASSUMPTION: here it is only a two tier structure. no recursive reading.
    if (dirp) {

	errno = 0;
	// read each file
	while (( dp = readdir(dirp)) != NULL) {
	    // if the file is . or .., skip it	
	    if(!strcmp(dp->d_name, ".") || !strcmp(dp->d_name, "..")){
		continue; // do nothing for "useless" directories
	    }
	    // if the file ends with ~, skip it also
	    if (dp->d_name[strlen(dp->d_name)-1]=='~') {
		continue;
	    }
	    // otherwise, use the fileName as the modelID
	    string modelID(dp->d_name);
	    string fullpath=dirName+"/"+dp->d_name;
	    // parse the single model file then
	    parseSingleModelFile(fullpath, modelID);
	}
    }
    // if the path is a file, use the file name as the modelID
    else if (fin) {
	bDir = false;
	int pos = dirName.find_last_of("/");
	if (pos!=-1) {
	    string modelID(dirName.substr(pos+1));
	    parseSingleModelFile(dirName,modelID);
	}
	else {
	    parseSingleModelFile(dirName,dirName);
	}
    }
    // otherwise, error msg
    else{
	cerr << "\n parseDirectoryModel: Can't open directory " << dirName << endl;
	exit (1);
    } 
#endif //USE_GCC

#ifdef _MSC_VER
    // I couldn't find good apis to read a directory. If this is the only
    // way, then I don't understand why microsoft doesn't provide a better
    // wrapper for these functions. -- Bing Bai.  

    long hfile = -1;
    int  status = 0;
    char fullpath[MAXPATHLEN] = {0};
    int curDrive;
    char curPath[MAXPATHLEN];
    
    // store the information of the file/directory
    struct stat st;    

    // preprocess the name, get rid the last char if it is '\'
    if(dirName.at(dirName.length()-1)=='\\') {
	dirName = dirName.substr(0, dirName.length()-1);
    }

    // if it is not empty
    if( stat(dirName.c_str(), &st)==0 ){
	// if is a directory
	if (st.st_mode & _S_IFDIR) {
	    // keep the current driver and path info
	    curDrive = _getdrive();
	    _getcwd(curPath, MAXPATHLEN);
	    // go into the directory
	    status = _chdir(dirName.c_str());
	    // check each file in the directory
	    SearchDirectory((char*)dirName.c_str());
	    // go back to the original place
	    _chdrive(curDrive);
	    _chdir(curPath);
	}
	// if it is a file
	else {
	    bDir = false;
	    int pos = dirName.find_last_of("\\");
	    if (pos!=-1) {
		string modelID(dirName.substr(pos+1));
		parseSingleModelFile(dirName,modelID);
	    }
	    else {
		parseSingleModelFile(dirName,dirName);
	    }
	}
    }

#endif //_MSC_VER
}

#ifdef _MSC_VER 
// added by shenzhi for MSVC	
void ModelDB::SearchDirectory(char *pathname)
{
    long handle;
    struct _finddata_t filestruct;  
    //info for the file(or directory)
    char path_search[500]; 
    
    // start the searching, find the first file or subdirectory under current path
    // "*" represents "search for everything", filestruct keeps the searching results
    handle = _findfirst("*", &filestruct);
    
    // if handle == -1, the directory is empty, stop search and return 
    // ERROR MESSAGE 12
    if((handle == -1)) 
	throw MyException("Wrong model file or path value");

    do{	
	// check whether the first object is a directory (filestruct.name is the pathname) 
	if(::GetFileAttributes(filestruct.name) & FILE_ATTRIBUTE_DIRECTORY) {

	    // if it is a directory, skip it; 

	    // the code commented below is to read recursively under the directory.
	    /*
	    // if it is a directory, enter it and recursively call search_directory
	    // note: skip "." or ".." files 
	    if(filestruct.name[0] != '.') {
	    _chdir(filestruct.name);
	    SearchDirectory(pathname);
	    // after searching, go back up a level
	    _chdir("..");
	    }
	    */
	}
	// if it is a file, and not ending with ~
	else if(filestruct.name[strlen(filestruct.name)-1]!='~') {
	    // get the full path
	    _getcwd(path_search, 500);  
	    // then get the pathname for the file (including the filename)
	    strcat(path_search,"\\");
	    strcat(path_search,filestruct.name);
	    // parse the file
	    string modelID(filestruct.name);
	    parseSingleModelFile(path_search, modelID);		
	} 
    } while (_findnext(handle, &filestruct)==0);

    _findclose(handle); 
    
}
#endif

/* ==============================================================
parseSingleModelFile(string fileName, string & modelID)

     Parameters: fileName is the whole file path for the file; while modelID is using the file name currently.

     Functionality: this method aims to read single BBRtrain model file and BMRtrain model file.
                    The most important goal is to get the feature:beta_value pairs, and store them in inverted index hmModel.
                    For BBRtrain model, the threshold value should be kept for class choosing.
                    The method will also check the correctness of the model file while parsing. 

     Usage: this method will be called by parseModelFileDir(string dirName), when -m directory is specified on command line.

     Background information:

     a. We produce a model from a BBRtrain model file as follows: 

         1). The MODELNAME for the corresponding model will be the filename.
         2). The classes for the model will be named "-1" and "+1". 
         3). The beta values in the model file are assumed to all be for class "+1". All beta values for class "-1" equal 0.

	 The format of the model file is to have these lines (some of them we will treat as optional):

	               Bayesian Binary Regression TOKEN TOKEN+
		       (tfMethod INTEGER)?
		       (idfMethod INTEGER)?
		       (cosineNormalize INTEGER)?
		       (featRestrict INTEGER+ )?
		       (endofheader)?
		       (topic <class>)?
		       (modelType INTEGER INTEGER INTEGER INTEGER)?
		       (design 0)? 
		       topicFeats INTEGER*
		       beta FLOAT+
		       (threshold FLOAT)?
		       endoftopic 


     b. The format of the BMR model file is to have these lines (some of them we will treat as optional):

		       FILETYPELINE
		       (classes CLASSNAME+)?
		       (tfMethod INTEGER)?
		       (idfMethod INTEGER)?
		       (cosineNormalize INTEGER)?
		       (featRestrict FEATUREID+)?
		       (nDesignParams INTEGER)?
		       (modelType INTEGER INTEGER INTEGER INTEGER INTEGER)?
		       (design INTEGER)?
		       (topicFeats FEATUREID+)?
		       endofheader 
	       followed by one or more lines of the form: 
		       betaClassSparse CLASSNAME (FEATUREID:COEFFICIENT)* 		       

================================================================*/

void ModelDB::parseSingleModelFile(string fileName, string & modelID){
  ifstream fin(fileName.c_str());
  string linebuf;
  int lineindex = 0;

  // the model file could be BBR model file or BMR model file
  string modeltype = "unknown";

  // <featureID, mid, cid, weight> will be added to datastructure hmModel
  string cateID, featureID;
  double weight;

  // store the feature ids for BBR models
  vector<string> vTopicFeats;

  // store the classes shown in "classes" line in BMR model file; 
  // when processing the betaClassSparse lines, remove the classname in the vector if it appears in betaClassSparse lines;
  // in the end, for each class remained in the vector, add a reference class.
  vector<string> vClasses;

  // keeps the thresh value in BBRtrain model. will be used in label choosing.
  string threshval;

  // flags for line existence
  bool hasBetaLine = false;
  bool hasTopicFeatsLine = false;
  bool hasThresholdLine = false;
  bool hasClsLine = false;
  bool hasCosineLine = false;
  bool hasTFLine = false;
  bool hasIDFLine = false;
  // the number of lines with "betaClassSparse" as beginning; use this when no "classes" line is present
  int nbetaClassSparseLines = 0;



  while(!fin.eof()){

    // read in a line; keep track of the line number, including empty lines & comment lines;
      char firstlinebuf[256];
      if(lineindex==0) {
	  fin.getline(firstlinebuf,256);
	  lineindex++;
	  linebuf = firstlinebuf;
      }
      else{
	  getline(fin, linebuf);
	  lineindex++;
      }

    // find out the model file type, BBR or BMR
    // !!!!!! ASSUMPTION: each model file has a header line to specify the model type; no empty line or comment lines ahead!!!!
    // !!!!!! Note: the header is defined in Multimodel.h, and is CASE SENSITIVE !!!!!!!!!
    if ( lineindex==1 && modeltype == "unknown") {
	// tell whether it is a BBR model file, which begins with "Bayesian Binary Regression".	
	if(linebuf.find(BBRheader,0)!=string::npos) {
	    modeltype = "BBRtrain";
	    cateID = "1";
	    continue;
	}
	// tell whether it is a BMR model file
	// could accept the following four:
	//    Bayesian Polytomous Regression TOKEN TOKEN+
	//    Bayesian Multinomial Regression TOKEN TOKEN+
	//    Multinomial logistic regression model format sparse VERSION TOKEN+
	//    Multinomial logistic regression model format sparse-symbolic VERSION TOKEN+
	else if (   linebuf.find(BPRheader,0)!=string::npos 
	    || linebuf.find(BMRheader,0)!=string::npos
	    || linebuf.find(MLRheader,0)!=string::npos) {
	    modeltype = "BMRtrain";
	    continue;
	}
	// tell whether it is a sparse model file
	else if ( linebuf.find(SPSheader,0)!=string::npos) {
	    modeltype = "Sparse";
	    fin.close();
	    parseSparseModel(fileName);
	}
	// otherwise, this is not a model file
	else {
	    fin.close();
	    // if the command line option specifies a model file instead of a directory, then throw error msg.
	    if(!bDir) {
		// ERROR MESSAGE 13
		string emsg = "Error: " + args.modelFile + " is not a valid model file -- No model file heading at line 1.";
		throw MyException(emsg);
	    }
	    // if a directory of model files is speciified, return, and get the next file. 
	    // the check will be done after reading the directory.
	    return;
	}
    }

    // processing other lines using istringstream, i.e., read word by word 
    istringstream sstr;
    sstr.str(linebuf);
    string token;
    // get the first word of the line
    sstr >> token;

    // if this is a comment line, skip it 
    if (token[0]=='#')
	continue;

    // if this is an optional line, skip it also
    // further checks could be inserted if necessary
    if (token=="featRestrict" || token=="endofheader" || token=="topic" || token =="modelType" || token =="nDesignParams" || token == "design")
	continue;
 
    // if this is tfMethod, idfMethod, cosNormalize line
    // NOTE: if no such lines, fine;
    //       if there is a line, the value must be 0; otherwise, quit
    //       if there are more than one lines beginning with the same key word like "tfMethod", error;
    if (token == "tfMethod" || token =="idfMethod" || token == "cosineNormalize") {
	// tell whether it is a first appearance
	// ERROR MESSAGE 14
	if ( (token == "tfMethod" && hasTFLine ) || (token=="idfMethod" && hasIDFLine) || (token == "cosineNormalize" && hasCosineLine)) {
	    string emsg = "Wrong data format: the " + token + " line appeared more than once in model file " + fileName;
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	if(token=="tfMethod") hasTFLine = true;
	if(token=="idfMethod") hasIDFLine = true;
	if(token=="cosineNormalize") hasCosineLine=true;
	// get the method value
	string methodval;
	sstr >> methodval;
	// if could not read successfully
	if (sstr.fail()) {
	    // ERROR MESSAGE 15
	    string emsg = "Wrong format in model file: " + fileName;
	    if(args.errorFile!="") (*errout)<<emsg<<" at line "<<lineindex<<endl;
	    throw MyException(emsg);
	}
	// if the value is not 0
	if (atoi(methodval.c_str())!=0) {
	    // ERROR MESSAGE 16
	    string emsg = token +" "+ methodval +" was specified in model file " + fileName 
		+ ". This version of BXRclassify does not support transforming feature values. Please transform data before using BMR.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	// skip the rest of the line, shall we? Or what might be right after it? :-)
	continue;
    }

    // if this is a modelname line, use this as the modelID
    if(token == "modelname") {
	string tempmodelname;
	sstr >> tempmodelname;
	if(sstr.fail()) {
	    string emsg = "Warning: No modelname specified";
	    if(args.errorFile!="") 
		(*errout)<<emsg<<endl;
	}
	modelID = tempmodelname;
    }

    // if this is the topicFeats line. This is required in BBRtrain model, while not important in BMRtrain 
    if(token == "topicFeats") {
	// set the flag, for BBRtrain model to check when it reads a "beta" line
	if(!hasTopicFeatsLine)
	    hasTopicFeatsLine = true;
	else{
	    // ERROR MESSAGE 17
	    string emsg = "Wrong format in model file " + fileName +": topicFeats line appeared more than once.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	// add the featureIDs into the vector
	while(!sstr.eof()) {
	    string topicfeat;
	    sstr >> topicfeat;
	    // if one reading fails, skip it to the next one; 
	    if(sstr.fail())
		continue;
	    vTopicFeats.push_back(topicfeat);
	}
	// also, add the "constant" feature
	vTopicFeats.push_back("@constant");
	continue;
    }

    if(token == "beta") {
	if(!hasBetaLine)
	    hasBetaLine = true;
	else{
	    // ERROR MESSAGE 18
	    string emsg = "Wrong format in model file " + fileName +": beta line appeared more than once.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}

	// if no topicFeats line before beta line in model file, quit with error msg
	if(vTopicFeats.size()==0) {
	    // ERROR MESSAGE 19
	    string emsg = "beta line was not preceeded by a topicFeats line in " + modeltype + " model file " + fileName;
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	// read in the beta value, get the corresponding featureID, add the entry; cateID is set to "1" before;
	int topicindex = 0;
	while(!sstr.eof()) {
	    string betaval;
	    sstr >> betaval;
	    if(sstr.fail()) {
		topicindex++;
		continue;
	    }
	    double weight = atof(betaval.c_str());
	    if(topicindex < vTopicFeats.size())
		addEntry(vTopicFeats.at(topicindex++),modelID,cateID,weight);
	    else {
		// ERROR MESSAGE 20
		string emsg = "Wrong model file format: more beta values required in " + fileName;
		if(args.errorFile!="") (*errout)<<emsg<<endl;
		throw MyException(emsg);
	    }
	} // finish reading the betas 

	continue;
    }

    // if the line is threshold line
    if (token == "threshold") {
	// if the line appears before, quit with error;
	if(!hasThresholdLine)
	    hasThresholdLine = true;
	else {
	    // ERROR MESSAGE 21
	    string emsg = "Wrong format in model file " + fileName +": threshold line appeared more than once.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	// get the threshold value
	sstr >> threshval;
	if (sstr.fail()) {
	    // ERROR MESSAGE 22
	    string emsg = "Error in reading the threshold value if file "+fileName;
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	continue;
    }

    // if the line is classes line in BMR model
    if (token == "classes") {
	// if the line appears before, quit with error;
	if(!hasClsLine)
	    hasClsLine = true;
	else {
	    // ERROR MESSAGE 23
	    string emsg = "Wrong format in model file " + fileName +": classes line appeared more than once.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	// read in the classes' names into the vector<string> vClasses
	while (!sstr.eof()) {
	    string classname;
	    sstr >> classname;
	    if(sstr.fail()) {
		continue;
	    }
	    vClasses.push_back(classname);
	}
	// if there are less than 2 classes specified, quit with error
	if (vClasses.size()<2) {
	    // ERROR MESSAGE 24
	    string emsg = "There must be a minimum of two classes in classes line of a BMRtrain files.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
	}
	continue;
    }

    // for BMR models, need to read in betaClassSparse line
    if (token == KW_betaSparse) {
	// the second word is the class ID
	sstr >> cateID;
	if(sstr.fail()){
	    // ERROR MESSAGE 25
	    string emsg = "Wrong format in model file: "+fileName;
	    if(args.errorFile!="") (*errout)<<emsg<<" at line "<<lineindex<<endl;
	    throw MyException(emsg);
	}

	// keep tract of the number of lines if there is a classID for this line;
	nbetaClassSparseLines++;

	// this boolean is used to see whether the class is a reference class or not;
	bool hasFeaturePairs = false;
	// read in the feature:weight pairs
	while(!sstr.eof()){
	    char featurePair[1024];
	    // read in one feature:weight pair
	    sstr >> featurePair;
	    if(sstr.fail()){
		continue;
	    }  
	    // if succeeds, this class is not a reference class
	    hasFeaturePairs = true;
	    // split the string into featureID and weight
	    char * pcol = strstr(featurePair, ":");
	    if(pcol == NULL){
		// ERROR MESSAGE 26
		string emsg = "Wrong format in model file: "+fileName +". Missing colon in " + featurePair + "\n";
		if(args.errorFile!="") (*errout)<<emsg<<endl;
		throw MyException(emsg);
	    }
	    *pcol = 0;
	    string featureID(featurePair);
	    double weight = atof(pcol+1);
	    // if the featureID is 0(in BBR/BMR model files) or -1(in -m vector option model files)
	    if(featureID=="0" || featureID=="-1")
		featureID = "@constant";
	    addEntry(featureID, modelID, cateID, weight);
	}


	// check whether the class is in "classes" line if there is a "classes" line
	// if yes, and the class is a reference class, do nothing;
	// if yes, and the class is not a reference class, remove the entry in vClasses; 
	// if no, error
	if (vClasses.size()>0) {
	    vector<string>::iterator vIt = find(vClasses.begin(),vClasses.end(),cateID);
	    if(vIt!=vClasses.end() && hasFeaturePairs)
		vClasses.erase(vIt);
	    else if(vIt==vClasses.end()) {
		// ERROR MESSAGE 27
		string emsg = "Class "+cateID+" is not in classes list, but is in a betaClassSparse line, in file " + fileName;
		if(args.errorFile!="") (*errout)<<emsg<<endl;
		throw MyException(emsg);
	    }
	}
    
	continue;
    }

    // for other vlaues of tokens, no check and no constraints


  }  // end of reading file

  fin.close();

  // post processing now
  if (modeltype == "BBRtrain") {
      // check whether the BBRtrain model has a beta line
      if (!hasBetaLine) {
	  // ERROR MESSAGE 28
	  string emsg = "beta line is missing from " + modeltype + " model file " + fileName;
	  if(args.errorFile!="") (*errout)<<emsg<<endl;
	  throw MyException(emsg);
      }
      // check whether the BBRtrain model has a threshold line
      if (!hasThresholdLine) {
	  // ERROR MESSAGE 29
	  string emsg = "threshold line is missing from " + modeltype + " model file " + fileName;
	  if(args.errorFile!="") (*errout)<<emsg<<endl;
	  throw MyException(emsg);
      }
      // check whether the BBRtrain model has a topicFeats line
      // if there is a beta line, then this error will be reported when reading the beta line already
      if (!hasTopicFeatsLine) {
	  // ERROR MESSAGE 30
	  string emsg = "topicFeats line is missing from " + modeltype + " model file " + fileName;
	  if(args.errorFile!="") (*errout)<<emsg<<endl;
	  throw MyException(emsg);
      }

      // if the BBRtrain model file has threshold line, keep the info in m_BBRModelThresh map
      // check whether m_BBRModelThresh already has an entry for the modelID, if yes, quit with error; otherwise, insert the entry.
      // to be consistent with BBRclassify, which uses 0.5 as default threshold value when model file did not give a value larger than 0.
      int mid = getIntID(modelID, modelKeys, reverseModelKeys, maxModelID);
      if (m_BBRModelThresh.find(mid)==m_BBRModelThresh.end()) {
	  m_BBRModelThresh[mid] = atof(threshval.c_str())>0.0 ? atof(threshval.c_str()) : 0.5;
      }
      else {
	  // ERROR MESSAGE 31
	  string emsg = "The model ID " + modelID + " has already existed.";
	  if(args.errorFile!="") (*errout)<<emsg<<endl;
	  throw MyException(emsg);
      }
  }

  if(modeltype == "BMRtrain") {
      // check when no classes line is available; the classes line check is done when reading the classes line
      if (nbetaClassSparseLines<2 && !hasClsLine) {
	  // ERROR MESSAGE 32
	  string emsg = "There must be a minimum of two betaClassSparse lines in BMRtrain files that omit the 'classes' line.";
	  if(args.errorFile!="") (*errout)<<emsg<<endl;
	  throw MyException(emsg);
      }
  }

  // add reference classes if there is any
  if (modeltype == "BBRtrain") {
      addEntry("@constant",modelID,"-1",0.0);
  }
  if(modeltype == "BMRtrain" && vClasses.size()>0) {
      for (int i=0; i<vClasses.size(); i++)
	  addEntry("@constant",modelID,vClasses.at(i),0.0);
  }

}



// ...
// Format 1 is what I'll call the "sparse" format.  It is an ASCII, line
// oriented format.  Here's pseudo-BNF for it:
//     FILE := (COMMENTLINE | DATALINE)*
//     COMMENTLINE := <whitespace>* # <any characters> <end-of-line>
//     DATALINE := MODELNAME CLASSNAME FEATURENAME COEFFICIENT <end-of-
// line>
//     MODELNAME:= NAME
//     CLASSNAME:= NAME
//     FEATURENAME:= (NAME | @constant)
//     NAME:= <any string of printable, non-whitespace, ASCII
// characters not starting with #>
//     COEFFICIENT:= <doubleing point number; exponential notation allowed>
// Using the notation of my 20-June-06 msg, a data line has the form
//      h k j beta[h,k,j]
// with h, k, and j allowed to be symbolic names.  (@constant is used
// for the constant term.)
// ...
// --Email from Dave Lewis 06/30/2006. Bing Bai  
void ModelDB::parseSparseModel(const string & fileName){
  
    vector<string> localModelNames;

  ifstream fin(fileName.c_str());
  if(!fin.is_open()){
      // ERROR MESSAGE 33
    string emsg = "Can't open model file "+fileName;
    if(args.errorFile!="") (*errout)<<emsg<<endl;
    throw MyException(emsg);
  }
  string buf;
  string modelID;
  string cateID;
  string featureID;
  double weight;
  string strRefclsName = "";
  int lineindex = 0;

  while(!fin.eof()){
    getline(fin, buf);
    istringstream sstr;
    sstr.str(buf.c_str());

    sstr >> modelID;
	
    // can't read a single string, empty line
    if(sstr.fail()){
      continue;
    }
    if(modelID[0] == '#'){ // comment line
      continue;
    }

    // empty line and comment lines are not counted
    lineindex ++;

    // if the line is reference line
    if(modelID=="@allmodels") {
	if(lineindex==1) {
	    sstr >> strRefclsName;
	    if(sstr.fail()) {
		// ERROR MESSAGE 34
		string emsg = "Can not read in Reference Class Name in " + fileName;
		if(args.errorFile!="") (*errout)<<emsg<<endl;
		throw MyException(emsg);
	    }
	    continue;
	}
	else {
	    // ERROR MESSAGE 35
	    string emsg = "The reference class line is not the first data line of the model file " + fileName;
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg); 
	}
    } 

    sstr >> cateID >> featureID >> weight;
    // error in readling lines
    if(sstr.fail()){
	// ERROR MESSAGE 36
      string emsg = "Wrong format in model file: "+ fileName;
      if(args.errorFile!="") (*errout)<<emsg<<endl;
      throw MyException(emsg);
    }
    // error check whether the class name is the same as the refclass name
    if (strRefclsName!="" && cateID==strRefclsName) {
	// ERROR MESSAGE 37
	    string emsg = "Model file " + fileName + " attempts to define a coefficient value for the combination (" + modelID + ", " + cateID + ", " + featureID +"). " + cateID + " was defined as the reference class, so this is not allowed.";
	    if(args.errorFile!="") (*errout)<<emsg<<endl;
	    throw MyException(emsg);
    }

    // cout << modelID << " " << cateID << " " << weight << endl;
    addEntry(featureID, modelID, cateID, weight);

    string upperMID = toUpperCase(modelID);
    if(find(localModelNames.begin(),localModelNames.end(),upperMID)==localModelNames.end())
	localModelNames.push_back(upperMID);
  } // end of reading file

  if (strRefclsName!="") {
      for(int i=0; i<localModelNames.size(); i++) {
	  addEntry("@constant",localModelNames.at(i),strRefclsName, 0.0);
      }
  }
}


void ModelDB::addEntry(const string & featureID, const string & modelID,
		       const string & cateID, const double weight){

    // if -k/--classic is specified, cateID should be ASCII strings for integers
    // if -I/--intergerclasses specified, cateID should be ASCII strings for integers also
    if ((args.classic || args.integerclasses)){
	// rough pre-check for the cateID string
	if(!intvalue(cateID)) {
	    // ERROR MESSAGE 38
	    string emsg = "The class id " + cateID + " for " + modelID + " is out of range.";
	    throw MyException(emsg);
	}
    }

    // featureID, modelID, cateID should not be case sensitive
    string upperFID = toUpperCase(featureID);
    string upperMID = toUpperCase(modelID);
    string upperCID = toUpperCase(cateID);

    THM::iterator hmit; // hash map
    hmit=hmModel.find(upperFID);
    // if the entry does not exist, add into hash map, and point the iterator to the newly added entry.

    if(hmit==hmModel.end()) {
	hmModel[upperFID] = map<pair<int,int>,double>();
	hmit = hmModel.find(upperFID);
    }
    map<pair<int,int>,double> &feEntries = hmit->second;

    int mid = getIntID(upperMID, modelKeys, reverseModelKeys, maxModelID);
    if(args.classic || args.integerclasses) {
	int tempcid = atoi(cateID.c_str());
	char buf[128];
	sprintf(buf,"%d",tempcid);
	upperCID = buf;
    }
    int cid = getIntID(upperCID, cateKeys, reverseCateKeys, maxCateID);

    pair<int,int> keys(mid,cid);
    if(feEntries.find(keys)!=feEntries.end()) {
	// ERROR MESSAGE 39
	string emsg = "The coefficient for the combination (" + upperMID + "," + upperCID + "," 
	    + upperFID +") is multiply defined in model file.";
	if(args.errorFile!="") (*errout)<<emsg<<endl;
	throw MyException(emsg);
    }
    else {
	feEntries.insert(pair<pair<int,int>,double>(keys,weight));
    }


    // build the accumulator structure. Basically it is a two-level
    // map or hash_map, the index for first level is the model id,
    // that for the second level is the class id
    TACM::iterator acmit = acm.find(mid);
    // if -k/--classic is specified, check the number of models in acm currently
    // otherwise, if the model is not in the table yet, create a new entry for it.
    if(acmit == acm.end()){

	if(args.classic && acm.size()>0 ) 
	    // ERROR MESSAGE 40
	    throw MyException("Only one model could be specified with '-k/--classic'.");

	acm[mid] = TCM();
	acmit = acm.find(mid);
    }
    TCM &acCates = acmit->second;
  
    TSM::iterator sit = sumIP.find(mid);
    if(sit == sumIP.end()){
	sumIP[mid] = 0;
	sit = sumIP.find(mid);
    }
    
    TCM::iterator accit = acCates.find(cid);
    if(accit == acCates.end()){
	TST ttt;
	ttt.mid = mid;
	ttt.cid = cid;
	ttt.innprod = 0;
	ttt.odds = 0;
	ttt.prob = 0;
	ttt.logodds = 0;
	ttt.logprob = 0;
	ttt.used = false;
	acCates[cid] = ttt;
    }
}
  
void ModelDB::addExtraClass(const string & extraClassName){
    THKL::iterator cateks=cateKeys.find(extraClassName);
    string upperExtraClassName = toUpperCase(extraClassName);

    TVKRL::iterator modelks;
    string constantKey("@constant");
    for(modelks = reverseModelKeys.begin(); modelks != reverseModelKeys.end();  modelks++){
	addEntry(constantKey, *modelks, upperExtraClassName, 0.0); 
    }	
}
