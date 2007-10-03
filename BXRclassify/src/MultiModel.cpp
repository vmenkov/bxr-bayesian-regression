// author: Bing Bai
// email: bbai@cs.rutgers.edu
// 2006/06/29

/*
changes made:
   1. SL 14Jan07. change #constant to @constant
   2. SL 15Jan07. change Results File Format 1 (sparse results format) to put doc id at the beginning, 
                  "true:" and "pred:" before true and predicted classes; 
   3. SL 16Jan07. add parseDirectoryModel2() && parseSingleModelFile2()
   4. SL 16Jan07. in processData(), allow the docID be read even no label file is specified
   5. SL 16Jan2007. in processData(), bug - error in empty line reading
   6. Jan 17,2007. Add parseDirectoryModel2(string) to replace parseDirectoryModel(string);
   7. Jan 17,2007. Add parseSingleModelFile2(string,string&) to replace parseSingleModelFile(string,string&)
   8. Jan 17,2007. in parseSingleModelFile2(), add error check when there are multiple lines about the same keywords which should only appear once;
   9. Jan 17,2007. in parseDirectoryModel2(), add check when the file ends with ~.
   10. Jan 18,2007. in MultiModel.h, add map<int,double> m_BBRModelThresh to keep track of thresh value for each BBRModel.
   11. Jan 18,2007. in parseSingleModelFile2(string,string&), keep <modelID,threshval> info in m_BBRModelThresh map.
   12. Jan 18,2007. in sparseOutput(), when outputing for printpred, add class choosing ability for BBRtrain model.
   13. Jan 18,2007. in parseSingleModelFile2(string,string&), when reading the threshold value, if it is 0.0, default it to 0.5. In order to be consistent with BBRclassify.
   14. Jan 18,2007. in parseSparseModel(string), add parse ability for REFCLASS line; 
   15. Jan 18,2007. change THM data structure in MultiModel.h; from <string,vector<tuple<int,int,double> > > to <string,map<pair<int,int>,double>;
   16. Jan 18,2007. modify addEntry method accordingly
   17. Jan 18,2007. modify getProbabilities method accordingly
   Jan 18,2007. add (mid,cid,fid) check, the error check 4 for sparse model file. 
   Jan 18,2007. no -R option
   SL 22Jan07. modify readLabelFile method; add refclsline; 
   SL 22Jan07. modify readLabelFile method; add error checks;
   SL 23Jan07. get rid of -m option; change parseDirectoryModel2(string) into parseModelFileDir(string); 
   SL 23Jan07. get rid of oldbbr code.		   
   SL 23Jan07. split code into several source files.
   SL 23~28Jan07.  add refclassline for readTestVectorFile; 
                   enable true class labels in testvector file; 
                   update command line options according to new specs;
		   update command line option check according to new specs;
		   modify result file format 1 output;
		   add tie situation handling in result file format output;
		   the default value for args.criterion was set during sanitycheck, instead of in cmd class;
		   result file format 2 was added;
		   result file format 2 error check (1 model file && numeric class name) in addEntry();
		   TOUTPUT data structure changed: add maxPCid, maxLinear (equivalent to maxPPred currently);
		   precision set to 8 in genResultFile();
		   modelname,classname,featurename,docIDname are not case sensitive any more;
		   read label file after reading model file; so that label file error check 2 could be done while reading;
   SL 30Jan07~31Jan07. Windows version modification and testing. it works on the lab machines with .net 2003 installed.
   SL 1Feb07. clean codes in readTestVector.cpp; modify the processData() for empty line reading; modify the getProbability()
              to get rid of char[]; add comments.
   SL 2Feb07. clean readLabelFile.cpp, add cntReadLabelFile() for different way of label file processing.	      
   SL 2Feb07. clean codes, add comments.
   SL 2Feb07. allow #ID: and # ID: for idline in testvector file reading
   SL 12-16Feb07. modify situations when no output for -linear, but still select one for -maxlinear;
                  add sparse header;
                  labelfile reading plus memory control. (not tested yes)
                  integerclass option (still add to string-int map)
   SL 20Feb07. add error handling for 0 model file reading.
               modify Windows directory reading, skip files ending with ~.
               bug: in parseSparseModel file, change the addExtraClasses to addEntry for the models only in this sparse 
                    model file.
   add stdin read for data file. ver0.10		    
		  
 */



#include "MultiModel.h"
#include <stdio.h>

string ModelDB::scales[] = {"prob", "logprob", "odds", "logodds", "linear", "all"};
string ModelDB::criteria[] = {"all", "prob", "logprob", "odds", "logodds",
			      "linear", "maxlinear", "maxlinearall", "maxproball", "predicted"};

void ModelDB::parseArgs(int argc, char ** argv, Args & args){
  // handle arguments
  // Reading arguments
  CmdLine cmd(argv[0], "Bayesian Multi-class Regression - classifying", VERSION );

  // X.A. The -k / --classic option
  // indicates that we should produce a result file in dense format (Format 2, section VIII)
  // instead of our default sparse results file (Format 1, described in section VII).
  SwitchArg classicResultFileArg("k", "classic", "Dense Result File Format (Format 2)", false);
  cmd.add(classicResultFileArg);
  
  // X.B. The -c / --criterion option
  // This option specified which (DOCID, MODELNAME, CLASSNAME) tuples to include in the results.
  // The possible arguments are: 
  // --criterion all: No restrictions, i.e., write one results line for each (DOCID,MODELNAME,CLASSNAME) combinations
  // --criterion prob: write one results line for each (DOCID,MODELNAME,CLASSNAME) combination such that the 
  //                   estimated probability pk of MODEL/CLASS for DOCID is >=VAL, where VAL is specified by 
  //                   --includethresh.
  // --criterion logprob: write one results line for each combination with ln(pk)>=VAL
  // --criterion odds: write out one results line for each combination where odds of CLASS (vs. union of all classes),
  //                   i.e., pk/(1-pk)>=VAL
  // --criterion logodds: write out one results line for each combination with ln(pk/(1-pk))>=VAL
  // --criterion linear: write one results line for each combination where CLASSNAME has a linear score >=VAL
  // --criterion maxlinear: for each (DOCID, MODELNAME) combination, exactly one dataline is written out. The line
  //                        is for the class with the highest linear score for this MODEL. In case of ties on linear
  // 									      score, we choose the predicted class.
  // --criterion maxlinearall: for each (DOCID, MODELNAME) combination, one dataline is written out for every CLASSNAME
  //                           that has a linear score equal to the highest linear score seen for this model.
  // --criterion maxproball: for each (DOCID, MODELNAME) combination one dataline is written out for every CLASS that
  //                         has a class probability equal to the highest class probability seen for this combination.
  //                         This differs from "--criterion maxlinearall" only when there is different linear scores
  //                         that saturate to the same probability.
  // --criterion predicted [DEFAULT]: for each (DOCID, MODELNAME) combination, exactly one dataline is written out,
  //                                  for the CLASSNAME that is predicted for that DOCID (see section II). under 
  //                                  current specifications this gives identical behavior to "--criterion maxlinear",
  //                                  however, that may change in the future (e.g. if we add utility specifications).
  string desc = "output by criterion. Output a line only if ";
  desc       += "the criterion >= VAL, which is specified by '-t";
  desc       += "/--includethresh' option. Possible values: ";
  for(int i = 0; i < sizeof(ModelDB::criteria)/sizeof(string); i++){
    desc = desc +"'"+ ModelDB::criteria[i] +"', ";
  }
  desc  += "Default is 'predicted'.";
  ValueArg<string> criterionArg("c", "criterion", desc, false, "", "output criterion");
  cmd.add(criterionArg);
  
  // X.C. The -t / --includethresh option
  // This option takes the form --includethresh VAL
  // where VAL (which should be numeric) is intended as a threshold to be used with scores of the type specified by 
  // --criterion.
  // If one of those five arguments is specified to --criterion, but --includethresh is not specified, then we use
  // a default value of the inclusion threshold as follows:
  //    --criterion prob:    --includethresh default VAL is 1/(number of classes)
  //    --criterion logprob: --includethresh default VAL is ln(1/(number of classes))
  //    --criterion odds:    --includethresh default VAL is 1/(number of classes - 1)
  //    --criterion logodds: --includethresh default VAL is ln(1/(number of classes - 1))
  //    --criterion linear:  --includethresh default VAL is 0.0
  desc  = "Threshold value for criterion. Works only if '-c/--criterion'";
  desc += "option is set to 'prob', 'logprob', 'odds', 'logodds' or 'linear'. ";
  desc += "Output a line whose criterion is greater than or equal to <threshold>";
  ValueArg<double> threshArg("t", "includethresh", desc, false, 0, 
			     "threshold value (only valid with '-c')");
  cmd.add(threshArg);
  
  // X.D. The -L / --labels option 
  // The option "--labels PATHNAME" specifies that the PATHNAME is a label file containing the true labels for the 
  // test vectors. If a label file is specified, then labels in the test vector file are ignored. 
  ValueArg<string> labelArg("L", "labels", "Label file", false, "", "label file");
  cmd.add(labelArg);


  // X.E. The -o / --omit option 
  // This option takes the form "--omit CLASSNAME"
  // It forbids any result line where the predicted class is a CLASSNAME from being written. 
  // The usual situation this will be used is when we have a large set of binary logistic regression models
  // where the classes are, for instances "0" and "1" for each model, e.g.
  //         Politics: 0, 1
  //         Sports  : 0, 1
  //         Economics:0, 1
  //         etc.
  // Then we might specify --omit 0
  // so that we only get the results lines for the classes (assumed to be "1" in all the above cases) that indicates
  // the test vector has certain content (as opposed to not having certain content). The --omit restriction is 
  // compatible with all --criterion restrictions, and should be viewed as being an additional filter applied 
  // in addition to whatever restriction is imposed by --criterion; tuples must pass both tests to be written out.
  ValueArg<string> omitArg("o", "omit", "Suppress a class label from being output", false, "", "class label");
  cmd.add(omitArg);

  
  // X.H. The -s / --scale option
  // The --scale option specifies which thransformation(s) of the scores we include in the sparse results file:
  // --scale prob: ALL RESULTS lines have the ONE_SCALE_LINE format. SCORE is the estimated probability, pk, of 
  //               membership of example DOCID in class CLASSNAME foe model MODELNAME.
  // --scale logprob: ALL RESULTS lines have ONE_SCALE_LINE format. SCORE is ln(pk), i.e., estimated log prob of 
  //                  membership.
  // --scale odds: ALL RESULTS line have ONE_SCALE_LINE format. SCORE is pk/(1-pk), i.e., estimated odds of membership
  // --scale logodds: ALL RESULTS lines have ONE_SCALE_LINE format. SCORE should is ln(pk/(1-pk)), ie, estimated 
  //                  logodds of membership.
  // --scale linear: ALL RESULTS lines have the ONE_SCALE_LINE format. SCORE is the linear score for the DOCID, MODEL,
  //                 CLASS triple, ie, the dot product of the beta vector for CLASS with the beta vector for DOCID. 
  //                 Note if there are exactly two classes, and the other classes has a beta vector which is all 0s,
  //                 then the linear score for CLASS is the same as the log odds for DOCID to belong to CLASS.
  // --scale all: ALL RESULTS lines have the ALL_SCALE_LINE format.
  //
  // If --scale is omitted, we look at --criterion and set the scale as if the following --scale commands had been used
  //      --criterion prob       ==> --scale prob
  //      --criterion logprob    ==> --scale logprb
  //      --criterion odds       ==> --scale odds
  //      --criterion logodds    ==> --scale logodds
  //      --criterion maxproball ==> --scale prob
  //      --criterion linear     ==> --scale linear
  //      --criterion maxlinear  ==> --scale linear
  //      --criterion maxlinearall ==> --scale linear
  //      --criterion predicted  ==> --scale linear
  //      --criterion not specified: defaults to --criterion predicted, so --scale linear.
  desc = "Choose the type of score(s) to output. Options are ";
  for(int i = 0; i < sizeof(ModelDB::scales)/sizeof(string); i++){
    desc = desc +"'"+ ModelDB::scales[i] +"', ";
  }
  desc += "Default is the same as 'criterion'. Only valid with ";
  desc += "sparse result file format. See manual for detailed restrictions. ";
  ValueArg<string> resultFileFmtArg("s", "scale", desc, false, "", "result file format");
  cmd.add(resultFileFmtArg);
  
  // X.F. The -j / --printlabel option
  // If --printlabel is specified we include the correct lable in each line in sparse format results files. 
  // If the correct label cannot be determined for a DOCID/MODELNAME combination, we use "@unknown" as the label.
  // If --classic is specified, then we always print labels.
  SwitchArg printTrueClassArg("j", "printtrue", "Print true class label", false);
  cmd.add(printTrueClassArg);


  // X.G. The -p / --printpred option
  // If --printpred is specified we include the predicted class for a model (section II) in each line in each sparse 
  // result file line for that model.
  // this option is redundant if --classic is specified.
  SwitchArg printPredClassArg("p", "printpred", "Print predicted class label", false);
  cmd.add(printPredClassArg);


  // X.I.  The -I / --integerclasses Option
  // This option changes the usual procedure for comparing NAMEs in certain limited cases.  
  // If this option is specified then we treat class names as being integers.  
  // We consider two strings to represent the same classname if they are numerically equal after being read by the standard string to integer conversion in C++.  
  // So +1 and 1 are the same name, as are +0, -0, and 0. (specs-20070212)
  SwitchArg intclsArg("I", "integerclasses", "classes labels are all integers", false);
  cmd.add(intclsArg);

  
  // results file
  ValueArg<string> outfileArg("r", "results", "results file", false, "", "results file");
  cmd.add(outfileArg);

  // error file
  ValueArg<string> errorfileArg("", "trace", "trace file", false, "", "trace file");
  cmd.add(errorfileArg);

  // test vector file
  UnlabeledValueArg<string>  datafileArg("datafile", "Data file; '-' for stdin","", "datafile"); 
  cmd.add( datafileArg );

  // model file / model directory
  UnlabeledValueArg<string>  modelfileArg("modelfile","Model file","", "modelfile"); 
  cmd.add( modelfileArg );


  // Parse command line
  cmd.parse(argc, argv);
  
  args.classic = classicResultFileArg.getValue();
  args.criterion = criterionArg.getValue();
  args.threshold = threshArg.getValue();
  args.scale = resultFileFmtArg.getValue();
  args.printTrue =  printTrueClassArg.getValue();
  args.printPred = printPredClassArg.getValue();
  args.omit = omitArg.getValue();
  args.labelFile = labelArg.getValue();
  args.modelFile = modelfileArg.getValue();
  args.dataFile = datafileArg.getValue();
  args.outFile = outfileArg.getValue();
  args.errorFile = errorfileArg.getValue();
  args.integerclasses = intclsArg.getValue();

}


void ModelDB::sanityCheck(Args & args){

    // Anomaly checking for -k / --classic option
    if (args.classic) {
	// it is an error to specify any of the following options with --classic:
	// --criterion, --includethresh, --omit, or --scale
	// ERROR MESSAGE 1
	string emsg = "";
	if(args.threshold != 0)
	    emsg += "'-t/--includethresh' is not allowed with '-k/--classic'.\n";
	if(args.scale != "")
	    emsg += "'-s/--scale' is not allowed with '-k/--classic'.\n";
	if(args.criterion != "")
	    emsg += "'-c/--criterion' is not allowed with '-k/--classic'.\n";
	if(args.omit != "")
	    emsg+= "'-o/--omit' is not allowed with '-k/--classic'.\n";
	if(emsg != "")
	    throw MyException(emsg);

	// it is redundant to specify --printlabel or --printpred with --classic.
	// ERROR MESSAGE 2
	if(args.printTrue || args.printPred)
	    cerr << "The options '-j/--printtrue' and '-p/--printpred' are redundant with '-k/--classic'" << endl;
    }

    // ERROR MESSAGE 3
    if(args.criterion != ""  &&  getIntCriterion(args.criterion) == -1)
	throw MyException("Value for '-c/--criterion' not valid. ");
  
    // make sure the range of specified value is valid with the selected criteria.
    // ERROR MESSAGE 4
    if(args.criterion == "prob" && (args.threshold < 0 || args.threshold > 1)) {
	char buf[300];
	sprintf(buf,"With prob class probability as an inclusion criterion, --includethresh must specify a value between 0 and 1 (inclusive), not %f",args.threshold);
	throw MyException(buf);
    }
  
    // ERROR MESSAGE 5
    if(args.criterion == "logprob" && args.threshold > 0) {
	char buf[300];
	sprintf(buf,"With log class probability as an inclusion criterion, --includethresh must specify a value <=0, not %f", args.threshold);
	throw MyException(buf);
    }

    // ERROR MESSAGE 6
    if(args.criterion == "odds" && args.threshold < 0 ) {
	char buf[300];
	sprintf(buf,"With odds as an inclusion criterion, --includethresh must specify a value >= 0, not %f", args.threshold);
	throw MyException(buf);
    }

    // check whether --criterion options are compatible with --includethresh
    // ERROR MESSAGE 7
    if(args.threshold!=0.0 && args.criterion!="prob" && args.criterion!="logprob" 
         && args.criterion!="odds" && args.criterion!="logodds" && args.criterion!="linear") {
	throw MyException("The '-t/--includethresh' option may be used if and only if one of the following option is also specified for '-c/--criterion': prob, logprob, odds, logodds, or linear");
    }

    // Anomaly checking for -l/--labels option
    // ERROR MESSAGE 8
    if(args.labelFile!="" && !args.classic && !args.printTrue) {
	cerr << "WARNING: '-l/--labels' was specified but you did not specify including the labels in the results." 
             << endl;
    }

    // check --scale value
    if(args.scale == "") {
	if( args.criterion == "prob" || args.criterion == "logprob" 
         || args.criterion == "odds" || args.criterion == "logodds" )
	    args.scale = args.criterion; 
	else if (args.criterion == "maxproball" || args.criterion == "all")
	    args.scale = "prob";
	else
	    args.scale = "linear";
    }
    else{
	// ERROR MESSAGE 9
	if(getIntScale(args.scale) == -1)
	    throw MyException("Invalid -s parameter.");	
    }

    // set the default value for args.criterion if not specified
    if (args.criterion == "" && !args.classic)
	args.criterion = "predicted";

   
}

// convert string modelID and classID into integer, to save some space.
int ModelDB::getIntID(const string & stringID, THKL & keys, TVKRL & reverseKeys, int & max_ID){
    int iIntID;
    THKL::iterator keyPair;
    keyPair = keys.find(stringID);
    if(keyPair != keys.end()){
	iIntID = keyPair->second;
    }else{
	iIntID = max_ID;
	reverseKeys.push_back(stringID);
	keys[stringID] = max_ID++;
    }
    return iIntID;
}


string & ModelDB::getModelStringID(int intID) {
    return reverseModelKeys[intID];
}

string & ModelDB::getCateStringID(int intID) {
    return reverseCateKeys[intID];
}


// constructor.
// read the model files, and do post check.
ModelDB::ModelDB(const Args & args){

    // set the error output
    if (args.errorFile != "") {
	errout = new ofstream(args.errorFile.c_str());
	if(!((ofstream *)errout)->is_open()) {
	    cerr<<"error in opening " << args.errorFile <<endl;
	    exit(1);
	}
    }
    else
	errout = & cout;

    // parsing file, building the hashtable. To be written after the 
    // file format is received. -- bing bai 06/26/2006.
    this->args = args;
    maxModelID = 0;
    maxCateID = 0;
    maxlabelnumber = 0;
    totalclasses = 0;
    lowmemory = false;
    emergencypool = new char[512000];
    bDir = true;

    // parse model files based on the path given in commandline;
    // if the path is a directory, read each model file under that directory;
    // if the path is a file, read the one model file only.

    // For each model file
    // if it begins with "Bayesian Binary Regression", then it is a BBR model file;
    // if it begins with "Bayesian Multinomial Regression" and other BMR header, then it is a BMR model file.
    // otherwise, it is a sparse model file.

    // NOTE: if a directory is specified, each file in the directory MUST be a model file.
    parseModelFileDir(args.modelFile);

    // check whether each model contains more than 1 classes
    postSanityCheck();

    // originally, the label file is read here, by calling readLabelFile(args.labelFile).
    // according to new specs, the label file is read with readTestVectorFile by calling cntReadLabelFile;
}


ModelDB::~ModelDB() {
    if(!lowmemory)
	delete[] emergencypool;

    if(args.errorFile!="") {
	((ofstream *)errout)->close();
	delete errout;
    }
}


// currently only one mission, see if there is a model with only one class.
void ModelDB::postSanityCheck(){
    if(acm.size()==0) {
	// ERROR MESSAGE 10
	string emsg = "Error! No model file existing in the directory: " + args.modelFile;
	throw MyException(emsg);
    }
  for(TACM::iterator acmit=acm.begin(); acmit != acm.end(); acmit++){
    TCM & t = acmit->second;
    if(t.size() == 1){
	// ERROR MESSAGE 11
      string emsg = "Error! Model has only one class: "+ getModelStringID(acmit->first);
      if(args.errorFile!="") (*errout)<<emsg<<endl;
      throw MyException(emsg);
    }      
  }
}



void ModelDB::clearEntryMarks(TMMM & v){
  TMMM::iterator vit;
  for(vit = v.begin(); vit != v.end(); vit++){
    TCM::iterator tit = *vit;
    TST & t = tit->second;
    t.innprod = 0;
    t.odds = 0;
    t.prob = 0;
    t.logodds = 0;
    t.logprob = 0;    
    t.used = false;
    sumIP[t.mid] =0;
  }  
}

// called only once. using the slow mode
int ModelDB::getIntCriterion(string & crit){
  for(unsigned int i = 0; i < sizeof(criteria)/sizeof(string); i++){

    if(crit == criteria[i])	return i;
  }
  return -1;
}

int ModelDB::getIntScale(string & sc){
  for(unsigned int i = 0; i < sizeof(scales)/sizeof(string); i++){
    if(sc == scales[i]) return i;
  }
  return -1;
}

// change the string into lowercase
string ModelDB::toUpperCase(const string &s) {
    string result;
    for (int i = 0; i < s.length(); i++) {
	result += toupper(s[i]);
    }
    return result;
}

// get the integer value for a string
// atoi will return 123 for "123ert"
// while we want error for "123ert".
bool ModelDB::intvalue(const string &s) {

    int x = 0;
    while(s[x]!='\0') {
	if((int(s[x])<48 || int(s[x])>57) && s[x]!='.' && s[x]!='+' && s[x]!='-' && s[x]!='e' && s[x]!='E')  
	    return false;
	x++;
    }
    return true;
}
