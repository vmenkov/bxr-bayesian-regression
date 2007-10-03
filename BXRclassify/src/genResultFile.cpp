
#include "MultiModel.h"
/*
  Generate the results file.

  for -k/--classic version, use the result file format 2:
       truelabel prob1 prob2 ... probn predictedlabel
    
  if no -k/--classic specfied, use the result file format 1.
     
           FILE:=(COMMENT_LINE | RESULTS_LINE)*
	   COMMENT_LINE:=<whitespace>*#<any character>*
	   RESULTS_LINE:=ONE_SCALE_LINE | ALL_SCALE_LINE
	   ONE_SCALE_LINE:=DOCID MODEL LINECLASS SCORE OTHERCLASSES
	   ALL_SCALE_LINE:=DOCID MODEL LINECLASS PROB LOGPROB ODDS LOGODDS LINEAR OTHERCLASSES
	   DOCID:=NAME
	   MODEL:=NAME
	   LINECLASS:=CLASS [the class that the score in this line is for.]
	   CLASS:=NAME
	   OTHERCLASSES:=PRED_CLASS_ENTRY? TRUE_CLASS_ENTRY?
	   PRED_CLASS_ENTRY:=pred CLASS
	   TRUE_CLASS_ENTRY:=true CLASS


  If --printpred is specified, show PRED_CLASS_ENTRY.
	   
  If --printlabel is specified, show TRUE_CLASS_ENTRY. The following is the logic for the true label reading: 

  /////////////////////////////////////////
  if --lable specified
  * FOR each test vector
      o IF test vector has no DOCID
          + Exit with error:
           	  "All vectors must have DOCIDs when label file is used. Vector ID missing for vector on line LINENUM of 
                   test vector file PATHNAME."
      o ELSE
          + DOCID = document ID of test vector
          + Read through label file until we reach first line for this  DOCID, i.e. "DOCID MODELNAME CLASSNAME" or until 
            we reach end of file.   
            Read and save in memory all triples (possibly none) for this DOCID.  
            Process any reference class line for DOCID in the context of the models used in this run. 
	  + FOR each MODELNAME for which we write results for DOCID
		* IF label file specifies CLASSNAME for pair (DOCID, MODELNAME)
                      o TRUE_CLASS = CLASSNAME
	        * ELSE
                      o TRUE_CLASS = "@unknown"
                * Write results


  if no --label specified		
  * IF there is more than one model loaded
     o Exit w/ error "Can't include true labels in result if there are multiple models but no label file"
  * ELSE
     o FOR each test vector
          + TRUE_CLASS = label that vector has in test vector file
          + Write results
  /////////////////////////////////////	  

  ******************************************************************************/


void ModelDB::sparseOutput(TMMM & v, string & docID, string & trueLabel, int crit, int sc, ostream * pout){
    // set the precision in output to 8
    pout->precision(8);
    cout.precision(8);

    map<int, TOUTPUT> vd;
    TMMM::iterator vit;
    int citPred;
    int omitCid;
    
    // get the integer index for the omit class name
    if(args.omit != "")
	omitCid=getIntID(args.omit, cateKeys, reverseCateKeys, maxCateID);


    // TACM: the two level map: hash_map on mid, then map on cid; keep different scores
    for(TACM::iterator acmit = acm.begin(); acmit != acm.end(); acmit++){
	
	TCM & cm = acmit->second;
	
	// calculate the default threshold for each model, if --includethresh not specified
	if (args.threshold == 0.0 && args.criterion!="") {
	    int nClasses = cm.size();
	    double probthresh = (double)1/(double)nClasses;
	    double oddsthresh = (double)1/(double)(nClasses-1);
	    switch(crit) {
		case CT_prob:       args.threshold = probthresh; break;
		case CT_logProb:    args.threshold = log(probthresh); break;
		case CT_odds:       args.threshold = oddsthresh; break; 
		case CT_logOdds:    args.threshold = log(oddsthresh); break;
		case CT_linear:     args.threshold = 0.0; break;
	    }
	} 


	TOUTPUT & ot = ((vd.insert(make_pair(acmit->first, TOUTPUT()))).first)->second;
	ot.used = false;
	ot.maxP = 0.0;

	// for each class of the model, calculate the different scores
	for(TCM::iterator cmit = cm.begin(); cmit != cm.end(); cmit++){
	    updateOutputVector(cmit, crit, vd); 
	}

    } // end of processing every model


    // generate the results line(s)

    map<int, TOUTPUT>::iterator toit;

    for(toit = vd.begin(); toit != vd.end(); toit++){
	// --common part first
	int mid = toit->first;
	TOUTPUT & to = toit->second;
	
	// if --classic is specifed, output the truelabel first
	if(args.classic) 
	    (*pout) << trueLabel << " ";
	
	vector<TCM::iterator>::iterator vcitit;
	for(vcitit = to.v.begin(); vcitit != to.v.end(); vcitit++){
	    TST & st = (*vcitit)->second;
	    
	    // omit the class;
	    if(args.omit!="" && st.cid == omitCid)
		continue;
	    
	    // if --classic is specified, output the prob of this class, then process the next class.
	    if(args.classic) {
		(*pout) << st.prob << " ";
		continue;
	    }

	    // if no --classic specified, 
            // FIRST output the "DOCID MODEL LINECLASS"
	    (*pout) << docID << " " << getModelStringID(mid) << " " << getCateStringID(st.cid) << " ";     
      
	    // ALL_SCALE_LINE: "PROB LOGPROB ODDS LOGODDS LINEAR"
	    if(sc == SC_all){	  
		(*pout) << st.prob << " " << st.logprob << " "
			<< st.odds << " " << st.logodds << " " << st.innprod << " ";
	    }
	    // ONE_SCALE_LINE: "SCORE"
	    else{   
		double score;
		switch(sc){
		    case SC_prob:	score = st.prob;	break;
		    case SC_logProb:	score = st.logprob;	break;
		    case SC_odds:	score = st.odds;	break;
		    case SC_logOdds:	score = st.logodds;	break;
		    case SC_linear:	score = st.innprod;	break;
		}
		(*pout) << score << " " ;
	    } 
      
	    // PRED_CLASS_ENTRY?
	    // if the model is a BBRtrain model, use thresh value for predict label.
	    if(args.printPred){
		// search the mid in the BBRtrain model id list
		map<int, double>::iterator mit_BBRModelThresh = m_BBRModelThresh.find(mid);
		// if this is a BBRtrain model, 
		// the predicted label is "1" if the score is larger than or equal to the threshold value; 
		// otherwise, the predicted label is "-1"
		if (mit_BBRModelThresh != m_BBRModelThresh.end()) {
		    double ppredCls1;
		    string cid = getCateStringID(to.maxPPredCid);
		    // get the ppred for class +1
		    ppredCls1 = (cid=="1"||cid=="+1") ? to.maxPPred : 1-to.maxPPred;
		    if (ppredCls1 >= mit_BBRModelThresh->second)
			(*pout) << "pred:" << "1" << " ";
		    else
			(*pout) << "pred:" << "-1" << " ";
		}
		// if this is not a BBRtrain model
		else
		    (*pout) << "pred:" << getCateStringID(to.maxPPredCid) << " ";

	    }
	    
	    // TRUE_CLASS_ENTRY?
	    if(args.printTrue){
		
		string trueCate = "@unknown";

		// if label file is specified, use the map<mid,cid> of the current test vector to get the true label
		if(args.labelFile!="") {
                    
		    THLL::iterator ltit = labelTable.find(toUpperCase(docID));
		    if(ltit != labelTable.end()){
			vector<TID> & vl = ltit->second;
			TID ttid(mid, 0);
			vector<TID>::iterator vtit = lower_bound(vl.begin(), vl.end(), ttid, less<TID>());
			if(vtit != vl.end() && vtit->modelID == mid){
			    trueCate = getCateStringID(vtit->cateID);
			}
		    }
		    else{
		    }
		    
                    /*
		    // search the map<mid,cid> for the true label for the mid of the current test vector 
		    map<int,int>::iterator mmit = m_trueLabelForOneVector.find(mid);
		    if(mmit!=m_trueLabelForOneVector.end()) {
			    trueCate = getCateStringID(mmit->second);
		    }
		    */
		}
		// if no label file is specified, and there is only one model, use the true label in test vector file.
		else {
		    if(acm.size()>1) {
			// ERROR MESSAGE 53
			string emsg = "Cannot include true labels in results if there are multiple models but no label file.";
			if(args.errorFile!="")
			    (*errout)<< emsg << endl;
			throw MyException(emsg);
		    }
		    else
			trueCate = trueLabel;
		}
		
		(*pout) << "true:" << trueCate<<" ";
	    }

	    (*pout) << endl;

	} // end of processing all the classes for this model
	
	// if -k/--classic specified, output predclassname
	if(args.classic)
	    (*pout) << getCateStringID(to.maxPCid) << endl;	
    } // for the function
}


void ModelDB::updateOutputVector(TCM::iterator & tit, int crit, map<int, TOUTPUT> & vd){

  TST & t = tit->second;
//  TOUTPUT & ot = 
//    ((vd.insert(make_pair(t.mid, TOUTPUT()))).first)->second;
//  ot.used = false;
  TOUTPUT & ot = (vd.find(t.mid))->second;

  double linear = t.innprod;
  double prob = t.prob;
  double logprob = t.logprob;
  double odds = t.odds;
  double logodds = t.logodds;

  // for classic output
  if(args.classic) {
      if(!ot.used) {
	  ot.maxP = prob; ot.maxPCid = t.cid; ot.used = true;
      }
      else if(prob > ot.maxP){
	  ot.maxP = prob;
	  ot.maxPCid = t.cid;
      }
      // if there is a tie, use lexicographic order
      else if (prob == ot.maxP) {
	  if (getCateStringID(ot.maxPCid) < getCateStringID(t.cid))
	      ot.maxPCid = t.cid;
      }

      if(prob < args.threshold) 
	  return;
      ot.v.push_back(tit); 

      return;
  }



  switch(crit){
      // Dave said: currently predicted is the maximum probability,
      // may be changed in the future -- Bing Bai
      // did change: predicted is the maximum linear now. -- shenzhi li (25Jan07)

      // please see MultiModel.cpp :: parseArgs for more details
      // --criterion all: no restrictions
      case CT_all:   
	  ot.v.push_back(tit);    
	  break;	  
      // --criterion prob: write out result line only when prob >=threshold	  
      case CT_prob:	
	  if(prob < args.threshold) 
	      return;
	  ot.v.push_back(tit); 
	  break;
      // --criterion logprob: write out result line only when ln(prob) >=threshold	  
      case CT_logProb:	
	  if(logprob < args.threshold)  
	      return;
	  ot.v.push_back(tit); 
	  break;
      // --criterion odds: write out result line only when odds >=threshold	  
      case CT_odds:	
	  if(odds < args.threshold)  
	      return;
	  ot.v.push_back(tit);	
	  break;
      // --criterion logodds: write out result line only when ln(odds) >=threshold	  
      case CT_logOdds:	
	  if(logodds < args.threshold)	  
	      return;
	  ot.v.push_back(tit);	
	  break;
      // --criterion linear: write out result line only when linear >=threshold	  
      case CT_linear:   
	  if(linear < args.threshold)    
	      return;
	  ot.v.push_back(tit); 
	  break;

      // --criterion maxlinear: write out one result line for each (docid,mid) pair. if there is a pair, choose the earliest one.

      case CT_maxlinear:
      case CT_predicted:
	  if(!ot.used) {
	      ot.used = true;
	      ot.maxLinear = linear; ot.v.clear(); ot.v.push_back(tit);
	  }
	  else if( linear > ot.maxLinear || (linear == ot.maxLinear 
                 && getCateStringID((tit->second).cid)<getCateStringID((ot.v.at(0))->second.cid) ))  { 
	      ot.maxLinear = linear; ot.v.clear(); ot.v.push_back(tit);
	  }

	  if(args.printPred) { ot.maxPPred = ot.maxLinear; ot.maxPPredCid = t.cid; }
	  break;
      // --criterion maxlinearall: write out all results lines with max score for each (docid,mid) pair. 
      case CT_maxlinearall:
	  if(!ot.used) {
	      ot.used = true;
	      ot.maxLinear = linear; ot.v.clear(); ot.v.push_back(tit);
	      if(args.printPred) { ot.maxPPred = ot.maxLinear; ot.maxPPredCid = t.cid; }
	  }
	  else if( linear > ot.maxLinear ){
	      ot.maxLinear = linear; ot.v.clear(); ot.v.push_back(tit);	  
	      if(args.printPred) { ot.maxPPred = ot.maxLinear; ot.maxPPredCid = t.cid; }
	  }
	  else if(linear == ot.maxLinear) {
	      ot.v.push_back(tit);
	      if(args.printPred) { 
		  string oldmaxpredcls = getCateStringID(ot.maxPPredCid);
		  string newmaxpredcls = getCateStringID(t.cid); 
		  if (newmaxpredcls<oldmaxpredcls)  ot.maxPPredCid = t.cid;
	      }
	  }
	  break;
      // --criterion maxproball: write out all results lines with max score for each (docid,mid) pair. 
      case CT_maxProbAll:
	  if(prob > ot.maxP){
	      ot.maxP = prob; ot.v.clear(); ot.v.push_back(tit); ot.maxPCid = t.cid;
	  }else if(prob == ot.maxP){
	      ot.v.push_back(tit);
	      string oldmaxpredcls = getCateStringID(ot.maxPCid);
	      string newmaxpredcls = getCateStringID(t.cid); 
	      if (newmaxpredcls<oldmaxpredcls)  ot.maxPCid = t.cid;

	  }
	  break;


  } // end of switch



  // use max linear currently for Pred instead of prob
  if(args.printPred){
      // if -c predicted or -c maxlinear or -c maxlinearall specified, the maxPPred and maxPPredCid have been taken care already.
      if(crit != CT_maxlinear && crit != CT_predicted && crit!=CT_maxlinearall) {
	  if(!ot.used) {
	      ot.used = true; ot.maxPPred = linear; ot.maxPPredCid = t.cid; }
	  else if(linear > ot.maxPPred){
	      ot.maxPPred = linear;
	      ot.maxPPredCid = t.cid;
	  }
	  // if there is a tie, use lexicographic order
	  else if (linear == ot.maxPPred) {
	      if (getCateStringID(ot.maxPPredCid) > getCateStringID(t.cid))
		  ot.maxPPredCid = t.cid;
	  }
	  
      }

  }
}

