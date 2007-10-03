#include "MultiModel.h"


int main(int argc, char ** argv){
  try{

    Args args;

    // Generate argument descriptions, parse arguments, perform
    // first-level argument checking.
    ModelDB::parseArgs(argc, argv, args);

    // sanity check, see format.txt for rules. Give higher-level
    // perform second-level argument checking.
    ModelDB::sanityCheck(args);

    ModelDB mdb(args);
    cerr << "index building finished" << endl;

    mdb.processData();

  } 
  catch (ArgException e)  {
      cerr << "***Command line error: " << e.error() 
	   << " for arg " << e.argId() << endl;
      return 1;
  }
  catch (MyException &e){
      cerr << "***Invalid options: " << e.error() << endl;
      return 1;
  }  
  catch( std::exception& e){
      //    Log(0)<<std::endl<<"***Exception: "<<e.what();
      cerr<<std::endl<<"***Exception: "<<e.what();
      return 1;
  }
  catch(...){
      //    Log(0)<<std::endl<<"***Unrecognized exception";
      cerr<<std::endl<<"***Unrecognized exception";
      return 1;
  }

  return 0;
}

