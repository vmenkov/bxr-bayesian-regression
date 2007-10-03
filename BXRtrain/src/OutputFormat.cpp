#include "PolyZO.h"
#include "OutputFormat.h"

void displayConfusionTable( ostream& o, const RowSetMetadata& metadata, const vector< vector<unsigned> >& CT )
{
    unsigned nAll=0, nCorrect=0;
    o<<endl<<setw(20)<<"Belongs\\Assigned ";
    for( unsigned k=0; k<metadata.getClassCount(); k++ )
        o<<setw(8)<<metadata.getClassId(k);

    for( unsigned k=0; k<metadata.getClassCount(); k++ ) {
        o<<endl<<setw(20)<<metadata.getClassId(k);
        for( unsigned r=0; r<metadata.getClassCount(); r++ ) {
            o<<setw(8)<<CT[k][r];
            nAll += CT[k][r];
            if( k==r )
                nCorrect += CT[k][k];
        }
    }
    o<<endl<<"Total errors: "<<nAll-nCorrect<<" out of "<<nAll<<", "<<100.0*(nAll-nCorrect)/nAll<<"%";
}
void makeConfusionTable( ostream& o, const RowSetMetadata& names, 
                           const vector<unsigned>& y, const vector<unsigned>& prediction )
{
    vector< vector<unsigned> > CT( names.getClassCount(), vector<unsigned>( names.getClassCount(), 0 ) );
    for( unsigned i=0; i<y.size(); i++ )
        CT [y[i]] [prediction[i]] ++;

    displayConfusionTable( o, names, CT );
}

void makeCT2by2( ostream& o, const RowSetMetadata& names, 
                           const vector<unsigned>& y, 
                           const vector< vector<double> >& allScores,
                           const vector<unsigned>& prediction)
//						   vector<double>& roc)
{
    vector<unsigned> TP(names.getClassCount(),0), FP(names.getClassCount(),0), FN(names.getClassCount(),0); //, TN(drs.c(),0);
    unsigned n=0;
    for( unsigned i=0; i<y.size(); i++ ) {
        if( prediction.at(i)==y.at(i) )  
            TP.at(y.at(i))++;
        else {
            FN.at(y.at(i))++;
            FP.at(prediction.at(i))++;
        }
        n++;
    }

    for( unsigned k=0; k<names.getClassCount(); k++ ) {
        o<<"\n\nOne-vs-All view: Class "<<names.getClassId(k);
		// TODO: Redo ROC
//        double roc = calcROC( allScores, prediction, k );
        displayCT2by2( o, TP.at(k), FP.at(k), FN.at(k), n-TP.at(k)-FP.at(k)-FN.at(k) /*TN*/, -1 );
    }
}


void displayCT2by2( ostream& o, int TP, int FP, int FN, int TN, double roc )
{
    o<<"\nConfusion matrix:   Relevant\tNot Relevant"
        <<"\n\tRetrieved    \t"<<TP<<"\t"<<FP
        <<"\n\tNot Retrieved\t"<<FN<<"\t"<<TN;
    double precision = (TP + FP)>0 ? 100.0*TP/(TP + FP) : 100;
    o<<"\nPrecision = "<< precision;
    double recall  = (TP + FN)>0 ? 100.0*TP/(TP + FN) : 100;
    o<<"\nRecall  = "<< recall;
    double F1  = (2*TP+FN+FP)>0 ? 100.0*2*TP/(2*TP+FN+FP) : 100;
    o<<"\nF1 = "<< F1;
    o<<"\nROC area under curve "<<roc;
}

