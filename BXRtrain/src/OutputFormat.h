#ifndef OUTPUT_FORMAT_H
#define OUTPUT_FORMAT_H

void displayConfusionTable( ostream& o, const RowSetMetadata& names, const vector< vector<unsigned> >& CT );
void makeConfusionTable( ostream& o, const RowSetMetadata& names, 
                           const vector<unsigned>& y, const vector<unsigned>& prediction );


// TODO: Why are there two such functions?
void displayCT2by2( ostream& o, int TP, int FP, int FN, int TN, double roc );
void makeCT2by2( ostream& o, const RowSetMetadata& names, 
                           const vector<unsigned>& y, 
                           const vector< vector<double> >& allScores,
                           const vector<unsigned>& prediction);


#endif // OUTPUT_FORMAT_H

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
