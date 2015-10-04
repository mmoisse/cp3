/**
 * Date: July, 2005
 * This program is used to generate 100 different random numbers and then 
 * get sequence IDs from the lines whose line numbers are those 100 random 
 * numbers in file 
 * /cise/research38/tamer/xli/LCR/data/swissprot/repeatPercentage/combinedRepPer.  
 * The output is a file that contains sequence IDs:  /cise/research38/tamer/xli/LCR/blastExp/input/seqId.
 * 
**/

package applications;

import java.io.*;
import java.util.*;

class SeqIdGettor {

  private static final int TOTALNUM = 475; 
  private static final int RANDOMNUM = 100;
  Integer[] sortedR = new Integer[ RANDOMNUM];
  RandomAccessFile rf;

  public SeqIdGettor() {
    try {
      File f = new File( "/cise/research38/tamer/xli/LCR/data/swissprot/repeatPercentage/combinedRepPer" );
      rf = new RandomAccessFile ( f, "r" );
    } 
    catch ( IOException ex ) {
    }
  }

  
  public void getRandomNums () {
    Random generator = new Random();
    int r, k = 0;
    Integer tmpIg;   
    Vector randoms = new Vector();
    while ( k <  RANDOMNUM ) {  
      r = generator.nextInt(TOTALNUM);
      if ( r != 0 ) {
        tmpIg = new Integer( r );
        if ( randoms.indexOf( tmpIg ) == -1 ) { // sampling without replacement
	  // System.out.print( r + " " );
	  randoms.add( tmpIg );
          sortedR[k] = tmpIg;
	  ++k;
        }
      }
    }
    //    System.out.println();
    //System.out.println();    
    MergeSort.mergeSort( sortedR ); 
    for ( int i = 0; i <  RANDOMNUM; i++ ) {
      r = sortedR[i].intValue();
      System.out.print( r + " " );
    }
    System.out.println();
  }
  

  public void getSeqIDs() {
    try {
      String seqId = new String();   
      int r, lineNum = 0;
      for ( int i = 0; i <  RANDOMNUM; i++ ) {
	r = sortedR[i].intValue();
	while ( lineNum != r ) {
	  seqId = rf.readLine();
	  ++ lineNum;
	}
	int index = seqId.indexOf( ":" );
	seqId = seqId.substring( 0, index );
	System.out.println( seqId );
      } 
      rf.close();
    }
    catch ( IOException ex ) {
    }  
  }


  public static void main ( String args[] ) {
    SeqIdGettor sig = new SeqIdGettor();
    sig.getRandomNums();
    sig.getSeqIDs();
  }
}
