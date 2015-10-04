/*************
 **Author: Xuehui Li
 **Date: May, 2005
 ** This program is used to sort 'gbaRecall'( or 'gbaPrecision' ), and change 'cardRecall'( or 'carPrecision') and 'segRecall'( or 'segPrecision') so that the i th entry in 'gbrRecallSorted' and the i th entry in 'cardRecallChanged' as well as 'segRecallChanged' correspond to the same sequence. 
 ** There three input parameters. After programm running, three output files are generated.
 ** An example of running the program: 
java applications.PrecisionRecallSorter /cise/research/tamer/xli/LCR/graphLCR/metrics/swissprot/wForgetRate/postProcessing/mim/sort/recall/gbaRecall  /cise/research/tamer/xli/LCR/graphLCR/metrics/swissprot/wForgetRate/postProcessing/mim/sort/recall/cardRecall /cise/research/tamer/xli/LCR/graphLCR/metrics/swissprot/wForgetRate/postProcessing/mim/sort/recall/segRecall
** Three files will be generated in the same folder: gbaRecallSorted, cardRecallChanged and segRecallChanged.
**************/

package applications;

import java.io.*;
import java.util.*;

class PrecisionRecallSorter {

  private RandomAccessFile gbaRfIn;
  private RandomAccessFile cardRfIn;
  private RandomAccessFile segRfIn;
  private RandomAccessFile gbaRfOut;
  private RandomAccessFile cardRfOut;
  private RandomAccessFile segRfOut;


  public PrecisionRecallSorter( String gba, String card, String seg ) {
    try {
      File f = new File( gba );
      gbaRfIn = new RandomAccessFile( f, "r" );
      f = new File( card );
      cardRfIn = new RandomAccessFile( f, "r" );
      f = new File( seg );
      segRfIn = new RandomAccessFile( f, "r" );
      f = new File( gba + "Sorted" );
      gbaRfOut = new RandomAccessFile( f, "rw" );
      f = new File( card + "Changed" );
      cardRfOut = new RandomAccessFile( f, "rw" );
      f = new File( seg + "Changed" );
      segRfOut = new RandomAccessFile( f, "rw" );
    }
    catch ( IOException ex ) {
    }
  }

  
  public Vector readInput( RandomAccessFile rf ) {
    Vector vc = new Vector();
    String line = new String();
    try {
      line = rf.readLine();
      while ( line != null ) {
	line = line.trim();
	vc.add( line );
	line = rf.readLine();
      } 
    }
    catch( IOException ex ) {
    }
    return vc;
  }


  public void sortGba() {
    Vector gbaVc = new Vector();
    Vector cardVc = new Vector();
    Vector segVc = new Vector();
    gbaVc = readInput( gbaRfIn );
    cardVc = readInput( cardRfIn );
    segVc = readInput( segRfIn );
    String valStr = new String();
    double valDb = 0, max = -1;
    int i = 0, maxIndex = 0;
    while ( !(gbaVc.isEmpty() )) {
      i = 0;
      max = -1;
      while ( i < gbaVc.size() ) {
	valStr = ( String ) gbaVc.elementAt( i );
	valDb = Double.parseDouble( valStr );
	if ( valDb > max ) {
	  max = valDb;
	  maxIndex = i;
	}
	++i;
      }
      try{
	gbaRfOut.writeBytes( (String)gbaVc.remove( maxIndex ) + "\r\n" );
	cardRfOut.writeBytes( (String)cardVc.remove( maxIndex ) + "\r\n" );
	segRfOut.writeBytes( (String)segVc.remove( maxIndex ) + "\r\n" );
      }
      catch ( IOException ex ) {
      } 
    }
  }
  

  public void closeAll( ) {
    try{
      gbaRfIn.close();
      cardRfIn.close();
      segRfIn.close();
      gbaRfOut.close();
      cardRfOut.close();
      segRfOut.close();
    }
    catch ( IOException ex ) {
    }
  }


  public static void main ( String args[] ) {
    PrecisionRecallSorter prs = new PrecisionRecallSorter( args[0], args[1], args[2] );
    prs.sortGba();
    prs.closeAll();
  }
  
}
