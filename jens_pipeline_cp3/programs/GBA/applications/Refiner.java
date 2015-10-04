/**
* Author: Xuehui Li
* Date: March, 2005
* This is program is used to refine the low-complexity regions of sequences identified by a LCR-identification algorithm.
* Given a masked LCR, if its length is less then 6 letters, we check its relationship with its neighbors. If the similarity with none of its neighbors is greater than 4 letters or its(the currently considered LCR) length, we remove it from the LCR set. 
* There are two input parameters. The first one is the LCR blocks file which tells which parts of sequences are LCRs identified by an algorithm. The second one is the sequence file.
* Suppose the algorithm is GBA and the sequence file is flybase. An example of running the program is:
* java applications.refiner swissprotLCRBlocks/wForgetRateMatrices/postProcess/cutting/meanStd/flyabseMeanStd ../data/swissprot/sequenceInfor/seqFromFlybase > o
 **/

package applications;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.formats.Pair;
import jaligner.matrix.MatrixLoader;
import jaligner.util.SequenceParser;
import java.io.*;

class Refiner {

  public static final int blockLenThreshold = 6;
  public static final int similarityLenThreshold = 4;

  private RandomAccessFile lcrRf, seqRf;
  private String[] tmpLCRs = new String[5000], LCRs = new String[5000];
  private int tmpIndexLCR = 0, indexLCR = 0;
  
  public Refiner ( String lcrFile, String seqFile ) {
    try {
      File f = new File( lcrFile ); // LCR file
      lcrRf = new RandomAccessFile( f, "r" );
      f = new File (seqFile); // sequence file
      seqRf = new RandomAccessFile (f, "r");
      
    }
    catch ( IOException ex ) {
    
    }  
  }


  /* mark is used to distinguish tmpLCRs (1) and repeats (0). tmp is used to distinguish uncombined repeats (false) and combined repeats (true).
   */  
  public void printLCRsOReps( int mark, boolean tmp ) {
    int limit = 0;
    if ( mark == 1 )
      if (tmp)
	limit = indexLCR;
      else 
	limit = tmpIndexLCR;
    else 
      if (tmp) {
	//limit = indexRepeat;
	System.out.print( "!!" );
      }
      else {
	//limit = tmpIndexRepeat;
	System.out.print( "**");
      }
    for ( int i = 0; i < limit; i++ ) {
      if ( mark == 1 )
	if (tmp) 
	  System.out.print( LCRs[i]+" " );	
	else 
	  System.out.print( tmpLCRs[i]+" " );
      else {
	/* 
	if (tmp)
	  System.out.print( repeats[i]+" " );
	else 
	  System.out.print( tmpRepeats[i]+" " );
	  */
      }      
    }
    System.out.println();
    if ( mark == 0 ) 
      if (tmp) 
	System.out.println( "**************************" ); 
  }

 
 //  mark: 1 for LCRs, 0 for repeats 
  public void getLCRsOReps (String line, int mark) {
   int i, j;
    String str = line;
    str = str.trim();
    if ( str.length() == 0 ) {
      str = null;
      System.out.println( "yeah" );
    }
    while ( str != null ) {
      i = str.indexOf(" ");
      j = str.indexOf("-");
      if ( i == -1 ) {
	if ( mark == 1) 
	  tmpLCRs[ tmpIndexLCR ] = str.substring( 0 );
	else { 
	  // tmpRepeats[ tmpIndexRepeat ] = str.substring( 0 );
	}
	str = null;
      }     
      else {
	if ( mark == 1 )
	  tmpLCRs[ tmpIndexLCR ] = str.substring( 0, i );
	else {
	  // tmpRepeats[ tmpIndexRepeat ] = str.substring( 0, i );
	
	}str = str.substring( i + 1 );     
      }
      if ( mark ==1 )
	++ tmpIndexLCR;
      else {
	//++ tmpIndexRepeat;
      }
    }
    /***********
    if ( mark == 1)
      System.out.println( "unprocessed LCRs:");
    else 
      System.out.println( "uncombined repeats:");
    printLCRsOReps(mark, false);
    *************/
    if (mark ==0) {
      //combineRepeats();
      System.out.println( "combined repeats:");
      printLCRsOReps(mark, true);
    }
  }


  public String getSeq() {
    String seq = new String(), line = new String();
    try {
      line = seqRf.readLine();
      while ( ( line != null ) &&  ( !( line.startsWith( ">" ))) ) {
	line = line.trim();
	seq = seq + line;
	line = seqRf.readLine();
      }
    }
    catch ( IOException ex ) {
    }
    return seq;
  }
  

  public String findAlignment( String seq1, String seq2 ) {
    String aliPos = new String();
    try {
      Sequence s1 = SequenceParser.parse( seq1 );
      Sequence s2 = SequenceParser.parse( seq2 );
      //System.out.println( "align sequences: " + seq1 + "???" + seq2 );
      Alignment alignment = SmithWatermanGotoh.align(s1, s2, MatrixLoader.load("BLOSUM62"), 10f, 0.5f);
      int similarLen = alignment.getSimilarity();  // get the length of the same and similar letters;
      if ((similarLen>similarityLenThreshold)||(similarLen == seq1.length())) { // only If the length of similar and same letters is greater than 4       
	aliPos = new Pair().format( alignment );
	//System.out.println( "the alignment: " +  aliPos + "     " + similarLen );
      } 
    }
    catch (Exception e) {
      //logger.log(Level.SEVERE, "Failed running example: " + e.getMessage(), e);
    }
    //////////////// System.out.println( seq1 + "!!!!"+ seq2 + "***" + aliPos);
    //System.out.println("aliPos: "+ aliPos);
    return aliPos;
  }
  
  
  public void refineLCRs( String seq ) {
    String block1 = new String(), block2 = new String(), lcr1 = new String(), lcr2 = new String();
    int i = 0, index = 0, start1 = 0, end1 = 0, start2 = 0, end2 = 0;
    while ( i < tmpIndexLCR ) {
      block1 = tmpLCRs[i];
      index = block1.indexOf( "-" );
      start1 = Integer.parseInt( block1.substring( 0, index ));
      end1 = Integer.parseInt( block1.substring( index + 1 ));
      if (( end1 - start1 + 1 ) > blockLenThreshold ) {
	LCRs[indexLCR] = tmpLCRs[i];
	++ indexLCR;
      }
      else { 
	lcr1 = seq.substring(start1 - 1, end1);
	if ( i != 0 ) {
	  block2 = tmpLCRs[i - 1];
	  index = block2.indexOf( "-" );
	  start2 = Integer.parseInt( block2.substring( 0, index ));
	  end2 = Integer.parseInt( block2.substring( index + 1 ));
	  lcr2 = seq.substring(start2 - 1, end2);
	  lcr2 = findAlignment( lcr1, lcr2 );
	  if ( lcr2.length() != 0 ) {	
	    LCRs[indexLCR] = tmpLCRs[i];
	    ++ indexLCR;
	  }
	  else if ( i != (tmpIndexLCR-1)) {
	    block2 = tmpLCRs[i + 1];
	    index = block2.indexOf( "-" );
	    start2 = Integer.parseInt( block2.substring( 0, index ));
	    end2 = Integer.parseInt( block2.substring( index + 1 ));
	    lcr2 = seq.substring(start2 - 1, end2);
	    lcr2 = findAlignment( lcr1, lcr2 );
	    if ( lcr2.length() != 0 ) {	
	      LCRs[indexLCR] = tmpLCRs[i];
	      ++ indexLCR;
	    } 
	  }
	}
	else {
	  if ( i != (tmpIndexLCR-1)) {
	    block2 = tmpLCRs[i + 1];
	    index = block2.indexOf( "-" );
	    start2 = Integer.parseInt( block2.substring( 0, index ));
	    end2 = Integer.parseInt( block2.substring( index + 1 ));
	    lcr2 = seq.substring(start2 - 1, end2);
	    lcr2 = findAlignment( lcr1, lcr2 );
	    if ( lcr2.length() != 0 ) {	
	      LCRs[indexLCR] = tmpLCRs[i];
	      ++ indexLCR;
	    } 
	  } 
	}
      }
      ++ i;
    }
    ///////////////System.out.println( "refined LCRs:");
    printLCRsOReps(1,true);
  }


  public void refine() {
    String lcrLine = new String();
    String seq = new String();
    try {
      lcrLine = lcrRf.readLine();
      seq = seqRf.readLine();
      while (lcrLine != null) {
	if (lcrLine.startsWith(">")) {
	  System.out.println(lcrLine);
	  lcrLine = lcrRf.readLine();
	  getLCRsOReps(lcrLine, 1); // 1 for LCRs
	  seq = getSeq();
	  //System.out.println(seq);	  
	  refineLCRs(seq);
	  tmpIndexLCR = 0;
	  indexLCR = 0;
	  lcrLine = lcrRf.readLine();
	  System.out.println();  
	}
	else 
	  lcrLine = lcrRf.readLine();
      }	
      lcrRf.close();
      seqRf.close();
    }
    catch (IOException ex) {
    }
  }
  

  public static void main ( String args[] ) {
     Refiner r = new Refiner(args[0], args[1] );
     r.refine();
  }
}
