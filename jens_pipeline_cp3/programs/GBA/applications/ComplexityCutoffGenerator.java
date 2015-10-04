/**
 * Author: Xuehui Li
 * Date: August, 2005
 *This program is used to calculate the complexity cutoff value in the longest path intervals extension of GBA.
 * There are four three parameters: three file names. 
* The first one specifies which regions are masked as LCRs/HCRs by the algorithm. 
* The second one is the sequence file.
* The third one is the output file which contains the generated cutoff value.
* An example to run the program when the normalized complexity measure is used is:
*  java applications.ComplexityCutoffGenerator swissprotLCRBlocks/wForgetRateMatrices/postProcess/cutting/meanStd/combinedMeanStd ../data/swissprot/sequenceInfor/combinedSeq > knowledge/comCutOfffValue
**/

package applications;

import java.io.*;
import java.util.*;

class ComplexityCutoffGenerator {

  private static final double C = 0.05;
  private static final int SEQNUM = 474;

  private RandomAccessFile lcrRf, seqRf;
  private Vector LCRs = new Vector();

  public ComplexityCutoffGenerator (String lcrFile, String seqFile) {
   try {
      File f = new File (lcrFile);
      lcrRf = new RandomAccessFile (f, "r");
      f = new File (seqFile);
      seqRf = new RandomAccessFile (f, "r");
   }
    catch (IOException ex) {
    }
  }


  public void printLCRs() {
    int len = LCRs.size();
    for ( int i = 0; i< len; i++ ) 
      System.out.print( (String) LCRs.elementAt(i) + " " );
    System.out.println();
  }


  public void getLCRs (String line ) {
    int i;
    String str = line;
    str = str.trim();
    if ( str.length() == 0 ) {
      str = null;
      System.out.println( "yeah" );
    }
    while ( str != null ) {
      i = str.indexOf(" ");
      if ( i == -1 ) {
	LCRs.add(str.substring( 0 ));
	str = null;
      }     
      else {
	LCRs.add(str.substring( 0, i ));
	str = str.substring( i + 1 );     
      }
    }
    printLCRs();
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


  public void getLcrSeqs( String seq) {
    System.out.println(seq);
    String str = new String();
    int index = 0, start, end;
    //int len = LCRs.size();
    // for ( int i = 0; i < len; i++) {
    int i = 0;
    while ( i < LCRs.size()) {
      str = (String) LCRs.elementAt(i);
      index = str.indexOf("-");
      start = Integer.parseInt(str.substring(0, index));
      end = Integer.parseInt(str.substring(index+1));
      //System.out.println(start +" "+end);
      if (start > end ) {
	LCRs.remove(i);
      }
      else {
	str = seq.substring( start-1, end);
	LCRs.remove(i);
	LCRs.add(i, str);
	++i;
      }
    }
    printLCRs();
    System.out.println("********************");
  }  


  public String getMeanStd(Vector v) {
    String meanStd = new String(); 
    double ele, mean = 0, std = 0;
    int num = v.size();
    for ( int i = 0; i < num; i++ ) { 
      ele = ((Double)v.elementAt(i)).doubleValue();
      System.out.print(ele+" ");
      mean = mean + ele;
      std = std + Math.pow(ele,2d);
    }
    mean = mean / (num);
    std = std / (num)   ;
    std = std - Math.pow(mean,2d);
    std= Math.sqrt(std);
    meanStd = mean + " " + std;
    System.out.println();
    System.out.println(meanStd);
    return meanStd;
  }

  public void generate () {
    ComplexityCalculator cc = new ComplexityCalculator();
    cc.initializeAlphabet();
    cc.initializeNewAlphabet();
    cc.computeNewScoringMatrix();
    cc.normalizeNewScoringMatrixPow();
    String lcrLine = new String();
    String seq = new String(), seqLine = new String();
    double com = 0, mean = 0, stanDev = 0;
    int i = 0, r = 0, preR = 0, k1 =0, k2 = 0;
    boolean less = false, first = true;;
    Random generator = new Random(); 
    Vector randoms = new Vector(), comV = new Vector();
    Integer tmpIg;
    try {
      while (( !less )&&(randoms.size() < SEQNUM)) {
	++i;
	r = generator.nextInt(SEQNUM);
	if ( r != 0 ) {
	  tmpIg = new Integer( r );
	  if ( randoms.indexOf( tmpIg ) == -1 ) { // sampling without replacement
	    randoms.add( tmpIg );
	    lcrRf.seek(0);
	    seqRf.seek(0);
	    k1 = 0;
	    k2 = 0;
	    while (k1 <r) {
	      lcrLine = lcrRf.readLine();
	      if (lcrLine.startsWith(">")) 
		++k1;
	    }
	    System.out.println(lcrLine);  
	    lcrLine = lcrRf.readLine();
	    getLCRs(lcrLine); 
	    while (k2<r) {
	      seqLine = seqRf.readLine();
	      //System.out.pritnln(seqLine);
	      if (seqLine.startsWith(">"))
		++k2;
	    }
	    seq = getSeq();
	    getLcrSeqs(seq);   
	    // com = cc.calculateSeqNormComplexity( LCRs); // used the normalized complexity measure
	    com = cc.calculateSeqUnNormComplexity( LCRs);// used the un-normalized complexity measure
	    comV.add(new Double(com));	   
	    String meanStd = getMeanStd(comV); ;
	    int index = meanStd.indexOf(" ");
	    mean = Double.parseDouble(meanStd.substring(0,index));
	    stanDev = Double.parseDouble(meanStd.substring(index+1));
	    if ((stanDev/Math.sqrt(i)) < C*mean) {
	      if (!first) {
		System.out.println("cutoff values is : " + mean);
		less = true;
	      }
	      else {
		first = false;
		System.out.println("first");
	      }
	    }
	    LCRs.clear();
	  }
	}
      }
      lcrRf.close();
      seqRf.close();
    }
    catch (IOException ex) {
    }  
  }
  

  public static void main (String args[]) {
    ComplexityCutoffGenerator cvpsc = new ComplexityCutoffGenerator(args[0], args[1]);
    cvpsc.generate();
  }
  
}
