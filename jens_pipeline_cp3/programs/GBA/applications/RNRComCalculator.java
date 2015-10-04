
/* 
** Author: Xuehui Li
** Date: Dec,  2005
** This program is used to compute the Shannon Entropy and our new complexity measure of all repeats and nonrepeats in a sequence file.
** There are two input parameters. The first one is a file that tells the repeats of a sequence file. The second one is the sequence file.
** The output is a file that lists all repeats and nonrepeats from all sequences, and their complexities in a non-decreasig order.
** An example to run the program:
**  java applications.RNRComCalculator ../data/swissprot/repeatsInfor/mgdRepeat ../data/swissprot/sequenceInfor/seqFromMgd > complexityMeasuresComparison/complexityNorMeasuresComparison
*/

package applications;

import java.io.*;

class RNRComCalculator {

  private static final int NUMOFSEQUENCES = 468;  // For mgd, this constant should be set as 128.

  private RandomAccessFile repeatRf;
  private RandomAccessFile seqRf;

  public RNRComCalculator (String repeatFile, String seqFile) {
    try {
      File f = new File(repeatFile);
      repeatRf = new RandomAccessFile(f, "r");
      f = new File(seqFile);
      seqRf = new RandomAccessFile(f, "r");
    }
    catch (IOException ex) {}
  }

  public void calculate() {
    double[] RSCom = new double[NUMOFSEQUENCES];
    double[] NRSCom = new double[NUMOFSEQUENCES];
    double[] RCom = new double[NUMOFSEQUENCES];
    double[] NRCom = new double[NUMOFSEQUENCES];
    Double[] RSComD = new Double[NUMOFSEQUENCES];
    Double[] NRSComD = new Double[NUMOFSEQUENCES];
    Double[] RComD = new Double[NUMOFSEQUENCES];
    Double[] NRComD = new Double[NUMOFSEQUENCES];
    int comIndex = 0;
    int index1 = 0, index2 = 0, start = 0, end = 0,  preEnd = 0;
    try {
      String repeatLine = repeatRf.readLine();
      String seqLine = seqRf.readLine();
      String seq = new String();
      String repeats = new String();
      String nonRepeats = new String();
      double com = 0;
      while (repeatLine != null) {
	if(repeatLine.startsWith(">")&&(repeatLine.indexOf("REPT_MOUSE")==-1)&&(repeatLine.indexOf("SON_MOUSE")==-1) && (repeatLine.indexOf("DERM_MOUSE")==-1) && (repeatLine.indexOf("HNRH1_MOUSE")==-1) && (repeatLine.indexOf("HNRH2_MOUSE")==-1)) {
	  repeatLine = repeatLine.trim();
	  while (!(seqLine.startsWith(repeatLine)))
	    seqLine = seqRf.readLine();
	  repeatLine = repeatRf.readLine();
	  repeatLine = repeatLine.trim();
	  seqLine = seqRf.readLine();
	  while ((seqLine!=null)&&(!(seqLine.startsWith(">")))) {
	    seqLine = seqLine.trim();
	    seq = seq + seqLine;
	    seqLine = seqRf.readLine();
	  }
	  //System.out.println(seqLine);
	  index1 = 0; index2 = 0; start = 0; end = 0;  preEnd = 0;
	  while (repeatLine!= " ") {
	    index1 = repeatLine.indexOf("-");
	    start = Integer.parseInt(repeatLine.substring(0,index1));
	    index2 = repeatLine.indexOf(" ");
	    if (index2 == -1) {
	      end = Integer.parseInt(repeatLine.substring(index1+1));
	      repeatLine = " ";
	    }
	    else {
	      end = Integer.parseInt(repeatLine.substring(index1+1, index2));
	      repeatLine = repeatLine.substring(index2+1);
	    }
	    repeats = repeats + seq.substring(start-1, end);
	    // System.out.println(start + " " + end +"       " + preEnd);
	    if (start !=0)
	      if (start != preEnd)
		nonRepeats = nonRepeats + seq.substring(preEnd, start-1);
	    preEnd = end;
	  }
	  nonRepeats = nonRepeats + seq.substring(preEnd);
	  //System.out.println(repeats);
	  //System.out.println(nonRepeats);
	  ComplexityCalculator cc = new ComplexityCalculator();
	  cc.initializeAlphabet();
	  com = cc.calculateEntropy(repeats);
	  RSComD[comIndex] = new Double(com);
	  RSCom[comIndex] = com;
	  com = cc.calculateEntropy(nonRepeats);
	  NRSComD[comIndex] = new Double(com);
	  NRSCom[comIndex] = com;
	  //System.out.println( RSCom[comIndex] + " "+ NRSCom[comIndex]);
	  /*
	  cc.initializeNewAlphabet();
	  cc.computeNewScoringMatrix();
	  cc.normalizeNewScoringMatrixPow(); 	
	  com = cc.calculateNor2LetterEntropyWScoMatrix(repeats);
	  */
	  com = cc.calculateModifiedEntropy(repeats);
	  //com = cc.calculate2LetterEntropyWScoMatrix(repeats);
	  RComD[comIndex] = new Double(com);
	  RCom[comIndex] = com;
	  //com = cc.calculateNor2LetterEntropyWScoMatrix(nonRepeats);
	  com = cc.calculateModifiedEntropy(nonRepeats);
	  //com = cc.calculate2LetterEntropyWScoMatrix(nonRepeats);
	  NRComD[comIndex] = new Double(com);
	  NRCom[comIndex] = com;
	  //System.out.println( RCom[comIndex] + " "+ NRCom[comIndex]);
	  ++comIndex;
	  repeats = new String();
	  nonRepeats = new String();
	  seq = new String();
	}
	else 
	  repeatLine = repeatRf.readLine(); 
      }
      repeatRf.close();
      seqRf.close();
    }
    catch (Exception ex) {
      System.out.println(preEnd + " " + start + " " + end);
      System.out.println(ex.getMessage());
      System.exit(0);
    }
    MergeSortDouble.mergeSort(RSComD);
    MergeSortDouble.mergeSort(NRSComD);
    MergeSortDouble.mergeSort(RComD);
    MergeSortDouble.mergeSort(NRComD);
    System.out.println("      SR                  SNR                  NewR                NewNR");
    double total = NUMOFSEQUENCES;
    for (int i =0; i<comIndex; i++ ) {
      System.out.println((i+1) + " " + RSComD[i].toString() + " "+ NRSComD[i].toString() + " "+ RComD[i].toString() + " "+ NRComD[i].toString());
      //System.out.println(RSComD[i].toString() + " " + NRSComD[i].toString() + " "+ RComD[i].toString() + " " + NRComD[i].toString());
    }
  }
  
  public static void main (String args[]) {
    RNRComCalculator rcc = new RNRComCalculator(args[0], args[1]); 
    rcc.calculate();
  }
}
