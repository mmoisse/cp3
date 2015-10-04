/* 
** Author: Xuehui Li
** Date: Jan,  2006
** This program is used to compute the Shannon Entropy and our new complexity me
asure of all repeats and nonrepeats in a sequence file.
** There are two input parameters. The first one is a file that tells the repeat
s of a sequence file. The second one is the sequence file.
** The output is a file that lists all repeats and nonrepeats from all sequences
, and their complexities in a non-decreasig order.
** An example to run the program:
** java applications.NorNRepCDFProcessor 200 complexityMeasuresComparison/combinedSet/analysis/ShannonRepeatsNormalized.txt complexityMeasuresComparison/combinedSet/analysis/ShannonNonRepeatsNormalized.txt > mplexityMeasuresComparison/combinedSet/analysis/ShannonRNR.txt
*/


package applications;

import java.io.*;
import java.util.Vector;

public class NorNRepCDFProcessor {

  double binNum;
  private Vector repeats;
  private Vector nonRepeats;
  
  
  public NorNRepCDFProcessor (String str) {
    binNum = Double.parseDouble(str);
    repeats = new Vector();
    nonRepeats = new Vector();
  }
  
  public Vector getRepeatPercentage(String file) {
    double com1 = 0;
    double com2 = 0;
    Vector com1s = new Vector();
    try{
      File f = new File(file);
      RandomAccessFile rf = new RandomAccessFile(f, "r");
      String line = rf.readLine();
      line = line.trim();
      int index = line.indexOf(" ");
      com2 = Double.parseDouble(line.substring(0,index));
      com1 = com2;
      double intervalLen = (1-com2)/binNum; 
      String percentage = line.substring(index+2);
      String prePercentage = percentage;
      while (com1 < 1) {
	  while ((line!=null) && (com2 <= com1)) {
	      System.out.println("shoot");
	      prePercentage = percentage;
	      line = rf.readLine();
	      if (line != null) {
		  //System.out.println(line);   
		  line = line.trim();
		  index = line.indexOf(" ");
		  com2 = Double.parseDouble(line.substring(0,index));
		  percentage = line.substring(index+2);	
	      }
	      // else 
	      //com1 = 1;
	  }
	  repeats.add(prePercentage);
	  com1s.add(Double.toString(com1));
	  while (com2 > com1) {
	      com1 = com1 + intervalLen;
	  }
      }
    }
    catch (IOException ex) {}
    System.out.println("???");   
 return com1s;
  }
  
  public void getNonRepeatPercentage(Vector v, String file) {
    double com1 = 0;
    double com2 = 0;
    int indexV = 1;
    try{
      File f = new File(file);
      RandomAccessFile rf = new RandomAccessFile(f, "r");
      String line = rf.readLine();
      line = line.trim();
      int index = line.indexOf(" ");
      com2 = Double.parseDouble(line.substring(0,index));
      com1 = Double.parseDouble((String)v.elementAt(0));
      String percentage = line.substring(index+1);
      String prePercentage = percentage;
      while ((nonRepeats.size()< repeats.size()) && (com1 <= 1)) {
	  while (com2 <= com1) {
	      //  System.out.println("************" + com2 + "      "+ com1);
	    prePercentage = percentage;
	  line = rf.readLine();
	  if (line!=null) {
	      line = line.trim();
	      index = line.indexOf(" ");
	      com2 = Double.parseDouble(line.substring(0,index));
	      percentage = line.substring(index+1);	
	  }
	  else 
	    com1 = 1;
	}
	nonRepeats.add(prePercentage);
	++ indexV;
	if (com1 != 1) {
	    //System.out.println(com1);
	    if (indexV < repeats.size())
		com1 = Double.parseDouble((String)v.elementAt(indexV));
	}
	else System.out.println("1111111");
      }
    }
    catch (Exception ex) {System.out.println(ex.getMessage());}
  }

  public void print() {
    for (int i=0; i<nonRepeats.size(); i++) 
      //System.out.println((String)nonRepeats.elementAt(i));
      System.out.println((String)repeats.elementAt(i) + " " + (String)nonRepeats.elementAt(i));
  }
  
  public static void main (String args[] ) {
    NorNRepCDFProcessor np = new NorNRepCDFProcessor(args[0]);
    Vector v = np.getRepeatPercentage(args[1]);
    System.out.println("OOOOOOOOOOOOOOOO");
    np.getNonRepeatPercentage(v, args[2]);
    np.print();
  }
}
