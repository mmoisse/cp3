/*
** Author: Xuehui Li
** Date: Dec, 2005
** This program is used to generate the data for a cdf (cumulative distribution function) plot.
** There are two input parameters.
** The first one is the file that contain all complexities. The second one is the number of bins all complexities will be put into.
** An example to run the program:
** java applications.CumulativeDFDataGenerator complexityMeasuresComparison/ShannonRepeats.txt 20 > complexityMeasuresComparison/ShannonRepeatsCDFData
*/

package applications;

import java.io.*;

public class CumulativeDFDataGenerator {

  private static final int TOTALNUM = 128;

  private double min, max;
  private double[] complexities;

  public CumulativeDFDataGenerator () {
    min = 0;
    max = 0;
    complexities = new double[TOTALNUM];
  }

  public void getComplexities(String comFile) {
    int index = 0;
    try {
      File f = new File(comFile);
      RandomAccessFile rf = new RandomAccessFile(f, "r");
      String line = new String();
      line = rf.readLine();
      while (line != null) {
	line = line.trim();
	complexities[index] = Double.parseDouble(line);
	++ index;
	line = rf.readLine();
      }
      rf.close();
    }
    catch (IOException ex) {}
    min = complexities[0];
    max = complexities[index-1];  
  }

  public void generate (String binNum) {
    double interval = (max-min)/Double.parseDouble(binNum);
    double[] percentages = new double[Integer.parseInt(binNum)+1];
    int index = 0;
    for (int i=0; i<TOTALNUM; i++) {
      if (complexities[i]>= (min+(index+1)*interval)) {
	//System.out.println((min+(index+1)*interval) +" "+ + percentages[index]);
	++ index;
	if (complexities[i] < (min+(index+1)*interval))
	  percentages[index] = percentages[index] +1;
	else {
	  //System.out.println((min+(index+1)*interval) +" "+ + percentages[index]);
	  ++ index;
	  percentages[index] = percentages[index] +1;
	}
      }
      else 
	percentages[index] = percentages[index] +1;
    }
    double sum = 0;
    for (int i=0; i<index; i++) {
      percentages[i] = percentages[i] / TOTALNUM;
      sum = sum + percentages[i];
      System.out.println((min+(i+1)*interval) +" "+ sum);
      //      System.out.println(percentages[i]);
    }
  }

  public static void main (String args[] ) {
    CumulativeDFDataGenerator cddfg = new CumulativeDFDataGenerator();
    cddfg.getComplexities(args[0]);
    cddfg.generate(args[1]);
  }

}
