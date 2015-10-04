/* 
** Author: Xuehui Li
** Date: Jan,  2006
** This program is used to normalize complexities in a file by the maximum complexity in this file. 
** There is one input parameter. It is the name of the file that lists the complexities and the percentage. Each line has the format of "complexity percentage". 
** The output file has the same format as the input file. However, the complexity is the normalized one of the corresponding complexity in the input file.
** An example to run the program:
**  java applications/Normalizer complexityMeasuresComparison/combinedSet/ShannonRepeats.txt > complexityMeasuresComparison/combinedSet/analysis/ShannonRepeatsNormalized.txt
*/



package applications;

import java.io.*;

public class Normalizer {

  public Normalizer () {
  }	

  public void normalize(String file) {
    try {
      File f = new File(file);
      RandomAccessFile rf = new RandomAccessFile(f,"r");
      String line = new String(), preLine = line;
      line = rf.readLine();
      while (line!=null) {
	preLine = line;	
	line = rf.readLine();
      }
      int index = preLine.indexOf(" ");
      preLine = preLine.substring(0,index);
      double max = Double.parseDouble(preLine);
      rf.seek(0);
      // System.out.println("max " + max);
      line = rf.readLine();
      while (line!=null) {
	index = line.indexOf(" ");
	preLine = line.substring(0,index);
	double com = Double.parseDouble(preLine) / max;
	//System.out.println();
	//System.out.println(line);
	System.out.println(com + line.substring(index));
	line = rf.readLine();
      }
      rf.close();
    }	
    catch (IOException ex) {}
  }	

  public static void main (String args[]) {
    Normalizer nn = new Normalizer();
    nn.normalize(args[0]);
  }

}
