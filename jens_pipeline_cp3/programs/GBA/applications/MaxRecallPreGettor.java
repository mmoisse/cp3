/**
 * Author: Xuehui Li
 * Date: July, 2005
 *This program is used to calculate the maximum recall and precision by combining the masked LCRs identified two LCR-identification algorithms.
* The maximum recall is TPmax/T. The maximum precison is TPmax/TPmax+FPmin; 
* We get TPmax by taking the intersection of all repeats and the union of the LCRs from the two algorithms.
*  We get FPmin by taking the intersection of all repeats and the intersection of the LCRs from the two algorithms.
* There are four input parameters: three file names. 
* The first one specifies which regions are masked as LCRs by the first algorithm.
* The second one specifies which regions are masked as LCRs by the second algorithm. 
* The third one is the metrics file from one of the algorithms and specifies the repeat regions and repeat length contained in each sequence.
* Suppose the two algorithms are GBA and SEG. An example to run the program is:
* java applications.MaxRecallPreGettor swissprotLCRBlocks/wForgetRateMatrices/postProcess/cutting/meanStd/flyabseMeanStd ../SEG/pseg/output/swissprot/LCRBlocks/flybaseLCRBlocks metrics/swissprot/wForgetRate/postProcessing/cutting/meanStd/flyabseMeanStdMetrics > /cise/research/tamer/xli/LCR/combinedAlgo/metrics/maxGS/flybaseMax
**/

package applications;

import java.io.*;
import java.util.Vector;

class MaxRecallPreGettor {

  private RandomAccessFile LCRBlocksRf1, LCRBlocksRf2, repeatRf;
  private String[] LCRs1 = new String[5000], LCRs2 = new String[5000];   
  private Vector repeatRegions = new Vector(), repeatLength = new Vector();
  private int indexLCR1 = 0, indexLCR2 = 0;
  
  public  MaxRecallPreGettor (String lcrBlockFile1, String lcrBlockFile2, String metricsFile) {
    try {
      File f = new File (lcrBlockFile1);
      LCRBlocksRf1 = new RandomAccessFile(f,"r");
      f = new File (lcrBlockFile2);
      LCRBlocksRf2 = new RandomAccessFile(f,"r");
      f = new File(metricsFile);
      repeatRf = new RandomAccessFile( f, "r");
    }
    catch (IOException ex) {
    }

  }


  public void printLCRs( int mark ) {
    int limit = 0;
    if ( mark == 1 )
      limit = indexLCR1;
    else 
      limit = indexLCR2;
    for ( int i = 0; i < limit; i++ ) {
      if ( mark == 1 ) 
	System.out.print( LCRs1[i]+" " );
      else 
	System.out.print( LCRs2[i]+" " );
    }
    System.out.println();
    if ( mark == 2 )
      System.out.println( "**************************" ); 
  }


  public void getLCRs ( String line, int mark ) {
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
          LCRs1[ indexLCR1 ] = str.substring( 0 );
        else 
          LCRs2[ indexLCR2 ] = str.substring( 0 );
        str = null;
      }     
      else {
        if ( mark == 1 )
          LCRs1[ indexLCR1 ] = str.substring( 0, i );
        else 
          LCRs2[ indexLCR2 ] = str.substring( 0, i );
        str = str.substring( i + 1 );     
      }
      if ( mark ==1 )
        ++ indexLCR1;
      else 
        ++ indexLCR2;
    }
    printLCRs( mark );
  }


  public void printVector ( Vector v) {
    String str = new String();
    for ( int j = 0; j < v.size(); j++ ) {
      str =  (String) v.elementAt( j ) + " ";
      System.out.print( str );
    }
   System.out.println();  
  }
  
  
  public Vector clean ( Vector tmpUnion ) {
    Vector union = tmpUnion;
    String str = new String();
    System.out.println( "before cleaning" );     
    printVector(union);
    int i = 0, len = union.size();
    int index, start1, end1, start2, end2;
    String lcr1 = new String(), lcr2 = new String();
    while ( i < (len-1) ) {
      lcr1 = (String) union.elementAt( i );
      lcr2 = (String) union.elementAt( i+1 );
      index = lcr1.indexOf( "-" );
      start1 = Integer.parseInt( lcr1.substring( 0, index ));
      end1 = Integer.parseInt( lcr1.substring( index+1 ));
      index = lcr2.indexOf( "-" );
      start2 = Integer.parseInt( lcr2.substring( 0, index ));
      end2 = Integer.parseInt( lcr2.substring( index+1 ));
      if ( start1 >= start2 ) {
	System.out.print( "begin1 " );
	System.out.print( lcr1 + " " + lcr2 );
	System.out.println( " end1" );
	if ( end2 < end1 )
	  end2 = end1;
	union.remove(i);
	union.remove(i);
	lcr2 = start2 + "-" + end2;
	union.add(i, lcr2);
      }   
      else {
	if ( start2 <= end1 ) {
	  System.out.print( "begin2 " );
	  System.out.print( lcr1 + " " + lcr2 );
	  System.out.println( " end2" );
	  if ( end1 < end2 )
	    end1 = end2;
	  union.remove(i);
	  union.remove(i);
	  lcr1 = start1 + "-" + end1;
	  union.add(i, lcr1);
	}
	else ++i;
      }
      len = union.size();
    }
    System.out.println( "after cleaning");
    printVector(union);
    return union;
  }


  public int getTotalLength( Vector union ) {
    int index, index2, start, end, length = 0;
    String str = new String();
    for ( int i = 0; i < union.size(); i++ ) {
      str =  (String) union.elementAt( i );
      index = str.indexOf( "-" );
      start = Integer.parseInt( str.substring( 0, index ));
      end = Integer.parseInt( str.substring( index + 1 ));
      length = length + end - start + 1;
      System.out.println( str + " " + start + " " + end);
    }
    System.out.println( "length " + length);
    return length;
  }
  

  public Vector unionLCRs() { // to get max TP
    int maxTp = 0, i = 0, j = 0, index, bound, start1=0, start2=0, end1=0, end2=0;
    Vector union = new Vector();
    while (( i < indexLCR1 ) && ( j < indexLCR2)) {
      //System.out.println( LCRs1[i] + "**" + LCRs2[j] );
      index = LCRs1[i].indexOf( "-" );
      start1 = Integer.parseInt( LCRs1[i].substring( 0, index ));
      end1 = Integer.parseInt( LCRs1[i].substring( index+1 ));
      index = LCRs2[j].indexOf( "-" );
      start2 = Integer.parseInt( LCRs2[j].substring( 0, index ));
      end2 = Integer.parseInt( LCRs2[j].substring( index+1 ));
    
      if ( end1 <= start2 ) { // no overlap
	//System.out.println( "no overlap1");
	union.add( LCRs1[i] );
	++i;     
      }
      else if ( end2 <= start1 ) { // no overlap
	//System.out.println( "no overlap2");
	union.add( LCRs2[j] );
	++j;
      }
      else if (( start1 <= start2 ) &&( end1 >= end2 )) { // full containment 
	//System.out.println( "full containment1");  
	union.add( start1 + "-" + end1 );
	bound = 0;
	if ( j == ( indexLCR2 - 1 ))
	  ++ j;
	while (( j < ( indexLCR2 - 1 ) ) && ( bound < end1 )) {
	  ++ j;
	  index = LCRs2[j].indexOf( "-" );
	  bound = Integer.parseInt( LCRs2[j].substring( index + 1 ));
	}
      }
      else if (( start2 <= start1 ) &&( end2 >= end1 )) { // full containment 
	//System.out.println( "full containment2");  
	union.add( start2 + "-" + end2 );
	bound = 0;
	if ( i == ( indexLCR1 -1 ) ) 
	  ++ i;
	while (( i < ( indexLCR1-1) ) && ( bound < end2 )) {
	  ++ i;
	  index = LCRs1[i].indexOf( "-" );
	  bound = Integer.parseInt( LCRs1[i].substring( index + 1 ));
	}
      } 
      else if (( start1 < start2 ) && ( end1 < end2 ) && ( end1 > start2)) {  
	// intersect, i.e., overlap, but no containment
	//System.out.println( "overlapping1" );
	union.add( start1 + "-" + end2 );
	bound = 0;
	if ( i == ( indexLCR1 - 1 ))	     
	  ++i;
	while (( i < ( indexLCR1 - 1 )) && ( bound < end2 )) {
	  ++i;	      
	  index = LCRs1[i].indexOf( "-" );
	  bound = Integer.parseInt( LCRs1[i].substring( index + 1 ) );   
	}
      }      
      else if (( start2 < start1 ) && ( end2 < end1 ) && ( end2 > start1)) {  
	// intersect, i.e., overlap, but no containment
	//System.out.println( "overlapping2" );
	union.add( start2 + "-" + end1 );
	bound = 0;
	if ( j == (indexLCR2-1))
	  ++ j;
	while (( j < (indexLCR2-1)) && ( bound < end1 )) {
	  ++j;			
	  index = LCRs2[j].indexOf( "-" );
	  bound = Integer.parseInt( LCRs2[j].substring( index + 1 ) );   
	}
      }      
      else 
	System.out.println( "wrong" );
    }

    if ( end2 <= start1 ) {  
      //System.out.println("oo1");
    }
    else
      ++ i;
    while ( i < indexLCR1 ) {
      union.add( LCRs1[i] );   
      ++ i;
    }
    if ( end1 <= start2 ) { 
      //  System.out.println( "oo2");   
    } 
    else  
      ++ j; 
    while ( j < indexLCR2 ) {
      union.add( LCRs2[j] );
      ++ j;
    }
    System.out.println( "union of LCRs1 and LCRs2");
    union = clean( union ); 
    return union;  
  }


  public Vector intersectLCRs() { 
    Vector intersection = new Vector();
    int minFp = 0;
    //    double minFp = 0; 
    int i = 0, j = 0, k = 0;
    int lcrStartInt1 = 0, lcrEndInt1 = 0, lcrStartInt2 = 0, lcrEndInt2 = 0;
    String lcrStaratStr1 = new String(), lcrEndStr1 = new String(); 
    String lcrStartStr2 = new String(), lcrEndStr2 = new String();
    boolean found = false;
    while ( i < indexLCR1 ) {
      k = LCRs1[i].indexOf( "-" );
      lcrStaratStr1 = LCRs1[i].substring( 0, k );
      lcrStartInt1 = Integer.valueOf( lcrStaratStr1 ).intValue();
      lcrEndStr1 = LCRs1[i].substring( k + 1 );
      lcrEndInt1 = Integer.valueOf( lcrEndStr1 ).intValue();
      while ( ( !found ) && ( j < indexLCR2 )){
	k = LCRs2[j].indexOf( "-" );
	lcrStartStr2 = LCRs2[j].substring( 0, k );
	lcrStartInt2 = Integer.valueOf( lcrStartStr2 ).intValue();
	lcrEndStr2 = LCRs2[j].substring( k + 1 );
	lcrEndInt2 = Integer.valueOf( lcrEndStr2 ).intValue();
	if (( lcrEndInt2 >= lcrStartInt1 ) && (lcrStartInt2 <= lcrStartInt1) && (lcrEndInt2 <= lcrEndInt1)) {
	  intersection.add(  lcrStartInt1 + "-" + lcrEndInt2 );
	  ++j;
	} 
	else if ((  lcrEndInt2 <= lcrEndInt1) && (  lcrStartInt2 >= lcrStartInt1 )) {
	  intersection.add( lcrStartInt2 + "-" + lcrEndInt2  );
	  ++j;
	}
	else if (( lcrEndInt2 > lcrEndInt1)&&(lcrStartInt2 <= lcrEndInt1)&&( lcrStartInt2 >= lcrStartInt1)) {
	  intersection.add( lcrStartInt2 + "-" +  lcrEndInt1);
	  found = true;	}
	else if ( lcrStartInt2 > lcrEndInt1 ) {	
	  found = true;
	}
	else if (( lcrStartInt2 <= lcrStartInt1) && ( lcrEndInt2 >= lcrEndInt1 )) {
	  intersection.add( lcrStartInt1 + "-" + lcrEndInt1);
	  found = true;
	}
	else {
	  ++j;
	}
      }
      ++i;
      found = false;
    }  
    System.out.println( "intersections between LCRs1 and LCRs2");
    printVector(intersection);
    return intersection;
  }

  
  public void getRepeats() {
    String line = new String();
    int i, j;
    try {
      line = repeatRf.readLine();
      while ( line != null ) {
	if ( line.startsWith("Repeat Infor.:")) {
	  line = repeatRf.readLine();
	  repeatRegions.add(line.trim());
	}
	if ( line.startsWith("Total length of all repeats:")) {
	  line = line.trim();
	  i = line.indexOf(":");
	  repeatLength.add( line.substring(i+1));
	}
	line = repeatRf.readLine();
      }
      repeatRf.close();
    }
    catch ( IOException ex ) {
    }
  }


  public Vector getRepeatVector(String line) {
    Vector v = new Vector();
    String str = line;
    int i, j;
    if ( str.length() == 0 ) {
      str = null;
      System.out.println( "no repeat" );
    }  
    while ( str != null ) {
      i = str.indexOf(" ");
      j = str.indexOf("-");
      if ( i == -1 ) {
	v.add( str.substring( 0 ));
	str = null;
      }     
      else {
	v.add(str.substring( 0, i ));
	str = str.substring( i + 1 );     
      }
    }
    System.out.println( "repeats");
    printVector(v);
    return v;
  }


  public Vector intersectVectors( Vector v1, Vector v2 ) {
    Vector v = new Vector();
    int minFp = 0;
    int len1 = v1.size(), len2 = v2.size();
    int i = 0, j = 0, k = 0;
    int startInt1 = 0, endInt1 = 0, startInt2 = 0, endInt2 = 0;
    String str1, str2;
    String startStr1 = new String(), endStr1 = new String(); 
    String startStr2 = new String(), endStr2 = new String();
    boolean found = false;
    while ( i < len1 ) {
      str1 = (String) v1.elementAt(i);
      k = str1.indexOf( "-" );
      startStr1 = str1.substring( 0, k );
      startInt1 = Integer.valueOf( startStr1 ).intValue();
      endStr1 = str1.substring( k + 1 );
      endInt1 = Integer.valueOf( endStr1 ).intValue();
      while ( ( !found ) && ( j < len2)){
	str2 = (String) v2.elementAt(j);
	k = str2.indexOf( "-" );
	startStr2 = str2.substring( 0, k );
	startInt2 = Integer.valueOf( startStr2 ).intValue();
	endStr2 = str2.substring( k + 1 );
	endInt2 = Integer.valueOf( endStr2 ).intValue();
	if (( endInt2 >= startInt1 ) && (startInt2 <= startInt1) && (endInt2 <= endInt1)) {
	  v.add( startInt1 + "-" + endInt2 );
	  ++j;
	} 
	else if ((  endInt2 <= endInt1) && (  startInt2 >= startInt1 )) {
	  v.add(  startInt2 + "-" +  endInt2 );
	  ++j;
	}
	else if (( endInt2 > endInt1)&&(startInt2 <= endInt1)&&( startInt2 >= startInt1)) {
	  v.add( startInt2+ "-" + endInt1 );
	  found = true;	
	}
	else if ( startInt2 > endInt1 ) {	
	  found = true;
	}
	else if (( startInt2 <= startInt1) && ( endInt2 >= endInt1 )) {
	  v.add( startInt1 + "-" + endInt1 );
	  found = true;
	}
	else {
	  ++j;
	}
      }
      ++i;
      found = false;
    }  
    return v;
  }

  
  public void getMaxRecallPre() {
    String line1 = new String(), line2 = new String();
    Vector union = new Vector(), intersection = new Vector(), repeatVector = new Vector();
    Vector recPre = new Vector();
    int i = 0;
    double tp, fp, repeat;
    double recall = 0, precision = 0, recSum = 0, preSum = 0;
    getRepeats();
    //System.out.println( "REPeatRegion:" + repeatRegions.size());
    try {
      line1 = LCRBlocksRf1.readLine();
      line2 = LCRBlocksRf2.readLine();
      while ( line1 != null ) {
        if (( line1.startsWith( ">")) || ( line1.startsWith("ID"))) {
          System.out.println( line1 );
          line1 = LCRBlocksRf1.readLine();
          getLCRs( line1, 1 );
          line2 = LCRBlocksRf2.readLine();
          getLCRs( line2, 2 );
	  union =  unionLCRs(); // get the union of 'LCRs1' and 'LCRs2'
	  repeatVector = getRepeatVector((String)repeatRegions.elementAt(i));
	  System.out.println( "intersection to get max TP"  + repeatVector.size());
	  union = intersectVectors(union, repeatVector); // get the inersection of the 'union' and 'repeats'
	  tp = getTotalLength(union);
	  intersection = intersectLCRs();
	  System.out.println( "intersection to get min FP:");
	  intersection = intersectVectors(intersection, repeatVector);
	  repeat = Double.parseDouble( (String)repeatLength.elementAt(i));
	  fp = repeat - getTotalLength(intersection);
	  recall = tp / repeat;
	  if ( (tp+fp)==0 ) 
	    System.out.println("nonsense");
	  else  
	    precision = tp / ( tp + fp);      
	  System.out.println( tp + " " + fp + " "+ repeat + " " + recall + " " + precision );
	  System.out.println();
	  recPre.add(recall+" " + precision );
	  recSum = recSum + recall;
	  preSum = preSum + precision;
	  indexLCR1 = 0;
          indexLCR2 = 0;
	  line1 = LCRBlocksRf1.readLine();
          line2 = LCRBlocksRf2.readLine();
	  ++ i;
	}
        else {
          line1 = LCRBlocksRf1.readLine();
          line2 = LCRBlocksRf2.readLine();
        }    
      }
      LCRBlocksRf1.close();
      LCRBlocksRf2.close();
    }
    catch ( IOException ex ) {
    }
    recall = recSum / i;
    precision = preSum / i;
    System.out.println();
    System.out.println("all recalls and pre");
    for ( int j = 0; j < recPre.size(); j++ ) 
      System.out.println( (String) recPre.elementAt(j));
    System.out.println();
    System.out.println( "max of ave recall and precison " + i + " " + recall + " " + precision );
  }


  public static void main (String args[]) {
    MaxRecallPreGettor mrpg = new MaxRecallPreGettor (args[0], args[1], args[2]);
    mrpg.getMaxRecallPre();
  }

}
