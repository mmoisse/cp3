/**
 * Date: June, 2005
 * 128 sample points (or sequences) are divided into 32 intervals 
 * Format of the output file: "knowledge/sampledLenRepPer":
 * shortestLength longestLength meanOfRepPer-stdOfRepPer meanOfRepPer+stdOfRepPer
**/

package applications;

import java.io.*;
import java.util.*;

class LenRepPerSampler {
  
  RandomAccessFile combinedRf;
  RandomAccessFile sampledRf;

  public LenRepPerSampler() {
    try {
      File f = new File( "/cise/research/tamer/xli/LCR/data/swissprot/repeatPercentage/combinedLenRepPer.txt" );
      combinedRf = new RandomAccessFile( f, "r" );
      f = new File( "knowledge/sampledLenRepPer");
      sampledRf = new RandomAccessFile( f, "rw" );
    }
    catch ( IOException ex ) {
    }
  }


  public void sample() {
    Random generator = new Random();
    int r, m;
    String line = new String(), len = new String(), repPer = new String();
    double aveRepPer = 0, expecVal = 0, expecSqu = 0, tmpD = 0;
    double shortestLen = 0, longestLen = 0;
    Integer[] sortedR = new Integer[128];
    Vector randoms = new Vector();
    Integer tmpIg;
    // the total number of sequences is 474
    int k = 0;
    while ( k < 128 ) {  // 128 sample sequences
      r = generator.nextInt(475);
      if ( r != 0 ) {
        tmpIg = new Integer( r );
        if ( randoms.indexOf( tmpIg ) == -1 ) { // sampling without replacement
          randoms.add( tmpIg );
          sortedR[k] = tmpIg;
	  ++k;
        }
	// System.out.println( r );
      }
    }
    MergeSort.mergeSort( sortedR );    
    try {
      for ( int i = 0; i < 128; i++ ) { 
      	r = sortedR[i].intValue();
	System.out.println( r );
	combinedRf.seek( 0 );
	for ( int j = 0; j < r; j++ )
	  line = combinedRf.readLine(); 
	line = line.trim();
	System.out.println( "line:"+ line );
	int index = line.indexOf( " " );
	len = line.substring( 0, index );
	repPer = line.substring( index + 1 );
	tmpD =  Double.parseDouble( repPer ); 
	expecVal = expecVal + tmpD / 8d;
	expecSqu = expecSqu + Math.pow( tmpD,2d ) / 8d;
	aveRepPer = aveRepPer + tmpD;	
	m = i % 8;
	if ( m == 0 )
	  sampledRf.writeBytes( len + " " ); // shortest length	  
	if ( m == 7 ) {
	  aveRepPer = aveRepPer / 8d;
	  expecSqu = Math.sqrt( expecSqu - Math.pow( expecVal, 2d ) ); 
	  expecVal = Math.abs( aveRepPer - expecSqu ); // mean - std
	  aveRepPer = aveRepPer + expecSqu; // mean + std
	  sampledRf.writeBytes(len+ " "+ expecVal+" " +aveRepPer+ "\r\n");//longest length
	  aveRepPer = 0;
	  expecVal = 0;
	  expecSqu = 0;
	  System.out.println();
	}
      }
      combinedRf.close();
      sampledRf.close(); 
    }
    catch ( IOException ex ) {
    }  
  }

  public Vector readSampledLenRepPer() {
    Vector v = new Vector();
    try {
      File f = new File( "knowledge/sampledLenRepPer" );
      RandomAccessFile rf = new RandomAccessFile( f, "r" );
      String line = rf.readLine();
      while ( line != null ) {
	v.add( line );
	line = rf.readLine();    
      }
      rf.close();
    }
    catch ( IOException ex ) {
    }
    return v;
  }
  

  public void useSample( int len ) {
    double limit = 0;
    Vector sample = readSampledLenRepPer();
    // find the right cut percentage
    boolean found = false;
    String range = new String();
    double per, prePer = 0;
    int i = 0, shortest, longest, preLongest = 0;
    while ( !found ) {
      range = (String)sample.elementAt( i );
      range = range .trim();
      int index = range.indexOf(",");
      shortest = Integer.parseInt( range.substring( 0, index ));
      int index2 = range.indexOf( " " );
      longest = Integer.parseInt( range.substring( index+1, index2 ) );
      index = range.indexOf(":");
      per =  Double.parseDouble( range.substring( index + 1));
      if ( ( len >= shortest ) && ( len <= longest)) {
	limit = len * ( 1 - per );
	found = true;
	System.out.println( "1: " + per);
      }
      else if ( len < shortest ) {
	int diff1 = shortest - len;
	int diff2 = len - preLongest;
	if ( ( diff2 > diff1 ) || ( preLongest == 0 )) {  
	  limit = len * ( 1- per );
	  found = true;
	  System.out.println( "2: " + per);
	}
	else {
	  limit = len * ( 1- prePer );
	  found = true;
	  System.out.println( "3: " + prePer);
	}
      }
      ++i;
      preLongest = longest;
      prePer = per;
    }
  }


  public static void main ( String args[] ) {
    LenRepPerSampler lrps = new LenRepPerSampler();
    //lrps.sample();
    lrps.useSample( Integer.parseInt( args[0] ) );  
  }

}
