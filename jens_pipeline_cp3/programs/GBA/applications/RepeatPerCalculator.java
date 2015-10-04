/**
Author: Xuehui Li
Date: April, 2005
This is program is used to compute the percentage of repeats in a sequence: total lenght of repeats / total length of the sequence.
The input is from  
       metrics/swissprot/wForgetRate/postProcessing/sgd/sim095Sgd2LetModEntrNor2LetModEntrMer63CutAdjBlockSmiWatNCheComMetrics
The output is at 
/cise/research/tamer/xli/LCR/data/swissprot/repeatPercentage

**/
package applications;

import java.io.*;
import java.util.Vector;

class RepeatPerCalculator {

  private RandomAccessFile in;
  
  public RepeatPerCalculator( String repeatFile ) {
    try {
      File f = new File( repeatFile );
      in = new RandomAccessFile( f, "r" );
    }
    catch ( IOException ex ) {
    }
  }


  public void calculatePer( ) {
    String line = new String(), id = new String(), lenStr = new String(), totalRepStr = new String();
    float per, len, totalRep;
    Vector lenVec = new Vector(), perVec = new Vector();
    try {
      int i, j;
      line = in.readLine();
      while( line != null){
	if ( line.startsWith( "ID" )) {
	  line = line.trim();
	  i = line.indexOf( " " );
	  line = line.substring( i + 1 );
	  i = line.indexOf( " " );
	  id = line.substring( 0, i );
	  i = line.lastIndexOf( " " );
	  lenStr = line.substring( i + 1 );
	  lenVec.add( lenStr );
	  len = Float.parseFloat( lenStr );
	  while ( !(line.startsWith( "Total length of all repeats"))) 
	    line = in.readLine();
	  line = line.trim();
	  i = line.indexOf( ":" );
	  totalRepStr = line.substring( i + 1 );
	  totalRep = Float.parseFloat( totalRepStr );
	  per = totalRep / len;
	  perVec.add( String.valueOf(per));
	  //System.out.println( id +":   repeatPer = " +  totalRep + " / " + len + " = " + per );	
  	}
	line = in.readLine();
      }       
      in.close();
      i = 0;
      while ( i < lenVec.size()){
	lenStr = ( String ) lenVec.remove( i );
	per = Float.parseFloat( ( String ) perVec.remove( i ) );
	//System.out.println( lenStr + "  " + per );
	j = 0;
	while ( j != -1 ) {
	  j = lenVec.indexOf( lenStr );
	  if ( j != -1 ) {
	    //System.out.println( lenVec.remove( j ) + "  " + perVec.remove( j ) );
	    lenVec.remove( j );
	    per = ( per + Float.parseFloat( ( String )( perVec.remove( j ))) ) / 2f;
	  }
	}
	System.out.println( lenStr + "  " + per );
      }
    }catch ( IOException ex ) {
    
    }
  }


  public static void main( String args[] ) {
    RepeatPerCalculator rpc = new RepeatPerCalculator( args[0] );
    rpc.calculatePer();
  }
  
}
