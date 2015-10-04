/*
** Author: Xuehui Li
** Date: April, 2005
** This program is used to learn the frequency patterns about the twenty amino acids in pfam-B repeats and non-pfam-B-repeats, respectively.
** The input consists of two parameters. The first one is a sequence file( say, /cise/research/tamer/xli/LCR/data/swissprot/sequenceInfor/seqFromFlybase ) and the second one is the corresponding file about repeats( say, /cise/research/tamer/xli/LCR/data/swissprot/repeatsInfor/flybaseRepeat ). The ith piece of information of the second file is the repeat information of the ith sequence in the first file.
The output is two matrices for repeats and non-repeats, repectively.
*/

package applications;

import java.io.*;
import java.util.*;

class SwissprotMatrixLearner {

  private Vector  alphabet;
  //private int[][] repeatMatrix;
  //private int[][] nonRepeatMatrix;
  private float[][] repeatMatrix;
  private float[][] nonRepeatMatrix;
  private  RandomAccessFile rfSeq, rfRep; 

  public SwissprotMatrixLearner ( String sequenceFile, String repeatFile ) {
    initializeAlphabet();
    //repeatMatrix = new int[20][20];
    //nonRepeatMatrix = new int[20][20];
    repeatMatrix = new float[20][20];
    nonRepeatMatrix = new float[20][20];
    try {
      File f = new File( sequenceFile );
      rfSeq = new RandomAccessFile( f, "r" );
      f = new File( repeatFile );
      rfRep = new RandomAccessFile( f, "r" );
    }
    catch ( IOException ex ) {
    }
  }


  public void initializeAlphabet() {
    alphabet = new Vector();
    alphabet.add( "A" );
    alphabet.add( "R" );
    alphabet.add( "N" );
    alphabet.add( "D" );
    alphabet.add( "C" );
    alphabet.add( "Q" );
    alphabet.add( "E" );
    alphabet.add( "G" );
    alphabet.add( "H" ); 
    alphabet.add( "I" ); 
    alphabet.add( "L" );  
    alphabet.add( "K" );  
    alphabet.add( "M" );  
    alphabet.add( "F" ); 
    alphabet.add( "P" );  
    alphabet.add( "S" );
    alphabet.add( "T" );  
    alphabet.add( "W" );  
    alphabet.add( "Y" );  
    alphabet.add( "V" );
  }


  public String generateConnectedSeq() {
    String str = new String(), seq = new String();
    try {
      str = rfSeq.readLine();
      while (( str != null ) && ( !( str.startsWith( ">" )))) {
	str = str.trim();
	seq = seq + str;
	str = rfSeq.readLine();
      }
    }
    catch ( IOException ex ){
    }
    return seq;
  }


  public String combineRepeats( String repeats ) {
    String str = repeats, repeat1 = new String(), repeat2 = new String();
    Vector subRepeats = new Vector();
    int k  = 0;
    str = repeats;
    while ( k != -1 ) {
      str = str.trim();
      k = str.indexOf( " " );
      if ( k != -1 )
	subRepeats.add( str.substring( 0, k ));
      else 
	subRepeats.add( str ); 
      if ( k != -1 )
	str = str.substring( k );
    }
    k = 0;
    int i = 0, j = 0, start = 0, end  = 0;
    while ( k < ( subRepeats.size() - 1 ) ) {
      repeat1 = (String)subRepeats.elementAt( k );
      repeat2 = (String)subRepeats.elementAt( k + 1 );
      i = repeat1.indexOf( "-" );
      end = Integer.parseInt( repeat1.substring( i + 1 ));
      j = repeat2.indexOf( "-" );
      start = Integer.parseInt( repeat2.substring( 0, j ));
      if ( start == ( end + 1 )) {
	str =  repeat1.substring( 0, i ) + "-" + repeat2.substring( j + 1 );
	subRepeats.remove( repeat1 );
	subRepeats.remove( repeat2 );
	subRepeats.add( k, str );
      }
      else 
	++k;
    }
    str = (String)subRepeats.elementAt( 0 );
    k = 1;
    while ( k < subRepeats.size()) {
      repeat1 = (String)subRepeats.elementAt( k );
      str = str + " " + repeat1;
      ++k;
    }
    return str;
  }


  //  adopt forget rate: ( 0.9 ) power k 
  public void workOnMatrix( String subSequence, int mark , float forgetRate ) {
    //int[][] matrix;
    float[][] matrix; 
    if ( mark == 0 )
      matrix = nonRepeatMatrix;
    else
      matrix = repeatMatrix;
    String subSeq = subSequence, tmpSeq = new String(), outsideLetter = new String(), insideLetter = new String();
    int row = 0, col = 0, len = 0;
    while ( subSeq.length() != 0 ) {
      tmpSeq = subSeq.substring( 1 );
      outsideLetter = subSeq.substring( 0, 1 );
      row = alphabet.indexOf( outsideLetter );
      //matrix[row][row] = matrix[row][row] + 1; //1f;
      int i = 1;
      while ( tmpSeq.length() != 0 ) {
	insideLetter = tmpSeq.substring( 0, 1 );
	col = alphabet.indexOf( insideLetter );
	Float f = new Float( Math.pow( forgetRate, i ) );
	float ff = f.floatValue();
	if  ( row == col ) 
	  matrix[row][col] = matrix[row][col] + ff;     
	else 
	  matrix[col][row] = matrix[col][row] + ff; 
	tmpSeq = tmpSeq.substring( 1 );
	//System.out.println( outsideLetter + " " + insideLetter + " " + row + " " + col + " " + matrix[row][col]);
	++i;
      }
      subSeq = subSeq.substring( 1 );      
    } 
    if ( mark == 0 )
      nonRepeatMatrix = matrix;
    else
      repeatMatrix = matrix;
  }


  /*  
// 0 for non-repeat, 1 for repeat. Time complexity: O(n*n).
  public void workOnMatrix( String subSequence, int mark ) {
    //int[][] matrix;
    float[][] matrix; 
    if ( mark == 0 )
      matrix = nonRepeatMatrix;
    else
      matrix = repeatMatrix;
    String subSeq = subSequence, tmpSeq = new String(), outsideLetter = new String(), insideLetter = new String();
    int row = 0, col = 0, len = 0;
    while ( subSeq.length() != 0 ) {
      tmpSeq = subSeq.substring( 1 );
      outsideLetter = subSeq.substring( 0, 1 );
      row = alphabet.indexOf( outsideLetter );
      matrix[row][row] = matrix[row][row] + 1; //1f;
      while ( tmpSeq.length() != 0 ) {
	insideLetter = tmpSeq.substring( 0, 1 );
	col = alphabet.indexOf( insideLetter );
	matrix[row][col] = matrix[row][col] + 1; //1f;    
	matrix[col][row] = matrix[col][row] + 1; //1f;
	tmpSeq = tmpSeq.substring( 1 );
      }
      subSeq = subSeq.substring( 1 );      
    } 
    if ( mark == 0 )
      nonRepeatMatrix = matrix;
    else
      repeatMatrix = matrix;
  }
  */


  public void workOnSequence( String seq, String rep, String len , float forgetRate) {
    String sequence = seq, subSequence = new String(), repeats = rep, subRepeat = new String();
    int nonRepeatStart = 0, nonRepeatEnd = 0, repeatStart = 0, repeatEnd = 0; 
    int index = 0, i = 0, j = 0;
    boolean beginWithOne = false;
    //System.out.println( seq );
    repeats = repeats.trim();
    if ( repeats.startsWith( "1-" )) 
      beginWithOne = true;
    else 
      nonRepeatStart = 1;
    while ( repeats != null ) {
      j = repeats.indexOf( "-" );
      i = repeats.indexOf( " " );
      if ( i == -1 ) {
	subRepeat = repeats;
	repeats = null;
      }
      else {
	subRepeat = repeats.substring( 0, i );
	repeats = repeats.substring( i + 1 ); 
	repeats = repeats.trim();    
      }
      repeatStart = Integer.parseInt( subRepeat.substring( 0, j ));
      if ( !beginWithOne ) {
	nonRepeatEnd =  repeatStart - 1;
	//System.out.println( "nonRepeat:  "+ nonRepeatStart + " " +  nonRepeatEnd );
	subSequence = sequence.substring( nonRepeatStart - 1, nonRepeatEnd );
	workOnMatrix( subSequence, 0, forgetRate ); // work on non-repeat
      }
      else 
	beginWithOne = false;
      repeatEnd = Integer.parseInt( subRepeat.substring( j + 1 ));
      //System.out.println( "repeat: " + repeatStart + " " + repeatEnd );
      subSequence = sequence.substring( repeatStart -1, repeatEnd );
      workOnMatrix( subSequence, 1, forgetRate ); // work on repeat
      nonRepeatStart = repeatEnd  + 1;     
    }
    if ( !(subRepeat.endsWith( len ))) { 
      nonRepeatEnd = Integer.parseInt( len );
      subSequence = sequence.substring( nonRepeatStart -1, nonRepeatEnd );
      workOnMatrix( subSequence, 0, forgetRate ); // work on non-repeat
    }
  }


  public void normalizeSingle( int mark ) {
    float sum = 0f;
    float[][] matrix;
    //int sum = 0;
    //int[][] matrix;
    if ( mark == 0 )
      matrix = nonRepeatMatrix;
    else
      matrix = repeatMatrix;
    for ( int i = 0; i < 20; i++ ) {
      sum = 0f;
      for ( int j = 0; j < 20; j++ ) 
	sum = sum + matrix[ i ][ j ];
      for ( int j = 0; j < 20; j++ )
	if ( sum != 0f )
	  matrix[ i ][ j ] = matrix[ i ][ j ] / sum; 
    }
    if ( mark == 0 )
      nonRepeatMatrix = matrix;
    else
      repeatMatrix = matrix;
  }
  

  public void normalizeBoth() {
    normalizeSingle( 0 );
    normalizeSingle( 1 ); 
  }


  public void learnPattern( float forgetRate ) {
    String id = new String(), seq =new String(), repeats = "dummy", len = new String();
    try { 
      rfSeq.readLine();
      seq = generateConnectedSeq();
      while (( seq != null ) && ( repeats != null )) {
	while ( repeats.indexOf( "-" ) == -1 ) { 
	  repeats = rfRep.readLine();
	}
	repeats = repeats.trim();
	len = rfRep.readLine().trim();
	repeats = combineRepeats( repeats );
	workOnSequence( seq, repeats, len, forgetRate );
	seq = generateConnectedSeq();
	repeats = rfRep.readLine();	
      }
    }
    catch ( IOException ex ) {
    }
    normalizeBoth();
  }
  

  public void closeBoth() {
    try {
      rfSeq.close();
      rfRep.close();
    }
    catch ( IOException ex ) {
    }  
  }


  public void printMatricesRowByRow() {
    System.out.println( "Non-Repeat matrix: " );
    for ( int i = 0; i < 20; i++ ) {
      for ( int j = 0; j < 20; j++ )
	System.out.print( nonRepeatMatrix[i][j] + "   " );
      System.out.println();
      System.out.println();
    }
    System.out.println();
    System.out.println();
    System.out.println( "Repeat matrix: " );
    for ( int i = 0; i < 20; i++ ) { 
      for ( int j = 0; j < 20; j++ )
	System.out.print( repeatMatrix[i][j] + "   " );
      System.out.println();
      System.out.println();
    }
    System.out.println();
  }


  public void printMatricesColByCol() {
    System.out.println( "Non-Repeat matrix: " );
    for ( int i = 0; i < 20; i++ ) {
      for ( int j = 0; j < 20; j++ )
	System.out.print( nonRepeatMatrix[j][i] + "   " );
      System.out.println();
      System.out.println();
    }
    System.out.println();
    System.out.println();
    System.out.println( "Repeat matrix: " );
    for ( int i = 0; i < 20; i++ ) { 
      for ( int j = 0; j < 20; j++ )
	System.out.print( repeatMatrix[j][i] + "   " );
      System.out.println();
      System.out.println();
    }
    System.out.println();
  }


  public static void main ( String args[] ) {
    
    SwissprotMatrixLearner swl = new SwissprotMatrixLearner( args[0],args[1] );
    float forgetRate = Float.parseFloat( args[2] );
    swl.learnPattern( forgetRate );
    swl.closeBoth();
    swl.printMatricesColByCol();
    
    // System.out.println( Math.pow( 2,3 ));
    /*
    float[][] m1, m2;
    m1 = new float[20][20];
    // m2 = new float[20][20];
    for ( int i = 0; i < 20; i++ ) 
      for ( int j = 0; j < 20; j++ )
	m1[i][j] = i + j;    

    m2 = m1;
    for ( int i = 0; i < 20; i++ ) { 
      for ( int j = 0; j < 20; j++ )
	System.out.print( m2[i][j] + "   " );
      System.out.println();
    }
    */
  }

}
