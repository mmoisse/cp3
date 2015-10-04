/*
** Author: Xuehui Li
** Date: April, 2005
** This program is used to learn the frequency patterns about the twenty amino acids in domains and non-domains, respectively.
** The input consists of two parameters. The first one is a sequence file( say, /cise/research/tamer/xli/LCR/data/Pfam/sequenceInfor/ABC\ tran.txt ) and the second one is the corresponding domain file( say, /cise/research/tamer/xli/LCR/data/Pfam/domainInfor/ABC\ tranDomain ). The ith piece of information of the second file is the domain information of the ith sequence in the first file.
The outpur is two matrices for domain regions and non-domain regions, repectively.
*/

package applications;

import java.io.*;
import java.util.*;

class PfamMatrixLearner {

  private Vector  alphabet;
  //private int[][] pfamAMatrix;
  //private int[][] nonDomainMatrix;
  private float[][] pfamAMatrix;
  private float[][] nonDomainMatrix;
  private  RandomAccessFile rfSeq, rfDm; 

  public PfamMatrixLearner ( String sequenceFile, String domainFile ) {
    initializeAlphabet();
    //pfamAMatrix = new int[20][20];
    //nonDomainMatrix = new int[20][20];
    pfamAMatrix = new float[20][20];
    nonDomainMatrix = new float[20][20];
    try {
      File f = new File( sequenceFile );
      rfSeq = new RandomAccessFile( f, "r" );
      f = new File( domainFile );
      rfDm = new RandomAccessFile( f, "r" );
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


  // 0 for non-domain, 1 for pfamA-domain. Time complexity: O(n*n).
  public void workOnMatrix( String subSequence, int mark ) {
    //int[][] matrix;
    float[][] matrix; 
    if ( mark == 0 )
      matrix = nonDomainMatrix;
    else
      matrix = pfamAMatrix;
    String subSeq = subSequence, tmpSeq = new String(), outsideLetter = new String(), insideLetter = new String();
    int row = 0, col = 0, len = 0;
    while ( subSeq.length() != 0 ) {
      tmpSeq = subSeq.substring( 1 );
      outsideLetter = subSeq.substring( 0, 1 );
      row = alphabet.indexOf( outsideLetter );
      if (( row < 0 ) || ( row >19 )) 
	  System.out.println( "outsideLetter: " + outsideLetter + " " + row );
      matrix[row][row] = matrix[row][row] + 1f;
      while ( tmpSeq.length() != 0 ) {
	insideLetter = tmpSeq.substring( 0, 1 );
	col = alphabet.indexOf( insideLetter );
	if (( col < 0 ) || ( col > 19 )) 
	  System.out.println( "insideLetter:" + insideLetter + " " + col );
	matrix[row][col] = matrix[row][col] + 1f;    
	matrix[col][row] = matrix[col][row] + 1f;
	tmpSeq = tmpSeq.substring( 1 );
      }
      subSeq = subSeq.substring( 1 );      
    } 
    if ( mark == 0 )
      nonDomainMatrix = matrix;
    else
      pfamAMatrix = matrix;
  }


  public void workOnSequence( String seq, String dom, String len ) {
    String sequence = seq, subSequence = new String(), domains = dom, subDomain = new String();
    int nonDomainStart = 0, nonDomainEnd = 0, domainStart = 0, domainEnd = 0; 
    int index = 0, i = 0, j = 0;
    boolean beginWithOne = false;
    //System.out.println( len ); 
    //System.out.println( len + " " + seq );
   domains = domains.trim();
    if ( domains.startsWith( "1-" )) 
      beginWithOne = true;
    else 
      nonDomainStart = 1;
    while ( domains != null ) {
      j = domains.indexOf( "-" );
      i = domains.indexOf( " " );
      if ( i == -1 ) {
	subDomain = domains;
	domains = null;
      }
      else {
	subDomain = domains.substring( 0, i );
	domains = domains.substring( i + 1 ); 
	domains = domains.trim();    
      }
      //System.out.println( domains + " " + j );
      domainStart = Integer.parseInt( subDomain.substring( 0, j ));
      if ( !beginWithOne ) {
	nonDomainEnd =  domainStart - 1;
	System.out.println( len + " " + nonDomainStart + " " + nonDomainEnd );
	System.out.println( sequence );
	subSequence = sequence.substring( nonDomainStart - 1, nonDomainEnd );
	workOnMatrix( subSequence, 0 ); // work on non-domain
      }
      else 
	beginWithOne = false;
      domainEnd = Integer.parseInt( subDomain.substring( j + 1 ));
      System.out.println( len + " " + domainStart + " " + domainEnd );
      System.out.println( sequence );
      subSequence = sequence.substring( domainStart -1, domainEnd );
      workOnMatrix( subSequence, 1 ); // work on domain
      nonDomainStart = domainEnd  + 1;     
    }
    if ( !(subDomain.endsWith( len ))) { 
      nonDomainEnd = Integer.parseInt( len );
      System.out.println( len + " " + nonDomainStart + " " + nonDomainEnd );
      System.out.println( sequence );
      subSequence = sequence.substring( nonDomainStart -1, nonDomainEnd );
      //System.out.println( subSequence );
      workOnMatrix( subSequence, 0 ); // work on non-domain
    }
  }


  public void normalizeSingle( int mark ) {
    float sum = 0f;
    float[][] matrix;
    //int sum = 0;
    //int[][] matrix;
    if ( mark == 0 )
      matrix = nonDomainMatrix;
    else
      matrix = pfamAMatrix;
    for ( int i = 0; i < 20; i++ ) {
      sum = 0f;
      for ( int j = 0; j < 20; j++ ) 
	sum = sum + matrix[ i ][ j ];
      for ( int j = 0; j < 20; j++ )
	if ( sum != 0f )
	  matrix[ i ][ j ] = matrix[ i ][ j ] / sum; 
    }
    if ( mark == 0 )
      nonDomainMatrix = matrix;
    else
      pfamAMatrix = matrix;
  }
  

  public void normalizeBoth() {
    normalizeSingle( 0 );
    normalizeSingle( 1 ); 
  }


  public void learnPattern() {
    String id = new String(), seq =new String(), domains = "dummy", len = new String();
    try { 
      rfSeq.readLine();
      seq = generateConnectedSeq();
      while (( seq != null ) && ( domains != null )) {
	while (( domains.length() == 0 ) || ( domains.indexOf( "-" ) == -1 )) { 
	  domains = rfDm.readLine();
	}
	if ( domains != null ) {
	  domains = domains.trim();
	  len = rfDm.readLine().trim();
	  workOnSequence( seq, domains, len );
	  seq = generateConnectedSeq();
	  domains = rfDm.readLine();	
	}
      }
    }
    catch ( IOException ex ) {
    }
    normalizeBoth();
  }
  

  public void closeBoth() {
    try {
      rfSeq.close();
      rfDm.close();
    }
    catch ( IOException ex ) {
    }  
  }


  public void printMatrices() {
    System.out.println( "Non-Domain matrix: " );
    for ( int i = 0; i < 20; i++ ) { 
      for ( int j = 0; j < 20; j++ )
	System.out.print( nonDomainMatrix[i][j] + "   " );
      System.out.println();
      System.out.println();
    }
    System.out.println();
    System.out.println();
    System.out.println( "Pfam-A matrix: " );
    for ( int i = 0; i < 20; i++ ) {
      for ( int j = 0; j < 20; j++ )
	System.out.print( pfamAMatrix[i][j] + "   " );
      System.out.println();
      System.out.println();
    }
    System.out.println();
  }


  public static void main ( String args[] ) {
    PfamMatrixLearner pml = new PfamMatrixLearner( args[0],args[1] );
    pml.learnPattern();
    pml.closeBoth();
    pml.printMatrices();
   
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
