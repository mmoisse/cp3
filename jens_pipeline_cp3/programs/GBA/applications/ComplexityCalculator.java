/* 
** Author: Xuehui Li
** Date: April, 2005
** This program is used to compute the different kinds of complexities of a given string: entropy, the reciprocal probability without a scoring matrix, or the reciprocal probability when using a scoring matrix,  of the string.
** There are two input parameters: the first one is a string whose complexity is to be computed and the second one is an indicator( 'en' means calculating the entropy, 'rec' means calculating the reciprocal complexity without using a scoring a matrix, and 'sco' means calculating the reciprocal complexity when using a scoring a matrix, and 'all' means calculating all complexities.).
*/
package applications;

import java.io.*;
import java.util.*;

class ComplexityCalculator {

  private Vector alphabet;
  private double[][] scoringMatrix = new double[20][20]; 
  private double[] prob = new double[20];
  private Vector newAlphabet = new Vector();
  private double[][] newScoringMatrix = new double[400][400];

  public ComplexityCalculator( ) {
   
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


  public void getProb( String str ) {
    int index = 0;
    int  len = str.length();
    String subStr = str, letter = new String();
    for ( int i = 0; i < len; i++ ) {
      letter = subStr.substring( 0, 1 );
      index = alphabet.indexOf( letter );
      prob[index] = prob[index] +1;
      subStr = subStr.substring( 1 );
    }
  }
 

  public void clearProb() {
    for ( int i = 0; i < 20; i++ )
      prob[i] = 0;
  }
  
  // expression (1) in DSR paper is used here, i.e., Shannon Entropy  
  public double calculateEntropy( String str ) {
    double entropy = 0f;
    clearProb();
    getProb( str );
    double  len = str.length();
    for ( int i = 0; i < 20; i++ ) 
      if ( prob[i] != 0 ) {
	prob[i] = prob[i] / len;
	entropy = entropy + ( 0 - prob[i] * Math.log( prob[i] ));
      } 
    return entropy;
  }


  // expression (5) in DSR paper is used here
  public double calculateReciprocalPro( String str ) {
    double recPro = 1f;
    clearProb();
    getProb( str );
    double len = str.length();
    for ( int i = 0; i < 20; i++ ) 
      if ( prob[i] != 0 ) {
	prob[i] =  Math.pow( ( len / prob[i] ), prob[i] );
	recPro = recPro * prob[i];
      } 
    recPro = recPro * Math.pow( ( 1f / 20f ), len );
    return recPro;
  }


  public void printScoringMatrix(){
    for ( int i = 0; i < 20; i++ ) {
      for ( int j = 0; j < 20; j++) {
	if ( scoringMatrix[i][j] < 0 )
	  System.out.print( scoringMatrix[i][j] + " " );
	else 
	  System.out.print( " " + scoringMatrix[i][j] + " " );
      }
      System.out.println();
    }
  }


  public void readScoringMatrix( String fileName ){
    String line = new String(), score = new String();
    try {
      File f = new File( fileName );
      RandomAccessFile rf = new RandomAccessFile( f, "r");
      for ( int i = 0; i < 20; i++ ) {
	line = rf.readLine();
	//System.out.println( line);
	int k = 0;
	for ( int j = 0; j < 20; j++) {
	  score = line.substring( k, k + 2 ).trim();
       	  scoringMatrix[i][j] = Integer.parseInt( score );
	  //System.out.println( score + "**" + scoringMatrix[i][j] );  
	  k = k +3;
	}
      }
      rf.close();
    }
    catch ( IOException ex ) {
    }
    //printScoringMatrix();
    //System.out.println();
  }


  // 
  public void  normalizeScoringMatrixPow() {
    for ( int i = 0; i < 20; i++ ) 
      for ( int j = 0; j < 20; j++) {
	scoringMatrix[i][j] =  Math.pow( 2, scoringMatrix[i][j] );
      } 
    //printScoringMatrix();
    //System.out.println();
  }


  public void  normalizeScoringMatrix() {
    int min = -4, max = 11; 
    for ( int i = 0; i < 20; i++ ) 
      for ( int j = 0; j < 20; j++) {
	scoringMatrix[i][j] = ( scoringMatrix[i][j] - min ) / ( max - min);
      } 
    //printScoringMatrix();
    //System.out.println();
  }


  public void normalizeMore() {
    double sum = 0;
    for ( int i = 0; i < 20; i++ ){
      sum = 0;
      for ( int j = 0; j < 20; j++ )
	sum = sum + scoringMatrix[i][j];
      for ( int j = 0; j < 20; j++ ) 
	scoringMatrix[i][j] = scoringMatrix[i][j] / sum;
    }
    //printScoringMatrix();
    //System.out.println();
  }


  //  expression (6) in DSR paper and Blosum62 Scoring Matrix are used 
  public double calculateRecProWScoringMatrix( String str ) {
    //System.out.println( str );
    double recProWScoringMatrix;
    readScoringMatrix( "knowledge/blosum62Matrix" );
    normalizeScoringMatrix();
    clearProb();
    getProb( str );    
    double len = str.length(), sum = 0;
    String letter = new String();
    int index = 0;
    for ( int i = 0; i < 20; i++ ) {
      double vi = 0;
      if ( prob[i] != 0 ) {
      	for( int j = 0; j < len; j++ ){
	  letter = str.substring( j, j + 1 );
	  //System.out.print( letter + " " );
	  index = alphabet.indexOf( letter );
	  if ( scoringMatrix[i][index] == 0 )
	    vi = vi + 0.0333333333;
	  else 
	    vi = vi + scoringMatrix[i][index];
	  //System.out.println( vi );
	}
	vi = vi * Math.log( vi );
	//System.out.println( "vi: " + vi );
	sum = sum + vi;
      }
    }
    
    recProWScoringMatrix = len * Math.log( len );
    recProWScoringMatrix = recProWScoringMatrix - len * Math.log( 20D );
    //System.out.println( recProWScoringMatrix + "  "+ sum );
    recProWScoringMatrix = recProWScoringMatrix - sum;
    return recProWScoringMatrix;
  }
  
  
  /** 
   *  scoring matrix is used
   * only those letters in the alphabet who appear in the sequence  are included in ps[] 
   **/  
  public double calculateModifiedEntropy( String str ) {
    double modEntr = 0, singleProb = 0, sum = 0;
    readScoringMatrix( "knowledge/blosum62Matrix" );
    normalizeScoringMatrixPow();
    normalizeMore();
    clearProb();
    getProb( str );
    double len = str.length();
    for ( int i = 0; i < 20; i++ )
      prob[i] = prob[i] / len;
    double[] ps = new double[20];
    int k = 0; // total # of different letters in the sequence
    for ( int i = 0 ; i < 20; i++ ) {
      singleProb = 0;
      if ( prob[i] != 0 ){
	for ( int j = 0; j < 20; j++ ) {
	  singleProb = singleProb + prob[j] * scoringMatrix[i][j];
	  //System.out.println( prob[j] + "  " + scoringMatrix[i][j]);
	}
	//System.out.println( singleProb );
	ps[k] = singleProb;
	sum = sum + singleProb;
	k = k + 1;
      }
    }
    for ( int i = 0; i < k; i++ ) {
      ps[i] = ps[i] / sum;
      //System.out.println( ps[i]);
      modEntr = modEntr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    return modEntr;
  }


  // version2 
  //// allletters in the alphabet who appear in the sequence  are included in ps[] 
  public double calculateModifiedEntropy2( String str ) {
    double modEntr = 0, singleProb = 0, sum = 0;
    readScoringMatrix( "knowledge/blosum62Matrix" );
    normalizeScoringMatrix();
    normalizeMore();
    clearProb();
    getProb( str );
    double len = str.length();
    for ( int i = 0; i < 20; i++ ){
      prob[i] = prob[i] / len;
      //      System.out.println( prob[i] + "  ");
    }
    double[] ps = new double[20];
    for ( int i = 0 ; i < 20; i++ ) {
      singleProb = 0;
      for ( int j = 0; j < 20; j++ ) {
	singleProb = singleProb + prob[j] * scoringMatrix[i][j];
	//System.out.println( prob[j] + "  " + scoringMatrix[i][j]);
      }
      //System.out.println( singleProb );
      ps[i] = singleProb;
    }
    for ( int i = 0; i < 20; i++ ) {
      //System.out.println( ps[i]);
      modEntr = modEntr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    return modEntr;
  }

  
  /**
   *  the complexity is normalized by the sequence length. 
   *  Only those letters in the alphabet who appear in the sequence  are included in ps[] 
   **/
  public double calculateNorModifiedEntropy( String str ) {
    double modEntr = 0, singleProb = 0, sum = 0;
    readScoringMatrix( "knowledge/blosum62Matrix" );
    normalizeScoringMatrixPow();
    normalizeMore();
    clearProb();
    getProb( str );
    double len = str.length();
    for ( int i = 0; i < 20; i++ )
      prob[i] = prob[i] / len;
    double[] ps = new double[20];
    int k = 0; // total # of different letters in the sequence
    for ( int i = 0 ; i < 20; i++ ) {
      singleProb = 0;
      if ( prob[i] != 0 ){
	for ( int j = 0; j < 20; j++ ) {
	  singleProb = singleProb + prob[j] * scoringMatrix[i][j];
	  //System.out.println( prob[j] + "  " + scoringMatrix[i][j]);
	}
	//System.out.println( singleProb );
	ps[k] = singleProb;
	sum = sum + singleProb;
	k = k + 1;
      }
    }
    for ( int i = 0; i < k; i++ ) {
      ps[i] = ps[i] / sum;
      //System.out.println( ps[i]);
      modEntr = modEntr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    double strLen = str.length();////////////////////////////////
    modEntr = modEntr / strLen;//////////////////////////////////
    return modEntr;
  }


  public void initializeNewAlphabet() {
    int si = alphabet.size();
    for ( int i = 0; i < 20; i++ )
      for ( int j = 0; j < 20; j++ ) { 
	newAlphabet.add( (String) alphabet.elementAt( i ) + (String) alphabet.elementAt( j ));  
      }
  }

  
  public double calculateNor2LetterEntropy( String str ) {
    double let2Entr = 0;
    initializeNewAlphabet();
    double[] newProb = new double[400];
    for ( int i = 0; i < 400; i++ )
      newProb[ i ] = 0; 
    int index = 0;
    double len = str.length();
    String subStr = str, letters = new String();
    len = len - 1;
    for ( int i = 0; i < len; i++ ) {
      //System.out.println( i + " " + subStr );
      letters = subStr.substring( 0, 2 );
      index = newAlphabet.indexOf( letters );
      newProb[index] = newProb[index] +1;
      subStr = subStr.substring( 1 );
    }
    for ( int i = 0; i < 400; i++ ) 
      if ( newProb[i] != 0 ) {
	newProb[i] = newProb[i] / len;
	//System.out.println( newProb[ i ] );
	let2Entr = let2Entr + ( 0 - newProb[i] * Math.log( newProb[i] ));
      } 
    let2Entr = let2Entr / len; // normalize
    return let2Entr;
  }
  

  public void computeNewScoringMatrix() {
    readScoringMatrix( "knowledge/blosum62Matrix" );
    String lets1 = new String(), lets2 = new String();
    String letter1 = new String(), letter2 = new String();
    double score = 0; 
    for ( int i = 0; i < 400; i++ ) 
      for ( int j = 0; j < 400; j++) {
	score = 0;
	lets1 = ( String ) newAlphabet.elementAt( i );
	lets2 = ( String ) newAlphabet.elementAt( j );
	//System.out.println( lets1 + " " + lets2 );
	letter1 = lets1.substring( 0, 1 );
	letter2 = lets2.substring( 0, 1 );
	int index1 = alphabet.indexOf( letter1 );
	int index2 = alphabet.indexOf( letter2 );
	//System.out.println( scoringMatrix[ index1 ][ index2 ] );
	score = score + scoringMatrix[ index1 ][ index2 ];
	letter1 = lets1.substring( 1, 2 );
	letter2 = lets2.substring( 1, 2 );
	index1 = alphabet.indexOf( letter1 );
	index2 = alphabet.indexOf( letter2 );
	//System.out.println( scoringMatrix[ index1 ][ index2 ] );
	score = score + scoringMatrix[ index1 ][ index2 ];
	newScoringMatrix[i][j] = score / 2D;
	//System.out.println( newScoringMatrix[i][j] );
      }
  }


  public void testNewScoringMatrix( String str1, String str2 ) {
    initializeNewAlphabet();
    computeNewScoringMatrix();
    int index1 = newAlphabet.indexOf( str1 );
    int  index2 = newAlphabet.indexOf( str2 );
    // System.out.println( index1 + " " + index2 + " " + newScoringMatrix[ index1 ][ index2 ] );
  }


  public void normalizeNewScoringMatrixPow() {
    for ( int i = 0; i < 400; i++ ) 
      for ( int j = 0; j < 400; j++) {
	newScoringMatrix[i][j] =  Math.pow( 2, newScoringMatrix[i][j] );
      }     
    double sum = 0;
    for ( int i = 0; i < 400; i++ ){
      sum = 0;
      for ( int j = 0; j < 400; j++ )
	sum = sum + newScoringMatrix[i][j];
      for ( int j = 0; j < 400; j++ ) 
	newScoringMatrix[i][j] = newScoringMatrix[i][j] / sum;
    }
  }  


  public double calculate2LetterEntropyWScoMatrix( String str ) {
    double let2Entr = 0;
    /*
    initializeNewAlphabet();
    computeNewScoringMatrix();
    normalizeNewScoringMatrixPow();
    */
    double[] newProb = new double[400];
    for ( int i = 0; i < 400; i++ )
      newProb[ i ] =0; 
    int index = 0;
    double len = str.length();
    String subStr = str, letters = new String();
    len = len - 1;
    for ( int i = 0; i < len; i++ ) {   // get newProb[]
      //System.out.println( i + " " + subStr );
      letters = subStr.substring( 0, 2 );
      index = newAlphabet.indexOf( letters );
      newProb[index] = newProb[index] +1;
      subStr = subStr.substring( 1 );
    }
    double[] ps = new double[400];
    double singleProb = 0, sum = 0;
    int k = 0; // total # of different 2-letter combination in the sequence
    for ( int i = 0 ; i < 400; i++ ) {
      singleProb = 0;
      if ( newProb[i] != 0 ){
	for ( int j = 0; j < 400; j++ ) {
	  singleProb = singleProb + newProb[j] * newScoringMatrix[i][j];
	}
	//System.out.println( singleProb );
	ps[k] = singleProb;
	sum = sum + singleProb;
	k = k + 1;
      }
    }
    // System.out.println( "k= " + k );
    for ( int i = 0; i < k; i++ ) {
      ps[i] = ps[i] / sum;
      //System.out.println( ps[i]);
      let2Entr = let2Entr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    return let2Entr;
  }


   public double calculateNor2LetterEntropyWScoMatrix( String str ) {
    double let2Entr = 0;
    /*
    initializeNewAlphabet();
    computeNewScoringMatrix();
    normalizeNewScoringMatrixPow();
    */
    double[] newProb = new double[400];
    for ( int i = 0; i < 400; i++ )
      newProb[ i ] =0; 
    int index = 0;
    double len = str.length();
    String subStr = str, letters = new String();
    len = len - 1;
    for ( int i = 0; i < len; i++ ) {   // get newProb[]
      //System.out.println( i + " " + subStr );
      letters = subStr.substring( 0, 2 );
      index = newAlphabet.indexOf( letters );
      newProb[index] = newProb[index] +1;
      subStr = subStr.substring( 1 );
    }
    double[] ps = new double[400];
    double singleProb = 0, sum = 0;
    int k = 0; // total # of different 2-letter combination in the sequence
    for ( int i = 0 ; i < 400; i++ ) {
      singleProb = 0;
      if ( newProb[i] != 0 ){
	for ( int j = 0; j < 400; j++ ) {
	  singleProb = singleProb + newProb[j] * newScoringMatrix[i][j];
	}
	//System.out.println( singleProb );
	ps[k] = singleProb;
	sum = sum + singleProb;
	k = k + 1;
      }
    }
    // System.out.println( "k= " + k );
    for ( int i = 0; i < k; i++ ) {
      ps[i] = ps[i] / sum;
      //System.out.println( ps[i]);
      let2Entr = let2Entr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    let2Entr = let2Entr / len; // normalize
    return let2Entr;
  }


  /* calculate the unNormalized complexity of a sequence(i.e., all LCRs of a sequence) by using Nor2LetterEntropyWScoMatrix */
  public double calculateSeqUnNormComplexity (Vector lcrs) {
    double let2Entr = 0;
    double[] newProb = new double[400];
    for ( int i = 0; i < 400; i++ )
      newProb[ i ] =0; 
    int index = 0;
    double sumLen = 0, len = 0;
    String str = new String(), letters = new String();
    for ( int j = 0; j < lcrs.size(); j++ ) {
      str = (String) lcrs.elementAt(j);     
      len = str.length();
      //sumLen = sumLen + len;
      len = len - 1;
      for ( int i = 0; i < len; i++ ) {   // get newProb[]
	//System.out.println( i + " " + str );
	letters = str.substring( 0, 2 );
	index = newAlphabet.indexOf( letters );
	newProb[index] = newProb[index] + 1;
	str = str.substring( 1 );
      }
      //System.out.println(sumLen);
    }
    double[] ps = new double[400];
    double singleProb = 0, sum = 0;
    int k = 0; // total # of different 2-letter combination in the sequence
    for ( int i = 0 ; i < 400; i++ ) {
      singleProb = 0;
      if ( newProb[i] != 0 ){
	//System.out.println( (String)newAlphabet.elementAt(i) +" " + newProb[i]);
	for ( int j = 0; j < 400; j++ ) {
	  singleProb = singleProb + newProb[j] * newScoringMatrix[i][j];
	}
	//System.out.println( singleProb );
	ps[k] = singleProb;
	sum = sum + singleProb;
	k = k + 1;
      }
    }
    // System.out.println( "k= " + k );
    for ( int i = 0; i < k; i++ ) {
      ps[i] = ps[i] / sum;
      //System.out.println( ps[i]);
      let2Entr = let2Entr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    //let2Entr = let2Entr / sumLen; // normalize
    return let2Entr; 
  }  


  /* calculate the normalized complexity of a sequence(i.e., all LCRs of a sequence) by using Nor2LetterEntropyWScoMatrix */
  public double calculateSeqNormComplexity (Vector lcrs) {
    double let2Entr = 0;
    double[] newProb = new double[400];
    for ( int i = 0; i < 400; i++ )
      newProb[ i ] =0; 
    int index = 0;
    double sumLen = 0, len = 0;
    String str = new String(), letters = new String();
    for ( int j = 0; j < lcrs.size(); j++ ) {
      str = (String) lcrs.elementAt(j);     
      len = str.length();
      sumLen = sumLen + len;
      len = len - 1;
      for ( int i = 0; i < len; i++ ) {   // get newProb[]
	//System.out.println( i + " " + str );
	letters = str.substring( 0, 2 );
	index = newAlphabet.indexOf( letters );
	newProb[index] = newProb[index] + 1;
	str = str.substring( 1 );
      }
      //System.out.println(sumLen);
    }
    double[] ps = new double[400];
    double singleProb = 0, sum = 0;
    int k = 0; // total # of different 2-letter combination in the sequence
    for ( int i = 0 ; i < 400; i++ ) {
      singleProb = 0;
      if ( newProb[i] != 0 ){
	//System.out.println( (String)newAlphabet.elementAt(i) +" " + newProb[i]);
	for ( int j = 0; j < 400; j++ ) {
	  singleProb = singleProb + newProb[j] * newScoringMatrix[i][j];
	}
	//System.out.println( singleProb );
	ps[k] = singleProb;
	sum = sum + singleProb;
	k = k + 1;
      }
    }
    // System.out.println( "k= " + k );
    for ( int i = 0; i < k; i++ ) {
      ps[i] = ps[i] / sum;
      //System.out.println( ps[i]);
      let2Entr = let2Entr + ( 0 - ps[i] * Math.log( ps[i] ));    
    }
    let2Entr = let2Entr / sumLen; // normalize
    return let2Entr; 
  }  


  public void calculateComplexity( String str ) {
    System.out.println( "entropy: " + calculateEntropy( str ));
    System.out.println( " reciprocal probablity: " + calculateReciprocalPro( str ) );
  }

  
  public static void main ( String args[] ) {
    ComplexityCalculator cc = new ComplexityCalculator();
    cc. initializeAlphabet();
    if( args[1].equals( "en" ) ) 
      System.out.println( "entropy: " + cc.calculateEntropy( args[0] ));
    else if ( args[1].equals( "rec") )
      System.out.println( "reciprocal probablity:" + cc.calculateReciprocalPro( args[0] ));
    else if ( args[1].equals( "sco" ) )
      System.out.println( "reciprocal probablity with a scoring matrix: " +  cc.calculateRecProWScoringMatrix( args[0] ));
    else if ( args[1].equals( "me" ) ) 
      System.out.println( "modified entropy: " + cc.calculateModifiedEntropy( args[0] ) );
    else if ( args[1].equals( "me2" ) )
      System.out.println( "modified entropy version2: " + cc.calculateModifiedEntropy2( args[0] ) );
    else if (  args[1].equals( "nme" ))
      System.out.println( "normalized modified entropy: " + cc.calculateNorModifiedEntropy( args[0] ) );
    else if ( args[1].equals( "n2let" ))
      System.out.println( "2-letter entropy: " + cc.calculateNor2LetterEntropy( args[0] ));
    else if ( args[1].equals( "2letsco" )) {
      cc.initializeNewAlphabet();
      cc.computeNewScoringMatrix();
      cc.normalizeNewScoringMatrixPow();
      System.out.println( "2-letter entropy with a scoring matrix: " + cc.calculate2LetterEntropyWScoMatrix( args[0] ));
    }    
    else if ( args[1].equals( "n2letsco" )) {
      cc.initializeNewAlphabet();
      cc.computeNewScoringMatrix();
      cc.normalizeNewScoringMatrixPow();
      System.out.println( "normalized 2-letter entropy with a scoring matrix: " + cc.calculateNor2LetterEntropyWScoMatrix( args[0] )); 
    }
    else if ( args[2].equals( "test" ) )
      cc.testNewScoringMatrix( args[0], args[1] );
    else if ( args[1].equals("seq")) { 
      Vector lcrs = new Vector();
      lcrs.add( "QEANQEYQEPVCSPVPEPEPEPEPEPEP");
      lcrs.add("PPPEPQPEPEPQPLPDPAPLPE");
      lcrs.add("EAEPEP");
      cc.initializeAlphabet();
      cc.initializeNewAlphabet();
      cc.computeNewScoringMatrix();
      cc.normalizeNewScoringMatrixPow();
      System.out.println(cc.calculateSeqUnNormComplexity( lcrs ));
    }
    else
      cc.calculateComplexity( args[0] );
  }
  
}
