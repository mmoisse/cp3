/*
** Author: Xuehui Li
** Date: April, 2005
** This program is used to compute five metrics: sensitivity, precision, recall, Jarccard Coefficient, and Minkwoski. 
* It collects information about LCRs, HCRs, repeat regions and non-repeat regions by comparing a file ( say, /cise/research/tamer/xli/LCR/data/swissprot/repeatInfor/flybaseRepeat ) that contains repeat information and a second file ( /cise/research/tamer/xli/LCR/graphLCR/swissprotLCRBlocks/wOMatrices/flybaseLCRBlocks) that contains LCR information of a sequence file.
* It assumes the ith sequence and the ith sequence refer to the same sequence( which is the function of Consistenter.java  ). 
* The first parameter is a repeat information file.
* The second is a LCR sequence file. 
** sensitivity = ( # of letters in intersections of LCRs and repeat regions in a sequence + # of letters in intersectioins of HCRs and non-repeat regions in the same sequence) / total length of the same sequence
** precision = # of letters in intersections of LCRs and repeat regions in a sequence / # of letters in all LCRs of the same sequence
** recall = # of lettejava applications.Gbm ../data/swissprot/sequenceInfor/seqFromSgd swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095 3 15 5 knowledge/lcrs 0 > tmp/samplingResults/sgdrs in intersections of LCRs and repeat regions in a sequence/ total # of letters in all repeat regions of  the same sequence
* There are three input parameters.
* The first two are the two files to be compared.
* The last one is the output file to be generated.
* An example to run the program:
* java applications.SwissprotMetricsCalculator ../data/swissprot/repeatsInfor/flybaseRepeat swissprotLCRBlocks/wForgetRateMatrices/postProcess/cutting/meanStd/flyabseMeanStd > metrics/swissprot/wForgetRate/postProcessing/cutting/meanStd/flybaseMeanStdMetrics 
*/

package applications;

import java.io.*;
import java.text.*;

class SwissprotMetricsCalculator {
  private String[] LCRs = new String[5000],  HCRs = new String[5000];
  private String[] domains = new String[100], nonDomains = new String[100], tmpNonDomains = new String[100];
  private String strDomain = new String();
  private int indexLCR = 0, indexHCR = 0, remainder = 1;
  private int indexDomain = 0, indexNonDomain = 0, indexTmpNonDomain = 0, indexMetric = 0;
  private double totalLCRs = 0, totalRepeats = 0;
  private double[] sensitivity = new double[200], precision = new double[200], recall= new double[200], jc = new double[200], mm = new double[200]; 
  private boolean head = true;
  private RandomAccessFile rf1, rf2;


  public SwissprotMetricsCalculator( String filename1, String filename2 ) {  
    try {
      File f1 = new File( filename1 ); // repeat file
      File f2 = new File( filename2 ); // LCR file
      rf1 = new RandomAccessFile( f1, "r" );
      rf2 = new RandomAccessFile( f2, "r" );
      for ( int i = 0; i < 55; i++ ) {
	sensitivity[i]= 5555.0d;
	precision[i] = 5555.0d;
	recall[i] = 5555.0d;
	jc[i] = 5555.0d;
	mm[i] = 5555.0d;
      }
    }
    catch ( IOException ex ) {
    
    }  
  }


  public void getLCRHCRs( String numOfLetters ) {
    int i = 0, j = 0;
    String str = new String(), strTmp = new String();
    String startLCRStr = new String(), endLCRStr = new String();
    String startHCRStr = new String(), endHCRStr = new String();
    int startHCRInt = 0, endHCRInt = 0; 
    boolean beginWithOne = false; // LCR begins with "1-"???
    try {
      str = rf2.readLine(); 
    }
    catch ( IOException ex ) {
      
    }
    strTmp = str.trim();
    if ( !(strTmp.startsWith( "1-" ))) 
      startHCRStr = "1";
    else 
      beginWithOne = true;
    if ( strTmp.length() == 0 )
      strTmp = null;
    //System.out.println( strTmp );
    while ( strTmp != null ) {
      i = strTmp.indexOf( " " );
      j = strTmp.indexOf( "-" );
      startLCRStr = strTmp.substring( 0, j );
      if ( i == -1 ) {
	LCRs[ indexLCR ] = strTmp.substring( 0 );
	endLCRStr = strTmp.substring( j + 1 );
	strTmp = null;
      }     
      else {
	LCRs[ indexLCR ] = strTmp.substring( 0, i );
	endLCRStr = strTmp.substring( j + 1, i );
	strTmp = strTmp.substring( i + 1 );     
      }
      totalLCRs = totalLCRs + Integer.parseInt( endLCRStr ) - Integer.parseInt( startLCRStr ) + 1;
      endHCRInt = Integer.parseInt( startLCRStr ) - 1;
      if ( !beginWithOne ) {
	HCRs[ indexHCR ] = startHCRStr + "-" + endHCRInt;
	++indexHCR;
      }      
      else 
	beginWithOne = false;
      startHCRStr = Integer.toString ( Integer.parseInt( endLCRStr ) + 1 );
      ++ indexLCR;
    }
    strTmp = str.trim();
    if ( (!strTmp.endsWith( "-" + numOfLetters ))) {
      HCRs[ indexHCR ] = startHCRStr + "-" + numOfLetters;
      ++indexHCR;
    }
  }


  public void printLCRs() {
    int i = 0;
    System.out.println( "LLLLLLLLLLLLLLLLLLLLLLLLLLL" );
    while ( i < indexLCR ) {
      System.out.print( LCRs[i] + " " );
      ++ i;
    }
    System.out.println();
    System.out.println( "LLLLLLLLLLLLLLLLLLLLLLLLLLL" );
  }


  public void printHCRs() {
    int i = 0;
    System.out.println( "HHHHHHHHHHHHHHHHHHHHHHHHHHH" );
    while ( i < indexHCR ) {
      System.out.print( HCRs[i] + " " );
      ++ i;
    }
    System.out.println();
    System.out.println( "HHHHHHHHHHHHHHHHHHHHHHHHHHH" );
  }


  public double getInterLCRRepeat() {
    double total = 0; 
    int i = 0, j = 0, k = 0;
    int lcrStartInt = 0, lcrEndInt = 0, nonDStartInt = 0, nonDEndInt = 0;
    String lcrStartStr = new String(), lcrEndStr = new String(); 
    String nonDStartStr = new String(), nonDEndStr = new String();
    boolean found = false;
    while ( i < indexLCR ) {
      k = LCRs[i].indexOf( "-" );
      lcrStartStr = LCRs[i].substring( 0, k );
      lcrStartInt = Integer.valueOf( lcrStartStr ).intValue();
      lcrEndStr = LCRs[i].substring( k + 1 );
      lcrEndInt = Integer.valueOf( lcrEndStr ).intValue();
      while ( ( !found ) && ( j < indexNonDomain )){
	k = nonDomains[j].indexOf( "-" );
	nonDStartStr = nonDomains[j].substring( 0, k );
	nonDStartInt = Integer.valueOf( nonDStartStr ).intValue();
	nonDEndStr = nonDomains[j].substring( k + 1 );
	nonDEndInt = Integer.valueOf( nonDEndStr ).intValue();
	if (( nonDEndInt >= lcrStartInt ) && (nonDStartInt <= lcrStartInt) && (nonDEndInt <= lcrEndInt)) {
	  total = total + nonDEndInt - lcrStartInt + 1;	
	  ++j;
	} 
	else if ((  nonDEndInt <= lcrEndInt) && (  nonDStartInt >= lcrStartInt )) {
	  total = total + nonDEndInt - nonDStartInt + 1;
	  ++j;
	}
	else if (( nonDEndInt > lcrEndInt)&&(nonDStartInt <= lcrEndInt)&&( nonDStartInt >= lcrStartInt)) {
	  total = total + lcrEndInt - nonDStartInt + 1;
	  found = true;	}
	else if ( nonDStartInt > lcrEndInt ) {	
	  found = true;
	  }
	else if (( nonDStartInt <= lcrStartInt) && ( nonDEndInt >= lcrEndInt )) {
	  total = total + lcrEndInt - lcrStartInt + 1;
	  found = true;
	}
	else {
	  ++j;
	}
      }
      ++i;
      found = false;
    }  
    System.out.println( "total length of intersections between LCRs and repeats:" + total );
    return total;
  }


  public double getInterHCRNonRepeat() {
    double total = 0;
    int i = 0, j = 0, k = 0;
    int hcrStartInt = 0, hcrEndInt = 0, domStartInt = 0, domEndInt = 0;
    String hcrStartStr = new String(), hcrEndStr = new String(); 
    String domStartStr = new String(), domEndStr = new String();
    boolean found = false;
    while ( i < indexHCR ) {
      k = HCRs[i].indexOf( "-" );
      hcrStartStr = HCRs[i].substring( 0, k );
      hcrStartInt = Integer.valueOf( hcrStartStr ).intValue();
      hcrEndStr = HCRs[i].substring( k + 1 );
      hcrEndInt = Integer.valueOf( hcrEndStr ).intValue();
      while ( ( !found ) && ( j < indexDomain )){
	k = domains[j].indexOf( "-" );
	domStartStr = domains[j].substring( 0, k );
	domStartInt = Integer.valueOf( domStartStr ).intValue();
	domEndStr = domains[j].substring( k + 1 );
	domEndInt = Integer.valueOf( domEndStr ).intValue();
	if (( domEndInt >= hcrStartInt ) && (domStartInt <= hcrStartInt) && (domEndInt <= hcrEndInt)) {
	  total = total + domEndInt - hcrStartInt + 1;	
	  ++j;
	} 
	else if ((  domEndInt <= hcrEndInt) && (  domStartInt >= hcrStartInt )) {
	  total = total + domEndInt - domStartInt + 1;
	  ++j;
	}
	else if (( domEndInt > hcrEndInt)&&(domStartInt <= hcrEndInt)&&( domStartInt >= hcrStartInt)) {
	  total = total + hcrEndInt - domStartInt + 1;
	  found = true;
	}
	else if ( domStartInt > hcrEndInt ) {	
	  found = true;
	  }
	else if (( domStartInt <= hcrStartInt) && ( domEndInt >= hcrEndInt )) {
	  total = total + hcrEndInt - hcrStartInt + 1;
	  found = true;
	} 
	else {
	  ++j;
	}
      }
      ++i;
      found = false;
    }
    System.out.println( "Total length of intersectsss between HCRs and nonRepteats:" + total );
    return total;
  }
  

  public void computeMetrics( double LCRRepeat, double HCRNonRepeat, double totalNum ) {
    sensitivity[ indexMetric ] = ( LCRRepeat + HCRNonRepeat ) / totalNum;
    if ( totalLCRs == 0 ) 
      precision[ indexMetric ] = 2222.0f;
    else 
      precision[ indexMetric] = LCRRepeat / totalLCRs;
    if ( totalRepeats == 0 ) 
      recall[indexMetric] = 2222.0f;
    else 
      recall[ indexMetric ] = LCRRepeat / totalRepeats;
    jc[ indexMetric ] = LCRRepeat / ( totalNum - HCRNonRepeat );
    mm[ indexMetric ] = Math.sqrt( ( totalLCRs - LCRRepeat + totalRepeats - LCRRepeat) / totalRepeats );
    System.out.println(sensitivity[ indexMetric ] + "  " + precision[ indexMetric] + "  " +  recall[ indexMetric ] + "  " + jc[indexMetric] + "  " + mm[indexMetric]);
    ++indexMetric;
  }



 public void getNonRepeats( int totalNum) {
    String startStr = new String(), endStr = new String(), str = new String();
    boolean beginWithOne = false; // the first repeat begins with one
    int startInt = 0, endInt = 0;
    int i = 0, j = 0;
    while ( i < indexNonDomain ) {
      if ( i == 0 ) {
 	if (!( nonDomains[i].startsWith( "1-" ))) {
	  startStr = "1";
	  startInt = 1;
	  int k = nonDomains[i].indexOf( "-" );
	  endStr = nonDomains[i];
	  endStr = endStr.substring( 0 , k );/////////////
	  endInt = Integer.valueOf( endStr ).intValue() - 1;
	}
	else 
	  beginWithOne = true;
      }
      else {
	beginWithOne = false;
	int k = nonDomains[i-1].indexOf( "-" );
	str = nonDomains[i-1].substring( k + 1 );
	startInt = Integer.valueOf( str ).intValue() + 1;
	endStr = nonDomains[i];
	k = endStr.indexOf( "-" );
	endStr = endStr.substring( 0, k );
	endInt = Integer.valueOf( endStr ).intValue() - 1;
	startStr = String.valueOf( startInt );
      }
      if ( !beginWithOne )
	if ( startInt <= endInt ) {
	  endStr = String.valueOf( endInt );
	  domains[j] = startStr + "-" + endStr;
	  ++j;
	}
      ++i;
    }
    indexDomain = j;
    endStr = Integer.toString( totalNum );
    str = nonDomains[ indexNonDomain - 1 ];
    if ( !( str.endsWith( endStr ))) {
      int k = str.indexOf( "-" );
      str = str.substring( k + 1 );
      startInt = Integer.valueOf( str ).intValue() + 1;
      startStr = String.valueOf( startInt );      
      domains[indexDomain] = startStr + "-" + endStr;
      ++indexDomain;
    }
  }



  public void printRepeats( int k) {
    int i = 0;
    if ( k == 0 ){
      System.out.println( "tmpRepeats" );
      while ( i < indexTmpNonDomain ) {
	System.out.print( tmpNonDomains[ i ] + " ");
	++i;
      }
      System.out.println();
    }
    else {
      System.out.println( "Repeat Infor.:" );
      while ( i < indexNonDomain ) {
	System.out.print( nonDomains[ i ] + " ");
	++i;
      }
      System.out.println();
      System.out.println( "Total length of all repeats:" + totalRepeats );
    }
  }


  public void printNonRepeats () {
    System.out.println( "Non-repeat Infor.:" );
    int i = 0;
    while ( i < indexDomain ) {
      System.out.print( domains[ i ] + " ");
      ++i;
    }
    System.out.println();
  }



 public void  combineRepeats() {
    String str = new String(), str1 = new String();
    int i = 0, k = 0, j = 0;
    while ( i < indexTmpNonDomain ) {
      str1 = tmpNonDomains[ i ];
      k = str1.indexOf( "-" );
      str1 = str1.substring( 0, k );
      k = j;
      while ( ( k == j ) && ( i < indexTmpNonDomain )) {
	str = tmpNonDomains[i];
	k = str.indexOf( "-" );
	str = str.substring( k + 1 );
	j = Integer.valueOf( str ).intValue() + 1; 
	if (( i + 1 ) < indexTmpNonDomain ) {	
	  str =  tmpNonDomains[i+1];
	  k = str.indexOf( "-" );
	  str = str.substring( 0, k );
	  k = Integer.valueOf( str ).intValue();
	}  
	++i;
      }
      if ( k != j ) {
	str =  tmpNonDomains[i-1];
	k = str.indexOf( "-" );
	str = str.substring( k + 1 );
	k = Integer.valueOf( str1 ).intValue();
	j = Integer.valueOf( str ).intValue();
	nonDomains[ indexNonDomain ] = str1 + "-" + str;
	totalRepeats = totalRepeats + j - k + 1;
	++indexNonDomain;
      }
    }
  }



  public String getRepeatsNonRepeats() {
    String id = new String(), numOfLetters = new String();
    indexNonDomain = 0;
    totalRepeats = 0; // NonDomains = 0;
     try{
       int k = 0; 
       while ( !(strDomain.startsWith( ">" )))
	 strDomain = rf1.readLine();	 
       indexNonDomain = 0;
        indexTmpNonDomain = 0;
       id = strDomain.substring( 1 );
       System.out.println();
       System.out.println( "************************************************" );
       strDomain = rf1.readLine();
       while ( k != -1 ) {
	 strDomain = strDomain.trim();	
	 k = strDomain.indexOf( " " );
	 if ( k != -1 )
	   tmpNonDomains[ indexTmpNonDomain ] = strDomain.substring( 0, k );
	 else 
	   tmpNonDomains[ indexTmpNonDomain ] = strDomain;
	 ++ indexTmpNonDomain;
	 if ( k != -1 )
	   strDomain = strDomain.substring( k );
       }
       numOfLetters = rf1.readLine();
       System.out.println(  "ID: " + id + "    TotalNumOfLetters: " + numOfLetters );
       combineRepeats();
       printRepeats( 1 );
       int totalNum = Integer.valueOf( numOfLetters ).intValue(); 
       getNonRepeats( totalNum );
       printNonRepeats(); 
     }
     catch (IOException e) {
      
     }
     return numOfLetters;
  }



  public void getLCRHCRRepeatNonRepeatMetrics(){
    try {
     String strMasked = rf2.readLine();
      double totalNum = 0, LCRRepeat = 0,  HCRNonRepeat = 0;
      while ( strMasked != null ) {
	if ( strMasked.startsWith( ">" )) {
	  String numOfLetters = getRepeatsNonRepeats();
	  getLCRHCRs( numOfLetters );
	  printLCRs();
	  System.out.println( "totalLCR:" + totalLCRs );
	  System.out.println( "<<<<<<<<<<" );
	  printHCRs();
	  LCRRepeat = getInterLCRRepeat();
	  HCRNonRepeat = getInterHCRNonRepeat();
	  totalNum = Integer.valueOf( numOfLetters ).intValue();	 
	  computeMetrics( LCRRepeat, HCRNonRepeat, totalNum );
	  indexLCR = 0;
	  indexHCR = 0;
	  totalLCRs = 0;
	  strMasked = rf2.readLine();
	}
	else 
	  strMasked = rf2.readLine();
      
      }
    }
    catch (IOException ex ) {
      
    }
  }
  

  public void printMetrics() {
    System.out.println();
    System.out.println( "      prec.     " + "recall      " + "Jac.Co.     " + "Min.Mea." );
    int i = 0;
    double senAve = 0d, preAve = 0d, recAve = 0d, jcAve = 0d, mmAve = 0d;
    NumberFormat nf = NumberFormat.getInstance();
    while ( i < indexMetric ) {
      System.out.println( ( i + 1)+ "     " + nf.format( precision[i] ) + "      " + nf.format( recall[i] ) + "      " + nf.format( jc[i] ) + "      " + nf.format(mm[i]) );
      senAve = senAve +  sensitivity[i];
      if ( precision[i] !=  2222.0f ) 
	preAve = preAve + precision[i];
      if ( recall[i] !=  2222.0f )
	recAve = recAve + recall[i]; 
      jcAve = jcAve + jc[i];
      mmAve = mmAve + mm[i];
      ++i;
    }
    senAve = senAve / indexMetric;
    preAve = preAve /indexMetric;
    recAve = recAve / indexMetric;
    jcAve = jcAve / indexMetric;
    mmAve = mmAve /indexMetric;
    System.out.println( "Ave." + "   " + nf.format( preAve )  + "        "+ nf.format( recAve ) + "        " + nf.format( jcAve ) + "        " + nf.format(mmAve) );
    System.out.println();
  }


  public void closeBoth() {   
    try{
      rf1.close();
      rf2.close();
    }
    catch ( IOException ex ) {
      
    }
  }


  public static void main( String[] args) {
    SwissprotMetricsCalculator mc = new SwissprotMetricsCalculator( args[0], args[1]
);
    mc.getLCRHCRRepeatNonRepeatMetrics();
    mc.closeBoth();
    mc.printMetrics();
  }

}
