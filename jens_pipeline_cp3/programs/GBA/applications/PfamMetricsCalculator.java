/*
** Author: Xuehui Li
** Date: April, 2005
** This program is used to compute three metrics: sensitivity, precision and recall. It collects information about LCRs, HCRs, non-domain regions and domain regions by comparing a file ( say, /cise/research/tamer/xli/LCR/data/Pfam/domainInfor/ABC\ tranDomain) that contains domain information and a second file ( /cise/research/tamer/xli/LCR/graphLCR/pfamLCRBlocks/ABC ) that contains LCR information of a sequence file. It assumes the ith sequence and the ith sequence refer to the same sequence( which is the function of Consistenter.java  ). The first parameter is a domain file, the second is a LCR sequence file. 
** sensitivity = ( # of letters in intersections of LCRs and non-domain regions in a sequence + # of letters in intersectioins of HCRs and domain regions in the same sequence) / total length of the same sequence
** precision = # of letters in intersections of LCRs and non-domain regions in a sequence / # of letters in all LCRs of the same sequence
** recall = # of letters in intersections of LCRs and non-domain regions in a sequence/ total # of letters in all non-domain regions of  the same sequence
*/
package applications;

import java.io.*;

class PfamMetricsCalculator {
  private String[] LCRs = new String[1000],  HCRs = new String[1000];
  private String[] domains = new String[100], nonDomains = new String[100];
  private String strDomain = new String();
  private int indexLCR = 0, indexHCR = 0, remainder = 1;
  private int indexDomain = 0, indexNonDomain = 0, indexMetric = 0;
  private float totalLCRs = 0, totalNonDomains = 0;
  private float[] sensitivity = new float[55], precision = new float[55], recall= new float[55]; 
  private boolean head = true;
  private RandomAccessFile rf1, rf2;
  private Integer[] domInt = new Integer[100];

  public PfamMetricsCalculator( String filename1, String filename2 ) {  
    try {
      File f1 = new File( filename1 ); // domain file
      File f2 = new File( filename2 ); // LCR file
      rf1 = new RandomAccessFile( f1, "r" );
      rf2 = new RandomAccessFile( f2, "r" );
      for ( int i = 0; i < 55; i++ ) {
	sensitivity[i]= 5555.0f;
	precision[i] = 5555.0f;
	recall[i] = 5555.0f;
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


  public float getInterLCRNonDomain() {
    float total = 0; 
    int i = 0, j = 0, k = 0;
    int lcrStartInt = 0, lcrEndInt = 0, nonDStartInt = 0, nonDEndInt = 0;
    String lcrStartStr = new String(), lcrEndStr = new String(); 
    String nonDStartStr = new String(), nonDEndStr = new String();
    boolean found = false;
    //System.out.println( "Intersections of LCRS and NonDomains:" );
    while ( i < indexLCR ) {
      k = LCRs[i].indexOf( "-" );
      lcrStartStr = LCRs[i].substring( 0, k );
      lcrStartInt = Integer.valueOf( lcrStartStr ).intValue();
      lcrEndStr = LCRs[i].substring( k + 1 );
      lcrEndInt = Integer.valueOf( lcrEndStr ).intValue();
      //System.out.println( lcrStartInt + "--" + lcrEndInt + "8888888888888888" );
      while ( ( !found ) && ( j < indexNonDomain )){
	k = nonDomains[j].indexOf( "-" );
	nonDStartStr = nonDomains[j].substring( 0, k );
	nonDStartInt = Integer.valueOf( nonDStartStr ).intValue();
	nonDEndStr = nonDomains[j].substring( k + 1 );
	nonDEndInt = Integer.valueOf( nonDEndStr ).intValue();
	//System.out.println( nonDEndInt + "  " + nonDStartInt + "  0000000000000000000" );
	if (( nonDEndInt >= lcrStartInt ) && (nonDStartInt <= lcrStartInt) && (nonDEndInt <= lcrEndInt)) {
	  //	  System.out.println( "11111111111111111111111111111" );
	  total = total + nonDEndInt - lcrStartInt + 1;	
	  //System.out.print( lcrStartStr + "-" + nonDEndStr + "  " );
	  ++j;
	} 
	else if ((  nonDEndInt <= lcrEndInt) && (  nonDStartInt >= lcrStartInt )) {
	  //System.out.println( "2222222222222222222222222222222222" );
	  total = total + nonDEndInt - nonDStartInt + 1;
	  //System.out.print( nonDStartStr + "-" + nonDEndStr + "  " );
	  ++j;
	}
	else if (( nonDEndInt > lcrEndInt)&&(nonDStartInt <= lcrEndInt)&&( nonDStartInt >= lcrStartInt)) {
	  //System.out.println( "3333333333333333333333" );
	  total = total + lcrEndInt - nonDStartInt + 1;
	  //System.out.print( nonDStartStr + "-" + lcrEndStr + "  " ); 
	  found = true;
	}
	else if ( nonDStartInt > lcrEndInt ) {	
	  //System.out.println( "4444444444444444444444" );
	  found = true;
	  }
	else if (( nonDStartInt <= lcrStartInt) && ( nonDEndInt >= lcrEndInt )) {
	  //System.out.println( "555555555555555555555555555" );
	  total = total + lcrEndInt - lcrStartInt + 1;
	  found = true;
	}
	else {

	  //System.out.println( "666666666666666666666nonDStart:" + nonDStartInt+" lcrStart:"+lcrStartInt+" nonDEnd:"+nonDEndInt+" lcrEnd:"+lcrEndInt); 
	  // found = true;
	  ++j;
	}
      }
      ++i;
      found = false;
    }
    System.out.println( "total length of intersections between LCRs and NonDomains:" + total );
    return total;
  }


  public float getInterHCRDomain() {
    float total = 0;
    int i = 0, j = 0, k = 0;
    int hcrStartInt = 0, hcrEndInt = 0, domStartInt = 0, domEndInt = 0;
    String hcrStartStr = new String(), hcrEndStr = new String(); 
    String domStartStr = new String(), domEndStr = new String();
    boolean found = false;
    while ( i < indexHCR ) {
      //System.out.println( " enter " );
      k = HCRs[i].indexOf( "-" );
      hcrStartStr = HCRs[i].substring( 0, k );
      hcrStartInt = Integer.valueOf( hcrStartStr ).intValue();
      hcrEndStr = HCRs[i].substring( k + 1 );
      hcrEndInt = Integer.valueOf( hcrEndStr ).intValue();
      //      System.out.println( hcrStartInt + "-" + hcrEndInt + "?88888888888888888" ); 
      while ( ( !found ) && ( j < indexDomain )){
	k = domains[j].indexOf( "-" );
	domStartStr = domains[j].substring( 0, k );
	domStartInt = Integer.valueOf( domStartStr ).intValue();
	domEndStr = domains[j].substring( k + 1 );
	domEndInt = Integer.valueOf( domEndStr ).intValue();
	if (( domEndInt >= hcrStartInt ) && (domStartInt <= hcrStartInt) && (domEndInt <= hcrEndInt)) {
	  //   System.out.println( domStartInt + "-" + domEndInt + "?1111111111111" );
	  total = total + domEndInt - hcrStartInt + 1;	
	  ++j;
	} 
	else if ((  domEndInt <= hcrEndInt) && (  domStartInt >= hcrStartInt )) {
	  //System.out.println( domStartInt + "-" + domEndInt + "?2222222222222222" );
	  total = total + domEndInt - domStartInt + 1;
	  ++j;
	}	else if (( domEndInt > hcrEndInt)&&(domStartInt <= hcrEndInt)&&( domStartInt >= hcrStartInt)) {
	  //System.out.println( domStartInt + "-" + domEndInt + "?3333333333333333" );
	  total = total + hcrEndInt - domStartInt + 1;
	  found = true;
	}
	else if ( domStartInt > hcrEndInt ) {	
	  //System.out.println( domStartInt + "-" + domEndInt + "?4444444444444444444" );
	  found = true;
	  }
	else if (( domStartInt <= hcrStartInt) && ( domEndInt >= hcrEndInt )) {
	  //System.out.println( domStartInt + "-" + domEndInt + "?5555555555555555555" );
	  total = total + hcrEndInt - hcrStartInt + 1;
	  found = true;
	} 
	else {
	  //System.out.println( domStartInt + "-" + domEndInt + "?66666666666666666666" );
	  ++j;
	}
      }
      ++i;
      //System.out.println( "i:" + i + "  indexHCR:" + indexHCR );
      found = false;
    }
    System.out.println( "total length of intersectsss between HCRs and Domains:" + total );
    return total;
  }
  

  public void computeMetrics( float LCRNonDomain, float HCRDomain, float totalNum ) {
    sensitivity[ indexMetric ] = ( LCRNonDomain + HCRDomain) / totalNum;
    if ( totalLCRs == 0 ) 
      precision[ indexMetric ] = 2222.0f;
    else 
      precision[ indexMetric] = LCRNonDomain / totalLCRs;
    // System.out.println( "totalNonDomains3: " + totalNonDomains + "   " + totalNum );
    if ( totalNonDomains == 0 ) 
      recall[indexMetric] = 2222.0f;
    else 
      recall[ indexMetric ] = LCRNonDomain / totalNonDomains;
    System.out.println(sensitivity[ indexMetric ] + "  " + precision[ indexMetric] + "  " +  recall[ indexMetric ]);
    ++indexMetric;
  }


  public void getLCRHCRDomainMetrics(){
    try {
     String strMasked = rf2.readLine();
      float totalNum = 0, LCRNonDomain = 0,  HCRDomain = 0;
      while ( strMasked != null ) {
	if ( strMasked.startsWith( ">" )) {
	  String numOfLetters = getSortedDomains();	  
	  System.out.println( "totalNonDomains1:" + totalNonDomains );
	  getLCRHCRs( numOfLetters );
	   
	  printLCRs();
	  System.out.println( "totalLCR:" + totalLCRs );
	  System.out.println( "<<<<<<<<<<" );
	  printHCRs();
	  LCRNonDomain = getInterLCRNonDomain();
	  HCRDomain = getInterHCRDomain();
	  totalNum = Integer.valueOf( numOfLetters ).intValue();	 
	  computeMetrics( LCRNonDomain, HCRDomain, totalNum );
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


  public void assignDomainInt() {
    for ( int i = 0; i < 100; i++ )
      domInt[i] = new Integer( 10000 );
  }


  public void printDomainInt() {
    for ( int i = 0; i < 3; i++ )
      System.out.println( domInt[i] );
    System.out.println("****************************");  
  }


  public void printDomains() {
    int i = 0;
    while ( i < indexDomain ) {
      System.out.print( domains[ i ] + " ");
      ++i;
    }
    System.out.println();
  }


  public void printNonDomains () {
    System.out.println( "Non-Domain Infor.:" );
    int i = 0;
    while ( i < indexNonDomain ) {
      System.out.print( nonDomains[ i ] + " ");
      ++i;
    }
    System.out.println();
  }


  public void getDomainStarts() {
    String str;
    int i = 0, j;
    while ( i < indexDomain ) {
      j = domains[i].indexOf( "-" );
      str = domains[i].substring( 0, j );
      //System.out.println( "start:" + str + "AA" );      
      domInt[i] = Integer.valueOf( str );
      ++i;
    }
  }


  public void sortDomains() {
    String str1, str2;
    String[] tmpDomains = new String[100];
    int i = 0, j, k;
    boolean found;
    while ( i < indexDomain ) {
      k = 0;
      found = false;
      str1 = Integer.toString( domInt[i] );
      //System.out.println( "start:" + str1 );
      while ( !found ) {
	//System.out.println( "domains:" + domains[k] );
	j = domains[k].indexOf( "-" );
	str2 = domains[k].substring( 0, j );
	//System.out.println( "str2:" + str2 );
	if ( str1.equals( str2 ) ){
	  //System.out.println( "found" );
	  tmpDomains[i] = domains[k];
	  found = true;
	}
	else
	  ++k;      
      }
      ++i;   
    }
    i = 0;
    while ( i < indexDomain ) {
      domains[i] = tmpDomains[i];
      ++i;
    } 
  }


  public void getNonDomains( int totalNum) {
    String startStr = new String(), endStr = new String(), str = new String();
    boolean beginWithOne = false; // the first domain begins with one
    int startInt = 0, endInt = 0;
    int i = 0, j = 0;
    while ( i < indexDomain ) {
      if ( i == 0 ) {
 	if (!( domains[i].startsWith( "1-" ))) {
	  startStr = "1";
	  startInt = 1;
	  endInt = domInt[i].intValue() - 1;
	}
	else 
	  beginWithOne = true;
      }
      else {
	beginWithOne = false;
	int k = domains[i-1].indexOf( "-" );
	str = domains[i-1].substring( k + 1 );
	startInt = Integer.valueOf( str ).intValue() + 1;
	endInt = domInt[i] - 1;
	startStr = String.valueOf( startInt );
	}
      if ( !beginWithOne )
	if ( startInt <= endInt ) {
	  endStr = String.valueOf( endInt );
	  nonDomains[j] = startStr + "-" + endStr;
	  totalNonDomains = totalNonDomains + endInt - startInt +1;
	  ++j;
	}
      ++i;
    }
    indexNonDomain = j;
    endStr = Integer.toString( totalNum );
    str = domains[ indexDomain - 1 ];
    if ( !( str.endsWith( endStr ))) {
      int k = str.indexOf( "-" );
      str = str.substring( k + 1 );
      startInt = Integer.valueOf( str ).intValue() + 1;
      startStr = String.valueOf( startInt );      
      nonDomains[indexNonDomain] = startStr + "-" + endStr;
      totalNonDomains = totalNonDomains + totalNum - startInt +1;
      ++indexNonDomain;
    }
  }


  
  public String getSortedDomains() {
    String id = new String(), numOfLetters = new String();
    indexDomain = 0;
     totalNonDomains = 0;
    try{
      if (  ( strDomain != null ) && (!strDomain.startsWith( ">" )))
	strDomain = rf1.readLine();    
      if ( strDomain != null )  {
	int k = 0; 
	int j = 0;
	if (strDomain.startsWith(">")) {
	  indexDomain = 0;
	  k = strDomain.lastIndexOf( "|" );
	  id = strDomain.substring( k+2, k+8 );
	  j = strDomain.indexOf( "a");
	  numOfLetters = strDomain.substring ( k+9, j-1 );
	  System.out.println();
	  System.out.println( "************************************************" );
	  System.out.println( id + "XXXX" + numOfLetters );
	  k = strDomain.indexOf( "|" );
	  strDomain = rf1.readLine();
	  //	  System.out.println(strDomain);
  while((strDomain!= null)&&(!strDomain.startsWith(">"))&&(strDomain.indexOf("-")!= -1)){
            strDomain = strDomain.substring( k-2); 
	    if ( strDomain.startsWith(" ")) 
	      strDomain = strDomain.substring( 1 ); 
	    j = strDomain.indexOf( " " );
	    String numOfDoms = strDomain.substring( 0, j );
	    //System.out.println( "numofdomains:" + numOfDoms + "AA");
	    if ( numOfDoms.equals( "1" ))  {
	      j = strDomain.lastIndexOf( " " );
	      domains[ indexDomain ] = strDomain.substring( j+1 );
	      //System.out.println( "domain:" + tmpDomain[ indexDomain ] );
	      ++ indexDomain;
	    }
	    else {
	      int numDom = Integer.valueOf( numOfDoms).intValue();
	      System.out.println( "numDom:" + numDom );
	      int i = 0;
	      while ( i < numDom ) {
		String tmpStr;
		j = strDomain.lastIndexOf( " " );
		domains[ indexDomain ] = strDomain.substring( j+1 );
		// System.out.println( "domain:" + tmpDomain[ indexDomain ] );
		++ indexDomain;
		strDomain = strDomain.substring( 0, j );
		++i;
	      }
	    }
	    strDomain = rf1.readLine();	 
	  }  
	}
	int i = 0;
	//System.out.println( "Domain Infor.:" );
	//printDomains();
	getDomainStarts();
	MergeSort.mergeSort( domInt );  
	sortDomains();
	System.out.println( "The sorted domains are:"  );
	printDomains();
	int totalNum = Integer.valueOf( numOfLetters ).intValue(); 
	getNonDomains( totalNum );
	printNonDomains();
	assignDomainInt();
      }
    }
    catch (IOException e) {
      
    }
    return numOfLetters;
  }

  
  public void printMetrics() {
    System.out.println();
    System.out.println( "     sensitivity       " + "precision       " + "recall" );
    int i = 0;
    float senAve = 0f, preAve = 0f, recAve = 0f;

    while ( i < indexMetric ) {
      System.out.println( ( i + 1 ) + "    " + sensitivity[i] + "          " + precision[i] + "          " + recall[i] );
      senAve = senAve +  sensitivity[i];
      if ( precision[i] !=  2222.0f ) 
	preAve = preAve + precision[i];
      if ( recall[i] !=  2222.0f )
	recAve = recAve + recall[i];
      ++i;
    }
    
    senAve = senAve / indexMetric;
    preAve = preAve /indexMetric;
    recAve = recAve / indexMetric;
    System.out.println( "Ave.  " + senAve  + "          " + preAve  + "          "+ recAve );
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
    PfamMetricsCalculator mc = new PfamMetricsCalculator( args[0], args[1]
);
    mc.assignDomainInt();
    mc.getLCRHCRDomainMetrics();
    mc.closeBoth();
    mc.printMetrics();
  }

}
