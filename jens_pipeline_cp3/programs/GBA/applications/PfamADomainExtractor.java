/*
** Author: Xuehui Li
** Date: April, 2005
** This program is used to extract pfam-A domain information of a sequence from the complete domain information of the sequence.
** There are two input parameters: the first one is the name of the file containing the complete information of sequences( say, /cise/research/tamer/xli/LCR/data/Pfam/domainInfor/COX1Domain) and the second one is the name of the file to be generated( say, /cise/research/tamer/xli/LCR/graphLCR/pfamADomainInfor/COX1PfamADomain), which only contains pfam-A domain information of sequences.
**The output is the newly generated file which contains pfam-A domain information.
*/


package applications;

import java.io.*;

class PfamADomainExtractor {

  private File f;
  private RandomAccessFile rfReader;
  private String[] domains = new String[100];
  private Integer[] domInt = new Integer[100];
  private int indexDomain = 0;

  public PfamADomainExtractor ( String completeDomainFile ) {
    try {
      f = new File ( completeDomainFile );
      rfReader = new RandomAccessFile ( f, "r" );
    }
    catch ( IOException ex ) {
    }
  }


  public void assignDomainInt() {
    for ( int i = 0; i < 100; i++ )
      domInt[i] = new Integer( 10000 );
  }


  public void printDomains() {
    int i = 0;
    while ( i < indexDomain ) {
      System.out.print( domains[ i ] + " ");
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


  public void sortDomains( ) {
    String str1, str2;
    String[] tmpDomains = new String[100];
    int i = 0, j, k;
    boolean found;
    while ( i < indexDomain ) {
      k = 0;
      found = false;
      str1 = Integer.toString( domInt[i] );
      while ( !found ) {
	j = domains[k].indexOf( "-" );
	str2 = domains[k].substring( 0, j );
	if ( str1.equals( str2 ) ){
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


  public void getSortedPfamADomains() {
    String id = new String(), numOfLetters = new String(), strDomain = new String();
    try{
      while ( strDomain != null ) {
	while (  ( strDomain != null ) && (!strDomain.startsWith( ">" )))
	  strDomain = rfReader.readLine();    
	if ( strDomain != null )  {
	  assignDomainInt();
	  int k = 0; 
	  int j = 0;
	  if (strDomain.startsWith(">")) {
	    indexDomain = 0;
	    k = strDomain.lastIndexOf( "|" );
	    id = strDomain.substring( k+2, k+8 );
	    j = strDomain.indexOf( "a");
	    numOfLetters = strDomain.substring ( k+9, j-1 );
	    k = strDomain.indexOf( "|" );
	    strDomain = rfReader.readLine();
	    while((strDomain!= null)&&(!strDomain.startsWith(">"))&&(strDomain.indexOf("-")!= -1)){
	      if (  !( strDomain.startsWith( "Pfam-B" ))) {  
		strDomain = strDomain.substring( k-2); 
		if ( strDomain.startsWith(" ")) 
		  strDomain = strDomain.substring( 1 ); 
		j = strDomain.indexOf( " " );
		String numOfDoms = strDomain.substring( 0, j );
		if ( numOfDoms.equals( "1" ))  {
		  j = strDomain.lastIndexOf( " " );
		  domains[ indexDomain ] = strDomain.substring( j+1 );
		  ++ indexDomain;
		}
		else {
		  int numDom = Integer.valueOf( numOfDoms).intValue();
		  int i = 0;
		  while ( i < numDom ) {
		    String tmpStr;
		    j = strDomain.lastIndexOf( " " );
		    domains[ indexDomain ] = strDomain.substring( j+1 );
		    ++ indexDomain;
		    strDomain = strDomain.substring( 0, j );
		    ++i;
		  }
		}
	      }
	      strDomain = rfReader.readLine();	 
	    }  
	  }
	  int i = 0;
	  if ( indexDomain != 0 ) {
	    System.out.println(">" +  id );
	    getDomainStarts();
	    MergeSort.mergeSort( domInt );  
	    sortDomains();
	    printDomains();
	    System.out.println( numOfLetters );
	    System.out.println();
	  }
	}
      }
      rfReader.close();
    }
    catch (IOException e) {
      
    }
  }


  public static void main ( String args[] ) {
    PfamADomainExtractor pade = new PfamADomainExtractor( args[ 0 ] );  
    pade.getSortedPfamADomains();  
  }

}
