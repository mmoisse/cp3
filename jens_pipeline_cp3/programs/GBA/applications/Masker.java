/**
* Author: Xuehui Li 
* Date: June 28, 2005
* This program is used to generate maksed sequence files. Low-complexity 
* regions (LCRs) are masked with "x". 
* When 'mask()' is called, the input is a sequence file with a sequence and 
* a second file that specifies which parts (positions) of the sequence are 
* LCRs. The output is the masked sequence file.
* When 'maskMutiple()' is called, the first input is changed to a sequence 
* file that can contain more than one sequence.
* 'mask()' is used by gbm.java to generate masked sequences once after lcr
* blocks are found.  
* 'maskMutiple()' is used in main(). The example is :   
* java applications.Masker /cise/research38/tamer/xli/LCR/blastExp/input/unmaskedSeqs /cise/research38/tamer/xli/LCR/blastExp/input/repeast > /cise/research38/tamer/xli/LCR/blastExp/input/maskedSeqs
**/
package applications;

import java.io.*;

class Masker{

  private RandomAccessFile srf;
  private RandomAccessFile prf;
 
  public Masker( String seqFile, String posFile ) {
    try {
      File f = new File ( seqFile );
      srf = new RandomAccessFile ( f, "r" );
      f = new File ( posFile );
      prf = new RandomAccessFile ( f, "r" );
    }
    catch ( IOException ex ) {
    }
  }


  public Masker () {

  }

  
  public String getSeqPosBlocks( ) {
    String seqLine = new String(), seq = new String(), posLine;
    try {
      seqLine = srf.readLine();
      System.out.println( seqLine ) ;      
      seqLine = srf.readLine();
      while ( seqLine != null ) {
	seq = seq + seqLine;
	seqLine = srf.readLine();
      }
      posLine = prf.readLine();
      posLine = prf.readLine();
      posLine = posLine.trim();
      seq = seq + "*" + posLine;
    }
    catch ( IOException ex ) {
    }
    return seq;
  }
  

  public void writeLCR( int start, int end ) {
    int i, j, div1, div2;
    div1 = start / 60;
    div2 = end /60;
    if ( div1 == div2 ) {
      if (( start % 60) != 0 ) 
	i =  end - start + 1;
      else {
	System.out.println( "x" );
  	i =  end - start; 
      }
      j = 0;
      while ( j < i ) {
	System.out.print( "x" );
	++ j;
      }
    }
    else {
      ++div1;
      j = 0;
      i = 60 - start % 60 + 1;
      //System.out.println( i );
      while ( j < i ) {
	System.out.print( "x" );
	++ j;	
      }
      System.out.println();
      while ( div1!= div2 ) {
	++ div1;
	j = 0;
	while ( j < 60 ) {
	  System.out.print( "x" );
	  ++ j;
	}
	System.out.println();
      }
      j = 0;
      i = end % 60 ;
      while ( j < i ) {
	System.out.print( "x" );
	++ j;
      }
    }
  }


  public void mask( String seq, String posBlocks  ) {
    String hcr;
    int indexSpace, indexDash, start, end, curPos, div1, div2;
    curPos = 1;
    indexSpace = posBlocks.indexOf( " " );
    indexDash = posBlocks.indexOf( "-" );
    start = Integer.parseInt( posBlocks.substring( 0, indexDash ) );
    if ( indexSpace == -1 ) {
      end = Integer.parseInt( posBlocks.substring( indexDash + 1 ));
    }
    else {
      end= Integer.parseInt(posBlocks.substring(indexDash+1,indexSpace));
      posBlocks = posBlocks.substring( indexSpace + 1 ); 
    }
    
    if ( start != 1 )
      curPos = 1;
    else { // write the lcr beginning at the first letter of the sequence
      writeLCR( start, end );
      curPos = end + 1;
      // get start and end
      indexSpace = posBlocks.indexOf( " " );
      indexDash = posBlocks.indexOf( "-" );
      start = Integer.parseInt( posBlocks.substring( 0, indexDash ) );
      if ( indexSpace == -1 ) {
	end = Integer.parseInt( posBlocks.substring( indexDash + 1 ));
      }
      else {
	end= Integer.parseInt(posBlocks.substring(indexDash+1,indexSpace));
	posBlocks = posBlocks.substring( indexSpace + 1 ); 
      }
    }
    //System.out.println( "?" + start + "?"+ end + "?" );
    boolean finished = false;
    while ( ! finished ) {
      if ( curPos > start ) // eg. 18-49 45 - 58
	start = curPos;
      if ( curPos < start ) { // write HCR
	div1 = curPos / 60;
	div2 = ( start - 1 ) / 60;
	if ( div1 == div2 ){
	  if ( ( curPos % 60 ) != 0 ) 
	    hcr = seq.substring( curPos-1 , start - 1 );
	  else {
	    System.out.println( seq.substring( curPos-1, curPos));
	    hcr = seq.substring( curPos, start - 1 );
	  }
	  System.out.print( hcr );
	}
	else {
	  if ((curPos % 60) == 0 ) {
	    System.out.println( seq.substring( curPos -1, curPos ));
	    ++curPos;
	  }
	  
	  while ( div1!= div2 ) {
	    ++ div1;
	    //System.out.println( curPos-1 + ":" + 60 * div1 );
	    hcr = seq.substring( curPos - 1, 60 * div1 );
	    System.out.println( hcr );
	    curPos = 60 * div1 + 1;
	  }
	  hcr = seq.substring( curPos-1, start - 1 );
	  System.out.print( hcr );
	}
      }
      //System.out.println(start +"?"+ end);
      writeLCR( start, end ); // write LCR
      //System.out.println( "????????????????" );
      curPos = end + 1;
      
      if ( indexSpace != -1 ) {
	// get start and end
	indexSpace = posBlocks.indexOf( " " );
	indexDash = posBlocks.indexOf( "-" );
	start = Integer.parseInt( posBlocks.substring( 0, indexDash ) );
	if ( indexSpace == -1 ) {
	  end = Integer.parseInt( posBlocks.substring( indexDash + 1 ));
	}
	else {
	  end= Integer.parseInt(posBlocks.substring(indexDash+1,indexSpace));
	  posBlocks = posBlocks.substring( indexSpace + 1 ); 
	}
      }
      else
	finished = true;
    }
    
    // write the last block of hcr that ends at the end of the sequence
    if ( end != seq.length() ) {
      start = seq.length() + 1;
      div1 = curPos / 60;
      div2 = ( start - 1 ) / 60;
      if ( div1 == div2 ){
	if ( ( curPos % 60 ) != 0 ) 
	  hcr = seq.substring( curPos-1 , start - 1 );
	else {
	  System.out.println( seq.substring( curPos-1, curPos));
	  hcr = seq.substring( curPos, start - 1 );
	}
	System.out.print( hcr );
	/*
	if ( ( start -1 ) % 60 == 0 )
	  System.out.println();
	*/
      }
      else {
	
	if ((curPos % 60) == 0 ) {
	  System.out.println( seq.substring( curPos -1, curPos ));
	  ++curPos;
	}
	++ div1;
	hcr = seq.substring( curPos - 1, 60 * div1 );
	System.out.println( hcr );
	curPos = 60 * div1;
	while ( div1!= div2 ) {
	  ++ div1;
	  hcr = seq.substring( curPos, 60 * div1 );
	  System.out.println( hcr );
	  curPos = 60 * div1;
	}
	hcr = seq.substring( curPos, start - 1 );
	System.out.print( hcr );
      }	
    }
  }


  public void maskMutiple() {
    String seqLine = new String(), seq = new String(), posLine;
    try {  
      seqLine = srf.readLine();
      seqLine = srf.readLine();
      while ( seqLine != null ) {
	seq = new String();
	while ( ( seqLine != null ) && ( !( seqLine.startsWith( ">")))) {
	  seq = seq + seqLine;
	  seqLine = srf.readLine();
	}
	//System.out.println( seq );
	posLine = prf.readLine();
	System.out.println( posLine );
	posLine = prf.readLine();
	//System.out.println( posLine );
	posLine = posLine.trim();	
	mask( seq, posLine );
	System.out.println();
	posLine = prf.readLine();
	seqLine = srf.readLine();
      }
    }
    catch ( IOException ex ) {
    }
  }


  public void closeAll() {
    try{
      srf.close();
      prf.close();
    }
    catch ( IOException ex ) {
    }
  }
  
  
  public static void main ( String args[] ) {
      
    // call maskMutiple()
    Masker ms = new Masker( args[0], args[1] );
    ms.maskMutiple();
    ms.closeAll();
    /*
    // call 'mask()
    masker ms = new masker( args[0], args[1] );  
    String seqPosBlocks = ms.getSeqPosBlocks();
    String seq, posBlocks;
    int index = seqPosBlocks.indexOf( "*" );
    seq = seqPosBlocks.substring( 0, index );
    posBlocks = seqPosBlocks.substring( index + 1 );
    ms.mask( seq, posBlocks );
    ms.closeAll();
    */
  }

}
