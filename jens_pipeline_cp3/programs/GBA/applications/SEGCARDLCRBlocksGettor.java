package applications;

import java.io.*;

class SEGCARDLCRBlocksGettor {

  private RandomAccessFile rf;

  public SEGCARDLCRBlocksGettor( String fileName) {
    try { 
      File f = new File( fileName );
      rf = new RandomAccessFile( f, "r" );      
    }
    catch ( IOException ex ) {
    }
  }


  public void get() {
    String line = new String();
    try { 
      line = rf.readLine();
      while ( line != null ) {
	while (( line != null ) && ( !(line.startsWith("ID:"))) ) {
	  line = rf.readLine();
	}
	if (( line != null) && ( line.startsWith( "ID"))) {
	  System.out.println( line );	
	  line = rf.readLine();
	  while ( !(line.startsWith("LLLLLL" ))) 
	    line = rf.readLine();
	  line = rf.readLine();
	  System.out.println( line );
	  System.out.println( );
	  line = rf.readLine();
	}
      }
      rf.close();
    }
    catch ( IOException ex ) {
    }
  }


  public static void main ( String args[] ) {

    SEGCARDLCRBlocksGettor sclg = new SEGCARDLCRBlocksGettor( args[0] );
    sclg.get();  
  }

}
