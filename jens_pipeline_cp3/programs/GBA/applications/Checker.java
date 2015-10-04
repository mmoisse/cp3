package applications;

import java.io.*;

class Checker {

  private RandomAccessFile rf;

  public Checker ( String file ) {
    try{
      File f = new File( file );
      rf = new RandomAccessFile( f, "r" );
    }
    catch ( IOException ex ) {
    }
  }
  

  public void check() {
    try {
      String line = rf.readLine();
      int len = 0, lineNum = 1;
      while ( line != null ) {
	len = line.length();
	if ( len > 60 ){
	  System.out.println( "longer than 60:" + lineNum ); 
	}
	else if ( len == 0 ) {
	  System.out.println( "empty line:" + lineNum );
	}
	else if ( line.indexOf(">") != -1 )
	  if ( !( line.startsWith( ">" )) )
	    System.out.println( "wrong at: " + lineNum );
	line = rf.readLine();
	++ lineNum;
      }
      rf.close();
    }
    catch ( IOException ex ) {
    }
  }


  public static void main ( String args[] ) {
    Checker ck = new Checker( args[0] );
    ck.check();
  }
}
