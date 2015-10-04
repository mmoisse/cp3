/**
 * Date: July, 2005
 * This program is used to get the repeat information based on seqIDs
 * The first input file is the combined repeat information file of the 5 
 * datasets in FAST format:
 * /cise/research38/tamer/xli/LCR/data/swissprot/repeatsInfor/combinedRepeats
 * The second input file specifies the IDs of those sequences whose repeat 
 * information to be extracted. 
 * This input file is generated from seqIdGettor.java: /cise/research38/tamer/xli/LCR/blastExp/input/seqId.
 * The output is a file that contains only repeat information of those sequences with ID in the second input file: 
 * /cise/research38/tamer/xli//LCR/blastExp/input/repeats.
 **/

package applications;

import java.io.*;

class RepeatGettor {

  RandomAccessFile repRf;
  RandomAccessFile idRf;  

  public RepeatGettor( String repFile, String idFile ) {
    try{
      File f = new File( repFile );
      repRf = new RandomAccessFile( f, "r" );
      f = new File( idFile ); 
      idRf = new RandomAccessFile( f, "r" );
    } 
    catch ( IOException ex ) {
    }
  }
  
  
  public void getRepeats() {
    String id = new String();
    String repeats = new String();
    try {
      id = idRf.readLine();
      while ( id != null ) {
	id = ">" + id;
	while ( !(repeats.startsWith( id )))
	  repeats = repRf.readLine();
	System.out.println( repeats );
	repeats = repRf.readLine();
	System.out.println( repeats );
	System.out.println();
	id = idRf.readLine();
      } 
      repRf.close();
      idRf.close();
    }
    catch ( IOException ex ) {
    }
  }


  public static void main ( String args[] ) {
    RepeatGettor rg = new RepeatGettor( args[0], args[1] );
    rg.getRepeats();
  }
}
