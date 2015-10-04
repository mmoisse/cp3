/**
 * Date: July, 2005
 * This program is used to get sequences from seqIDs
 * The first input file is the combined sequence file of the 5 datasets in FAST
 * format: /cise/research38/tamer/xli/LCR/data/swissprot/sequenceInfor/combinedSeq.
 * The second input file specifies the IDs of those sequences to be extracted. 
 * This input file is generated from seqIdGettor.java: /cise/research38/tamer/xli/LCR/blastExp/input/seqId.
 * The output is a sequence file that contains only sequences with ID in the 
 * second input file: /cise/research38/tamer/xli//LCR/blastExp/input/unmaskedSeqs. 
 **/

package applications;

import java.io.*;

class SeqGettor {

  RandomAccessFile seqRf;
  RandomAccessFile idRf;  

  public SeqGettor( String seqFile, String idFile ) {
    try{
      File f = new File( seqFile );
      seqRf = new RandomAccessFile( f, "r" );
      f = new File( idFile ); 
      idRf = new RandomAccessFile( f, "r" );
    } 
    catch ( IOException ex ) {
    }
  }
  
  
  public void getSequences() {
    String id = new String();
    String seq = new String();
    try {
      id = idRf.readLine();
      while ( id != null ) {
	id = ">" + id;
	while ( !(seq.startsWith( id )))
	  seq = seqRf.readLine();
	System.out.println( seq );
	seq = seqRf.readLine();
	while ( !(seq.startsWith( ">" ))) { 
	  System.out.println( seq );
	  seq = seqRf.readLine();
	}
	id = idRf.readLine();
      } 
      seqRf.close();
      idRf.close();
    }
    catch ( IOException ex ) {
    }
  }


  public static void main ( String args[] ) {
    SeqGettor gs = new SeqGettor( args[0], args[1] );
    gs.getSequences();
  }
}
