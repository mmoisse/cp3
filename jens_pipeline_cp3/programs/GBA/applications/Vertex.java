package applications;

public class Vertex {

  private String startLetter;
  private String endLetter;
  private int indegree;
  private float weight;
  private float lp;  // the logest path from a source vertex to this vertex 
  private Vertex predecessor;

  public Vertex () {
    startLetter = new String();
    endLetter = new String();
    indegree = 0;
    weight = 0.0f;
    lp = -8888.0f;
    predecessor = null;
  }


  public Vertex ( String start, String end, int d, float w ) {
    startLetter = start;
    endLetter = end;
    indegree = d;
    weight = w;
    lp = 0.0f; // modified?????
    predecessor = null;
  }


  public String getStartLetter() {
    return startLetter;
  }


  public String getEndLetter() {
    return endLetter;
  }
  

  public int getIndegree() {
    return indegree;
  }


  public float getWeight() {
    return weight;
  }


  public float getLP() {
    return lp;
  }


  public void setLP( float l) {
    lp = l;
  }


  public Vertex getPredecessor() {
    return predecessor;
  }


  public void setPredecessor( Vertex pre ) {
    predecessor = pre;
  }


  public void incIndegree() {
    ++indegree;
  }


  public void decIndegree() {
    --indegree;
  }

  public boolean equals( Vertex v ) {
    boolean b = false;
    String startV = v.getStartLetter();
    String endV = v.getEndLetter();
    if ( startLetter.equals( startV ) && ( endLetter.equals( endV )))
      b = true;  
    return b;
  }


  public static void main( String args[] ) {
    String str = "89";
    Integer ig = new Integer( str );
    System.out.println( ig.toString());
    Vertex v1 = new Vertex( args[0], args[1], 0, 0 );
    Vertex v2 = new Vertex( "2", "34", 0, 0 );
    int i = 2 + Integer.parseInt ( "5" );
    System.out.println( "i = " + i );
    if ( 3<=3 )
      System.out.println( "OK" );
    
    if ( v1.equals( v2 ))
      System.out.println( "equal ");
    else 
      System.out.println( "NOT equal" );
  }
  
}

