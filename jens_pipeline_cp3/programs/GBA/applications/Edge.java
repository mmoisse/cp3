/*
** This class represents the class of edges of a graph
 */

package applications;

class Edge {
  
  private Vertex source;
  private Vertex sink;
  private float weight;

  public Edge() {
    source = new Vertex();
    sink = new Vertex();  
    weight = 0;
  }
   

  public Edge( Vertex so, Vertex si, float w) {
    source = so;
    sink = si;
    weight = w;
  }


  public Vertex getSource() {
    return source;
  }


  public Vertex getSink() {
    return sink;
  }

  
  public float getWeight() {
    return weight;  
  }

  public void setSource( Vertex so ) {
    source = so;
  }

  public void setSink( Vertex si ) {
    sink = si;
  }
  
}
