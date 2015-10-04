/* 
 ** Author: Xuehui Li
 ** Date: March, 2005
 ** "gbm" or "GBM" is the abbrevaition of "A graph-based method for detecting low-complexity regions".
 ** This is program is used to find low-complexity regions in sequences
 ** NOTE: all vertices and edges have topological orders
 ** There are totally six input parameters. The first one is the sequence file 
 ** name. The second one the learned matrix file (/cise/research/tamer/xli/LCR
 ** /graphLCR/swissprotLearnedMatrices). The nex three are the threshold 1, 
 ** threshold 2 and threshold 3, respectively. At this time, all LCR Blocks 
 ** (or masked sequences) generated in both /cise/research/tamer/xli/LCR/
 ** graphLCR/swissprotLCRBlocks/ (/cise/research/tamer/xli/LCR/graphLCR/
 ** swissprotMaskedSeqs) and /cise/research/tamer/xli/LCR/graphLCR/
 ** pfamLCRBlocks/ are based on the three thresholds: " 3 15 5". The last
 ** parameter is used to choose the output format. "0" means LCR blocks will 
 ** be generated. "1" means masked sequences will be generated.
 ** How to run the program? 
 ** One example (to generate LCR blocks instead of masked sequences) with fixed sampling,
 * i.e., 27 LCRs from five sequences masked by SEG:
 ** java applications.Gbm ../data/swissprot/sequenceInfor/seqFromFlybase 
 *swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095 3 15 5 
 *knowledge/lcrs 0 > swissprotLCRBlocks/wForgetRateMatrices/postProcess/
 *flybase/sim095flybase2LetModEntrNor2LetModEntrMer63CutAdjBlockSmiWatNCheCom
 ** if conditional-sampling (i.e., G_t / sqrt(t) <= C * u_t)is used, the "knowledge/lcrs" 
 *will not make sense any more. 
 */
package applications;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.formats.Pair;
import jaligner.matrix.MatrixLoader;
import jaligner.util.SequenceParser;

import java.io.*;
import java.util.*;

class Gbm {

	// vertices and edges are vectors used to keep all the vertices, edges in a
	// graph generated from a sequence, respectively. vertexQueue is a vector
	// used to keep all vertices whose indegree is zero. lps is a vector used to
	// keep all the longest paths in all connected subgraphs of a sequence.
	// Every longest path in lps is a vector of vertices, excluding the dummy
	// source.
	// subVertices and subEdges are vectors used to keep all the vertices, edges
	// in a connected-graph which is a subgraph of the graph generated from a
	// sequence, respectively.

	private static final double comCut = 3.7491631225622157; // the
	// complexity
	// cutOff value
	// for extending
	// the
	// longest-path
	// intervals.

	private File f;
	private RandomAccessFile rf;
	private Vector vertices, subVertices, edges, subEdges, vertexQueue, lps;
	private float[][] repeatMatrix, nonRepeatMatrix;
	private float[] fVecNor, fVecUnNor;
	private double[][] scoringMatrix;
	private Vector alphabet;

	// private double comCut = 0; // the complexity cut-off value
	// private ComplexityCalculator cc;

	public Gbm(String fileName) {
		// the graph is given in a file where every line represents an edge and
		// has the fromat of "source sink weight"
		initializeAlphabet();
		try {
			f = new File(fileName);
			rf = new RandomAccessFile(f, "r");
		} catch (IOException ex) {

		}
		vertices = new Vector();
		subVertices = new Vector();
		edges = new Vector();
		subEdges = new Vector();
		vertexQueue = new Vector();
		lps = new Vector();
		repeatMatrix = new float[20][20];
		nonRepeatMatrix = new float[20][20];
		fVecNor = new float[20];
		fVecUnNor = new float[20];
		scoringMatrix = new double[20][20];
		for (int i = 0; i < 20; i++)
			fVecUnNor[i] = 0f;
	}

	public void initializeAlphabet() {
		alphabet = new Vector();
		alphabet.add("A");
		alphabet.add("R");
		alphabet.add("N");
		alphabet.add("D");
		alphabet.add("C");
		alphabet.add("Q");
		alphabet.add("E");
		alphabet.add("G");
		alphabet.add("H");
		alphabet.add("I");
		alphabet.add("L");
		alphabet.add("K");
		alphabet.add("M");
		alphabet.add("F");
		alphabet.add("P");
		alphabet.add("S");
		alphabet.add("T");
		alphabet.add("W");
		alphabet.add("Y");
		alphabet.add("V");
	}

	// get repeat/non-repeat matrices
	public void readRNRMatrices(String matricesFile) {
		try {
			File f = new File(matricesFile);
			RandomAccessFile rfm = new RandomAccessFile(f, "r");
			String row = new String();
			for (int i = 0; i < 20; i++) {
				row = rfm.readLine();
				row = rfm.readLine();
				row = row.trim();
				for (int j = 0; j < 20; j++) {
					int index = row.indexOf("   ");
					if (index != -1) {
						nonRepeatMatrix[i][j] = Float.parseFloat(row.substring(
								0, index));
					} else {
						nonRepeatMatrix[i][j] = Float.parseFloat(row);
					}
					row = row.substring(index + 3);
				}
			}
			rfm.readLine();
			rfm.readLine();
			rfm.readLine();
			for (int i = 0; i < 20; i++) {
				row = rfm.readLine();
				row = rfm.readLine();
				row = row.trim();
				for (int j = 0; j < 20; j++) {
					int index = row.indexOf("   ");
					if (index != -1) {
						repeatMatrix[i][j] = Float.parseFloat(row.substring(0,
								index));
					} else {
						repeatMatrix[i][j] = Float.parseFloat(row);
					}
					row = row.substring(index + 3);
				}
			}
			rfm.close();
		} catch (IOException ex) {
		}
	}

	public void readScoringMatrix(String fileName) {
		// System.out.println(" From readScoringMatrix: " + fileName);
		String line = new String(), score = new String();
		try {
			File f = new File(fileName);
			RandomAccessFile rf = new RandomAccessFile(f, "r");
			// System.out.println("rf: " + rf);
			for (int i = 0; i < 20; i++) {
				line = rf.readLine();
				int k = 0;
				for (int j = 0; j < 20; j++) {
					score = line.substring(k, k + 2).trim();
					scoringMatrix[i][j] = Integer.parseInt(score);
					k = k + 3;
				}
			}
			rf.close();
		} catch (IOException ex) {
			System.out.println("Exception Happened from readScoringmatrix");
		}
	}

	public void printMatricesRowByRow() {
		System.out.println("Non-Repeat matrix: ");
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++)
				System.out.print(nonRepeatMatrix[i][j] + "   ");
			System.out.println();
			System.out.println();
		}
		System.out.println();
		System.out.println();
		System.out.println("Repeat matrix: ");
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++)
				System.out.print(repeatMatrix[i][j] + "   ");
			System.out.println();
			System.out.println();
		}
		System.out.println();
	}

	public void createFirstVector(String window) {
		int len = window.length();
		String tmpWindow = window, letter = new String();
		for (int i = 0; i < len; i++) {
			letter = tmpWindow.substring(0, 1);
			int index = alphabet.indexOf(letter);
			fVecUnNor[index] = fVecUnNor[index] + 1f;
			tmpWindow = tmpWindow.substring(1);
		}
		for (int i = 0; i < 20; i++)
			fVecNor[i] = fVecUnNor[i];
		for (int i = 0; i < 20; i++)
			fVecNor[i] = fVecNor[i] / len;
	}

	// o for unNormalized, 1 for normalized
	public void printVector(int mark) {
		float[] tmpVector = new float[20];
		if (mark == 0)
			tmpVector = fVecUnNor;
		else
			tmpVector = fVecNor;
		for (int i = 0; i < 20; i++)
			System.out.print(tmpVector[i] + "  ");
		System.out.println();

	}

	public void constructSingleVertex(int start, int end) {
		String startLetter = Integer.toString(start + 1);
		String endLetter = Integer.toString(end + 1);
		Vertex v = new Vertex(startLetter, endLetter, 1, 1.0f);
		Vertex dummySource = new Vertex("0", "0", 0, 1.0f);
		v.setPredecessor(dummySource);
		vertices.add(v);
	}

	public void workOnFirstWindow(String window) {
		String letter1, letter2;
		createFirstVector(window);
		int len = window.length(), j = 0, row = 0, col = 0;
		for (int i = 0; i < len - 1; i++) {
			j = 1;
			letter1 = window.substring(i, i + 1);
			while ((i + j) < len) {
				letter2 = window.substring(i + j, i + j + 1);
				row = alphabet.indexOf(letter1);
				col = alphabet.indexOf(letter2);
				if (scoringMatrix[row][col] > 1)
					if (checkProbablity(row, col))
						constructSingleVertex(i, i + j);
				++j;
			}
		}
	}

	// at this time, row == col, since we only consider same letters
	public boolean checkProbablity(int row, int col) {
		boolean construct = false;
		float difference1 = Math.abs(repeatMatrix[row][col] - fVecNor[row]);
		float difference2 = Math.abs(nonRepeatMatrix[row][col] - fVecNor[row]);
		if (difference2 > difference1) {
			construct = true;
		}
		return construct;
	}

	public void addVertices(String window, int index, int startPos) {
		// System.out.println("From addVertices: window: " + window + " index: "
		// + index + " startPos: " + startPos );
		int len = window.length();
		// System.out.println("Window length: " + len);
		String letter = new String();
		int tmpIndex = 0;
		for (int i = 0; i < len - 1; i++) {
			// System.out.println("Entered the for loop");
			letter = window.substring(i, i + 1);
			tmpIndex = alphabet.indexOf(letter);
			// System.out.println("index: " + index + "tmpIndex: " + tmpIndex +
			// " scoringMatrix: " + scoringMatrix[index][tmpIndex]);
			if (scoringMatrix[index][tmpIndex] > 1) {
				// System.out.println("Index: " + index + " tmpIndex: " +
				// tmpIndex);
				if (checkProbablity(index, tmpIndex)) {
					// System.out.println("Constructing single vertex");
					constructSingleVertex(startPos + i, startPos + len - 1);
				}
			}
		}
	}

	public void constructVertices(String str, int th4) {
		char c1, c2;
		int len = str.length();
		boolean similar = false;
		Vertex v = new Vertex();
		float[] previousVecUnNOr = new float[20];
		String window = str.substring(0, th4);
		String oldLetter = new String(), newLetter = new String();
		workOnFirstWindow(window);
		int startPos = 1;
		while (startPos <= (len - th4)) {
			// System.out.println("Adding vertices");
			oldLetter = window.substring(0, 1);
			window = str.substring(startPos, startPos + th4);
			newLetter = window.substring(th4 - 1, th4);
			int index = alphabet.indexOf(oldLetter);
			fVecUnNor[index] = fVecUnNor[index] - 1;
			fVecNor[index] = fVecUnNor[index] / th4;
			index = alphabet.indexOf(newLetter);
			fVecUnNor[index] = fVecUnNor[index] + 1;
			fVecNor[index] = fVecUnNor[index] / th4;
			addVertices(window, index, startPos);
			++startPos;
		}
		// System.out.println("Vertices size: " + vertices.size());
	}

	// return the actual weight of the vertex ???
	// / to be modified later to include similar cases
	public boolean checkSimilarity(String c1, String c2) {
		boolean similar = false;
		if (c1.equals(c2))
			similar = true;
		return similar;
	}

	// find the percentage of letters appearing in vertices
	public void findLetterPercentageVer(float lF) {
		int len = vertices.size();
		Vector appearedLetters = new Vector();
		Vertex v = new Vertex();
		for (int i = 0; i < len; i++) {
			v = (Vertex) vertices.elementAt(i);
			String str = v.getStartLetter();
			int index = appearedLetters.indexOf(str);
			if (index == -1)
				appearedLetters.add(str);
			str = v.getEndLetter();
			index = appearedLetters.indexOf(str);
			if (index == -1)
				appearedLetters.add(str);
		}
		// /// Sort first
		len = appearedLetters.size();
		float per = len / lF;
		// System.out.println( "The letter percentage after vertex construction
		// is: " + len + " / " + lF + " = " + per );
		Integer[] posInt = new Integer[len];
		for (int i = 0; i < len; i++) {
			String str = (String) appearedLetters.get(i);
			posInt[i] = new Integer(str);
		}
		MergeSort.mergeSort(posInt);
		for (int i = 0; i < len; i++) {
			String str = posInt[i].toString();
			System.out.print(str + "  ");
		}
		System.out.println();

		// computePercentage( lF, appearedLetters );
	}

	public void computePercentage(float l, Vector al) {
		int len = al.size();
		float per = len / l;
		System.out.println("The percentage is: " + len + " / " + l + " = "
				+ per);
		for (int i = 0; i < len; i++) {
			String str = (String) al.elementAt(i);
			System.out.print(str + "  ");
		}
		System.out.println();
	}

	public void constructEdges(int th1, int th2, int th3) {
		// System.out.println("From constructEdges: th1: " + th1 + " th2: " +
		// th2
		// + " th3: " + th3);
		String startLetter = new String(), endLetter = new String();
		int i = 0, j = 0, p = 0, q = 0, l = vertices.size();
		Vertex v1 = new Vertex(), v2 = new Vertex();
		Edge e = new Edge();
		// System.out.println("Before entering for loop");
		// System.out.println("Value of l: " + l);
		for (int k = 0; k < l - 1; k++) {
			v1 = (Vertex) vertices.get(k);
			startLetter = v1.getStartLetter();
			i = Integer.parseInt(startLetter);
			endLetter = v1.getEndLetter();
			j = Integer.parseInt(endLetter);
			boolean end = false;
			int m = k + 1;
			// System.out.println("Before adding edges");
			while ((!end) && (m < l)) {
				v2 = (Vertex) vertices.get(m);
				startLetter = v2.getStartLetter();
				p = Integer.parseInt(startLetter);
				endLetter = v2.getEndLetter();
				q = Integer.parseInt(endLetter);

				if (checkConditions(i, j, p, q, th1, th2, th3)) {
					// System.out.println("Adding edges: ");

					v2.incIndegree();
					// modify the vertex in vertices
					vertices.setElementAt(v2, m);
					e = new Edge(v1, v2, 1.0f);
					edges.add(e);
				} else if ((p - i) > th2) // k2 = 25
					end = false;
				++m;
			}
			// System.out.println("After adding edges");
		}
		// System.out.println("After entering for loop");
	}

	// return the actual weight of the edge ???
	public boolean checkConditions(int i, int j, int p, int q, int th1,
			int th2, int th3) {
		boolean satisfied = false;
		int foo = (j - i) - (q - p);
		foo = Math.abs(foo);
		if (foo <= th1) // condition #1, k1 = 5
			if ((p - i) <= th2) // condition #2, k2 = 26
				if (((i <= p) && (p <= j)) && (j <= q)) // condition #3
					if ((i == p) || (p == j) || (j == q)) { // condition #4
						if (((j - i) <= th3) && ((q - p) <= th3))
							satisfied = true;
					} else
						satisfied = true;
		return satisfied;
	}

	public void modifyVertexQueue(Vector tmpQueue, Vector tmpVertices) {
		Vertex dummySource = new Vertex("0", "0", 0, 1.0f);
		vertexQueue.clear();
		vertexQueue.add(dummySource);
	}

	public void printVertexQueue() {
		System.out.println("All vertices in the queue:");
		Vertex v = new Vertex();
		Vertex previous = v;
		int l = vertexQueue.size();
		for (int i = 0; i < l; i++) {
			v = (Vertex) vertexQueue.get(i);
			previous = v.getPredecessor();
			if (previous != null)
				System.out.println(v.getStartLetter() + "A" + v.getEndLetter()
						+ "  " + v.getIndegree() + " " + v.getWeight() + " "
						+ v.getLP() + " " + previous.getStartLetter() + "A"
						+ previous.getEndLetter());
			else
				System.out.println(v.getStartLetter() + "A" + v.getEndLetter()
						+ "  " + v.getIndegree() + " " + v.getWeight() + " "
						+ v.getLP());
		}
	}

	// m = 0 for subgraph, m = 1 for the whole sequence
	public void printVertices(int m) {
		Vector tmpVertices = new Vector();
		if (m == 0)
			tmpVertices = subVertices;
		else
			tmpVertices = vertices;
		Vertex v = new Vertex();
		Vertex previous = v;
		int l = tmpVertices.size();
		System.out.println("All vertices in the graph: " + l);
		for (int i = 0; i < l; i++) {
			v = (Vertex) tmpVertices.get(i);
			previous = v.getPredecessor();
			if (previous != null)
				System.out.println(v.getStartLetter() + "A" + v.getEndLetter()
						+ "  " + v.getIndegree() + " " + v.getWeight() + " "
						+ v.getLP() + " " + previous.getStartLetter() + "A"
						+ previous.getEndLetter());
			else
				System.out.println(v.getStartLetter() + "A" + v.getEndLetter()
						+ "  " + v.getIndegree() + " " + v.getWeight() + " "
						+ v.getLP());
		}
	}

	// m = 0 for subgraph, m = 1 for the whole sequence
	public void printEdges(int m) {
		Vector tmpEdges = new Vector();
		if (m == 0)
			tmpEdges = subEdges;
		else
			tmpEdges = edges;
		Edge e = new Edge();
		Vertex sourceVer = new Vertex(), sinkVer = new Vertex();
		String str = new String();
		int length = tmpEdges.size();
		System.out.println("All edges in the graph: " + length);
		for (int i = 0; i < length; i++) {
			e = (Edge) tmpEdges.get(i);
			sourceVer = e.getSource();
			sinkVer = e.getSink();
			System.out.println(sourceVer.getStartLetter() + "A"
					+ sourceVer.getEndLetter() + "   lp: " + sourceVer.getLP()
					+ "   indegree: " + sourceVer.getIndegree() + "     "
					+ sinkVer.getStartLetter() + "A" + sinkVer.getEndLetter()
					+ "   lp: " + sinkVer.getLP() + "   indegree: "
					+ sinkVer.getIndegree() + "    weight: " + e.getWeight());
		}
	}

	public void addDummySource() {
		Vertex dummySource = new Vertex("0", "0", 0, 1.0f);
		addToSubVerticesVertexQueue(dummySource);
		addToSubEdges(dummySource);
	}

	public void addToSubVerticesVertexQueue(Vertex dummySource) {
		subVertices.add(0, dummySource);
		vertexQueue.clear();
		vertexQueue.add(dummySource);
	}

	public void addToSubEdges(Vertex dummySource) {
		Vertex v = new Vertex();
		Edge e = new Edge();
		int l = subVertices.size();
		for (int i = l - 1; i > 0; i--) {
			v = (Vertex) subVertices.get(i);
			e = new Edge(dummySource, v, 0.0f);
			subEdges.add(0, e);
		}
	}

	public Vector findLongestPath() {
		Vertex v = new Vertex();
		while (!(vertexQueue.isEmpty())) {
			v = (Vertex) vertexQueue.remove(0);
			traverseSubEdges(v);
		}
		// find the vertex to which the path from the source is the longest
		Vector lp = traverseVertices();
		return lp;
	}

	public void traverseSubEdges(Vertex ver) {
		boolean end = false, first = true;
		Edge e = new Edge();
		Vertex sourceVer = new Vertex();
		Vertex sinkVer = new Vertex();
		String str = new String();
		float w = 0, sourceLP = 0, sinkLP = 0;
		int len = subEdges.size(), j = 0;
		if (!(subEdges.isEmpty()))
			j = findEdges(ver);
		if (j == -1)
			end = true;
		else if (j != 0) {
		}
		while ((!end) && (j < len) && (!(subEdges.isEmpty()))) {
			e = (Edge) subEdges.elementAt(j);
			sourceVer = e.getSource();
			if (sourceVer.equals(ver)) {
				sourceLP = sourceVer.getLP();
				sinkVer = e.getSink();
				int i = subVertices.indexOf(sinkVer);
				sinkLP = sinkVer.getLP();
				w = e.getWeight();
				if ((sourceLP + w) > sinkLP) {
					sinkLP = sourceLP + w;
					sinkVer.setLP(sinkLP);
					sinkVer.setPredecessor(sourceVer);
				}
				sinkVer.decIndegree();
				subVertices.setElementAt(sinkVer, i);
				modifyVertexInEdges(sinkVer);
				i = sinkVer.getIndegree();
				if (i == 0)
					vertexQueue.add(sinkVer);
				subEdges.remove(j);
				len = subEdges.size();
				first = false;
			} else if (first) {
				System.out.println("This is a vertex with outdegree zero");
				end = true;
			} else
				end = true;
		}
	}

	public int findEdges(Vertex ver) {
		int i = 0, l = subEdges.size();
		boolean find = false;
		Vertex v = new Vertex();
		Edge e = new Edge();
		while ((!find) && (i < l)) {
			e = (Edge) subEdges.get(i);
			v = e.getSource();
			if (v.equals(ver))
				find = true;
			else
				++i;
		}
		if (!find)
			i = -1;
		return i;
	}

	public void modifyVertexInEdges(Vertex sinkVer) {
		Vertex v = new Vertex();
		Edge e = new Edge();
		int l = subEdges.size();
		for (int i = 0; i < l; i++) {
			e = (Edge) subEdges.get(i);
			v = e.getSource();
			if (v.equals(sinkVer))
				e.setSource(sinkVer);
			else {
				v = e.getSink();
				if (v.equals(sinkVer))
					e.setSink(sinkVer);
			}
			subEdges.setElementAt(e, i);
		}
	}

	public Vector traverseVertices() {
		Vertex v = new Vertex(), maxVer = new Vertex();
		int l = subVertices.size();
		float length = 0f, maxLp = -2222.0f;
		for (int i = 0; i < l; i++) {
			v = (Vertex) subVertices.get(i);
			length = v.getLP();
			if (length > maxLp) {
				maxLp = length;
				maxVer = v;
			}
		}
		Vector lp = constructLongestPath(maxVer);
		return lp;
	}

	public Vector constructLongestPath(Vertex maxVer) {
		Vector lp = new Vector();
		Vertex dummySource = new Vertex("0", "0", 0, 1.0f);
		Vertex v = maxVer;
		while (!(v.equals(dummySource))) {
			lp.add(0, v);
			v = v.getPredecessor();
		}
		return lp;
	}

	public void printLongestPath(Vector lp) {
		int i = 0;
		Vertex v = new Vertex();
		String str = new String();
		int length = lp.size();
		for (i = 0; i < length; i++) {
			v = (Vertex) lp.get(i);
			System.out
					.print(v.getStartLetter() + "A" + v.getEndLetter() + "  ");
		}
		System.out.println();
	}

	public Vector identifyLCRs() {
		addDummySource();
		Vector lp = findLongestPath();
		return lp;
	}

	public boolean checkExistence(Vertex v) {
		int index = subVertices.indexOf(v);
		if (index == -1)
			return false;
		else
			return true;
	}

	public void copy(Vector vt1, Vector vt2) {
		int l = vt1.size();
		Edge e = new Edge();
		for (int i = 0; i < l; i++) {
			e = (Edge) vt1.get(i);
			vt2.add(e);
		}
	}

	// make all edges beginning with the same vertex stay together
	public void clusterSubEdges() {
		Vector tmpSubEdges = new Vector();
		Edge e = new Edge();
		Vertex v = new Vertex(), ver = new Vertex();
		while ((!subEdges.isEmpty())) {
			e = (Edge) subEdges.remove(0);
			tmpSubEdges.add(e);
			v = e.getSource();
			int m = 0;
			int len = subEdges.size();
			while ((m < len) && (!(subEdges.isEmpty()))) {
				e = (Edge) subEdges.elementAt(m);
				ver = e.getSource();
				if (v.equals(ver)) {
					tmpSubEdges.add(e);
					subEdges.remove(m);
				} else
					++m;
				len = subEdges.size();
			}
		}
		copy(tmpSubEdges, subEdges);
	}

	// assign values to subVertices and subEdges ( BFS )
	public void extractConnectedGraph() {
		boolean first = true;
		Vertex v = new Vertex();
		Edge e = new Edge();
		Vector tmpQueue = new Vector();
		while ((first) || (!(tmpQueue.isEmpty()))) {
			if (first) { // start the first edge of a new connected subgraph
				e = (Edge) edges.remove(0);
				subEdges.add(e);
				v = e.getSource();
				subVertices.add(v);
				if (!(vertices.remove(v)))
					System.out.println("wrong1");
				v = e.getSink();
				subVertices.add(v);
				tmpQueue.add(v);
				if (!(vertices.remove(v)))
					System.out.println("wrong2");
				v = e.getSource();
				first = false;
				boolean same = true;
				int m = 0;
				int len = edges.size();
				while ((same) && (m < len)) { // remove all those edges having
					// the same source vertex as the
					// first edge
					e = (Edge) edges.elementAt(m);
					Vertex ver = e.getSource();
					if (v.equals(ver)) {
						edges.remove(m); // remove the edge who starts with v
						subEdges.add(e);
						ver = e.getSink();
						subVertices.add(ver);
						tmpQueue.add(ver); // put ver ( the sink of the edge )
						// into tmpQueue;
						if (!(vertices.remove(ver)))
							System.out.println("wrong3");
					} else
						same = false;
					len = edges.size();
				}
			} else {
				v = (Vertex) tmpQueue.remove(0);
				int m = 0;
				boolean found = false;
				int len = edges.size();
				// find the starting positon of those edges who start at the
				// first vertex from tmpQueue
				while ((!found) && (m < len)) { // skip all edges starting with
					// the vertex from the tmpQueue
					e = (Edge) edges.elementAt(m);
					Vertex sr = e.getSource();
					if (v.equals(sr))
						found = true;
					else {
						Vertex si = e.getSink();
						if (v.equals(si)) {
							boolean exist = checkExistence(sr);
							if (!exist) {
								subVertices.add(sr);
								tmpQueue.add(sr); // work on edges whose sink
								// vertex is the same as the
								// vertex from tempQueue
								if (!(vertices.remove(sr)))
									System.out.println("wrong4");
							}
							e = (Edge) edges.remove(m);
							subEdges.add(e);
							len = edges.size();

						} else
							++m;
					}
				}
				// System.out.println( "m = " + m );
				boolean same = true;
				while (same) { // remove all those edges starting with the
					// vertex from tmpQueue to subEdges
					len = edges.size();
					if (m < len) {
						e = (Edge) edges.elementAt(m);
						Vertex ver = e.getSource();
						if (v.equals(ver)) {
							edges.remove(m); // remove the edge who starts
							// with v
							subEdges.add(e);
							ver = e.getSink();
							boolean exist = checkExistence(ver); // check
							// whether
							// the
							// sinkVertex
							// is in the
							// subVertices
							// or not
							if (!exist) {
								subVertices.add(ver);
								tmpQueue.add(ver); // put ver ( the sink of the
								// edge ) into tmpQueue;
								if (!(vertices.remove(ver)))
									System.out.println("wrong5");
							}
						} else {
							same = false;
						}
					} else
						same = false;
				}
			}
		}
		clusterSubEdges();
	}

	public int workOnSequence(String str, int th1, int th2, int th3, int th4) {
		constructVertices(str, th4);
		// printVertices( 1 );//////////////////////////////
		int i = 0;
		constructEdges(th1, th2, th3);
		// printEdges(1);/////////////////////////
		boolean find = true;
		while (!(edges.isEmpty())) {
			// System.out.println("Enter the while loop, i:" + i);
			subVertices.clear();
			subEdges.clear();
			extractConnectedGraph();
			Vector lp = identifyLCRs(); // the longest path in a connected
			// subgraph
			// System.out.println("longest path"); /////////////////
			// printLongestPath(lp);//////////////
			lps.add(lp);
			++i;
		}
		// System.out.println("Within workOnSequence, i: " + i);
		return i;
	}

	// combine all letters from a sequence on different lines( stings) into a
	// single line ( string )
	public String generateSequence(String str) {
		String sequence = new String();
		String strTmp = str;
		boolean lastSeq = false;
		try {
			if (strTmp == null)
				sequence = null;
			if ((strTmp != null) && (strTmp.startsWith(">"))) {
				// System.out.println(
				// "*******************************************" );
				System.out.println(str);
				strTmp = rf.readLine();
			}
			while ((strTmp != null) && (!(strTmp.startsWith(">")))) {
				strTmp = strTmp.trim();
				sequence = sequence + strTmp;
				strTmp = rf.readLine();
			}
			if (strTmp != null)
				if (strTmp.startsWith(">")) {
					sequence = strTmp + "!" + sequence;
				}
		} catch (IOException ex) {
		}

		return sequence;
	}

	public void printPositions(Vector pos) {
		int l = pos.size();
		String str = new String();
		for (int i = 0; i < l; i++) {
			str = (String) pos.get(i);
			System.out.print(str + "  ");
		}
		System.out.println();
	}

	public void printLCRBlocks(Vector LCRBlocks) {
		String str = new String();
		int l = LCRBlocks.size();
		// System.out.println( "LCR Blocks: " );
		for (int i = 0; i < l; i++) {
			str = (String) LCRBlocks.get(i);
			int index = str.indexOf("-");
			String start = str.substring(0, index);
			String end = str.substring(index + 1);
			;
			int difference = Integer.parseInt(end) - Integer.parseInt(start);
			if (difference > 1)
				System.out.print(str + " ");
		}
		System.out.println();
	}

	public Vector getPositions(int k) {
		Vector pos = new Vector();
		// Vector posSingleVertexOnly = new Vector();
		String str = new String();
		int l = lps.size(), index = 0;
		Vector lp = new Vector();
		Vertex v = new Vertex();
		// longest path
		for (int i = 0; i < k; i++) { // get positions from those vertices in
			// lps
			lp = (Vector) lps.get(i);
			int len = lp.size();
			for (int j = 0; j < len; j++) {
				v = (Vertex) lp.get(j);
				str = v.getStartLetter();
				index = pos.indexOf(str);
				if (index == -1)
					pos.add(str);
				str = v.getEndLetter();
				index = pos.indexOf(str);
				if (index == -1)
					pos.add(str);
			}
		}
		return pos;
	}

	public Vector sortPositions(Vector pos) {
		Vector positions = new Vector(); // used to keep the sorted positions
		String str = new String();
		int len = pos.size(), current = 0, previous = 0;
		Integer[] posInt = new Integer[len];
		for (int i = 0; i < len; i++) {
			str = (String) pos.get(i);
			posInt[i] = new Integer(str);
		}
		MergeSort.mergeSort(posInt);
		Vector tmpLCRBlocks = new Vector();
		String start = new String();
		previous = posInt[0].intValue() - 1;
		start = Integer.toString((previous + 1));
		for (int i = 0; i < len; i++) {
			str = posInt[i].toString();
			current = posInt[i].intValue();
			// generate blocks of continuous positions. Say, the sorted integer
			// array is 3,4,5,6 8,9,10,11,12,13,29,30,31. It can be represented
			// as a vector of three strings( blocks ): 3-6, 8-13, 29-31.
			if (current != (previous + 1)) {
				tmpLCRBlocks.add(start + "-" + Integer.toString(previous));
				start = str;
			}
			previous = current;
			positions.add(str);
		}
		tmpLCRBlocks.add(start + "-" + Integer.toString(previous));
		len = tmpLCRBlocks.size();
		Vector LCRBlocks = new Vector();
		for (int i = 0; i < len; i++) {
			str = (String) tmpLCRBlocks.get(i);
			int index = str.indexOf("-");
			start = str.substring(0, index);
			String end = str.substring(index + 1);
			;
			int difference = Integer.parseInt(end) - Integer.parseInt(start);
			if (difference > 1) // the interval is at least 3-letter long
				LCRBlocks.add(str);
		}
		// printLCRs( LCRBlocks );
		return LCRBlocks;
	}

	public Vector extend(int startPos, int endPos, int limit, String direction,
			String seq) {
		Vector decRegs = new Vector();
		int pointer = 0, startDecPos = 0, endDecPos = 0;
		double com1 = 0, com2 = 0;
		String extReg = seq.substring(startPos - 1, endPos);
		ComplexityCalculator cc = new ComplexityCalculator();
		cc.initializeAlphabet();
		cc.initializeNewAlphabet();
		cc.computeNewScoringMatrix();
		cc.normalizeNewScoringMatrixPow();
		if (direction.equals("left")) { // extend to the left( front )
			boolean dec = false;
			pointer = startPos - 2;
			while ((pointer > limit) && (pointer > (startPos - 17))) {
				// System.out.println( "111111111111extReg:" + extReg );
				com1 = cc.calculate2LetterEntropyWScoMatrix(extReg);
				// com1 = cc.calculateModifiedEntropy( extReg );
				extReg = seq.substring(pointer, endPos);
				// com2 = cc.calculateModifiedEntropy( extReg );
				com2 = cc.calculate2LetterEntropyWScoMatrix(extReg);
				// com2 = cc.calculateReciprocalPro( extReg );
				// com2 = cc.calculateRecProWScoringMatrix( extReg );
				// System.out.println( com1 + " 11111 " + com2);
				if (com1 > com2) {
					if (!dec) {
						dec = true;
						// System.out.println( "from false to true111111111" );
						endDecPos = pointer + 2;
					}
				} else if (dec) {
					dec = false;
					// System.out.println( "from true to false1111111111" );
					startDecPos = pointer + 2;
					if (com1 < comCut) {
						decRegs.add(0, (startDecPos) + "-" + (endDecPos));
						// System.out.println( "decRegs added 11111111 " +
						// startDecPos + "-" + endDecPos );
					}
				}
				--pointer;
			}
			if ((dec) && (pointer == (startPos - 17))) {
				// System.out.println( "keeping decreasing1111111111" );
				while ((pointer > limit) && (dec)) {
					// System.out.println( "22222222222222extReg:" + extReg );
					// com1 = cc.calculateModifiedEntropy( extReg );
					com1 = cc.calculate2LetterEntropyWScoMatrix(extReg);
					// com1 = cc.calculateReciprocalPro( extReg );
					// com1 = cc.calculateRecProWScoringMatrix( extReg );
					extReg = seq.substring(pointer, endPos);
					// System.out.println("keep: " + extReg );
					// com2 = cc.calculateModifiedEntropy( extReg );
					com2 = cc.calculate2LetterEntropyWScoMatrix(extReg);
					// com2 = cc.calculateReciprocalPro( extReg );
					// com2 = cc.calculateRecProWScoringMatrix( extReg );
					// System.out.println( com1 + " 22222222222 " + com2);
					if (com1 < com2) {
						// System.out.println( "from true to false2222222222" );
						dec = false;
						startDecPos = pointer + 2;
						if (com1 < comCut) {
							decRegs.add(0, (startDecPos) + "-" + (endDecPos));
							// System.out.println( "decRegs added 222222222222 "
							// + startDecPos + "-" + endDecPos );
						}
					}
					--pointer;
				}
			}
			// the left extension touches the end of the last block of the
			// current lcr blocks
			if ((pointer == limit) && (dec)) {
				startDecPos = pointer + 2;
				if (com2 < comCut) {
					decRegs.add(0, (startDecPos) + "-" + (endDecPos));
					// System.out.println( "decRegs added 333333333333 " +
					// startDecPos + "-" + endDecPos );
				}
			}
			if (decRegs.size() == 0) {
				// System.out.println( "left: Empty" );
			} else {
				// System.out.print( "left: ");
				// printLCRs( decRegs );
			}
		} else {
			boolean dec = false;
			pointer = endPos + 1;// //////
			// extend to the right( back )
			while ((pointer < limit) && (pointer < (endPos + 15))) {
				// System.out.println( "333333333333333333333extReg:" + extReg
				// );
				// com1 = cc.calculateModifiedEntropy( extReg );
				com1 = cc.calculate2LetterEntropyWScoMatrix(extReg);
				// com1 = cc.calculateReciprocalPro( extReg );
				// com1 = cc.calculateRecProWScoringMatrix( extReg );
				extReg = seq.substring(startPos - 1, pointer);
				// com2 = cc.calculateModifiedEntropy( extReg );
				com2 = cc.calculate2LetterEntropyWScoMatrix(extReg);
				// com2 = cc.calculateReciprocalPro( extReg );
				// com2 = cc.calculateRecProWScoringMatrix( extReg );
				// System.out.println( com1 + " 33333333 " + com2);
				if (com1 > com2) {
					if (!dec) {
						dec = true;
						// System.out.println( "from false to
						// true33333333333333333" );
						startDecPos = pointer - 1;
					}
				} else if (dec) {
					dec = false;
					// System.out.println( "from true to false3333333333333" );
					endDecPos = pointer - 1;
					if (com1 < comCut) {
						decRegs.add((startDecPos) + "-" + (endDecPos));
						// System.out.println( "decRegs added 4444444 " +
						// startDecPos + "-" + endDecPos );
					}
				}
				++pointer;
			}
			if ((dec) && (pointer == (endPos + 15))) {
				// keep extending until the complexity starts increasing, which
				// means that several blocks generated from the longest path can
				// be included into lcrs during one call of the 'extend()' based
				// on a block
				while ((dec) && (pointer < limit)) {
					// System.out.println( "444444444444extReg:" + extReg );
					com1 = cc.calculate2LetterEntropyWScoMatrix(extReg);
					// com1 = cc.calculateModifiedEntropy( extReg );
					// com1 = cc.calculateReciprocalPro( extReg );
					// com1 = cc.calculateRecProWScoringMatrix( extReg );
					extReg = seq.substring(startPos - 1, pointer);
					// com2 = cc.calculateModifiedEntropy( extReg );
					com2 = cc.calculate2LetterEntropyWScoMatrix(extReg);
					// com2 = cc.calculateReciprocalPro( extReg );
					// com2 = cc.calculateRecProWScoringMatrix( extReg );
					// System.out.println( com1 + " 4444444 " + com2);
					if (com1 < com2) {
						dec = false;
						// System.out.println( "from true to false444444444" );
						endDecPos = pointer - 1;
						if (com1 < comCut) {
							decRegs.add((startDecPos) + "-" + endDecPos);
							// System.out.println( "decRegs added
							// 5555555555555555555 " + startDecPos + "-" +
							// endDecPos );
						}
					}
					++pointer;
				}
				if ((pointer == limit) && (dec)) {
					endDecPos = limit - 1;
					// System.out.println( "decRegs added 66666666666 " +
					// startDecPos + "-" + endDecPos );
					decRegs.add(0, (startDecPos) + "-" + (endDecPos));
				}
			}
			if (decRegs.size() == 0) {
				// System.out.println( "right: Empty" );
			} else {
				// System.out.print( "right: ");
				// printLCRs( decRegs );
			}
		}
		return decRegs;
	}

	public boolean shareLetter(String str1, String str2) {
		boolean shared = false;
		String str = str1, letter = new String();
		while ((str.length() != 0) && (!shared)) {
			letter = str.substring(0, 1);
			int index = str2.indexOf(letter);
			if (index != -1)
				shared = true;
			else if (str.length() != 0)
				str = str.substring(1);
		}
		return shared;
	}

	public boolean checkContribution(String currentBlock, Vector decRegs,
			String seq) {
		boolean contributed = false;
		Vector regs = decRegs;
		String block = new String();
		int i = 0, len = regs.size();
		while ((i < len) && (!contributed)) {
			block = (String) regs.elementAt(i);
			int index = block.indexOf("-");
			int start = Integer.parseInt(block.substring(0, index));
			int end = Integer.parseInt(block.substring(index + 1));
			block = seq.substring(start - 1, end);
			// System.out.println( "block: " + block+ " currentBlock: "+
			// currentBlock );
			contributed = shareLetter(currentBlock, block);
			++i;
		}
		return contributed;
	}

	public Vector appendLcrs(Vector lcrs, Vector appendedLcrs) {
		Vector lowComRegs = lcrs, tmpLcrs = appendedLcrs;
		while (!(tmpLcrs.isEmpty()))
			lowComRegs.add((String) tmpLcrs.remove(0));
		return lowComRegs;
	}

	public Vector pickUpDrop(Vector blocks, String seq) {
		Vector frontLcrs = new Vector(), backLcrs = new Vector(), lcrs = new Vector(), tmpBLOCKS = blocks;
		String currentBlock = new String(), tmpBlock = new String();
		boolean isFirstBlock = true;
		int limit = 0, index = 0, startPos = 0, endPos = 0;
		ComplexityCalculator cc = new ComplexityCalculator();
		cc.initializeAlphabet();
		cc.initializeNewAlphabet();
		cc.computeNewScoringMatrix();
		cc.normalizeNewScoringMatrixPow();
		while ((!tmpBLOCKS.isEmpty())) {
			frontLcrs.clear();
			backLcrs.clear();
			int lcrBlockStart = 0, lcrBlockEnd = 0;
			boolean extendToLeft = true; // whether to extend towards the
			// left
			boolean find = false; // find the current extending block
			// currentBlock can start in the middle of a block, or has the same
			// starting position as a block and it doesn't have to be the block
			// after the previous currentBlock
			// get the end position of the last block in lcrs
			if (!(lcrs.isEmpty())) {
				tmpBlock = (String) lcrs.lastElement();
				index = tmpBlock.indexOf("-");
				lcrBlockEnd = Integer.parseInt(tmpBlock.substring(index + 1));
				// System.out.println( "lcrBlockEnd: " + lcrBlockEnd );
				while ((!find) && (!(tmpBLOCKS.isEmpty()))) { // get the
					// current block
					currentBlock = (String) tmpBLOCKS.remove(0);
					index = currentBlock.indexOf("-");
					startPos = Integer.parseInt(currentBlock
							.substring(0, index));
					endPos = Integer
							.parseInt(currentBlock.substring(index + 1));
					// System.out.println( "find currentBlock: "+ startPos + " "
					// + endPos);
					if (startPos < lcrBlockEnd) {
						if (endPos > lcrBlockEnd)
							if ((endPos - lcrBlockEnd) >= 3) {
								startPos = lcrBlockEnd + 1;
								extendToLeft = false;
								find = true;
							}
					} else
						find = true;
				}
			} else {
				currentBlock = (String) tmpBLOCKS.remove(0);
				index = currentBlock.indexOf("-");
				startPos = Integer.parseInt(currentBlock.substring(0, index));
				endPos = Integer.parseInt(currentBlock.substring(index + 1));
				find = true;
			}
			if (find) {
				// System.out.println( "currentBlock:" + currentBlock );
				if (isFirstBlock) {
					limit = -1;
					isFirstBlock = false;
					// extend to the left( front )
					frontLcrs = extend(startPos, endPos, limit, "left", seq);
				} else if (extendToLeft) {
					limit = lcrBlockEnd - 1;
					// extend to the left( front )
					frontLcrs = extend(startPos, endPos, limit, "left", seq);
				}
				limit = seq.length() + 1;
				// extend to the right( back )
				backLcrs = extend(startPos, endPos, limit, "right", seq);
				double com = 0;
				index = currentBlock.indexOf("-");
				int cbStart = Integer
						.parseInt(currentBlock.substring(0, index)) - 1;
				int cbEnd = Integer.parseInt(currentBlock.substring(index + 1));
				// System.out.println("current block String:" + seq.substring(
				// cbStart,cbEnd));
				// com = cc.calculateModifiedEntropy( seq.substring( cbStart,
				// cbEnd ) );
				com = cc.calculate2LetterEntropyWScoMatrix(seq.substring(
						cbStart, cbEnd));
				// com = cc.calculateModifiedEntropy( seq.substring( cbStart,
				// cbEnd ) );
				// com = cc.calculateReciprocalPro( seq.substring( cbStart,
				// cbEnd ) );
				// com = cc.calculateRecProWScoringMatrix( seq.substring(
				// cbStart, cbEnd ) );
				boolean contributed = false;
				if (frontLcrs.size() != 0) {
					// get the start position of the first block in frontLcrs as
					// the start position of the block to be added into lcrs
					tmpBlock = (String) frontLcrs.elementAt(0);
					index = tmpBlock.indexOf("-");
					lcrBlockStart = Integer.parseInt(tmpBlock.substring(0,
							index));
					if (com > comCut) {
						// check whether the current block contributes to the
						// complexity-decreasing regions or not
						contributed = checkContribution(seq.substring(cbStart,
								cbEnd), frontLcrs, seq);
						if (!contributed)
							lcrBlockEnd = startPos - 1;
						else {
							lcrBlockEnd = endPos;
						}
					} else {
						lcrBlockEnd = endPos;
						// System.out.println( "com of currentBlock: " + com );
					}
					lcrs.add(lcrBlockStart + "-" + lcrBlockEnd);
				}

				boolean combine = false; // whether to combine the last block
				// in lcrs from frontLcrs and the
				// block to be added into lcrs from
				// backLcrs
				if ((!contributed) && (com > comCut)) {
					contributed = checkContribution(seq.substring(cbStart,
							cbEnd), backLcrs, seq);
					if (!contributed) {
						lcrBlockStart = endPos + 1;
					} else {
						if (frontLcrs.size() != 0)
							combine = true;
					}
				} else if (frontLcrs.size() != 0)
					combine = true;
				// get the end position of the last block in backLcrs as the end
				// position of the block to be added into lcrs
				if (!(backLcrs.isEmpty())) {
					tmpBlock = (String) backLcrs.lastElement();
					index = tmpBlock.indexOf("-");
					lcrBlockEnd = Integer.parseInt(tmpBlock
							.substring(index + 1));
					if (combine) {
						// System.out.println( "combine" );
						limit = lcrs.size();
						tmpBlock = (String) lcrs.remove(limit - 1);
						index = tmpBlock.indexOf("-");
						lcrBlockStart = Integer.parseInt(tmpBlock.substring(0,
								index));
					} else {
						if (com < comCut) {
							lcrBlockStart = startPos;
							// System.out.println( "Here, com" );
						} else if (contributed) {
							// System.out.println( "contributed to the back, com
							// > comCut " );
							lcrBlockStart = startPos;
						} else {
							// System.out.println( "OOOOOOOOOOOOOOOOOOOOOOOO" );
							lcrBlockStart = endPos + 1;
						}
					}
					lcrs.add(lcrBlockStart + "-" + lcrBlockEnd);
				} else {
					if ((frontLcrs.size() == 0) && (!contributed)
							&& (com < comCut)) {
						lcrs.add(currentBlock);
					}
				}
				// check whether to combine the last two blocks in the current
				// lcrs
				// len = lcrs.length();
				/*
				 * System.out.print( "current lcrs: " ); printLCRs( lcrs );
				 */
			}
		}
		return lcrs;
	}

	public Vector mergePurge(Vector lcrs) {
		Vector tmpLcrs = lcrs;
		String currentBlock = new String(), nextBlock = new String();
		int len = tmpLcrs.size(), i = 0;
		while (i < len) {
			if ((i + 1) < len) {
				currentBlock = (String) tmpLcrs.elementAt(i);
				int endIndex = currentBlock.indexOf("-");
				int end = Integer
						.parseInt(currentBlock.substring(endIndex + 1));
				nextBlock = (String) tmpLcrs.elementAt(i + 1);
				int startIndex = nextBlock.indexOf("-");
				int start = Integer
						.parseInt(nextBlock.substring(0, startIndex));
				if ((end == (start - 1)) || (end == start)) {
					// System.out.println( currentBlock + " " + nextBlock );
					currentBlock = currentBlock.substring(0, endIndex) + "-"
							+ nextBlock.substring(startIndex + 1);
					tmpLcrs.remove(i);
					tmpLcrs.remove(i);
					tmpLcrs.add(i, currentBlock);
				} else
					++i;
				len = tmpLcrs.size();
			} else
				++i;
		}
		i = 0;
		len = tmpLcrs.size();
		// System.out.println("After the extension:"); ////////////////////////
		// printLCRs( tmpLcrs," ", "0" ); //////////////////
		/*
		 * while ( i < len ) { currentBlock = (String) tmpLcrs.elementAt( i );
		 * int index = currentBlock.indexOf( "-" ); int start =
		 * Integer.parseInt( currentBlock.substring( 0, index )); int end =
		 * Integer.parseInt( currentBlock.substring( index + 1 )); if (( end -
		 * start ) < 7 ) tmpLcrs.remove( i ); else ++i; len = tmpLcrs.size(); }
		 */
		return tmpLcrs;
	}

	public boolean checkCombinedSubBlock(String seq1, String seq2, double cCut) {
		boolean delete = true;
		ComplexityCalculator cc = new ComplexityCalculator();
		cc.initializeAlphabet();
		cc.initializeNewAlphabet();
		cc.computeNewScoringMatrix();
		cc.normalizeNewScoringMatrixPow();
		// double com = cc.calculateReciprocalPro( seq1 + seq2 );
		double com = cc.calculateNor2LetterEntropyWScoMatrix(seq1 + seq2);
		// double com = cc.calculateNorModifiedEntropy( seq1 + seq2 );
		// System.out.println( "combined:" + seq1 + seq2 + " " + com );
		if (com > cCut)
			delete = false;
		return delete;
	}

	public String findAlignment(String seq1, String seq2) {
		String aliPos = new String();
		try {
			Sequence s1 = SequenceParser.parse(seq1);
			Sequence s2 = SequenceParser.parse(seq2);
			// System.out.println( "alignment sequences: " + seq1 + "???" + seq2
			// );
			Alignment alignment = SmithWatermanGotoh.align(s1, s2, MatrixLoader
					.load("BLOSUM62"), 10f, 0.5f);
			int similarLen = alignment.getSimilarity(); // get the length of the
			// same and similar
			// letters;
			if (similarLen > 4) { // only if the length of similar and same
				// letters is greater than 4
				aliPos = new Pair().format(alignment);
				// System.out.println( "the alignment: " + aliPos + " " +
				// similarLen );
			}
		} catch (Exception e) {
			// logger.log(Level.SEVERE, "Failed running example: " +
			// e.getMessage(), e);
		}
		// System.out.println(aliPos);
		return aliPos;
	}

	public Vector checkLeftRegs(int aliStart, int aliEnd, int start, int end,
			String seq, double cCut) {
		Vector left = new Vector();
		double com = 0;
		ComplexityCalculator cc = new ComplexityCalculator();
		cc.initializeAlphabet();
		cc.initializeNewAlphabet();
		cc.computeNewScoringMatrix();
		cc.normalizeNewScoringMatrixPow();
		if (aliStart > 7) { // the length of the left region must be longer than
			// 7
			com = cc.calculateNor2LetterEntropyWScoMatrix(seq.substring(
					start - 1, start + aliStart - 2));
			// com = cc.calculateNorModifiedEntropy( seq.substring( start - 1,
			// start + aliStart - 2 ));
			// System.out.println( "left1: " + seq.substring( start - 1, start +
			// aliStart - 2 ) + " " + com + " " + start + "-" + ( start +
			// aliStart - 2 ));
			if (com <= cCut)
				left.add(0, start + "-" + (start + aliStart - 2));
		}
		if ((end - start + 1 - aliEnd) > 7) {
			com = cc.calculateNor2LetterEntropyWScoMatrix(seq.substring(start
					+ aliEnd - 1, end));
			// com = cc.calculateNorModifiedEntropy( seq.substring( start +
			// aliEnd - 1, end ));
			// System.out.println( "left2:" + seq.substring( start + aliEnd - 1,
			// end ) + " " + com + " " + ( start + aliEnd ) + "-" + end );
			if (com <= cCut)
				left.add((start + aliEnd) + "-" + end);
		}
		return left;
	}

	public Vector addToResult(Vector result, Vector left) {
		Vector tmpResult = result;
		String str1 = new String(), str2 = new String();
		int j = 0;
		for (int i = 0; i < left.size(); i++) {
			str1 = (String) left.elementAt(i);
			int index = str1.indexOf("-");
			int endLeft = Integer.parseInt(str1.substring(index + 1));
			boolean found = false;
			while (!found) {
				if (j < result.size()) {
					str2 = (String) result.elementAt(j);
					index = str2.indexOf("-");
					int startResult = Integer
							.parseInt(str2.substring(0, index));
					if (endLeft < startResult) {
						found = true;
						// System.out.println( "Insert left into result: " +
						// str1 + " " + str2 );
						result.add(j, str1);
						j = j + 2;
					} else
						j++;

				} else {
					result.add(str1);
					// System.out.println( "append to the end of result" );
					found = true;
				}
			}
		}
		return result;
	}

	public Vector checkAdjBlock(int start1, int end1, String adjBlock,
			String seq, double cCut, String mark) {
		Vector result = new Vector();
		double com = 0;
		String seq1 = seq.substring(start1 - 1, end1);
		// System.out.println( "current block: " +start1 + " " + end1 + " " +
		// seq1 );
		int index1 = adjBlock.indexOf("-");
		int start2 = Integer.parseInt(adjBlock.substring(0, index1));
		int end2 = Integer.parseInt(adjBlock.substring(index1 + 1));
		String seq2 = seq.substring(start2 - 1, end2);
		String aliPos = new String();
		if (mark.equals("front"))
			aliPos = findAlignment(seq2, seq1);
		else
			aliPos = findAlignment(seq1, seq2);
		// System.out.println("ASJ " +adjBlock);
		// System.out.println("1 " + aliPos);
		if (aliPos.length() != 0) { // format of aliPos: 'a1-a2 b1-b2'
			index1 = aliPos.indexOf("-");
			int index2 = aliPos.indexOf(" ");
			int aliStart2 = 0;
			
				aliStart2 = Integer.parseInt(aliPos.substring(0, index1));
			
			int aliEnd2 = Integer
					.parseInt(aliPos.substring(index1 + 1, index2));
			String aliSeq1 = new String(), aliSeq2 = new String();
			if (mark.equals("front")) {
				aliSeq2 = seq.substring(start2 + aliStart2 - 2, start2
						+ aliEnd2 - 1);
				// System.out.println( "the first aligned subSeq: " + aliStart2
				// + " " + aliEnd2 + " " + aliSeq2 );
			} else {
				aliSeq1 = seq.substring(start1 + aliStart2 - 2, start1
						+ aliEnd2 - 1);
				// System.out.println( "the first aligned subSeq: " + aliStart2
				// + " " + aliEnd2 + " " + aliSeq1 );
			}
			aliPos = aliPos.substring(index2 + 1);

			// System.out.println("2 " + aliPos);
			index1 = aliPos.indexOf("-");
			int aliStart1 = Integer.parseInt(aliPos.substring(0, index1));
			int aliEnd1 = Integer.parseInt(aliPos.substring(index1 + 1));
			if (mark.equals("front")) {
				aliSeq1 = seq.substring(start1 + aliStart1 - 2, start1
						+ aliEnd1 - 1);
				// System.out.println( "the second subSeq: " + aliStart1 + " " +
				// aliEnd1 + " " + aliSeq1 );
			} else {
				aliSeq2 = seq.substring(start2 + aliStart1 - 2, start2
						+ aliEnd1 - 1);
				// System.out.println( "the second subSeq: " + aliStart1 + " " +
				// aliEnd1 + " " + aliSeq2 );
			}

			boolean decOrNot = true;
			if (mark.equals("front")) {
				decOrNot = true;
				// decOrNot = checkCombinedSubBlock( aliSeq2,aliSeq1, cCut );
				if (decOrNot) {
					result.add((start1 + aliStart1 - 1) + "-"
							+ (start1 + aliEnd1 - 1));
					// System.out.println( "added to result1: " + ( start1 +
					// aliStart1 - 1 )+ "-" + ( start1 + aliEnd1 - 1 ));
					Vector left = checkLeftRegs(aliStart1, aliEnd1, start1,
							end1, seq, cCut);
					result = addToResult(result, left);
				}
			} else {
				// decOrNot = checkCombinedSubBlock( aliSeq1,aliSeq2, cCut );
				decOrNot = true;
				if (decOrNot) {
					result.add((start1 + aliStart2 - 1) + "-"
							+ (start1 + aliEnd2 - 1));
					// System.out.println( "added to redult2: " + ( start1 +
					// aliStart2 - 1) + "-" + ( start1 + aliEnd2 - 1));
					Vector left = checkLeftRegs(aliStart2, aliEnd2, start1,
							end1, seq, cCut);
					result = addToResult(result, left);
					result.add((start2 + aliStart1 - 1) + "-"
							+ (start2 + aliEnd1 - 1));
					// System.out.println( "added to redult3: " + ( start2 +
					// aliStart1 - 1 ) + "-" + ( start2 + aliEnd1 - 1 ));
					left = checkLeftRegs(aliStart1, aliEnd1, start2, end2, seq,
							cCut);
					result = addToResult(result, left);
				}
			}
		}
		return result;
	}

	public Vector checkDeletability(Vector lcrs, int maxIndex, String seq,
			double cCut) {
		// System.out.println(lcrs);
		Vector result = new Vector();
		String block = new String();
		int start1 = 0, end1 = 0;
		block = (String) lcrs.elementAt(maxIndex);
		int index = block.indexOf("-");
		start1 = Integer.parseInt(block.substring(0, index));
		end1 = Integer.parseInt(block.substring(index + 1));
		if (maxIndex != 0) {
			block = (String) lcrs.elementAt(maxIndex - 1);
			//System.out.println( "front adjacent block: " + block );
			try {
			    result = checkAdjBlock(start1, end1, block, seq, cCut, "front");
			} catch (Exception e) {
				//System.out.println("start1: " + start1 + " end1: " + end1
					//	+ " block: " + block + " seq: " + seq + " cCut: "
						//+ cCut);
			}
		}
		if (result.size() == 0) {
			if (maxIndex != (lcrs.size() - 1)) {
				block = (String) lcrs.elementAt(maxIndex + 1);
				//System.out.println( "back adjacent block: " + block );
				try {
					result = checkAdjBlock(start1, end1, block, seq, cCut,
							"back");
				} catch (Exception e) {
					//System.out.println("start1: " + start1 + " end1: " + end1
						//+ " block: " + block + " seq: " + seq + " cCut: "
							//+ cCut);
				}
				if (result.size() != 0)
					result.add("back");
			}
		}
		return result;
	}

	public Vector readSampledLenRepPer() {
		Vector v = new Vector();
		try {
			File f = new File("knowledge/sampledLenRepPer");
			RandomAccessFile rf = new RandomAccessFile(f, "r");
			String line = rf.readLine();
			while (line != null) {
				v.add(line);
				line = rf.readLine();
			}
			rf.close();
		} catch (IOException ex) {
		}
		return v;
	}

	public Vector filter(Vector lcrs, String seq) {
		Vector tmpLcrs = lcrs, com = new Vector();
		float len = tmpLcrs.size();
		double singleCom = 0, max = -222222222;
		String str = new String(), block = new String();
		if (tmpLcrs.size() != 1) {
			ComplexityCalculator cc = new ComplexityCalculator();
			cc.initializeAlphabet();
			cc.initializeNewAlphabet();
			cc.computeNewScoringMatrix();
			cc.normalizeNewScoringMatrixPow();
			for (int i = 0; i < len; i++) {
				str = (String) tmpLcrs.elementAt(i);
				int index = str.indexOf("-");
				str = seq.substring(
						Integer.parseInt(str.substring(0, index)) - 1, Integer
								.parseInt(str.substring(index + 1)));
				singleCom = cc.calculateNor2LetterEntropyWScoMatrix(str);
				// singleCom = cc.calculateReciprocalPro( str );
				// singleCom= cc.calculateNorModifiedEntropy( str );
				com.add(Double.toString(singleCom));
			}
			int i = 0, j = 0, maxIndex = 0;
			double limit = 0, cCut = 0;
			/*
			 * if ( seq.length() > 500 ) limit = len * 0.6; // 0.6 else limit =
			 * len * 0.3; //0.3
			 */
			/*
			 * if ( len < 350 ) len = 350; double p = 1 - 350/len; limit = len *
			 * p;
			 */
			Vector sample = readSampledLenRepPer();
			int seqLen = seq.length();
			// find the right cut percentage
			boolean found = false;
			String range = new String();
			double per = 0, prePer = 0;
			int shortest, longest, preLongest = 0;
			while ((!found) && (i < sample.size())) {
				range = (String) sample.elementAt(i);
				range = range.trim();
				int index = range.indexOf(",");
				shortest = Integer.parseInt(range.substring(0, index));
				int index2 = range.indexOf(" ");
				longest = Integer.parseInt(range.substring(index + 1, index2));
				index = range.indexOf(":");
				per = Double.parseDouble(range.substring(index + 1));
				if ((seqLen >= shortest) && (seqLen <= longest)) {
					limit = len * (1 - per);
					found = true;
				} else if (seqLen < shortest) {
					int diff1 = shortest - seqLen;
					int diff2 = seqLen - preLongest;
					if ((diff2 > diff1) || (preLongest == 0)) {
						limit = len * (1 - per);
						found = true;
					} else {
						limit = len * (1 - prePer);
						found = true;
					}
				}
				++i;
				preLongest = longest;
				prePer = per;
			}
			if (!found)
				limit = len * (1 - per);
			i = 0;
			while (i < limit) {
				int l = com.size();
				j = 0;
				max = -222222222;
				while (j < l) {
					str = (String) com.elementAt(j);
					singleCom = Double.parseDouble(str);
					if (singleCom > max) {
						max = singleCom;
						maxIndex = j;
					}
					++j;
				}
				// System.out.println( "com: " + com.elementAt( maxIndex ));
				cCut = Double.parseDouble((String) com.remove(maxIndex));
				++i;
			}
			// System.out.println( "cCut: " + cCut );
			j = 0;
			i = 0;
			len = tmpLcrs.size();
			Vector result = new Vector();
			while ((i < limit) && (len != 1) && (j < len)) {
				str = (String) tmpLcrs.elementAt(j);
				int index = str.indexOf("-");
				// System.out.println( "current block?" + str.substring( 0,
				// index ) + "?" + str.substring( index + 1 ) );
				str = seq.substring(
						Integer.parseInt(str.substring(0, index)) - 1, Integer
								.parseInt(str.substring(index + 1)));
				// System.out.println( "cur subseq???" + str );
				// singleCom = cc.calculateReciprocalPro( str );
				// singleCom= cc.calculateNorModifiedEntropy( str );
				singleCom = cc.calculateNor2LetterEntropyWScoMatrix(str);
				// System.out.println( "singlecom:"+ singleCom + " " + cCut );
				if (singleCom >= cCut) {
					result = new Vector();
					result = checkDeletability(tmpLcrs, j, seq, cCut);
					int rSize = result.size();
					boolean fromBack = false;
					if (rSize != 0) {
						str = (String) result.elementAt(rSize - 1);
						if (str.equals("back")) {
							--rSize;
							fromBack = true;
						}
						// System.out.println( "remove: "+
						// (String)tmpLcrs.remove( j ));
						tmpLcrs.remove(j);
						for (int k = 0; k < rSize; k++) { // add 'result'
							// Vector into
							// tmpLcrs in order
							str = (String) result.elementAt(k);
							// System.out.print( "*"+ str + "*" );
							tmpLcrs.add(j, str);
							++j;
						}
						// System.out.println();
						// if ( j < tmpLcrs.size())
						// System.out.println( "the next one:" +
						// tmpLcrs.elementAt( j ));
						if (fromBack) {
							// System.out.println( "yes, from back " );
							// System.out.println( "removed: "+ tmpLcrs.remove(
							// j ));
							tmpLcrs.remove(j);
						}
						len = tmpLcrs.size();
					} else {
						tmpLcrs.remove(j);
						// System.out.println( "removed coz of high complexity"
						// );
						++i;
						len = tmpLcrs.size();
					}
				} else
					++j;
			}
			// System.out.println(i + " " + j );
		}
		return tmpLcrs;
	}

	public Vector postProcess(int i, String seq, String mark) {
		Vector lcrs = new Vector();
		Vector blocks = new Vector();
		Vector pos = getPositions(i);
		if (pos.size() != 0) {
			blocks = sortPositions(pos);
			// blocks = filter( blocks, seq );
			// printLCRs( blocks );
		}
		lcrs = pickUpDrop(blocks, seq);
		// printLCRs( lcrs );
		lcrs = mergePurge(lcrs);
		// printLCRs( lcrs );
		lcrs = filter(lcrs, seq);
		lcrs = filter(lcrs, seq); // /////////// filter for a second time
		// System.out.println( "AFTER************" );
		// System.out.println("i: " + i + " Seq: " + seq + " Mark: " + mark);
		printLCRs(lcrs, seq, mark);
		return lcrs;
	}

	public void computeLCRPercentage(Vector lcrs, String str) {

	}

	public void printLCRs(Vector LCRs, String seq, String mark) {
		// System.out.println("Seq: " + seq + " mark: " + mark);
		int len = LCRs.size();
		// System.out.println("Len: " + len);
		String posBlocks = new String();
		for (int i = 0; i < len; i++) {
			String str = (String) LCRs.get(i);
			posBlocks = posBlocks + " " + str;
		}
		// System.out.println("posBlocks before: " + posBlocks);
		posBlocks = posBlocks.trim();

		if (mark.equals("0")) // generate LCR blocks
			System.out.println(posBlocks);
		else { // generate masked sequences
			Masker ms = new Masker();
			// System.out.println("seq: " + seq + " posBlocks: " + posBlocks + "
			// End");
			ms.mask(seq, posBlocks);
		}
	}

	public void startt(int th1, int th2, int th3, int th4, String mark) {
		// System.out.println("From startt:");
		// System.out.println("th1: " + th1 + " th2: " + th2 + " th3: " + th3
		// + " th4: " + th4 + " mark: " + mark);
		String str = new String(), id = new String(), nextId = new String();
		int index = 0;
		boolean first = true;
		try {
			while (str != null) {
				str = rf.readLine();
				str = generateSequence(str);
				if (str != null) {
					if (str.indexOf(">") != -1) {
						index = str.indexOf("!");
						id = nextId;
						nextId = str.substring(0, index);
						str = str.substring(index + 1);
					} else
						id = nextId;
					if (!first) {
						System.out.println();
						System.out.println(id);
					} else {
						id = nextId;
						first = false;

					}
					vertices.clear();
					edges.clear();
					lps.clear();
					for (int i = 0; i < 20; i++) {
						fVecNor[i] = 0f;
						fVecUnNor[i] = 0f;
					}
					str = str.trim();
					// System.out.println(str);////////////////

					int i = workOnSequence(str, th1, th2, th3, th4);
					// System.out.println("Before work on sequence, i: " + i
					// + " th1: " + th1 + " th2: " + th2 + " th3: " + th3
					// + " th4: " + th4);
					// System.out.println(" Before work on sequence : str : "
					// + str);
					Vector lcrs = postProcess(i, str, mark);// process all
					// longest paths
					// from every
					// connected
					// subgraph

					/*
					 * computeLCRPercentage( lcrs, str ); printLCRs( lcrs, str );
					 */
				}
			}
			rf.close();
		} catch (IOException ex) {
		}
	}

	public static void main(String args[]) {
		for (int i = 0; i < 7; i++) {

			//System.out.println("Args[" + i + "]: " + args[i]);
		}

		int th1 = Integer.parseInt(args[2]);
		int th2 = Integer.parseInt(args[3]);
		int th3 = Integer.parseInt(args[4]);
		int th4 = th2;
		Gbm g = new Gbm(args[0]);
		g.readRNRMatrices(args[1]);

		/**
		 * Fixed by Nirmalya :11-05-2008 The location of the blosum matrix was
		 * hardcoaded. Now it needs to be given through args[5].
		 */
		// g.readScoringMatrix("knowledge/blosum62Matrix");
		g.readScoringMatrix(args[5]);
		g.startt(th1, th2, th3, th4, args[6]);

	}

}

/*
 * 
 * 1. Format of the output: SACACPQTSOP......( 60 letters) XXXX
 * TPQSKAQ..........( 60 letters)
 * 
 */
