/**
 * @author: Vikas Dayananda
**/

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.ArrayList;

/**
 *Used to signal violations of preconditions for
 *various shortest path algorithms.
 */
class GraphException extends RuntimeException {
	
	public GraphException(String name) 
	{
		super(name);
	}
}

/**
 * Edge class has information about vertex source name,destination name, cost and status (availability)  
 */

class Edge {
	public Vertex source, dest;
	public double cost;
	boolean status;

	public Edge(Vertex from, Vertex to, double cost, boolean status) {
		source = from;
		dest = to;
		this.cost = cost;
		this.status = status;
	}

	public Edge() {
	}

	public String getsource() {
		return source.name;
	}

	public String getdest() {
		return dest.name;
	}

	public double getCost() {
		return cost;
	}

	public boolean getstatus() {
		return status;
	}

	public void setstatus(boolean status) {
		this.status = status;
	}

	public void setCost(double cost) {
		this.cost = cost;
	}

}

/**Vertex class has information about vertex name,status,previous vertex, distance and
 *List of edges stores the neighbour edges information of a vertex.
 */
class Vertex {

	public Vertex() {

	}

	public String name;
	public boolean status;
	public double dist;
	public Vertex prev;
	public ArrayList<Edge> neighbour;

	public Vertex(String nm, boolean status) {
		name = nm;
		this.status = status;
		neighbour = new ArrayList<Edge>();
		reset();
	}

	public boolean getstatus() {
		return status;
	}

	public void setstatus(boolean status) {
		this.status = status;
	}

	public double getDist() {
		return dist;
	}

	public void setDist(double dist) {
		this.dist = dist;
	}

	public Vertex getPrev() {
		return prev;
	}

	public void setPrev(Vertex prev) {
		this.prev = prev;
	}

	public void reset() {
		dist = Double.POSITIVE_INFINITY;
		prev = null;

	}

}

/**
 * Graph Class initilizes and handles all the operations related to the Graph.
 * It has main method which reads the input file and executes queries.
// ******************PUBLIC OPERATIONS**********************
 * boolean processQuery(Scanner in,Graph g)  --> Executes queries.
 * void initializeGraph(String s, String d, double c, boolean st) --> Build an undirected graph.
 * void addEdge( String v, String w )		 --> Add additional edge
 * void deleteEdge( String v, String w )	 --> Delete the edge
 * void edgedown (String v, String w )       --> Make edge unavailable
 * void edgeup (String v, String w )         --> Make edge available	
 * void printPath( String w )  				 --> Print path after alg is run
 * void unweighted( String s ) 				 --> Single-source unweighted
 * booelan isEdgeDown(String source, String dest) --> Find if edge is down
 * void vertexUp ( String v ) 				 --> Make vertex available
 * void vertexDown ( String v )				 --> Make vertex unavailable
 * void Dijkstra ( String src )  			 --> Find shortest path from source
 * void printPath ( String dest)             --> Print path from source to destination
 * void DFS ()                               --> Visit each vertex and call DFS_Visit\
 * void DFS_Visit(String node)               --> Iterate through all adjacent edges
 */

public class Graph {
 
	private static TreeMap<String, Vertex> vertices = new TreeMap<String, Vertex>(); 	// Store vertices information																				 
	private LinkedList<Pair> downEdges = new LinkedList<Pair>(); 						// Store list of down edges.																// 
	public static final double INFINITY = Double.POSITIVE_INFINITY; 					// Define Infinity
	private TreeSet<String> reachableNodes = new TreeSet<String>(); 					// Stores set of reachable vertex.																
	public static File f; 																// Type file to read input file
	
	/**
     * A main routine that:
     * 1. Reads a file containing edges (supplied as a command-line parameter);
     * 2. Forms the graph;
     * 3. Repeatedly prompts for queries and
     *    runs the corresponding functions.
     * The data file is a sequence of lines of the format
     *    source destination cost
     */
	public static void main(String[] args) {
		Graph g = new Graph();											//Graph object is created to call functions
		g.downEdges.add(new Pair(" ", " ")); 							// Add empty pair of down edges ( initialization )
												 
		boolean status = true; 											// Initially all vertices are available
		try {
			FileReader fin = new FileReader(args[0]);					// Read the input file
			Scanner graphFile = new Scanner(fin);

			String line;
			System.out.println("Reading file " + args[0] + ".....");

			while (graphFile.hasNextLine()) {
				line = graphFile.nextLine();
				StringTokenizer st = new StringTokenizer(line);

				try {
					if (st.countTokens() != 3) {
						System.err.println("Skipping ill-formatted line " + line);
						continue;
					}
					String source = st.nextToken(); 							// get Source vertex
					String dest = st.nextToken(); 								// get destination vertex
					double distance = Double.parseDouble(st.nextToken()); 		// get  distance
																			
					g.initializeGraph(source, dest, distance, status);			 // Build an undirected graph
													 
				} catch (NumberFormatException e) {
					System.err.println("File content format not correct " + line);
				}
			}
			graphFile.close();
		} catch (IOException e) {
			System.err.println(e);
		}

		Scanner in = new Scanner(System.in);

		while (processQuery(in, g))
			;
	}
	
	/**
	 * The function will ask for queries from the user repeatedly and calls corresponding functions
	 * The program will exit when the command "quit" is entered.
	 */
	
	public static boolean processQuery(Scanner in, Graph g) {
		try {
			boolean status = true;
			String qry = in.nextLine();

			String[] arr = qry.split(" ");

			switch (arr[0]) {

			case "print":
				System.out.println("");
				g.print();
				break;

			case "quit":
				return false;

			case "reachable":
				System.out.println("");
				g.DFS();

				break;

			case "path":
				if (arr.length == 3) {
					String source = arr[1];
					String dest = arr[2];
					if ((vertices.get(source) == null) || (vertices.get(dest) == null)) {
						System.out.println("Vertex does not Exist");
						return true;
					}
					g.dijkstraAlgo(source);
					g.printPath(dest);
				} else
					System.out.println("Invalid Path Query");
				break;

			case "edgeup":
				if (arr.length == 3) {
					String source = arr[1];
					String dest = arr[2];
					if ((vertices.get(source) == null) || (vertices.get(dest) == null)) {
						System.out.println("Edge does not Exist");
						return true;
					}
					g.edgeUp(source, dest);
				} else
					System.out.println("Invalid EdgeUp Query");
				break;

			case "edgedown":
				if (arr.length == 3) {
					String source = arr[1];
					String dest = arr[2];
					if ((vertices.get(source) == null) || (vertices.get(dest) == null)) {
						System.out.println("Edge does not Exist");
						return true;

					}
					g.edgeDown(source,dest);
				} else
					System.out.println("Invalid EdgeDown Query");
				break;

			case "vertexup":
				if (arr.length == 2) {
					String vertex = arr[1];
					if ((vertices.get(vertex) == null)) {
						System.out.println("Vertex" + vertex + " does not Exist");
						return true;
					}
					g.vertexUp(vertex);
				} else
					System.out.println("Invalid VertexUp Query");
				break;

			case "vertexdown":
				if (arr.length == 2) {
					String vertex = arr[1];
					if ((vertices.get(vertex) == null)) {
						System.out.println("Vertex" + vertex + " does not Exist");
						return true;
					}
					g.vertexDown(vertex);
				} else
					System.out.println("Invalid VertexDown Query");
				break;

			case "deleteedge":
				if (arr.length == 3) {
					String source = arr[1];
					String dest = arr[2];
					if ((vertices.get(source) == null) || (vertices.get(dest) == null)) {
						System.out.println("Edge does not Exist");
						return true;

					}
					g.deleteEdge(source,dest);
				} else
					System.out.println("Invalid DeleteEdge Query");
				break;

			case "addedge":
				if (arr.length == 4) {
					Double dist = Double.parseDouble(arr[3]);
					if (dist < 0) {
						System.out.println("Negative Edges not accepted");
						return true;
					}
					g.addEdge(arr[1], arr[2], Double.parseDouble(arr[3]), status);
				} else
					System.out.println("Invalid AddEdge Query");
				break;

			default:
				System.out.println("Exit !!");
				break;
			}
		} catch (NoSuchElementException e) {
			System.out.println(" ");
			return false;
		} catch (GraphException e) {
			System.err.println(e);
		}
		return true;
	}
	
	 /**
     * Build a initial graph with undirected edges. 
     * If vertex name is not present in vertices map, add it.
     */

	public void initializeGraph(String sourceName, String destName, double cost, boolean status) {

		Vertex src = vertices.get(sourceName);
		Vertex dest = vertices.get(destName);
		if (src == null) {
			src = new Vertex(sourceName, status);
			vertices.put(sourceName, src);
		}
		if (dest == null) {
			dest = new Vertex(destName, status);
			vertices.put(destName, dest);
		}

		src.neighbour.add(new Edge(src, dest, cost, status));
		dest.neighbour.add(new Edge(dest, src, cost, status));
	}
	
	 /**
     * Add a new edge to the graph.
     * command - addeedge < source > < dest > < cost >
     */
	
	public void addEdge(String sourceName, String destName, double cost, boolean status) {
		Vertex src = vertices.get(sourceName);
		Vertex dest = vertices.get(destName);
		if (src == null) {
			src = new Vertex(sourceName, status);
			vertices.put(sourceName, src);
		}
		if (dest == null) {
			dest = new Vertex(destName, status);
			vertices.put(destName, dest);
		}
		ListIterator<Edge> it = src.neighbour.listIterator();
		while (it.hasNext()) {
			Edge edge = (Edge) it.next();
			if (destName.equals(edge.getdest())) {
				it.remove();
			}
		}
		src.neighbour.add(new Edge(src, dest, cost, status));
	}

	/**
     * Delete the edge entered by user. 
     * command - deleteedge < source > < dest >
     */
	
	public void deleteEdge(String sourceName, String destName) {
		Vertex src = vertices.get(sourceName);
		ListIterator<Edge> it = src.neighbour.listIterator();
		while (it.hasNext()) {
			Edge edge = (Edge) it.next();
			if (destName.equals(edge.getdest())) {
				it.remove();
			}
		}
	}
	
	/**
     * Make Edge available to use. 
     * command - edgeup < source > < dest >
     */

	public void edgeUp(String source, String dest) {

		Iterator<Pair> itr = downEdges.iterator();
		int exist = 0;
		Pair temp = null;
		if (downEdges.size() != 0) {
			while (itr.hasNext()) {
				temp = itr.next();
				if (source.equals(temp.key()) && dest.equals(temp.value())) {
					exist = 1;
					break;
				}
			}
		}
		if (exist == 1)
			downEdges.remove(temp);
	}
	
	/**
     * Make Edge unavailable to use. 
     * command - edgedown < source > < dest >
     */

	public void edgeDown(String source, String dest) {
		Iterator<Pair> itr = downEdges.iterator();
		int add = 1;

		if (downEdges.size() != 0) {
			while (itr.hasNext()) {
				Pair temp = itr.next();
				if (temp.key() == source && temp.value() == dest) {
					add = 0;
					break;
				}
			}
		}

		if (add == 1)
			downEdges.add(new Pair(source, dest));
	}
	
	/**
	 * This function will return if the edge is down/disabled else it will return false.
	 */

	public Boolean isEdgeDown(String source, String dest) {
		Iterator<Pair> itr = downEdges.iterator();
		int down = 0;

		if (downEdges.size() != 0) {
			while (itr.hasNext()) {
				Pair temp = itr.next();
				if (temp.key().equals(source) && temp.value().equals(dest)) {
					down = 1;
					break;
				}
			}
		}

		if (down == 0)
			return false;
		else
			return true;
	}

	/**
	 * Makes the Vertex Up by setting status as true
	 * command - vertexup < name > 
	 */
	
	public void vertexUp(String vertex) {
		// TODO Auto-generated method stub
		Vertex v = vertices.get(vertex);
		v.setstatus(true);
	}

	/**
	 * Makes the Vertex Down by setting status as false
	 * command - vertexdown < name > 
	 */
	
	public void vertexDown(String source) {
		// TODO Auto-generated method stub
		Vertex v = vertices.get(source);
		v.setstatus(false);
	}
	
	/**
	 * Sets all distances of vertices as Infinity and prev vertex as null prior to running shortest path
	 * algorithm.
	 */
	private void clearAll() {
		for (Vertex v : vertices.values())
			v.reset();
	}
	
	/**
	 * Calculates the shortest path between to vertex using Dijkstra's Algorithm and Min Heap.
	 */
	
	public void dijkstraAlgo(String startName) {
		
		clearAll();
		int size= vertices.size();
		MinHeap pQueue = new MinHeap(size);
		Vertex start = vertices.get(startName);
		if (start.status == false) {
			System.out.println("Vertex is Down");
			return;
		}
		start.dist = 0.0;
		pQueue.insert(new Path(start.name, start.dist));
		while (pQueue.size() != 0) {
			Vertex v = vertices.get((pQueue.remove().src));

			for (Edge edge : v.neighbour) {
				Vertex w = edge.dest;
				if (w.status == false)
					continue;
				if (isEdgeDown(v.name, w.name))
					continue;
				if (edge.cost < 0)
					throw new GraphException("Graph has negative edges");
				if (w.dist > v.dist + edge.cost ) {
					w.dist = v.dist + edge.cost ;
					w.prev = v;
					pQueue.insert(new Path(w.name, w.dist));
				}
			}
		}
	}
	
	/**
	 * prints shortest path from source to destination vertex after running shortest path algorithm.
	 */

	private void printPath(Vertex dest) {
		if (dest.prev != null) {
			printPath(dest.prev);
			System.out.print(" ");
		}
		System.out.print(dest.name);
	}

	public void printPath(String destName) {
		Vertex w = vertices.get(destName);
		if (w == null) {
			System.out.println("Destination vertex not found");
			return;
		} else if (w.dist == INFINITY)
			System.out.println(destName + " is unreachable");
		else {
			System.out.println();
			printPath(w);
			System.out.printf(" %.2f", w.dist);
			System.out.println();
		}
	}
	
	
	/**
	 * This function will print all the contents of the graph in Alphabetical Order.
	 * If a vertex or an edge is down then it will append the string "DOWN" near it.
	 */

	public void print() {
		for (Vertex v : vertices.values()) {
			Collections.sort(v.neighbour, new Comparator<Edge>() {
				public int compare(Edge e1, Edge e2) {

					return e1.getdest().compareTo(e2.getdest()) < 0 ? -1 : 1;
				}

			});
		
			if (v.status == false) {
				System.out.println(v.name + " " + "DOWN");
			} else {
				System.out.println(v.name);
			}
			ListIterator<Edge> it = v.neighbour.listIterator();
			while (it.hasNext()) {
				Edge edge = new Edge();
				edge = (Edge) it.next();
				String dest = edge.getdest();
				double cos = edge.getCost();
				if (edge.status == false) {
					System.out.println(" " + dest + " " + cos + " " + "DOWN");
				} else {
					System.out.println(" " + dest + " " + cos);
				}
			}

		}
		System.out.println("");
	}

	/**
	 * REACHABLE VERTICES ALGORITHM ( DFS )
	 * This algorithm will display all the vertices that are reachable from each of the Vertex in alphabetical order.
	 * Instead of Maintaing color , i created a structure called reachableNodes of Type TreeSet to maintain list of nodes visited.
	 * Source - Lecture Slide 11.
	 *--------------------------------------------------------------------------------------------
	 * 		ALGORITHM DFS( )
	 * 		1. Sort all the vertices in alphabetical order
	 * 		2. foreach Vertex which is UP
	 * 		3. 		print Vertex_Name
	 * 		4.		add Vertex to visited_nodes list
	 * 		5.		call REACHABLE(Vertex) algorithm
	 * 		6.	endforeach;
	 * 		7.	print all the Vertices from visited_nodes in alphabetical order
	 * 
	 * 		ALGORITHM DFS_Visit(Vertex)
	 * 		1. forach adjacent vertices
	 * 		2.		if(edge is down or vertex is down) then
	 * 		3.			continue to next iteration
	 * 		4.		endif;
	 * 		5.		if(Vertex exists in visited_node list) then
	 * 		6.			continue
	 * 		7.		else
	 * 		8.			add vertex to visited nodes list
	 * 		9.			call REACHABLE(Vertex)
	 * 		10.		endif;
	 * 		11.	endfor;
	 * ------------------------------------------------------------------------------------------------------
	 * TIME COMPLEXITY ->
	 * ------------------------------------------------------------------------------------------------------
	 * Here the algorithm DFS( ) will be called for each of the V vertices. DFS_Visit is called once for each Vertex. 
	 * So total cost of executing DFS_Visit is |E|.
	 * The Total Running time of DFS is O( |V| + |E| )											
	 * ------------------------------------------------------------------------------------------------------
	 */

	private void DFS() {

		TreeSet<String> tree = new TreeSet<String>();
		for (String key : vertices.keySet()) {
			if (vertices.get(key).status == false)
				continue;
			tree.add(key);
		}
		Iterator<String> itr = tree.iterator();
		String adjV;

		while (itr.hasNext()) {
			adjV = itr.next();
			System.out.print(adjV);
			System.out.println();

			reachableNodes.clear();
			reachableNodes.add(adjV);

			DFS_Visit(adjV);

			Iterator<String> adjItr = reachableNodes.iterator();
			while (adjItr.hasNext()) {
				String adjItrN = adjItr.next();
				if (adjItrN != adjV)
					System.out.println("  " + adjItrN);
			}
		}
		System.out.println("");
	}

	/**
	 * This recursive algorithm will be called by the dfsAlgo() function. This will iterate through all the adjacent edges.
	 */
	private void DFS_Visit(String node) {

		TreeSet<String> neighbour = new TreeSet<String>();
		for (Edge key : vertices.get(node).neighbour) {
			if (key.dest.status == false)
				continue;
			if (isEdgeDown(node, key.dest.name))
				continue;

			neighbour.add(key.dest.name);
		}
		for (String str : neighbour) {
			if (reachableNodes.contains(str))
				continue;
			else {
				reachableNodes.add(str);   
				DFS_Visit(str);
			}
		}
	}
}

/**
 * Path Class which stores the name of a Path.
 */
class Path {
	public String src;
	public double cost;

	public Path(String s, Double d) {
		src = s;
		cost = d;
	}

}

/**
 * This Class is a User-defined datatype which stores a pair of Strings.
 */
class Pair {
	private final String str1;
	private final String str2;

	public Pair(String aStr1, String aStr2) {
		str1 = aStr1;
		str2 = aStr2;
	}

	public String key() {
		return str1;
	}

	public String value() {
		return str2;
	}
}

/**
 * MinHeap Class which is used by Dijkstra's Algorithm. This class will
 * calculate Min-Heap of the tree, acting as a Priority Queue.
 */
class MinHeap {
	private Path[] Heap; 				// Heap variable declared of Path datatype
	private int size; 					// Size of the heap
	private int maxsize; 				// Maxsize of the heap (total number of vertices)

	private static final int TOP = 1; // Constant storing the TOP value

	/**
	 * Constructor which initializes the Min-Heap.
	 */
	public MinHeap(int maxsize) {
		this.maxsize = maxsize;
		this.size = 0;
		Heap = new Path[this.maxsize + 1];
		Heap[0] = new Path("null", (double) Integer.MIN_VALUE);
	}

	/**
	 * Returns the size of the Heap
	 */
	public int size() {
		return size;
	}
 
	/**
	 * Returns the Parent Node of a Vertex
	 */
	private int parent(int p) {
		return p / 2;
	}

	/**
	 * Returns the left child of a node
	 */
	private int leftNode(int p) {
		return (2 * p);
	}
 
	/**
	 * Returns the right child of a node
	 */
	private int rightNode(int p) {
		return (2 * p) + 1;
	}

	/**
	 * Checks if the node is a leaf node or not
	 */
	private boolean isLeaf(int p) {
		if (p - 1 >= (size / 2) && p <= size) {
			return true;
		}
		return false;
	}
	/**
	 * Swaps the position in the heap if the child node is smaller than the parent.
	 */
	private void swap(int p, int pos) {
		Path tmp;
		tmp = Heap[p];
		Heap[p] = Heap[pos];
		Heap[pos] = tmp;
	}

	private void minHeapify(int pos) {
		if (!isLeaf(pos) && this.size > 0) {
			if (Heap[pos].cost > Heap[leftNode(pos)].cost || Heap[pos].cost > Heap[rightNode(pos)].cost) {
				if (Heap[leftNode(pos)].cost < Heap[rightNode(pos)].cost) {
					
					swap(pos, leftNode(pos));
					minHeapify(leftNode(pos));
				} else {
					swap(pos, rightNode(pos));
					minHeapify(rightNode(pos));
				}
			}
		}
	}
	
	/**
	 * Inserts a new element in the heap
	 */
	public void insert(Path element) {
		Heap[++size] = element;
		int current = size;
		while (Heap[current].cost < Heap[parent(current)].cost) {
			swap(current, parent(current));
			current = parent(current);
		}
	}
	
	/**
	 * Perform minHeapify
	 */
	public void minHeap() {
		for (int pos = (size / 2); pos >= 1; pos--) {
			minHeapify(pos);
		}
	}
	
	/**
	 * Returns the smallest node i.e root node from heap and the perform minHeapify
	 * to determine new root.
	 */
	public Path remove() {
		Path popp = Heap[TOP];
		Heap[TOP] = Heap[size--];
		minHeapify(TOP);
		return popp;
	}

}
