/**
 * @author UCSD MOOC development team and Stacy Lynn Day
 * 
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
package roadgraph;


import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Consumer;

import geography.GeographicPoint;
import util.GraphLoader;

/**
 * @author UCSD MOOC development team and Stacy Lynn Day
 * 
 * A class which represents a graph of geographic locations
 * Nodes in the graph are intersections between 
 *
 */
public class MapGraph {
	private HashMap<GeographicPoint, MapNode> vertices;
	private HashSet<MapEdge> edges;
		
	/** 
	 * Create a new empty MapGraph 
	 */
	public MapGraph()
	{
		vertices = new HashMap<GeographicPoint, MapNode>();
		edges = new HashSet<MapEdge>();
	}
	
	/**
	 * Get the number of vertices (road intersections) in the graph
	 * @return The number of vertices in the graph.
	 */
	public int getNumVertices()
	{
		return vertices.values().size();
	}
	
	/**
	 * Return the intersections, which are the vertices in this graph.
	 * @return The vertices in this graph as GeographicPoints
	 */
	public Set<GeographicPoint> getVertices()
	{
		HashSet<GeographicPoint> intersections = (HashSet<GeographicPoint>) vertices.keySet();
		return intersections;
	}
	
	/**
	 * Get the number of road segments in the graph
	 * @return The number of edges in the graph.
	 */
	public int getNumEdges()
	{
		return edges.size();
	}
	
	
	/** Add a node corresponding to an intersection at a Geographic Point
	 * If the location is already in the graph or null, this method does 
	 * not change the graph.
	 * @param location  The location of the intersection
	 * @return true if a node was added, false if it was not (the node
	 * was already in the graph, or the parameter is null).
	 */
	public boolean addVertex(GeographicPoint location)
	{
		if (vertices.containsKey(location) || location == null)  {
			return false;
		}
		else {
			MapNode node = new MapNode(location);
			this.vertices.put(location, node);
			return true;
		}
	}
	
	/**
	 * Adds a directed edge to the graph from pt1 to pt2.  
	 * Precondition: Both GeographicPoints have already been added to the graph
	 * @param from The starting point of the edge
	 * @param to The ending point of the edge
	 * @param roadName The name of the road
	 * @param roadType The type of the road
	 * @param length The length of the road, in km
	 * @throws IllegalArgumentException If the points have not already been
	 *   added as nodes to the graph, if any of the arguments is null,
	 *   or if the length is less than 0.
	 */
	
	public void addEdge(GeographicPoint from, GeographicPoint to, String roadName, 
			String roadType, double length) throws NullPointerException {
		
	    // check to make sure there exist MapNodes for these points
		MapNode fromNode = vertices.get(from);
		MapNode toNode = vertices.get(to);

		// check these nodes are valid
		if (fromNode == null) {
			throw new NullPointerException("addEdge: fromNode:" + from + "is not in MapGraph");
		}	
		if (toNode == null) {
			throw new NullPointerException("addEdge: toNode:" + to + "is not in MapGraph");
		}	

		MapEdge e = new MapEdge(from, to, roadName, roadType, length);
		this.edges.add(e);
		
		// since all vertices are added first, we need to add this edge to the node
		// edge it's from point originates
		fromNode.addEdge(e);
	}
	

	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return bfs(start, goal, temp);
	}
	
	/** Find the path from start to goal using breadth first search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest (unweighted)
	 *   path from start to goal (including both start and goal).
	 */
	public List<GeographicPoint> bfs(GeographicPoint start, 
			 					     GeographicPoint goal, Consumer<GeographicPoint> nodeSearched) {
		
		// Every time you explore from a node (which is associated with a GeographicPoint), 
		// report that node to the consumer using the above instruction (assuming that 
		// next is the GeographicPoint associated with the node currently being explored):
		
		if (!arePreconditionsFulfilled(start, goal)) {
			return null;
		}
		
		//System.out.println("Beginning bfs");
		
		// Initialize start and goal nodes, queue, visited HashSet, parent HashMap
		MapNode startNode = vertices.get(start);
		MapNode goalNode = vertices.get(goal);
		Queue<MapNode> nodequeue = new LinkedList<MapNode>();
		HashSet<MapNode> visited = new HashSet<MapNode>();
		HashMap<MapNode, MapNode> parentMap = new HashMap<MapNode, MapNode>();
		
		//System.out.println("Initialization done");
		
		// Enqueue start onto the queue and add to visited
		nodequeue.add(startNode);
		visited.add(startNode);
		
		// while queue is not empty: 
		while (!nodequeue.isEmpty()) {
			// dequeue node curr from front of queue
			MapNode curr = nodequeue.remove();
			
			//System.out.println("At node: " + curr.toString());
				
			// if curr equals goal return parent map
			if (curr == goalNode) {
				// should I add to nodeSearched?
				nodeSearched.accept(curr.getLocation());
				break;
			}
			// get curr's neighbors
			Set<MapNode> neighbors = new HashSet<MapNode>();
			for (MapEdge e : curr.getEdges()) {
				neighbors.add(vertices.get(e.getToPoint()));
			}
				
			// for each of curr's neighbors, next, not visited set:
			for (MapNode next : neighbors) {
				if (!(visited.contains(next))) {
					// add next to visited set
					visited.add(next);
					// add curr as next's parent in parent map
					parentMap.put(next, curr);
					// enqueue n onto the queue
					nodequeue.add(next);
					nodeSearched.accept(curr.getLocation());
				}
			}	
		}

		// return the reconstructed path
		return reconstructPath(parentMap, startNode, goalNode);		
	}
	
	private boolean arePreconditionsFulfilled(GeographicPoint start, GeographicPoint goal) {
		if (start == null || goal == null) {
			throw new NullPointerException("Start or goal location is null!  No path exists.");
		}
		if (vertices.get(start) == null) {
			System.err.println("Start node " + start + " does not exist.");
			return false;
		}
		if (vertices.get(start) == null) {
			System.err.println("Goal node " + goal + " does not exist.");
			return false;
		}
		return true;
	}

	private List<GeographicPoint> reconstructPath(HashMap<MapNode, MapNode> parentMap, 
			MapNode startNode, MapNode goalNode) {
		List<GeographicPoint> path = new ArrayList<GeographicPoint>();
		MapNode curr = goalNode;
		while (curr != startNode) {
			path.add(curr.getLocation());
			curr = parentMap.get(curr);
		}
		path.add(startNode.getLocation());
		Collections.reverse(path);
		return path;
	}

	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
		Consumer<GeographicPoint> temp = (x) -> {};
        return dijkstra(start, goal, temp);
	}
	
	/** Find the path from start to goal using Dijkstra's algorithm
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> dijkstra(GeographicPoint start, 
										  GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		System.out.println("Beginning Dijkstra");
		
		return searchOnWeightedGraph(start, goal, nodeSearched, (a, b) -> 0.0);
	}	
		
	private List<GeographicPoint> searchOnWeightedGraph(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched, BiFunction<MapNode, MapNode, Double> f) {

		if (!arePreconditionsFulfilled(start, goal)) {
			return null;
		}
		
		// Initialize start and goal nodes,
		// Initialize PriorityQueue pq, visited HashSet, parent HashMap and
		// distances to infinity
		MapNode startNode = vertices.get(start);
		MapNode goalNode = vertices.get(goal);
		PriorityQueue<MapNode> nodequeue = new PriorityQueue<MapNode>();
		HashSet<MapNode> visited = new HashSet<MapNode>();
		HashMap<MapNode, MapNode> parentMap = new HashMap<MapNode, MapNode>();
		
		for (GeographicPoint gp : vertices.keySet()) {
			MapNode thisNode = vertices.get(gp);
			thisNode.setDistance(Double.POSITIVE_INFINITY);
			thisNode.setActualDistance(Double.POSITIVE_INFINITY);
		}
		
		//System.out.println("Initialization done");
		
		// Enqueue {s, 0} onto priority queue
		startNode.setDistance(0.0);
		startNode.setActualDistance(0.0);
		
		nodequeue.add(startNode);

		int numnodesvisited = 0;
		
		// while priority gueue (i.e. nodequeue) is not empty:
		while (!nodequeue.isEmpty()) {
			// dequeue node curr from front of queue
			MapNode curr = nodequeue.remove();
			numnodesvisited++;
			
			//System.out.println("At node: " + curr.toString());
		
			// if curr is not visited
			if (!(visited.contains(curr))) {
				// add curr to visited set
				visited.add(curr);
				// if curr == goal return parent map
				if (curr == goalNode) {
					nodeSearched.accept(curr.getLocation());
					break;
				}
				
				// create a hashmap of curr's neighbor's and edges distances
				HashMap<MapNode, Double> neighbors = new HashMap<MapNode, Double>();
				for (MapEdge e : curr.getEdges()) {
					neighbors.put(vertices.get(e.getToPoint()), e.getLength());
				}	
				
				// for each of curr's neighbors, next, not in visited set:
				for (MapNode next : neighbors.keySet()) {
					if (!(visited.contains(next))) {
						// if path through curr to next is shorter
						// update curr as next's parent in parent map
						// and enqueue {next, distance} into the pq
						double distance = curr.getActualDistance() + neighbors.get(next);
						if (distance < next.getActualDistance()) {
							//parentMap.remove(next);???
							next.setActualDistance(distance);
							distance += f.apply(next, goalNode);
							next.setDistance(distance);
							parentMap.put(next, curr);
							nodequeue.add(next);
						}
					}
				}	
			}
		}	
		
		System.out.println("Num nodes visited is : " + numnodesvisited);
		
		// return the reconstructed path
		return reconstructPath(parentMap, startNode, goalNode);
	}

	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return aStarSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using A-Star search
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> aStarSearch(GeographicPoint start, 
											 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());
		
		System.out.println("Beginning aStarSearch");
		
		return searchOnWeightedGraph(start, goal, nodeSearched, (a, b) -> a.getLocation().distance(b.getLocation()));
	}
	
	//--------------------------------------------------------------------------
	
	/** Find the path from start to goal using timed search using speed limits
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> timedDijkstraSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return timedDijkstraSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using timed search using speed limits
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */	
	public List<GeographicPoint> timedDijkstraSearch(GeographicPoint start, 
			 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());

		System.out.println("Beginning Timed Dijkstra Search");

		return searchOnWeightedGraphWithTime(start, goal, nodeSearched, (a, b) -> 0.0);
	}	
	
	/** Find the path from start to goal using timed search using speed limits
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> timedAStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return timedAStarSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using timed search using speed limits
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */	
	public List<GeographicPoint> timedAStarSearch(GeographicPoint start, 
			 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());

		System.out.println("Beginning Timed aStar Search");

		return searchOnWeightedGraphWithTime(start, goal, nodeSearched, (a, b) -> 10*(a.getLocation().distance(b.getLocation())));
	}	
	
	private List<GeographicPoint> searchOnWeightedGraphWithTime(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched, BiFunction<MapNode, MapNode, Double> f) {

		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());
		
		if (!arePreconditionsFulfilled(start, goal)) {
			return null;
		}
		
		// Initialize start and goal nodes,
		// Initialize PriorityQueue pq, visited HashSet, parent HashMap and
		// distances to infinity
		MapNode startNode = vertices.get(start);
		MapNode goalNode = vertices.get(goal);
		PriorityQueue<MapNode> nodequeue = new PriorityQueue<MapNode>();
		HashSet<MapNode> visited = new HashSet<MapNode>();
		HashMap<MapNode, MapNode> parentMap = new HashMap<MapNode, MapNode>();
		
		for (GeographicPoint gp : vertices.keySet()) {
			MapNode thisNode = vertices.get(gp);
			thisNode.setDistance(Double.POSITIVE_INFINITY);
			thisNode.setActualDistance(Double.POSITIVE_INFINITY);
		}
		
		//System.out.println("Initialization done");
		
		// Enqueue {s, 0} onto priority queue
		startNode.setDistance(0.0);
		startNode.setActualDistance(0.0);
		
		nodequeue.add(startNode);

		int numnodesvisited = 0;
		
		// while priority gueue (i.e. nodequeue) is not empty:
		while (!nodequeue.isEmpty()) {
			// dequeue node curr from front of queue
			MapNode curr = nodequeue.remove();
			numnodesvisited++;
			
			//System.out.println("At node: " + curr.toString());
		
			// if curr is not visited
			if (!(visited.contains(curr))) {
				// add curr to visited set
				visited.add(curr);
				// if curr == goal return parent map
				if (curr == goalNode) {
					nodeSearched.accept(curr.getLocation());
					break;
				}
				
				// create a hashmap of curr's neighbor's and edges time to drive edge
				HashMap<MapNode, Double> neighbors = new HashMap<MapNode, Double>();
				for (MapEdge e : curr.getEdges()) {
					neighbors.put(vertices.get(e.getToPoint()), e.getTime());
				}	
				
				// for each of curr's neighbors, next, not in visited set:
				for (MapNode next : neighbors.keySet()) {
					if (!(visited.contains(next))) {
						// if path through curr to next is shorter
						// update curr as next's parent in parent map
						// and enqueue {next, distance} into the pq
						double distance = curr.getActualDistance() + neighbors.get(next);
						if (distance < next.getActualDistance()) {
							//parentMap.remove(next);???
							next.setActualDistance(distance);
							distance += f.apply(next, goalNode);
							next.setDistance(distance);
							parentMap.put(next, curr);
							nodequeue.add(next);
						}
					}
				}	
			}
		}	
		
		System.out.println("Num nodes visited is : " + numnodesvisited);
		
		// return the reconstructed path
		return reconstructPath(parentMap, startNode, goalNode);
	}

	//--------------------------------------------------------------------------
	
	/** Find the path from start to goal using timed search using speed limits and time of day
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> timedWithTrafficDijkstraSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return timedWithTrafficDijkstraSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using timed search using speed limits
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */	
	public List<GeographicPoint> timedWithTrafficDijkstraSearch(GeographicPoint start, 
			 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());

		System.out.println("Beginning Timed With Traffic Dijkstra Search");

		return searchOnWeightedGraphWithTimeAndTraffic(start, goal, nodeSearched, (a, b) -> 0.0);
	}	
	
	/** Find the path from start to goal using timed search using speed limits and traffic
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */
	public List<GeographicPoint> timedWithTrafficAStarSearch(GeographicPoint start, GeographicPoint goal) {
		// Dummy variable for calling the search algorithms
        Consumer<GeographicPoint> temp = (x) -> {};
        return timedWithTrafficAStarSearch(start, goal, temp);
	}
	
	/** Find the path from start to goal using timed search using speed limits
	 * 
	 * @param start The starting location
	 * @param goal The goal location
	 * @param nodeSearched A hook for visualization.  See assignment instructions for how to use it.
	 * @return The list of intersections that form the shortest path from 
	 *   start to goal (including both start and goal).
	 */	
	public List<GeographicPoint> timedWithTrafficAStarSearch(GeographicPoint start, 
			 GeographicPoint goal, Consumer<GeographicPoint> nodeSearched)
	{
		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());

		System.out.println("Beginning Timed With Traffic aStar Search");

		return searchOnWeightedGraphWithTimeAndTraffic(start, goal, nodeSearched, (a, b) -> 10*(a.getLocation().distance(b.getLocation())));
	}	
	
	private List<GeographicPoint> searchOnWeightedGraphWithTimeAndTraffic(GeographicPoint start, GeographicPoint goal,
			Consumer<GeographicPoint> nodeSearched, BiFunction<MapNode, MapNode, Double> f) {

		// Hook for visualization.  See writeup.
		//nodeSearched.accept(next.getLocation());
		
		if (!arePreconditionsFulfilled(start, goal)) {
			return null;
		}
		
		// Initialize start and goal nodes,
		// Initialize PriorityQueue pq, visited HashSet, parent HashMap and
		// distances to infinity
		MapNode startNode = vertices.get(start);
		MapNode goalNode = vertices.get(goal);
		PriorityQueue<MapNode> nodequeue = new PriorityQueue<MapNode>();
		HashSet<MapNode> visited = new HashSet<MapNode>();
		HashMap<MapNode, MapNode> parentMap = new HashMap<MapNode, MapNode>();
		
		for (GeographicPoint gp : vertices.keySet()) {
			MapNode thisNode = vertices.get(gp);
			thisNode.setDistance(Double.POSITIVE_INFINITY);
			thisNode.setActualDistance(Double.POSITIVE_INFINITY);
		}
		
		//System.out.println("Initialization done");
		
		// Enqueue {s, 0} onto priority queue
		startNode.setDistance(0.0);
		startNode.setActualDistance(0.0);
		
		nodequeue.add(startNode);

		int numnodesvisited = 0;
		
		// while priority gueue (i.e. nodequeue) is not empty:
		while (!nodequeue.isEmpty()) {
			// dequeue node curr from front of queue
			MapNode curr = nodequeue.remove();
			numnodesvisited++;
			
			//System.out.println("At node: " + curr.toString());
		
			// if curr is not visited
			if (!(visited.contains(curr))) {
				// add curr to visited set
				visited.add(curr);
				// if curr == goal return parent map
				if (curr == goalNode) {
					nodeSearched.accept(curr.getLocation());
					break;
				}
				
				// create a hashmap of curr's neighbor's and edges time to drive edge
				HashMap<MapNode, Double> neighbors = new HashMap<MapNode, Double>();
				for (MapEdge e : curr.getEdges()) {
					SimpleDateFormat dateFormat = new SimpleDateFormat("hh aa");
					String formattedDate = dateFormat.format(new Date()).toString();
					//System.out.println(formattedDate);
					
					if (e.getRoadType().equals("primary") || e.getRoadType().equals("secondary") ||
						e.getRoadType().equals("tertiary") || e.getRoadType().equals("motorway") &&
						formattedDate.equals("04 PM") || formattedDate.equals("05 PM") ||
						formattedDate.equals("06 PM") || formattedDate.equals("07 AM") ||
						formattedDate.equals("08 AM") || formattedDate.equals("09 AM")) {
						neighbors.put(vertices.get(e.getToPoint()), 2*(e.getTime()));
					}
					else {
						neighbors.put(vertices.get(e.getToPoint()), e.getTime());
					}
				}	
				
				// for each of curr's neighbors, next, not in visited set:
				for (MapNode next : neighbors.keySet()) {
					if (!(visited.contains(next))) {
						// if path through curr to next is shorter
						// update curr as next's parent in parent map
						// and enqueue {next, distance} into the pq
						double distance = curr.getActualDistance() + neighbors.get(next);
						if (distance < next.getActualDistance()) {
							//parentMap.remove(next);???
							next.setActualDistance(distance);
							distance += f.apply(next, goalNode);
							next.setDistance(distance);
							parentMap.put(next, curr);
							nodequeue.add(next);
						}
					}
				}	
			}
		}	
		
		System.out.println("Num nodes visited is : " + numnodesvisited);
		
		// return the reconstructed path
		return reconstructPath(parentMap, startNode, goalNode);
	}

	
	public static void main(String[] args)
	{
/*		System.out.print("Making a new map...");
		MapGraph firstMap = new MapGraph();
		System.out.print("DONE. \nLoading the map...");
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", firstMap);
		System.out.println("DONE.");
		
		System.out.println("NumVertices: " + firstMap.getNumVertices());
		System.out.println("NumEdges: " + firstMap.getNumEdges());
		
		GeographicPoint testStart = new GeographicPoint(1.0, 1.0);
		GeographicPoint testEnd = new GeographicPoint(8.0, -1.0);
		
		System.out.println("Test 1 using simpletest: BFS should be 3");
		//System.out.println(firstMap.vertices.toString());
		
		List<GeographicPoint> testroute = firstMap.bfs(testStart,testEnd);
		if (testroute == null) {
			System.out.println("testroute is null");
		}
		System.out.println("Test 1 using simpletest is done: BFS should be 4 nodes, 3 hops");
		
		if (!(testroute == null)) {
			System.out.println("BFS takes " + (testroute.size()-1) + " hops.");
		}
		System.out.println("BFS takes " + (testroute.size()-1) + " hops.");
		/*for (GeographicPoint gp : testroute) {
			System.out.println(gp.getX() + gp.getY());
		}*/
		
		// You can use this method for testing.  
		
		
		/* Here are some test cases you should try before you attempt 
		 * the Week 3 End of Week Quiz, EVEN IF you score 100% on the 
		 * programming assignment.
		 */
		
/*		MapGraph simpleTestMap = new MapGraph();
		GraphLoader.loadRoadMap("data/testdata/simpletest.map", simpleTestMap);
		
		GeographicPoint testStart1 = new GeographicPoint(1.0, 1.0);
		GeographicPoint testEnd1 = new GeographicPoint(8.0, -1.0);
		
		System.out.println("Test 1 using simpletest: Dijkstra should be 9 and AStar should be 5");
		List<GeographicPoint> testroute1 = simpleTestMap.dijkstra(testStart1,testEnd1);
		List<GeographicPoint> testroute2 = simpleTestMap.aStarSearch(testStart1,testEnd1);
		
		System.out.println("Dijkstra takes " + (testroute1.size()-1) + " hops.");
		System.out.println("aStarSearch takes " + (testroute2.size()-1) + " hops.");
		
		
		MapGraph testMap = new MapGraph();
		GraphLoader.loadRoadMap("data/maps/utc.map", testMap);
		
		// A very simple test using real data
		GeographicPoint testStart2 = new GeographicPoint(32.869423, -117.220917);
		GeographicPoint testEnd2 = new GeographicPoint(32.869255, -117.216927);
		System.out.println("Test 2 using utc: Dijkstra should be 13 and AStar should be 5");
		List<GeographicPoint> testroute3 = testMap.dijkstra(testStart2,testEnd2);
		List<GeographicPoint> testroute4 = testMap.aStarSearch(testStart2,testEnd2);
		
		System.out.println("Dijkstra takes " + (testroute3.size()-1) + " hops.");
		System.out.println("aStarSearch takes " + (testroute4.size()-1) + " hops.");
		
		
		// A slightly more complex test using real data
		GeographicPoint testStart3 = new GeographicPoint(32.8674388, -117.2190213);
		GeographicPoint testEnd3 = new GeographicPoint(32.8697828, -117.2244506);
		System.out.println("Test 3 using utc: Dijkstra should be 37 and AStar should be 10");
		List<GeographicPoint> testroute5 = testMap.dijkstra(testStart3,testEnd3);
		List<GeographicPoint> testroute6 = testMap.aStarSearch(testStart3,testEnd3);
		System.out.println("Dijkstra takes " + (testroute5.size()-1) + " hops.");
		System.out.println("aStarSearch takes " + (testroute6.size()-1) + " hops.");  */
		
		
		/* Use this code in Week 3 End of Week Quiz */
		MapGraph theMap = new MapGraph();
		System.out.print("Loading the map...");
		GraphLoader.loadRoadMap("data/maps/utc.map", theMap);
		System.out.println("DONE.");

		GeographicPoint start = new GeographicPoint(32.8648772, -117.2254046);
		GeographicPoint end = new GeographicPoint(32.8660691, -117.217393);
		
		
		List<GeographicPoint> route1 = theMap.dijkstra(start,end);
		List<GeographicPoint> route2 = theMap.aStarSearch(start,end);
		List<GeographicPoint> route3 = theMap.timedDijkstraSearch(start,end);
		List<GeographicPoint> route4 = theMap.timedAStarSearch(start,end);
		List<GeographicPoint> route5 = theMap.timedWithTrafficDijkstraSearch(start,end);
		List<GeographicPoint> route6 = theMap.timedWithTrafficAStarSearch(start,end);
		
		System.out.println("Dijkstra takes " + (route1.size()-1) + " hops.");
		System.out.println("aStarSearch takes " + (route2.size()-1) + " hops.");
		System.out.println("timedDijkstraSearch takes " + (route3.size()-1) + " hops.");
		System.out.println("timedAStarSearch takes " + (route4.size()-1) + " hops.");
		System.out.println("timedWithTrafficDijkstraSearch takes " + (route5.size()-1) + " hops.");
		System.out.println("timedWithTrafficAStarSearch takes " + (route6.size()-1) + " hops.");
		
		//SimpleDateFormat dateFormat = new SimpleDateFormat("hh aa");
		//String formattedDate = dateFormat.format(new Date()).toString();
		//System.out.println(formattedDate);
		
	}
	
}
