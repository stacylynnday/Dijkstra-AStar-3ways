package roadgraph;

import java.util.HashSet;
import java.util.Set;

import geography.GeographicPoint;

/**
 * MapNode.java
 * 
 * @author Stacy Day
 * 
 * A class to represent a vertice (or Node or Intersection) in a graph 
 * which is a geographical Map that a car navigates.
 */

public class MapNode implements Comparable<MapNode> {
	// Since our map is always a geographical map, nodes keep track of 
	// location as GeographicPoint
	private GeographicPoint location;
	private HashSet<MapEdge> nodeedges;
	private double distance;
	private double predicteddistance;
	
	// default MapNode constructor
	public MapNode(GeographicPoint location)
	{
		this.location = location;
		nodeedges = new HashSet<MapEdge>();
		distance = Double.POSITIVE_INFINITY;
		predicteddistance = Double.POSITIVE_INFINITY;
	}

	// public void??
	public boolean addEdge(MapEdge edge) 
	{
		if (nodeedges.contains(edge) || edge == null) {
			return false;
		}
		else {
			nodeedges.add(edge);
			return true;
		}
	}
	
	// get the edges out of this node
	public Set<MapEdge> getEdges()
	{
		return nodeedges;
	}
	
	// get the location of a node
	public GeographicPoint getLocation()
	{
		return location;
	}
		
	// get the neighbors of this MapNode
/*	public Set<MapNode> getNeighbors() {
		Set<MapNode> neighbors = new HashSet<MapNode>();
		for (MapEdge e : nodeedges) {
			MapNode neighbor = vertices.;
			neighbors.add(neighbor);
		}
		return neighbors;
	}*/

	// predicted distance (infinity for dijkstra, straight line
	// distance for AStar
	public double getDistance() {
		return this.distance;
	}
	
	// predicted distance (infinity for dijkstra, straight line
	// distance for AStar
	public boolean setDistance(double distance) {
		this.distance = distance;
		return true;
	}
	
	public double getActualDistance() {
		return this.predicteddistance;
	}
	
	public boolean setActualDistance(double distance) {
		this.predicteddistance = distance;
		return true;
	}
	
	public String toString() {
		String roadnames = "";
		for (MapEdge e : nodeedges) {
			roadnames += e.getRoadName();
			roadnames += ", ";
		}
		return ("NODE: " + getLocation() + "\n is connected to " + nodeedges.size() 
			+ " streets: " + roadnames + "\n");
	}

	@Override
	public int compareTo(MapNode node) {
		// TODO Auto-generated method stub
		return ((Double)this.getDistance()).compareTo((Double)node.getDistance());
	}

}