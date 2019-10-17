package roadgraph;

import geography.GeographicPoint;

/**
 * @author Stacy Day
 * 
 * Class representing a MapEdge (or road) in our MapGraph
 *
 */

public class MapEdge {

	//private MapNode from;
	//private MapNode to;
	private GeographicPoint from;
	private GeographicPoint to;
	private String roadName;
	private String roadType;
	private double length;
	private double time;
	private int speedLimit;
	
	static final double DEFAULT_LENGTH = 0.00; // can't be < 0
	
	
	/** Create a new MapEdge object - 3 different constructors
	 * 
	 * @param from  The point at the beginning of the segment
	 * @param to  The point at the end of the segment
	 * @param roadName
	 * @param roadType
	 * @param length - DEFAULT_LENGTH = 0
	 * @param time - the time it take to drive this edge
	 * @param speedLimit
	 * 
	 */
	
	public MapEdge(GeographicPoint from, GeographicPoint to, String roadName, String roadType,
			double length) {
		this.from = from;
		this.to= to;
		this.roadName = roadName;
		this.roadType = roadType;
		//this.length = length;
		//double x = from.x - to.x;
		//double y = from.y - to.y;
		//this.length = Math.abs(Math.sqrt(x*x + y*y));
		this.length = to.distance(from);
		switch (roadType) {
			case "residential":
				this.speedLimit = 25;
				break;
			case "motorway_link":
				this.speedLimit = 35;
				break;
			case "unclassified":
				this.speedLimit = 10;
				break;
			case "primary":
				this.speedLimit = 65;
				break;
			case "secondary":
				this.speedLimit = 55;
				break;
			case "tertiary":
				this.speedLimit = 45;
				break;
			case "living_street":
				this.speedLimit = 25;
				break;
			case "motorway":
				this.speedLimit = 75;
				break;	
			default:
				this.speedLimit = 10;
			}
		this.time = speedLimit * length; 
	}
	
	// return the location of the from point
	public GeographicPoint getFromPoint()
	{
		return from;
	}
	
	// return the location of the end point
	public GeographicPoint getToPoint()
	{
		return to;
	}
	
	// return the length
	public double getTime()
	{
		return time;
	}
	
	// return the length
	public double getLength()
	{
		return length;
	}
	
	// return road name
	public String getRoadName() {
		return roadName;
	}
	
	// return road type
	public String getRoadType() {
		return roadType;
	}
	
	// print out Edge info
	public String toString() {
		return ("EDGE starts at: " + getFromPoint() + "\n ends at: " + getToPoint() + "\n is named: " +
	            getRoadName() + "\n and is of type: " + getRoadType() + "\n and is of length " + 
				getLength() + "\n");
	}


}
