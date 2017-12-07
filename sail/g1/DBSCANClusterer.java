package sail.g1;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

import sail.sim.Point;

/**
 * Implementation of density-based clustering algorithm DBSCAN.
 * 
 * Publication: 
 * Ester, Martin; Kriegel, Hans-Peter; Sander, Jörg; Xu, Xiaowei (1996). 
 * Simoudis, Evangelos; Han, Jiawei; Fayyad, Usama M., eds. 
 * A density-based algorithm for discovering clusters in large spatial 
 * databases with noise. Proceedings of the Second International Conference 
 * on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231
 * 
 * Usage:
 * - Identify type of input values.
 * - Implement metric for input value type using DistanceMetric interface.
 * - Instantiate using {#DBSCANClusterer(Collection, int, double, DistanceMetric)}.
 * - Invoke {#performClustering()}.
 * 
 * See tests and metrics for example implementation and use.
 * 
 * "mailto:cf@christopherfrantz.org Christopher Frantz
 *
 * param V: Input value type
 */
public class DBSCANClusterer {

	/** maximum distance of values to be considered as cluster */
	private double epsilon = 1f;

	/** minimum number of members to consider cluster */
	private int minimumNumberOfClusterMembers = 2;

	/** distance metric applied for clustering **/
	//private DistanceMetric<V> metric = null;

	/** internal list of input values to be clustered */
	private ArrayList<Point> inputValues = null;

	/** index maintaining visited points */
	private HashSet<Point> visitedPoints = new HashSet<Point>();

	/**
	 * Creates a DBSCAN clusterer instance. 
	 * Upon instantiation, call {link #performClustering()} 
	 * to perform the actual clustering.
	 * 
	 * param inputValues Input values to be clustered
	 * param minNumElements Minimum number of elements to constitute cluster
	 * param maxDistance Maximum distance of elements to consider clustered
	 * param metric Metric implementation to determine distance
	 * throws DBSCANClusteringException 
	 */
	public DBSCANClusterer(final ArrayList<Point> inputValues, int minNumElements, double maxDistance) {
		setInputValues(inputValues);
		setMinimalNumberOfMembersForCluster(minNumElements);
		setMaximalDistanceOfClusterMembers(maxDistance);
		//setDistanceMetric(metric);
	}

	/**
	 * Sets the distance metric
	 * 
	 * param metric
	 * throws DBSCANClusteringException 
	 */
	// public void setDistanceMetric(final DistanceMetric<V> metric) throws DBSCANClusteringException {
	// 	if (metric == null) {
	// 		throw new DBSCANClusteringException("DBSCAN: Distance metric has not been specified (null).");
	// 	}
	// 	this.metric = metric;
	// }

	/**
	 * Sets a collection of input values to be clustered.
	 * Repeated call overwrite the original input values.
	 * 
	 * param collection
	 * throws DBSCANClusteringException 
	 */
	public void setInputValues(final Collection<Point> collection){
		// if (collection == null) {
		// 	throw new DBSCANClusteringException("DBSCAN: List of input values is null.");
		// }
		this.inputValues = new ArrayList<Point>(collection);
	}

	/**
	 * Sets the minimal number of members to consider points of close proximity
	 * clustered.
	 * 
	 * param minimalNumberOfMembers
	 */
	public void setMinimalNumberOfMembersForCluster(final int minimalNumberOfMembers) {
		this.minimumNumberOfClusterMembers = minimalNumberOfMembers;
	}

	/**
	 * Sets the maximal distance members of the same cluster can have while
	 * still be considered in the same cluster.
	 * 
	 * param maximalDistance
	 */
	public void setMaximalDistanceOfClusterMembers(final double maximalDistance) {
		this.epsilon = maximalDistance;
	}

	/**
	 * Determines the neighbours of a given input value.
	 * 
	 * param inputValue Input value for which neighbours are to be determined
	 * return list of neighbours
	 * throws DBSCANClusteringException 
	 */
	private ArrayList<Point> getNeighbours(final Point inputValue) {
		ArrayList<Point> neighbours = new ArrayList<Point>();
		for(int i=0; i<inputValues.size(); i++) {
			Point candidate = inputValues.get(i);
			if (Point.getDistance(candidate,inputValue) <= epsilon) {
				neighbours.add(candidate);
			}
		}
		return neighbours;
	}

	/**
	 * Merges the elements of the right collection to the left one and returns
	 * the combination.
	 * 
	 * param neighbours1 left collection
	 * param neighbours2 right collection
	 * return Modified left collection
	 */
	private ArrayList<Point> mergeRightToLeftCollection(final ArrayList<Point> neighbours1,
			final ArrayList<Point> neighbours2) {
		for (int i = 0; i < neighbours2.size(); i++) {
			Point tempPt = neighbours2.get(i);
			if (!neighbours1.contains(tempPt)) {
				neighbours1.add(tempPt);
			}
		}
		return neighbours1;
	}

	/**
	 * Applies the clustering and returns a collection of clusters (i.e. a list
	 * of lists of the respective cluster members).
	 * 
	 * return
	 * throws DBSCANClusteringException 
	 */
	public ArrayList<ArrayList<Point>> performClustering(){

		// if (inputValues == null) {
		// 	throw new DBSCANClusteringException("DBSCAN: List of input values is null.");
		// }

		// if (inputValues.isEmpty()) {
		// 	throw new DBSCANClusteringException("DBSCAN: List of input values is empty.");
		// }

		// if (inputValues.size() < 2) {
		// 	throw new DBSCANClusteringException("DBSCAN: Less than two input values cannot be clustered. Number of input values: " + inputValues.size());
		// }

		// if (epsilon < 0) {
		// 	throw new DBSCANClusteringException("DBSCAN: Maximum distance of input values cannot be negative. Current value: " + epsilon);
		// }

		// if (minimumNumberOfClusterMembers < 2) {
		// 	throw new DBSCANClusteringException("DBSCAN: Clusters with less than 2 members don't make sense. Current value: " + minimumNumberOfClusterMembers);
		// }

		ArrayList<ArrayList<Point>> resultList = new ArrayList<ArrayList<Point>>();
		visitedPoints.clear();

		ArrayList<Point> neighbours;
		int index = 0;

		while (inputValues.size() > index) {
			Point p = inputValues.get(index);
			if (!visitedPoints.contains(p)) {
				visitedPoints.add(p);
				neighbours = getNeighbours(p);

				if (neighbours.size() >= minimumNumberOfClusterMembers) {
					int ind = 0;
					while (neighbours.size() > ind) {
						Point r = neighbours.get(ind);
						if (!visitedPoints.contains(r)) {
							visitedPoints.add(r);
							ArrayList<Point> individualNeighbours = getNeighbours(r);
							if (individualNeighbours.size() >= minimumNumberOfClusterMembers) {
								neighbours = mergeRightToLeftCollection(
										neighbours,
										individualNeighbours);
							}
						}
						ind++;
					}
					resultList.add(neighbours);
				}
			}
			index++;
		}
		return resultList;
	}

}