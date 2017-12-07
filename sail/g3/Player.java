package sail.g3;

import sail.sim.Point;
import sail.sim.Simulator;
import java.util.concurrent.TimeUnit;
import java.util.*;

public class Player extends sail.sim.Player {
    List<Point> targets;
    Map<Integer, Set<Integer>> visited_set;
    Random gen;
    int id;
    int numPlayers;
    Point initial;
    Point wind_direction;
    List<Point> prevGroupLocations;
    List<Point> groupMoves;
    List<Point> groupLocations;
    int currentTargetIdx;
    HashMap<Point, List<Point>> nnMap;
    HashMap<Point, Integer> targetNumMap;
    
    Point q1, q2, q3, q4;
    double qdiff = Math.pow(2.5, 0.5);
    
    int maxValue;
    
    Point toQuadrant = null;
	
	Point currentLocation;
    double BESTANGLE = 0.511337; 
    double BESTANGLE_UPWIND = 0.68867265359;
    Point bestDirection1;
    Point bestDirection2;
    Point bestDirection1_upwind;
    Point bestDirection2_upwind;
    double timeOnBestDirection1;
    double timeOnBestDirection2;
    double timeSpent;
    boolean upwind;
    boolean lastMoveIs1 = true;
    
    private int valueOfTarget(Point target) {
        int value = numPlayers;
        for (int p = 0; p < numPlayers; p++) {
            if (visited_set.get(p).contains(numOfTarget(target))) {
                value--;
            }
        }
        return value;
    }
    
    private void buildTargetNumMap(List<Point> targets) {
        for (int i = 0; i < targets.size(); i++) {
            this.targetNumMap.put(targets.get(i), i);
        }
    }
    
    public int numOfTarget(Point target) {
        return targetNumMap.get(target);
    }

    public void calculateNearestNeighbors() {
        for (int i = 0; i < this.targets.size(); i++) {
            Point target = targets.get(i);
            PriorityQueue<Point> neighbors = new PriorityQueue<Point>(this.targets.size() - 1, new Comparator<Point>() {
                @Override
                public int compare(Point o1, Point o2) {
                    return new Double(approximateTimeToTarget(target, o1)).compareTo(approximateTimeToTarget(target, o2));
                }
            });
            for (int j = 0; j < this.targets.size(); j++) {
                if (i == j) {
                    continue;
                }
                Point neighbor = targets.get(j);
                neighbors.add(neighbor);
            }
            List<Point> neighborsList = new ArrayList<Point>(neighbors.size());
            while (neighbors.size() > 0) {
                neighborsList.add(neighbors.poll());
            }
            this.nnMap.put(target, neighborsList);
        }
    }

    public boolean playerWillReachTargetFirst(int p, Point target) {
        int targetNum = numOfTarget(target);
        if (visited_set.get(p).contains(targetNum)) {
            return true;
        }
        if (!playerMovingToTarget(p, target)) {
            return false;
        }
        if (approximateTimeToTarget(groupLocations.get(p), target) < approximateTimeToTarget(groupLocations.get(id), target)) {
            return true;
        }
        return false;
    }

    private double approximateTimeToTarget(Point pos, Point target) {
        Point towardTarget = Point.getDirection(pos, target);
        double speed = Simulator.getSpeed(towardTarget, wind_direction);

        return Point.getDistance(pos, target) / speed;
    }

    private boolean playerMovingToTarget(int player, Point target) {
        Point playerMove = groupMoves.get(player);
        Point prevPlayerLoc = prevGroupLocations.get(player);
        Point towardTarget = Point.getDirection(prevPlayerLoc, target);

        return Math.abs(Point.angleBetweenVectors(playerMove, towardTarget)) < Math.PI / 6.0;
    }

    private int estimatedValue(Point target) {
        if (visited_set.get(id).contains(numOfTarget(target))) {
            return 0;
        }
        int value = numPlayers;
        for (int p = 0; p < numPlayers; p++) {
            if (playerWillReachTargetFirst(p, target)) {
                value--;
            }
        }
        return value;
    }
    
    private double estimatedAmortizedValue(Point from, Point target, int order, List<Point> excluding, double distance, int value) {
        if (order == 1) {
            return value / Math.pow(distance, 2.0);
        }
        List<Point> neighbors = nnMap.get(target);
        double max = 0.0;
        for (Point neighbor : neighbors) {
            if (excluding.contains(neighbor)) {
                continue;
            }
            excluding.add(neighbor);
            double neighborVal = estimatedAmortizedValue(target, neighbor, order - 1, excluding, distance + approximateTimeToTarget(target, neighbor), value + estimatedValue(neighbor));
            excluding.remove(excluding.size() - 1);
            if (neighborVal > max) {
                max = neighborVal;
            }
        }
        double UNCERTAINTY_FACTOR = 1.0;
        return max;
    }

    private double nnAdjustment(Point target, int targetNum) {
        List<Point> neighbors = nnMap.get(target);
        double total = 0.0;
        for (Point neighbor : neighbors) {
            //total += estimatedValue(neighbor) / Math.pow(approximateTimeToTarget(target, neighbor), 1.0);
        }
        return total / neighbors.size();
    }

    private List<Double> getTargetWeights() {
        List<Double> weights = new ArrayList<Double>();
        for (int i = 0; i < targets.size(); i++) {
            if (visited_set.get(id).contains(i)) {
                weights.add(-1.0);
            } else {
                Point target = targets.get(i);
                double value = estimatedAmortizedValue(groupLocations.get(id), target, 2, new ArrayList<Point>(), approximateTimeToTarget(groupLocations.get(id), target), estimatedValue(target));
                weights.add(value);
            }
        }
        return weights;
    }
    
    private boolean playerInQuadrant(int player, Point quadrant) {
        Point loc = groupLocations.get(player);
        
        double left = quadrant.x - qdiff;
        double right = quadrant.x + qdiff;
        double top = quadrant.y - qdiff;
        double bottom = quadrant.y + qdiff;
        
        return (loc.x >= left && loc.x <= right && loc.y >= top && loc.y <= bottom);
    }
    
    private int playersInQuadrant(Point quadrant) {
        int count = 0;
        for (int p = 0; p < numPlayers; p++) {
            if (playerInQuadrant(p, quadrant)) {
                count++;
            }
        }
        return count;
    }
    
    private boolean pointTowardQuadrant(Point p, Point quadrant) {
        Point loc = groupLocations.get(id);
        Point d1 = Point.getDirection(p, loc);
        Point d2 = Point.getDirection(quadrant, loc);
        double angle = Point.angleBetweenVectors(d1, d2);
        return Math.cos(angle) > Math.sqrt(0.5);
    }
    
    private Point getNextTarget() {
        if (toQuadrant != null && visited_set.get(id).contains(numOfTarget(toQuadrant))) {
            toQuadrant = null;
        }
        if (!visited_set.get(id).contains(currentTargetIdx)) {
		   return targets.get(currentTargetIdx);
			
        } else {
            int maxAvailableValue = 0;
            Point maxAvailableTarget = null;
            
            List<Point> neighbors = nnMap.get(targets.get(currentTargetIdx));
//            for (int i = 0; i < neighbors.size(); i++) {
//                if (i == 0 && toQuadrant != null) {
//                    break;
//                }
//                Point neighbor = neighbors.get(i);
//                int neighborNum = numOfTarget(neighbor);
//                int value = valueOfTarget(neighbor);
//                if (visited_set.get(id).contains(neighborNum)) {
//                    neighbors.remove(i);
//                    i--;
//                } else {
//                    if (value == numPlayers) {
//                        toQuadrant = neighbor;
//                        break;
//                    }
//                }
//            }
            
            //if (toQuadrant != null) System.out.println("Headed to (" + toQuadrant.x + ", " + toQuadrant.y + ")");
            Point closest = null;
            Point closestMax = null;
            int maxValSeen = numPlayers - 1;
            for (int i = 0; i < neighbors.size(); i++) {
                Point neighbor = neighbors.get(i);
                int neighborNum = numOfTarget(neighbor);
                int value = valueOfTarget(neighbor);
                if (visited_set.get(id).contains(neighborNum)) {
                    neighbors.remove(i);
                    i--;
                } else {
                    if (toQuadrant == null) {
                        if (maxValue <= 1) {
                            currentTargetIdx = numOfTarget(neighbor);
                            return neighbor;
                        } else {
                            if (closest == null) {
                                closest = neighbor;
                            }
                            if (valueOfTarget(neighbor) > maxValSeen) {
                                closestMax = neighbor;
                            }
                        }
                    } else {
                        if (pointTowardQuadrant(neighbor, toQuadrant)) {
                            currentTargetIdx = numOfTarget(neighbor);
                            return neighbor;
                        }
                    }
                }
            }
            if (closest == null) {
                currentTargetIdx = -1;
                return initial;
            }
            if (closestMax == null) {
                maxValue = 0;
                currentTargetIdx = numOfTarget(closest);
                return closest;
            } else if (closest != null) {
                if (approximateTimeToTarget(groupLocations.get(id), closest) < approximateTimeToTarget(groupLocations.get(id), closestMax) / 2) {
                    currentTargetIdx = numOfTarget(closest);
                    return closest;
                }   else {
                    currentTargetIdx = numOfTarget(closestMax);
                    return closestMax;
                }
            }
            currentTargetIdx = -1;
            return initial;
        }
    }

    private List<Point> calculateMoves(List<Point> initialLocations, List<Point> finalLocations) {
        ArrayList<Point> moves = new ArrayList<Point>();
        for (int i = 0; i < initialLocations.size(); i++) {
            Point init = initialLocations.get(i);
            Point fin = finalLocations.get(i);
            Point move = new Point(fin.x - init.x, fin.y - init.y);
            moves.add(move);
        }
        return moves;
    }

	private int getClosestTarget(){
		int smallestIdx = 0;

		for(int i = 1; i < targets.size(); i++){
			double time = approximateTimeToTarget(initial, targets.get(i));
			if(time < approximateTimeToTarget(initial, targets.get(smallestIdx))){
				//smallest = time;
				smallestIdx = i;
			}

		}
		return smallestIdx;
	}

    private boolean isEarlyGame() {
        int total = 0;
        for (Map.Entry<Integer, Set<Integer>> entry : visited_set.entrySet()) {
            total += entry.getValue().size();
        }
        double average = ((double) total) / numPlayers;
        return average < targets.size() * 0.25;
    }
    
    private Point getQuadrant() {
        Point loc = groupLocations.get(id);
        if (loc.x < 5.0 && loc.y < 5.0) {
            return q1;
        } else if (loc.x >= 5.0 && loc.y < 5.0) {
            return q2;
        } else if (loc.x < 5.0 && loc.y >= 5.0) {
            return q3;
        } else {
            return q4;
        }
    }
    
    private Point bestQuadrant() {
        Point best = q1;
        int num = playersInQuadrant(q1);
        
        int nq2 = playersInQuadrant(q2);
        if (nq2 < num) {
            num = nq2;
            best = q2;
        }
        int nq3 = playersInQuadrant(q3);
        if (nq3 < num) {
            num = nq3;
            best = q3;
        }
        int nq4 = playersInQuadrant(q2);
        if (nq4 < num) {
            num = nq4;
            best = q4;
        }
        
        return best;
    }
    
    private boolean inCrowdedArea() {
        int count = 0;
        for (int p = 0; p < numPlayers; p++) {
            if (p == id) {
                continue;
            }
            Point loc1 = groupLocations.get(id);
            Point loc2 = groupLocations.get(p);
            double dist = Point.getDistance(loc1, loc2);
            if (dist <= qdiff) {
                count++;
            }
        }
        return count >= 3;
    }

  @Override
  public Point chooseStartingLocation(Point wind_direction, Long seed, int t) {
      // you don't have to use seed unless you want it to
      // be deterministic (wrt input randomness)
      this.wind_direction = wind_direction;
      gen = new Random(seed);

      String temp = "speed_off_center";
      switch (temp) {
          case "geo_center" :
              initial = new Point((double) 5,(double) 5);
              break;
          case "geo_off_center" :
              initial = new Point(5.0 + gen.nextDouble(), 5.0 + gen.nextDouble());
              break;
          case "corner" :
              if(wind_direction.x < 0){
                  if(wind_direction.y < 0){
                      initial = new Point(7.5,7.5);
                  } else {
                      initial = new Point(7.5,2.5);
                  }
              } else {
                if(wind_direction.y < 0){
                    initial = new Point(2.5,7.5);
                } else {
                    initial = new Point(2.5,2.5);
                }
              }
              break;
          case "speed_off_center" :
              initial = new Point(5.0 + 2 * wind_direction.x, 5.0 -  2 * wind_direction.y);
              break;
          default :
              initial = new Point(gen.nextDouble()*10, gen.nextDouble()*10);
              break;
        }

        double speed = Simulator.getSpeed(initial, wind_direction);
        return initial;
    }

    @Override
    public void init(List<Point> group_locations, List<Point> targets, int id) {
        this.targets = targets;
        this.targetNumMap = new HashMap<Point, Integer>();
        buildTargetNumMap(targets);
        this.id = id;
        this.numPlayers = group_locations.size();
        this.maxValue = numPlayers;
        this.prevGroupLocations = group_locations;
        this.visited_set = new HashMap<>();
        for (int i = 0; i < numPlayers; i++) {
            visited_set.put(i, new HashSet<Integer>());
        }

		    this.currentTargetIdx = getClosestTarget();
        this.nnMap = new HashMap<Point, List<Point>>();
        calculateNearestNeighbors();
        this.q1 = new Point(5.0 - qdiff, 5.0 - qdiff);
        this.q2 = new Point(5.0 + qdiff, 5.0 - qdiff);
        this.q3 = new Point(5.0 - qdiff, 5.0 + qdiff);
        this.q4 = new Point(5.0 + qdiff, 5.0 + qdiff);
		this.currentLocation = group_locations.get(id);
		initializeBestSpeedVector();
    }

    @Override
    public Point move(List<Point> group_locations, int id, double dt, long time_remaining_ms) {
        this.groupMoves = calculateMoves(this.prevGroupLocations, group_locations);
        this.groupLocations = group_locations;
		this.currentLocation = group_locations.get(id);
        // testing timeouts...
        // try {
        //     TimeUnit.MILLISECONDS.sleep(1);
        // } catch(Exception ex) {
        //     ;
        // }
        // just for first turn
        prevGroupLocations = group_locations;
        if(visited_set.get(id).size() == targets.size()) {
            //this is if finished visiting all
            return Point.getDirection(group_locations.get(id), initial);
        } else{
//  				List<Double> weights = getTargetWeights();
//  				int highestIdx = 0;
//  				for(int i = 0; i < weights.size(); i++){
//  					if(weights.get(i) > weights.get(highestIdx)){
//  						highestIdx = i;
//            }
//				  }
//				  currentTargetIdx = highestIdx;
//
            if (toQuadrant != null) {
                if (!isEarlyGame()) {
                    toQuadrant = null;
                } else if (Point.getDistance(groupLocations.get(id), toQuadrant) < qdiff / 2) {
                    toQuadrant = null;
                }
            } else if (isEarlyGame()) {
                Point quadrant = getQuadrant();
                if (inCrowdedArea()) {//playersInQuadrant(quadrant) > numPlayers / 4) {
                    //System.out.println("Its crowded!");
                    //toQuadrant = bestQuadrant();
                }
            }
          return moveToTarget(group_locations, dt);
        }
    }

	private Point moveToTarget(List<Point> group_locations, double dt){
		Point pos;
		if(groupLocations == null){
			pos = initial;
		}else{
			pos = groupLocations.get(id);
    }
        Point target = targets.get(currentTargetIdx);
        if (visited_set.get(id).contains(currentTargetIdx)) {
            target = getNextTarget();
        }
		
		return computeNextDirection(target, dt);
	}

    /**
    * visited_set.get(i) is a set of targets that the ith player has visited.
    */
    @Override
    public void onMoveFinished(List<Point> group_locations, Map<Integer, Set<Integer>> visited_set) {
        this.visited_set = visited_set;
    }

	/*below are methods to aid in movement borrowed from g1*/
		
	
	private double getSpeedRelativeToWind(double angle) {
        double x = 2.5 * Math.cos(angle + Math.PI) - 0.5;
        double y = 5 * Math.sin(angle + Math.PI);
        return Math.sqrt((x)*(x) + (y)*(y));
    }
	
	private void initializeBestSpeedVector() {
        double windAngle = Math.atan2(this.wind_direction.y, this.wind_direction.x);
        
        double speed = getSpeedRelativeToWind(BESTANGLE);
        this.bestDirection1 = new Point(speed * Math.cos(windAngle + BESTANGLE), speed * Math.sin(windAngle + BESTANGLE));
        this.bestDirection2 = new Point(speed * Math.cos(windAngle - BESTANGLE), speed * Math.sin(windAngle - BESTANGLE));
        
        speed = getSpeedRelativeToWind(Math.PI + BESTANGLE_UPWIND);
        this.bestDirection1_upwind = new Point(speed * Math.cos(windAngle + Math.PI + BESTANGLE_UPWIND), speed * Math.sin(windAngle + Math.PI + BESTANGLE_UPWIND));
        this.bestDirection2_upwind = new Point(speed * Math.cos(windAngle + Math.PI - BESTANGLE_UPWIND), speed * Math.sin(windAngle + Math.PI - BESTANGLE_UPWIND));
    }

	private Point alternateBetween1And2(List<Point> group_locations, double dt) {
        this.currentLocation = group_locations.get(id);
        Point move;
        Point newLocation = new Point(bestDirection1.x * dt + group_locations.get(id).x, bestDirection1.y * dt + group_locations.get(id).y);
        Point newLocation2 = new Point(bestDirection2.x * dt + group_locations.get(id).x, bestDirection2.y * dt + group_locations.get(id).y);

        if (lastMoveIs1 && upwind) {
            timeOnBestDirection2 -= dt;
            move = bestDirection2_upwind;
        }
        else if (lastMoveIs1 && !upwind) {
            timeOnBestDirection2 -= dt;
            move = bestDirection2;
        }
        else if (!lastMoveIs1 && upwind) {
            timeOnBestDirection1 -= dt;
            move = bestDirection1_upwind;
        }
        else {
            timeOnBestDirection1 -= dt;
            move = bestDirection1;
        }
        lastMoveIs1 = !lastMoveIs1;
        return move;
    }
	
	private Point computeExpectedPosition(Point moveDirection, double dt) {
        Point unitMoveDirection = Point.getUnitVector(moveDirection);
        double speed = Simulator.getSpeed(unitMoveDirection, this.wind_direction);
        Point distanceMoved = new Point(
                unitMoveDirection.x * speed * dt,
                unitMoveDirection.y * speed * dt
        );
        Point nextLocation = Point.sum(this.currentLocation, distanceMoved);
        if (nextLocation.x < 0 || nextLocation.y > 10 || nextLocation.y < 0 || nextLocation.x > 10) {
            return this.currentLocation;
        }
        return nextLocation;
    }

	private Point computeNextDirection(Point target, double dt) {
		//System.out.println(target.x+" "+target.y);
        Point directionToTarget = Point.getDirection(this.currentLocation, target);
        Point perpendicularLeftDirection = Point.rotateCounterClockwise(directionToTarget, Math.PI/2.0);
        Point perpendicularRightDirection = Point.rotateCounterClockwise(directionToTarget, -Math.PI/2.0);
        return findBestDirection(perpendicularLeftDirection, perpendicularRightDirection, target, 100, dt);
    }
	
	private Point findBestDirection(Point leftDirection, Point rightDirection, Point target, int numSteps,
                                    double dt) {
        // First, check if the target is reachable in one time step from our current position.
        double currentDistanceToTarget = Point.getDistance(this.currentLocation, target);
        Point directionToTarget = Point.getDirection(this.currentLocation, target);
        double speedToTarget = Simulator.getSpeed(directionToTarget, this.wind_direction);
        double distanceWeCanTraverse = speedToTarget * dt;
        double distanceTo10MeterAroundTarget = currentDistanceToTarget-0.01;
        if (distanceWeCanTraverse > distanceTo10MeterAroundTarget) {
            return directionToTarget;
        }

        // If that is not the case, choose the direction that will get us to a point where the time to reach
        // the target is minimal, if going directly to the target.
        double minTimeToTarget = Double.POSITIVE_INFINITY;
        Point minTimeToTargetDirection = null;

        double totalRadians = Point.angleBetweenVectors(leftDirection, rightDirection);
        double radiansStep = totalRadians / (double) numSteps;
        for (double i = 0.0; i < totalRadians; i+=radiansStep) {
            Point direction = Point.rotateCounterClockwise(rightDirection, i);
            Point expectedPosition = computeExpectedPosition(direction, dt);
            Point nextDirection = Point.getDirection(expectedPosition, target);
            double distance = Point.getDistance(expectedPosition, target);
            double speed = Simulator.getSpeed(nextDirection, this.wind_direction);
            double timeToTarget = distance / speed;

            if (timeToTarget < minTimeToTarget) {
                minTimeToTargetDirection = direction;
                minTimeToTarget = timeToTarget;
            }
        }

        return minTimeToTargetDirection;
    }
}