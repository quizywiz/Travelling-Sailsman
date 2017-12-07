package sail.g1;

import sail.sim.Point;
import sail.sim.Simulator;
import java.util.concurrent.TimeUnit;
import java.util.*;

public class Player extends sail.sim.Player {
    private int K;

    int t;
    List<Point> targets;
    Map<Integer, Set<Integer>> visited_set;
    int[] targetScore;
    double[] estimatedTargetScores;


    Random gen;
    int id;
    Point wind_direction;
    Point currentLocation;
    Point initial;

    int curIndex;
    List<Integer> nextTargetIndexes;
    int nextTargetIndex;
    int skip;

    int turn_num;

    List<Point> groupLocations;

    double simThreshold;
    int simLookahead;

    double BESTANGLE = 0.511337; 
    double BESTANGLE_UPWIND = 0.45378265359;
    Point bestDirection1;
    Point bestDirection2;
    Point bestDirection1_upwind;
    Point bestDirection2_upwind;
    double timeOnBestDirection1;
    double timeOnBestDirection2;
    double timeSpent;
    boolean upwind;
    boolean lastMoveIs1 = true;

    HashMap<Integer, Point[]> playerMoves = new HashMap<Integer, Point[]>();
    int curSeqenceLength;

    Map<Integer, List<Integer>> bestPathFromTarget;

    ArrayList<ArrayList<Point>> clusters;
    HashMap<Point, Integer> pointToTargetIndex = new HashMap<>();
    int currentCluster = -1;
    int nextCluster = -1;

    @Override
    public Point chooseStartingLocation(Point wind_direction, Long seed, int t) {
        // you don't have to use seed unless you want it to 
        // be deterministic (wrt input randomness)
        this.t = t;
        gen = new Random(seed);
        if(t <= 5) {
            initial = initStartPoint(wind_direction, 2.0, 5.0, 5.0);
        } else if(t <= 100) {
            initial = initStartPoint(wind_direction, 4.0, 5.0, 5.0);
        } else {
            initial = initStartPoint(wind_direction, 2.0, 5.0, 5.0);
        }
        this.wind_direction = wind_direction;
        double speed = Simulator.getSpeed(initial, wind_direction);
        return initial;
    }

    private Point initStartPoint(Point wind_direction, Double shift, Double xstart, Double ystart){
        Point origin = new Point(1, 0);
        Double angle = Point.angleBetweenVectors(origin, wind_direction);
        Double x, y;
        x = -shift * Math.cos(angle) + xstart;
        y = -shift * Math.sin(angle) + ystart;
        Point start = new Point(x, y);
        return start;
    }

    @Override
    public void init(List<Point> group_locations, List<Point> targets, int id) {
        this.targets = targets;
        this.id = id;
        this.nextTargetIndex = -1;
        if(t <= 20) {
            simThreshold = 1.0;
            this.K = Math.min(1, t);
        } else if(t <= 50) {
            this.K = Math.min(1, t);
            simThreshold = 0.7;
            simLookahead = 4;
        } else if(t <= 200) {
            simThreshold = 1.0;
            simLookahead = 4;
            this.K = Math.min(1, t);
        } else {
            simThreshold = 0.7;
            simLookahead = 5;
            this.K = Math.min(1, t);
        }

        this.groupLocations = group_locations;

        targetScore = new int[targets.size()];
        estimatedTargetScores = new double[targets.size()];
        for(int i = 0; i < targetScore.length; i++) {
            targetScore[i] = group_locations.size();
        }

        this.bestPathFromTarget = new HashMap<>();
        curIndex = -1;
        skip = -1;
        this.currentLocation = group_locations.get(id);
        initializeBestSpeedVector();
        for (int i=0; i<group_locations.size(); i++){
            printPoint(group_locations.get(i));
            playerMoves.put(i, new Point[] {null, group_locations.get(i)});
        }


        for (int i=0; i<targets.size();i++) {
            pointToTargetIndex.put(targets.get(i),i);
        }
        
        ArrayList<Point> temp = new ArrayList<>();
        temp.addAll(targets);
        DBSCANClusterer clusterer = new DBSCANClusterer(temp, 1, .3);
        clusters = clusterer.performClustering();
        for(ArrayList<Point> pts : clusters) {
            System.out.println(pts.size());
        }
       /* printMap();
        System.out.println("targets: ");
        for (Point t : targets){
            printPoint(t);
        }*/
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

    private double getSpeedRelativeToWind(double angle) {
        double x = 2.5 * Math.cos(angle + Math.PI) - 0.5;
        double y = 5 * Math.sin(angle + Math.PI);
        return Math.sqrt((x)*(x) + (y)*(y));
    }

    @Override
    public Point move(List<Point> group_locations, int id, double dt, long time_remaining_ms) {
        this.currentLocation = group_locations.get(id);
        if (this.t < 10) {
            return christofideMove(group_locations, id, dt, time_remaining_ms);
        }
        this.groupLocations = group_locations;
        for (int i=0; i<group_locations.size(); i++){
            playerMoves.put(i, new Point[] {playerMoves.get(i)[1], group_locations.get(i)});
            if (i == id) continue;
            
            ArrayList<Point> predictList = predictNextTarget(i, playerMoves.get(i)[0], playerMoves.get(i)[1]);
            if (nextTargetIndexes != null && curIndex != -1 && curIndex < curSeqenceLength) {
                if (nextTargetIndexes.get(curIndex) == targets.size()) break;
                Point next = targets.get(nextTargetIndexes.get(curIndex));
                double oppTime = getTimeToTravel(next, group_locations.get(i));
                double ourTime = getTimeToTravel(next, group_locations.get(id)); 
                if (predictList.contains(next) && oppTime < ourTime) {
                    skip = nextTargetIndexes.get(curIndex);
                } else if (predictList.contains(next) && oppTime >= ourTime) {
                    skip = -1;
                }
            } else {
                for(int j=0; j<predictList.size(); j++) {
                    if (visited_set != null && !visited_set.get(this.id).contains(predictList.get(j))) {
                        double oppTime = getTimeToTravel(predictList.get(j), group_locations.get(i));
                        double ourTime = getTimeToTravel(predictList.get(j), group_locations.get(id));
                        if (ourTime < oppTime) {
                            nextTargetIndexes = selectKBest(predictList.get(j));
                            curIndex = 0;
                        }
                    } 
                }
            }
        }

        if (skip != -1) {
            nextTargetIndexes = selectKBest(group_locations.get(id));
            curIndex = 0;
            skip = -1;
        }
        
      //  System.out.println("---------------");
        return nearestNeighborMove(group_locations, id, dt, time_remaining_ms);
    }

    public Point christofideMove(List<Point> group_locations, int id, double dt, long time_remaining_ms) {
        boolean justHitATarget = false;
        if (visited_set != null && visited_set.get(id).size() == targets.size()) {
            return moveHelper(group_locations, initial, dt, justHitATarget);
        }

        if (currentCluster != -1 && curIndex != -1 && visited_set != null && 
                visited_set.get(id).contains(pointToTargetIndex.get(clusters.get(currentCluster).get(nextTargetIndexes.get(curIndex))))) {
            curIndex++;
            justHitATarget = true;
        }

        if (currentCluster == -1 || curIndex >= nextTargetIndexes.size()) {
            if (currentCluster != -1) {
                clusters.remove(currentCluster);
            }

            currentCluster = selectCluster(group_locations.get(id));
            if(clusters.size() > 1) {
                nextCluster = selectCluster(clusters.get(currentCluster).get(0));
            } else {
                nextCluster = -1;
            }
            nextTargetIndexes = dijkstra(clusters.get(currentCluster), nextCluster);
            curIndex = 0;
        }

       return moveHelper(group_locations, clusters.get(currentCluster).get(nextTargetIndexes.get(curIndex)), dt, justHitATarget);
    }

    private ArrayList<Integer> dijkstra(ArrayList<Point> cluster, int nextCluster) {
        if(cluster.size() == 1) {
            ArrayList<Integer> temp = new ArrayList<>();
            temp.add(0);
            return temp;
        }
        
            cluster.add(0,currentLocation);
        double[][] weight = new double[cluster.size()][cluster.size()];
        for(int i=0; i<cluster.size(); i++) {
            for(int j=0; j<cluster.size(); j++) {
                if(i == j) {
                    weight[i][j] = -1;
                } else {
                    if(i!=0 && j!=0) {
                    weight[i][j] = (targetScore[pointToTargetIndex.get(cluster.get(i))]+targetScore[pointToTargetIndex.get(cluster.get(j))])
                        /getTimeToTravel(cluster.get(i), cluster.get(j));
                    } else {
                        weight[i][j] = 1/getTimeToTravel(cluster.get(i), cluster.get(j));
                    }
                }
            }
        }
        Christofides ch = new Christofides(false, this.wind_direction);
        int [] shortestPath = ch.solve(weight);
        ArrayList<Integer> temp= new ArrayList<>();
        for(int test : shortestPath) {
            temp.add(test-1);
        }
        temp.remove(0);
        cluster.remove(currentLocation);
        return temp;
        /*  Point target = nextCluster==-1 ? initial : clusters.get(nextCluster).get(0);
        cluster.add(target);
        cluster.add(0, currentLocation);

        double[][] weight = new double[cluster.size()][cluster.size()];
        for(int i=0; i<cluster.size(); i++) {
            for(int j=0; j<cluster.size(); j++) {
                if(i == j) {
                    weight[i][j] = -1;
                } else {
                    weight[i][j] = getTimeToTravel(cluster.get(i), cluster.get(j));
                }
            }
        }

        ArrayList<Integer> possible = new ArrayList<>();
        double[] distances = new double[cluster.size()];
        int[] previous = new int[cluster.size()];

        for(int i=0; i<cluster.size(); i++) {
            possible.add(i);
            distances[i] = weight[0][i];
            if(distances[i] != -1) {
                previous[i] = 0;
            }
        }

        for(int i=0; i<cluster.size()-1; i++) {
            double l = Integer.MAX_VALUE;
            int p = 0;
            for(int j=0; j<possible.size(); j++) {
                if(distances[possible.get(j)] < l) {
                    p = j;
                    l = distances[possible.get(j)];
                }
            }
            possible.remove(p);

            for(int j=0; j<cluster.size()-1; j++) {
                if(distances[p] != -1 && weight[p][j] != -1 && distances[p]+weight[p][j] < distances[j]) {
                    distances[j] = distances[p] + weight[p][j];
                    previous[j] = p;
                }
            }
        }
        possible.clear();

        int loc = cluster.size()-1;
        while(cluster.get(previous[loc]) != cluster.get(0)) {
            possible.add(0, previous[loc]);
            loc = previous[loc];
        }
        cluster.remove(loc);
        cluster.remove(0);

        return possible;*/
    }

    private ArrayList<Point> predictNextTarget(int playerId, Point prev, Point current){
        ArrayList<Point> predictList = new ArrayList<>();
        if (prev != current){
            double slope = (current.y-prev.y) / (current.x-prev.x);
            double y, x;
            y = current.y;
            x = current.x;
            for (Point a : targets){
                if (visited_set.get(playerId).contains(a)) continue;
                double y_delta = a.y-y;
                double x_delta = a.x-x;
                double delta = y_delta/x_delta;
                double slope_variation = delta/slope;
                if (slope_variation >=0.9 && slope_variation <= 1.1){
                    double prev_y_delta = current.y-prev.y;
                    if ((prev_y_delta < 0 && y_delta < 0) || (prev_y_delta > 0 && y_delta > 0))
                        predictList.add(a);
                }
            }
        }

     /*   for (Point a : predictList){
            printPoint(a);
        }*/

        return predictList;
    }

    public int selectCluster(Point curLoc) {
        double maxScore = -1;
        int selectedCluster = 0;
        for(int i = 0; i < clusters.size(); i++) {
            double totalScore = 0;
            for(int j = 0; j < clusters.get(i).size(); j++) {
                totalScore += 1/getTimeToTravel(currentLocation, clusters.get(i).get(j));
                if(clusters.get(i).get(j).equals(curLoc)) {
                    totalScore = -1;
                    break;
                }
            } 
            if(clusters.get(i).size() > 0) totalScore /= clusters.get(i).size();
            if(totalScore > maxScore) {
                maxScore = totalScore;
                selectedCluster = i;
            }
        }

        return selectedCluster;
    }

    public Point nearestNeighborMove(List<Point> group_locations, int id, double dt, long time_remaining_ms) {
        turn_num++;
        boolean justHitATarget = false;
        while(curIndex < curSeqenceLength && visited_set != null && nextTargetIndexes != null && visited_set.get(id).contains(nextTargetIndexes.get(curIndex))) {
            curIndex++;
            justHitATarget = true;
        }

        if(nextTargetIndexes == null || curIndex == curSeqenceLength || curIndex == -1) {
            //System.out.println("Updating");
            for(int i = 0; i < targetScore.length; i++) {
                estimatedTargetScores[i] = targetScore[i];
                adjustScore(i);
                //System.out.println(estimatedTargetScores[i]);
            }
            nextTargetIndexes = selectKBest(group_locations.get(id));
            curIndex = 0;
            justHitATarget = true;
        }

        if(nextTargetIndexes.get(curIndex) == targets.size()) {
            return moveHelper(group_locations, initial, dt, justHitATarget);
        } else {
            return moveHelper(group_locations, targets.get(nextTargetIndexes.get(curIndex)), dt, justHitATarget);
        }
    }

    private Point moveHelper(List<Point> group_locations, Point target, double dt, boolean justHitATarget) {
        if (dt <= 0.004) {
            if (justHitATarget) {
                double x1 = group_locations.get(id).x;
                double y1 = group_locations.get(id).y; 
                double x = target.x;
                double y = target.y;
                Point newTarget = new Point(x,y);
                Point direction = Point.getDirection(this.currentLocation, newTarget);

                double theta = Point.angleBetweenVectors(direction, wind_direction);
    
                if (theta <= BESTANGLE ||
                    theta >= -BESTANGLE + 2*Math.PI) {
                    timeOnBestDirection1 = (direction.y - bestDirection2.y/bestDirection2.x * direction.x)/(bestDirection1.y - bestDirection2.y/bestDirection2.x * bestDirection1.x);
                    timeOnBestDirection2 = (direction.y - bestDirection1.y/bestDirection1.x * direction.x)/(bestDirection2.y - bestDirection1.y/bestDirection1.x * bestDirection2.x);
                    this.upwind = false;
                } else if (theta >= Math.PI - BESTANGLE_UPWIND &&
                    theta <= Math.PI + BESTANGLE_UPWIND) {
                    timeOnBestDirection1 = (direction.y - bestDirection2_upwind.y/bestDirection2_upwind.x * direction.x)/(bestDirection1_upwind.y - bestDirection2_upwind.y/bestDirection2_upwind.x * bestDirection1_upwind.x);
                    timeOnBestDirection2 = (direction.y - bestDirection1_upwind.y/bestDirection1_upwind.x * direction.x)/(bestDirection2_upwind.y - bestDirection1_upwind.y/bestDirection1_upwind.x * bestDirection2_upwind.x);
                    this.upwind = true;
                } 
                if (upwind) {
                    return bestDirection1_upwind;
                }
                else {
                    return bestDirection1;
                }
            }
            else {
                if (timeOnBestDirection1 > dt && timeOnBestDirection2 > dt) {
                    return alternateBetween1And2(group_locations, dt);
                }
                else if (timeOnBestDirection1 > dt) {
                    timeOnBestDirection1 -= dt;
                    return this.upwind ? bestDirection1_upwind : bestDirection1;
                }
                else if (timeOnBestDirection2 > dt) {
                    timeOnBestDirection2 -= dt;
                    return this.upwind ? bestDirection2_upwind : bestDirection2;
                }
                else {
                    return computeNextDirection(target, dt);
                }
            }
        }
        else {
            return computeNextDirection(target, dt);
        }
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

    private Point computeNextDirection(Point target, double dt) {
        Point directionToTarget = Point.getDirection(this.currentLocation, target);
        Point perpendicularLeftDirection = Point.rotateCounterClockwise(directionToTarget, Math.PI/2.0);
        Point perpendicularRightDirection = Point.rotateCounterClockwise(directionToTarget, -Math.PI/2.0);
        return findBestDirection(perpendicularLeftDirection, perpendicularRightDirection, target, 100, dt);
    }

    private Point findBestDirection(Point leftDirection, Point rightDirection, Point target, int numSteps,
                                    double dt) {
        double currentDistanceToTarget = Point.getDistance(this.currentLocation, target);
        Point directionToTarget = Point.getDirection(this.currentLocation, target);
        double speedToTarget = Simulator.getSpeed(directionToTarget, this.wind_direction);
        double distanceWeCanTraverse = speedToTarget * dt;
        double distanceTo10MeterAroundTarget = currentDistanceToTarget-0.01;
        if (distanceWeCanTraverse > distanceTo10MeterAroundTarget) {
            return directionToTarget;
        }
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

    public List<Integer> selectKBest(Point curLoc) {
        double maxScore = -1;
        List<Integer> bestSeq = new ArrayList<>();
        for(int i = 1; i <= K; i++) {
            List<Integer> seq = selectKBest(curLoc, i);
            double score = getTargetsScore(curLoc, seq);
            //System.out.println(i + ": " + score);
            if(score > maxScore) {
                bestSeq = seq;
                curSeqenceLength = i;
                maxScore = score;
            }
        }

        //System.out.println(curSeqenceLength);
        return bestSeq;
    }

    public List<Integer> selectKBest(Point curLoc, int L) {
        double maxScore = -1;
        int targetIndex = -1;
        List<Integer> p = new ArrayList<Integer>();
        List<Integer> maxP = new ArrayList<Integer>();
        if(visited_set != null && visited_set.get(this.id).size() >= targets.size() - K) {
            int nn = findNearestNeighbor(curLoc);
            maxP.add(nn);
            /*
            for(int j = 0; j < L; j++) {
                maxP.add(nn);
            }
            */
            return maxP;
        }

        for(int i = 0; i < targets.size(); i++) {
            if(visited_set != null && visited_set.get(this.id).contains(i) || i == skip) continue;

            p = findNearestNeighbors(i, L);

            // Compute time to get to target
            double score = getTargetsScore(curLoc, p);
            //System.out.println(score);
            if(score > maxScore) {
                maxScore = score;
                targetIndex = i;
                maxP = p;
            } else if(score > simThreshold*maxScore) {
                /*
                int k = 0;
                if(visited_set != null) {
                    k = (t - visited_set.get(this.id).size())/(t/4) + 1;
                } else {
                    k = 4;
                }
                */
                //int k = (t - visited_set.get(this.id).size())/100 + 1;
                //System.out.println(k);
                List<Integer> p1 = findNearestNeighbors(i, simLookahead);
                List<Integer> p2 = findNearestNeighbors(targetIndex, simLookahead);
                double s1 = getTargetsScore(curLoc, p1);
                double s2 = getTargetsScore(curLoc, p2);
                if(s1 > s2) {
                    maxScore = score;
                    targetIndex = i;
                    maxP = p;
                }
            }
        }


        if(targetIndex == -1) {
            targetIndex = targets.size();
        }

        //System.out.println(p.size());
        return maxP;
    }

    public List<Integer> findNearestNeighbors(int targetIndex, int L) {
        boolean cached = true;
        if(bestPathFromTarget.containsKey(targetIndex) && bestPathFromTarget.get(targetIndex).size() >= L) {
            for(int r : bestPathFromTarget.get(targetIndex)) {
                if(visited_set != null && visited_set.get(this.id).contains(r)) {
                    cached = false;
                }
            }
        } else {
            cached = false;
        }


        if(!cached) {
            List<Integer> nextPoints = new ArrayList<>();
            nextPoints.add(targetIndex);
            Set<Integer> excluded = new HashSet<Integer>();
            excluded.add(targetIndex);
            for(int i = 1; i < Math.max(K, L); i++) {
                int nxt = findNearestTarget(nextPoints.get(i-1), excluded);
                nextPoints.add(nxt);
                excluded.add(nxt);
            }

            bestPathFromTarget.put(targetIndex, nextPoints);
        }

        return bestPathFromTarget.get(targetIndex).subList(0, L);

    }

     private int findNearestTarget(int targetIndex, Set<Integer> excluded) {
        int nearestTargetIndex = -1;
        double maxScore = -1;
        Point targetLoc = targets.get(targetIndex);
        for(int i = 0; i < targets.size(); i++) {
            if(excluded.contains(i) || visited_set != null && visited_set.get(this.id).contains(i)) continue;

            // Compute time to get to target
            double score = getTargetScore(targetLoc, i);
            if(score > maxScore) {
                maxScore = score;
                nearestTargetIndex = i;
            }
        }

        // No more neighbors so return to starting point
        if(nearestTargetIndex == -1) {
            nearestTargetIndex = targets.size();
        }

        return nearestTargetIndex;
    }



    private int findNearestNeighbor(Point curLoc) {
        Point nearestTarget = null;
        int targetIndex = -1;
        double maxScore = -1;
        for(int i = 0; i < targets.size(); i++) {
            if(visited_set != null && visited_set.get(this.id).contains(i)) continue;

            // Compute time to get to target
            double score = getTargetScore(curLoc, i);

            if(score > maxScore) {
                nearestTarget = targets.get(i);
                maxScore = score;
                targetIndex = i;
            }
        }

        // No more neighbors so return to starting point
        if(nearestTarget == null) {
            nearestTarget = initial;
            targetIndex = targets.size();
        }
        
        return targetIndex;
    }

    private double getTimeToTravel(Point from, Point to) {
        return Point.getDistance(from, to)/Simulator.getSpeed(Point.getDirection(from, to), wind_direction);
    }

    private double getTargetScore(Point from, int targetIndex) {
        double travelTime = getTimeToTravel(from, targets.get(targetIndex));
        //System.out.println("time: " + travelTime);
        //System.out.println("score: " + targetScore[targetIndex]);
        return estimatedTargetScores[targetIndex]/Math.pow(travelTime, 1);
    }

    private double getTargetScore(Point from, Point to) {
        double travelTime = getTimeToTravel(from, to);
        return targetScore[pointToTargetIndex.get(to)]/travelTime;
    }

    private double getTargetsScore(Point from, List<Integer> targetIndexes) {
        if(targetIndexes.get(0) == t) return 0;
        double travelTime = getTimeToTravel(from, targets.get(targetIndexes.get(0)));
        double score = targetScore[targetIndexes.get(0)];

        for(int i = 1; i < targetIndexes.size(); i++) {
            travelTime += getTimeToTravel(targets.get(targetIndexes.get(i-1)), targets.get(targetIndexes.get(i)));
            score += targetIndexes.get(i) == t ? 0 : estimatedTargetScores[i];
        }

        //System.out.println("time: " + travelTime);
        //System.out.println("score: " + score);
        return score/(Math.pow(travelTime, 1));
    }

    /**
    * visited_set.get(i) is a set of targets that the ith player has visited.
    */
    @Override
    public void onMoveFinished(List<Point> group_locations, Map<Integer, Set<Integer>> visited_set) {
        this.groupLocations = group_locations;
        this.visited_set = visited_set;
        for(int i = 0; i < targetScore.length; i++) {
            targetScore[i] = group_locations.size();
        }

        for(int k : visited_set.keySet()) {
            Set<Integer> visitedTargets = visited_set.get(k);
            for(int l : visitedTargets) {
                targetScore[l]--;
            }
        }
    }

    private void adjustScore(int targetIndex) {
        for(int pid = 0; pid < groupLocations.size(); pid++) {
            double ourTime = getTimeToTravel(groupLocations.get(this.id), targets.get(targetIndex));
            double oppTime = getTimeToTravel(groupLocations.get(pid), targets.get(targetIndex));
            if(oppTime < ourTime) {
                if(this.visited_set == null || !this.visited_set.get(pid).contains(targetIndex)) {
                    estimatedTargetScores[targetIndex] -= (ourTime - oppTime)/ourTime;
                }
            }
        }
    }

    private double computeUnvisitedPlayersTimeTo(int targetIndex) {
        double time = 1.0;
        Point target = this.targets.get(targetIndex);
        for (int pid = 0; pid < this.groupLocations.size(); pid++) {
            if (pid == this.id) continue; // Skip our own.
            if (!this.visited_set.get(pid).contains(targetIndex)) {
                // This means that this player hasn't visited this target yet, so
                // compute her time to target.
                Point playerLoc = groupLocations.get(pid);
                time += getTimeToTravel(playerLoc, target);
            }
        }
        return time == 0 ? 1 : time;
    }

    private double getAngleWithTarget(int targetIndex) {
        return Math.min(Point.angleBetweenVectors(Point.getDirection(playerMoves.get(this.id)[1], playerMoves.get(this.id)[0]), Point.getDirection(targets.get(targetIndex), playerMoves.get(this.id)[1])),
            Point.angleBetweenVectors(Point.getDirection(targets.get(targetIndex), playerMoves.get(this.id)[1]), Point.getDirection(playerMoves.get(this.id)[1], playerMoves.get(this.id)[0])));
    }

    private void printPoint(Point a){
        System.out.println("(" + a.x + "," + a.y + ")");
    }

    private void printPointArray(Point prev, Point current){
        if (prev != null && current != null)
            System.out.println("[(" + prev.x + "," + prev.y + ") , (" + current.x + "," + current.y + ")");
        else if (prev==null)
            System.out.println("[(" + "null) , (" + current.x + "," + current.y + ")");
    }

    private void printMap(){
        for (Map.Entry<Integer, Point[]> entry : playerMoves.entrySet()) {
            Integer key = entry.getKey();
            Point[] value = entry.getValue();
            System.out.print(key + ": ");
            printPointArray(value[0], value[1]);
        }
    }
}
