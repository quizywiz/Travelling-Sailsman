package sail.g4;

import sail.sim.Point;
import sail.sim.Simulator;
import java.util.concurrent.TimeUnit;
import java.util.*;

public class Player extends sail.sim.Player {
    
    public final double WINDSPEED = 50.0;
    List<Point> targets;
    Map<Integer, Set<Integer>> visited_set;
    Random gen;
    int id;
    Point initial;
    int numPlayers = 0;
    Point wind_direction;
    double wind_angle;
    boolean inStraightMode = false;
    SailingHelper sHelper;
    List<Point> nextKpoints;
    int k;
    int pointsVisited = 0;

    public enum InitMode { 
        CENTER, WIND 
    }
    InitMode initMode = InitMode.WIND;

    public enum PlayMode { 
        ONE_STEP, K_STEPS, TRAVEL
    }
    PlayMode playMode = PlayMode.ONE_STEP;

    @Override
    public Point chooseStartingLocation(Point wind_direction, Long seed, int t) {

        // seed to be deterministic (wrt input randomness)
        this.gen = new Random(seed);
        this.wind_direction = wind_direction;
        this.wind_angle = Math.atan2(wind_direction.y, wind_direction.x); 
        this.k = Math.max(3, t);
        this.nextKpoints = new ArrayList<Point>();

        switch(initMode)   {
            case CENTER:
                this.initial = new Point(5.0, 5.0);
                break;

            case WIND:

                Double angle = Point.angleBetweenVectors(new Point(1, 0), wind_direction);
                Double x, y;
                x = -3.0 * Math.cos(angle) + 5.0;
                y = -3.0 * Math.sin(angle) + 5.0;
                this.initial = new Point(x, y);
                
                break;
        }
        
        return initial;
    }

    @Override
    public void init(List<Point> group_locations, List<Point> targets, int id) {
        
        this.id = id;
        this.targets = targets;
        this.numPlayers = group_locations.size();
        this.sHelper = new SailingHelper(targets, id, initial, wind_direction, wind_angle);        
    }

    @Override
    public Point move(List<Point> group_locations, int id, double dt, long time_remaining_ms) {
        // testing timeouts... 
        // try {
        //     TimeUnit.MILLISECONDS.sleep(1);
        // } catch(Exception ex) {
        //     ;
        // }
        // just for first turn

        Point goalLoc = null;
        Point prevTarget = null;
        Point currentLoc = group_locations.get(id);
        ArrayList<Integer> availableTargetIndices = new ArrayList<Integer>();
        

        if(visited_set != null && visited_set.get(id).size() == targets.size()) {
            goalLoc = initial;
        } else {
            

            for(int i = 0; i < targets.size(); i++) {
                if (visited_set == null || !visited_set.get(id).contains(i)) {
                    availableTargetIndices.add(i);
                }
            }

            switch(playMode)   {
                case K_STEPS:

                    if(nextKpoints.size() == 0) {
                        nextKpoints = sHelper.getkOptimalTargets(availableTargetIndices, currentLoc, group_locations, visited_set, this.id, k);
                    }
                    goalLoc = nextKpoints.get(0);
                    break;

                case ONE_STEP:                

                    int indexNextTarget = sHelper.getHeuristicDistance(availableTargetIndices, currentLoc, group_locations, visited_set, this.id).targetIndex;
                    goalLoc = targets.get(indexNextTarget);
                    break;
            }


        }
        double relAngleWouldNeedToGoIn = sHelper.getRelAngle(currentLoc, goalLoc);        
        if (sHelper.shouldCurve(relAngleWouldNeedToGoIn) && !sHelper.closeToBoundary(currentLoc, dt)/*&& isInsideBox(newPointCurvingWouldGive)*/) {
            return sHelper.getBestAbsoluteDirection(currentLoc, goalLoc);
        } else {
            return Point.getDirection(currentLoc, goalLoc);
        }
    }

    /**
    * visited_set.get(i) is a set of targets that the ith player has visited.
    */
    @Override
    public void onMoveFinished(List<Point> group_locations, Map<Integer, Set<Integer>> visited_set) {

        switch(playMode)   {
            case K_STEPS:
                if(this.visited_set != null && this.visited_set.get(id).size() > pointsVisited) {
                    pointsVisited++;                    
                    nextKpoints.remove(0);
                }
                break;
        }

        this.visited_set = visited_set;
    }
}
