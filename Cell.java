/******************************************************************************
  *  Name:         Robert Freeman
  *  Project:      Terraformer
  *
  *  Description:  Datatype that holds the location and stats of an arbitrarily
  *                small portion of land in an instance of Planet.java. Location,
  *                ID, and base color are immutable, but everything else is
  *                mutable.
  *******************************************************************************/

import java.awt.Color;

public class Cell {
    
    // height of the land in the cell
    private double elevation;
    // distance from center. only 1 or -1
    private final double radius;
    // angle from the x-axis of the cell. in pi-radians
    private final double theta;
    // angle from the z-axis of the cell. in pi-radians
    private final double phi;
    // temporary color of the cell
    private final Color color;
    // ID number of the cell
    private final long ID;
    
    // constructs a cell at sea level
    public Cell(double radius, double theta, double phi, long ID) {
        
        elevation = 0.0;
        this.radius = radius;
        this.theta = theta;
        this.phi = phi;
        this.ID = ID;
        color = new Color((int) (255 * Math.random()), (int) (255 * Math.random()),
                          (int) (255 * Math.random()));
        
    }
    
    // gets the base color of the tile
    public Color getColor() {
        
        return color;
        
    }
    
    // gets the base color of the tile
    public long getID() {
        
        return ID;
        
    }
    
    // set the elevation of the cell to something else
    public void setElevation(double newElevation) {
        
        elevation = newElevation;
        
    }
    
    // get the elevation of the cell
    public double getElevation() {
        
        return elevation;
        
    }
    
    // get the top left corner of the cell and return as a double array
    // (rad, theta, phi)
    public double[] getLocation() {
        
        double[] sphCoords = new double[3];
        sphCoords[0] = radius;
        sphCoords[1] = theta;
        sphCoords[2] = phi;
        return sphCoords;
        
    }
    
    // unit testing
    public static void main(String[] args) {
        
        
        
    }
    
}