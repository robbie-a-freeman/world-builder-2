/******************************************************************************
  *  Name:         Robert Freeman
  *  Project:      Terraformer
  *
  *  Description:  Creates a 2D projection of a 3D planet using the Princeton
  *                StdDraw library. The planet is supposed to be earth-like, in
  *                that the planet has the same interstellar properties as the
  *                Earth and is the same size/density, but the planet's tectonic
  *                activity is simulated and assembled differently. The idea is
  *                that the landmasses/ocean bodies will be completely different,
  *                resulting in varying weather patterns, terrain/biome
  *                formations, and large-scale human behaviors. Eventually,
  *                civilizations will be inserted into the program to coordinate
  *                with Planet.java. The planet is represented as a perfect
  *                sphere with cells representing plots of land. For ease of
  *                some calculations, the cells are projected onto the sphere
  *                but are actually located on a cube with the side length of 2r,
  *                where r is the radius of the sphere.
  *******************************************************************************/

import edu.princeton.cs.algs4.StdDraw;
import edu.princeton.cs.algs4.Stack;
import java.awt.Color;

public class Planet {
    
    // an Iterable that contains all of the gradients affecting the terrain of
    // the Planet
    private Stack<LandGradient> landGradients;
    // an Iterable that contains all of the cells representing the terrain of
    // the Planet
    private Cell[] cells;
    // number of faces in cube
    private final int SIDES_OF_CUBE;
    // length of side of a cell
    private final double CELL_SIDE_LENGTH;
    // radius of the sphere representing the Planet, not yet implemented
    // private final double PLANET_RADIUS;
    
    // generate an empty planet divided into cells. has a given number of
    // subdivisions. subdivisions must be a multiple of a perfect square and 6
    // (6, 24, 54, ...)
    public Planet(int subdivisions) {
        
        SIDES_OF_CUBE = 6;
        landGradients = new Stack<LandGradient>();
        final int cellRowCol = (int) Math.sqrt((double) subdivisions / (double) SIDES_OF_CUBE);
        CELL_SIDE_LENGTH = 2.0 / cellRowCol;
        
        // create cells, put them in the array
        cells = new Cell[subdivisions];
        int numberOfCells = 0;
        // top
        for (int i = 0; i < cellRowCol; i++) {
            for (int j = 0; j < cellRowCol; j++) {
                double[] sphCoords = toSpherical(2.0 * (double) i / (double) cellRowCol - 1.0,
                                                 -2.0 * (double) j / (double) cellRowCol + 1.0,
                                                 1.0);
                cells[numberOfCells] = new Cell(sphCoords[0], sphCoords[1], sphCoords[2], numberOfCells);
                numberOfCells++;
            }
        }
        
        // front
        for (int i = 0; i < cellRowCol; i++) {
            for (int j = 0; j < cellRowCol; j++) {
                double[] sphCoords = toSpherical(1.0,
                                                 2.0 * (double) i / (double) cellRowCol - 1.0,
                                                 -2.0 * (double) j / (double) cellRowCol + 1.0);
                
                cells[numberOfCells] = new Cell(sphCoords[0], sphCoords[1], sphCoords[2], numberOfCells);
                numberOfCells++;
            }
        }
        
        // right
        for (int i = 0; i < cellRowCol; i++) {
            for (int j = 0; j < cellRowCol; j++) {
                double[] sphCoords = toSpherical(2.0 * (double) i / (double) cellRowCol - 1.0,
                                                 1.0,
                                                 -1.0 * (double) j / (double) cellRowCol + 1.0);
                cells[numberOfCells] = new Cell(sphCoords[0], sphCoords[1], sphCoords[2], numberOfCells);
                numberOfCells++;
            }
        }
        
        // back
        for (int i = 0; i < cellRowCol; i++) {
            for (int j = 0; j < cellRowCol; j++) {
                double[] sphCoords = toSpherical(-1.0,
                                                 2.0 * (double) i / (double) cellRowCol - 1.0,
                                                 -2.0 * (double) j / (double) cellRowCol + 1.0);
                cells[numberOfCells] = new Cell(sphCoords[0], sphCoords[1], sphCoords[2], numberOfCells);
                numberOfCells++;
            }
        }
        
        // left
        for (int i = 0; i < cellRowCol; i++) {
            for (int j = 0; j < cellRowCol; j++) {
                double[] sphCoords = toSpherical(2.0 * (double) i / (double) cellRowCol - 1.0,
                                                 -1.0,
                                                 -2.0 * (double) j / (double) cellRowCol + 1.0);
                cells[numberOfCells] = new Cell(sphCoords[0], sphCoords[1], sphCoords[2], numberOfCells);
                numberOfCells++;
            }
        }
        
        // bottom
        for (int i = 0; i < cellRowCol; i++) {
            for (int j = 0; j < cellRowCol; j++) {
                double[] sphCoords = toSpherical(-2.0 * (double) i / (double) cellRowCol + 1.0,
                                                 2.0 * (double) j / (double) cellRowCol - 1.0,
                                                 -1.0);
                cells[numberOfCells] = new Cell(sphCoords[0], sphCoords[1], sphCoords[2], numberOfCells);
                numberOfCells++;
            }
        }
        
    }
    
    
// convert from cartesian to spherical coordinates assuming r = 1
    private static double[] toSpherical(double x, double y, double z) {
        
        double theta, phi;
        
        // theta = tan^-1(y / x)
        theta = Math.atan(y / x);
        if (x < 0.) {
            theta = Math.PI + theta;
        }
        if (Double.compare(x, 0.) == 0.) {
            if (y > 0.) {
                theta = Math.PI / 2.0;
            } else if (y < 0.) {
                theta = -Math.PI / 2.0;
            } else {
                theta = 0.;
            }
        }
        
        // phi = tan^-1(sqrt(x^2 + y^2) / z)
        phi = Math.atan(Math.sqrt(x * x + y * y) / z);
        if (z < 0.) {
            phi = Math.PI + phi;
        } else if (Double.compare(z, 0.) == 0) {
            phi = Math.PI / 2.0;
        }
        
        // if theta or phi are suuuuper small, just make them 0
        if (Math.abs(theta) < Math.pow(10.0, -15.0))
            theta = 0.;
        if (Math.abs(phi) < Math.pow(10.0, -15.0))
            phi = 0.;
        
        // find the radius, return the spherical coordinates as an array
        double[] sphCoords = new double[3];
        sphCoords[0] = 1.;
        sphCoords[1] = theta;
        sphCoords[2] = phi;
        
        return sphCoords;
        
    }
    
// convert from spherical to cartesian coordinates assuming abs(r) = 1
    private static double[] toCartesian(double r, double theta, double phi) {
        
        // find the line that runs through the necessary point
        double[] cartesianCoords = new double[3];
        // x = r * sin(phi) * cos(theta)
        cartesianCoords[0] = r * Math.sin(phi) * Math.cos(theta);
        // y = r * sin(phi) * sin(theta)
        cartesianCoords[1] = r * Math.sin(phi) * Math.sin(theta);
        // z = r * cos(phi)
        cartesianCoords[2] = r * Math.cos(phi);
        
        // scale the line so that it reaches the point on the cube
        double maxCoord = Math.max(Math.max(Math.abs(cartesianCoords[0]),
                                            Math.abs(cartesianCoords[1])),
                                   Math.abs(cartesianCoords[2]));
        
        // return the correct coordinates
        if (Math.abs(cartesianCoords[0]) < Math.pow(10.0, -15.0))
            cartesianCoords[0] = 0.;
        cartesianCoords[0] /= maxCoord;
        if (Math.abs(cartesianCoords[1]) < Math.pow(10.0, -15.0))
            cartesianCoords[1] = 0.;
        cartesianCoords[1] /= maxCoord;
        if (Math.abs(cartesianCoords[2]) < Math.pow(10.0, -15.0))
            cartesianCoords[2] = 0.;
        cartesianCoords[2] /= maxCoord;
        
        return cartesianCoords;
        
    }
    
// adds a LandGradient object to the list of gradient objects. changes the
// elevation of the Planet surface according to the gradient formula.
    public void sculpt(double third, double second, double first,
                       double constant, char var) {
        
        landGradients.push(new LandGradient(third, second, first, constant,
                                            landGradients.size(), var));
        
    }
    
// draw each cell from the top left corner with StdDraw
    private void drawCell(double x, double y, double CELL_SIDE_LENGTH, Color cellColor) {
        
        StdDraw.setPenColor(Color.BLACK);
        StdDraw.line(x, y, x + CELL_SIDE_LENGTH, y);
        StdDraw.line(x, y, x, y - CELL_SIDE_LENGTH);
        StdDraw.line(x + CELL_SIDE_LENGTH, y,
                     x + CELL_SIDE_LENGTH, y - CELL_SIDE_LENGTH);
        StdDraw.line(x, y - CELL_SIDE_LENGTH,
                     x + CELL_SIDE_LENGTH, y - CELL_SIDE_LENGTH);
        //StdDraw.setPenColor(cellColor);
        //for (int 
        
    }
    
    // test the given point and print the output to Sysout
    private void testPointConversion(double x, double y, double z) {
        
        double[] testSph = toSpherical(x, y, z);
        System.out.printf("\n(%5.2f, %5.2f, %5.2f) to sph - radius: %5.2f theta: %5.2f phi: %5.2f\n",
                          x, y, z, testSph[0], testSph[1], testSph[2]);
        double[] testCart = toCartesian(testSph[0], testSph[1], testSph[2]);
        System.out.printf("(%5.2f, %5.2f, %5.2f) to cart - x: %5.2f y: %5.2f z: %5.2f\n",
                           testSph[0], testSph[1], testSph[2], testCart[0], testCart[1], testCart[2]);
        
    }
    
// unit testing
    public static void main(String[] args) {
        
        Planet p = new Planet(Integer.parseInt(args[0]));
        // p.sculpt(1.0, 3.0, 4.0, -2.0, 'x');
        
        System.out.println("===================POINT TESTING===================");
        p.testPointConversion(1.0, 0., 0.);
        p.testPointConversion(0., 1.0, 0.);
        p.testPointConversion(-0.5, 0.5, 0.);
        p.testPointConversion(-1.0, 1.0, 0.);
        p.testPointConversion(-1.0, 0., 0.);
        p.testPointConversion(0., 0., 1.0);
        p.testPointConversion(0., 0., -1.0);
        p.testPointConversion(1.0, -1.0, 1.0);
        System.out.println();
        
        
        // update the planet cell-by-cell once. assumes elevation starts from 0.
        for (Cell c : p.cells) {
            double elevation = c.getElevation();
            double coords[] = toCartesian(c.getLocation()[0], c.getLocation()[1], c.getLocation()[2]);
            for (LandGradient lg : p.landGradients) {
                c.setElevation(elevation + lg.calculate(coords[0], coords[1], coords[2]));
            }
        }
        
        // create the canvas/base picture
        StdDraw.setCanvasSize(1024, 1024);
        StdDraw.setXscale(-2.0, 2.0);
        StdDraw.setYscale(-2.0, 2.0);
        StdDraw.setPenColor(Color.BLUE);
        StdDraw.filledCircle(0, 0, 1.0);
        
        // draw the cells
        System.out.println("===================CELL DRAWING===================");
        for (Cell c : p.cells) {
            double[] sphCoords = c.getLocation();
            double[] cartCoords = toCartesian(sphCoords[0], sphCoords[1], sphCoords[2]);
            if (cartCoords[0] > 0.999 && cartCoords[0] < 1.001) {
                p.drawCell(cartCoords[1], cartCoords[2], p.CELL_SIDE_LENGTH, c.getColor());
                System.out.printf("Cell number %d -   x: %5.4f   y: %5.4f   z: %5.4f   elevation: %5.4f\n",
                                  c.getID(), cartCoords[0], cartCoords[1], cartCoords[2], c.getElevation());
            }
        }
        System.out.println();
        
        System.out.println("===================ALL CELLS===================");
        for (Cell c : p.cells) {
            double[] sphCoords = c.getLocation();
            double[] cartCoords = toCartesian(sphCoords[0], sphCoords[1], sphCoords[2]);
            System.out.printf("\nCell number %d -   x: %5.4f   y: %5.4f   z: %5.4f   elevation: %5.4f\n",
                              c.getID(), cartCoords[0], cartCoords[1], cartCoords[2], c.getElevation());
        }
        System.out.println();
        
        StdDraw.show();
        
    }
    
}