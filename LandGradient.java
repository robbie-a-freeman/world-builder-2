/******************************************************************************
  *  Name:         Robert Freeman
  *  Project:      Terraformer
  *
  *  Description:  Immutable data type that constitutes a function that can be
  *                constant, linear, quadratic, or cubic. The function represents
  *                the change in elevation at a given point in a given Planet.
  *                The Planet is described in the class Planet.java. LandGradients
  *                can be dependent on x, y, or z. LandGradients can be merged
  *                together.
  *******************************************************************************/

public class LandGradient {
    
    // large ID number of land gradient, used for finding this object
    private final long ID;   
    // constant of function
    private final double constant;
    // linear coefficient of function
    private final double firstDegreeX;
    // quadratic coefficient of function
    private final double secondDegreeX;
    // cubic coefficient of function
    private final double thirdDegreeX;
    // linear coefficient of function
    private final double firstDegreeY;
    // quadratic coefficient of function
    private final double secondDegreeY;
    // cubic coefficient of function
    private final double thirdDegreeY;
    // linear coefficient of function
    private final double firstDegreeZ;
    // quadratic coefficient of function
    private final double secondDegreeZ;
    // cubic coefficient of function
    private final double thirdDegreeZ;
    
    // generate a land gradient with the given coefficients
    // only dependent on one variable var
    public LandGradient(double third, double second, double first,
                        double constant, long ID, char var) {
        
        this.constant = constant;
        this.ID = ID;
        
        switch(var) {
            case 'x' :
                this.firstDegreeX = first;
                this.secondDegreeX = second;
                this.thirdDegreeX = third;
                this.firstDegreeY = 0.;
                this.secondDegreeY = 0.;
                this.thirdDegreeY = 0.;
                this.firstDegreeZ = 0.;
                this.secondDegreeZ = 0.;
                this.thirdDegreeZ = 0.;
                break;
            case 'y' :
                this.firstDegreeX = 0.;
                this.secondDegreeX = 0.;
                this.thirdDegreeX = 0.;
                this.firstDegreeY = first;
                this.secondDegreeY = second;
                this.thirdDegreeY = third;
                this.firstDegreeZ = 0.;
                this.secondDegreeZ = 0.;
                this.thirdDegreeZ = 0.;
                break;
            case 'z' :
                this.firstDegreeX = 0.;
                this.secondDegreeX = 0.;
                this.thirdDegreeX = 0.;
                this.firstDegreeY = 0.;
                this.secondDegreeY = 0.;
                this.thirdDegreeY = 0.;
                this.firstDegreeZ = first;
                this.secondDegreeZ = second;
                this.thirdDegreeZ = third;
                break;
            default :
                throw new IllegalArgumentException();
        }
        
    }
    
    // string representation of a given LandGradient
    public String toString(LandGradient lg) {
        
        return "ID: " + lg.ID;
        
    }
    
    // calculates the scalar quantity at the given 3D point (x, y, z)
    public double calculate(double x, double y, double z) {
        
        double result = constant;
        result += thirdDegreeX * x * x * x + secondDegreeX * x * x + firstDegreeX * x;
        result += thirdDegreeY * y * y * y + secondDegreeY * y * y + firstDegreeY * y;
        result += thirdDegreeZ * z * z * z + secondDegreeZ * z * z + firstDegreeZ * z;
        
        return result;
        
    }
    
    // unit testing
    public static void main(String[] args) {
        
        LandGradient test = new LandGradient(1.0, 3.0, 4.0, -2.0, 1, 'x');
        System.out.println(test.toString());
        System.out.println(test.calculate(3.0, 3.0, 3.0));
        
    }
    
}