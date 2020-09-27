
public class Numeric
{  /**
      Tests whether two floating-point numbers are
      equal, except for a roundoff error
      x a floating-point number
      y a floating-point number
      returns true if x and y are approximately equal
   */
	
	// This is done by Raj starting
	public static void main(String[] args) {
		
		System.out.println(approxEqual(10.5555555555, 10.5555555555));
		System.out.println(approxEqual(10.15, 10.55));
		System.out.println(approxEqual(13, 12));
	}
	// This is done by Raj ending
	
	
   public static boolean approxEqual(double x, double y)
   {  
      final double EPSILON = 1E-10;
      if (Math.abs(x-y)<EPSILON)
      {
    	  return(true);
      }
      return(false);
   }
   // more numeric methods can be added here
}
