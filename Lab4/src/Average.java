
public class Average {

	public static void main(String[] args) {
		
		double [] x = {1,11,5,7,10}; // declaring a DOUBLE ARRAY called x and inputting numbers
		
		double [] A = PrefixAverages1(x); // declaring an ARRAY called A which is equal to calling the PREFIXAVERAGES1 METHOD with the input of x array
		
		for(int i = 0; i < A.length; i++) { 
			System.out.println(A[i]); //printing out all the elements in Array A 
		}
		
		
	}
	
	
	
	public static double [] PrefixAverages1 (double x []) {
		
		double A [] = new double [x.length];// declaring a double array called "A" which has a size (number of elements) equal to those of array x which is being passed into the method
		
		for(int i = 0; i < A.length; i++) { // going through all elements in array A and declaring s equal to the first element in x (technically it would have been fine 
			double s = x[0]; // s is equal to the first element in x
		
		for(int j = 1; j < A.length; j++) { 
			
			if(j <=  i) { // if J is less than or equal to I
				 s = s + x[j]; // s is equal to S added with element j in array x
			}
		}
		
		A[i] = s/(i+1); // Array A in element i is equal to s divided by i added with 1
		}
		return A;// returning A
	}
	
		
	public static double [] PrefixAverages2 (double x [] ) {
		
		double[] A = new double [x.length];
		double s = 0;
	
		for(int i = 0; i < x.length; i++) {
			s = s + x[i];
			A[i] = s/(i+1);
		}
		
		
		return A;
	}
	

}
	
	

