
public class ArrayMaxExercise {
	
	
	public static void main(String[] args) {
		
		double Array [] = {1.9999,3.323,5.34,6.5,9.43};
		
		System.out.println(ArrayMax(Array));
		
	}
	
	public static double ArrayMax (double [] Array){
		
		double CurrentMax = 0;
		
		for(int i = 1; i < Array.length; i++) {
			if(Array[i] > CurrentMax) {
				CurrentMax = Array[i];
			}
			
		}

		return CurrentMax;
		
	}
	

}
