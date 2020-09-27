import java.util.Random;

public class Lab7 {

	public static void main(String[] args) {
	
		double g [][] = {{0,1,2},{1,0,3},{2,3,0}};
		
		double mst [][] = MST.PrimsMST(g);
				
				
		double g2[][] = { {0,1,2,3,0},{1,0,6,0,5},{2,6,0,4,1},{3,0,4,0,2},{0,5,1,2,0}};
		
		double mst2 [][] = MST.PrimsMST(g2);
		
		//printArray(g2, mst2);

		//printArray(g, mst);
		
		double [][] g3 = RandomArray(100);
		
		double [][] mst3;
		
		
		long startTime = System.nanoTime();
		mst3 = MST.PrimsMST(g3);
		long endTime = System.nanoTime();
		long duration = (endTime - startTime);
		System.out.println("it took" + duration);
		
		//printArray(g3, mst3);
		
		
		
	}
	
	public static void printArray(double [][] G, double [][] MST) {
		
		System.out.println("Start of G");
	
	for(int i = 0; i < G.length; i++ ) {
		for(int j = 0; j < G[0].length; j++) {
			
			System.out.println(G[i][j]);
			
		}
		
		
	}
	
	System.out.println("Start of MST");
		
		for(int i = 0; i < MST.length; i++ ) {
			for(int j = 0; j < MST[0].length; j++) {
				
				System.out.println("I is:" + i);
				System.out.println(MST[i][j]);
				
			}
			
			
		}
	}
		
		public static double [][] RandomArray(int n){
			
			double Array[][] = new double[n][n];
			
			
			Random r = new Random();
			
			for(int i = 0; i < Array.length; i++) {
				for(int j = 0; j < Array[0].length; j++) {
					
					Array[i][j] = r.nextInt(101);
				
				}
				
			}
			
			return Array;
			
		}
	

	
		
		
	

	}

