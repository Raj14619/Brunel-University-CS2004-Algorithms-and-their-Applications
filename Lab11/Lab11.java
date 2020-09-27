import java.util.ArrayList;

public class Lab11 
{
	public static void main(String args[])
	{
		//Read in the weights
		//Make sure you change the filename as appropriate!
		//This will only work if "c:\temp\1000 Primes.txt" exists!
		ArrayList<Double> w = CS2004.ReadNumberFile("C:\\Users\\Raj\\eclipse-workspace\\Lab11\\src\\1000 Primes");
		//Set the weights
		ScalesChrome.SetWeights(w);
		//Run 10 repeats
		for(int i=0;i<10;++i)
		{
			//Reset the fitness count
			ScalesChrome.ClearFC();
			//The following parameters are not very good!
			//These are the ones you should try and optimise!
			int popsize = 100; // originally 10
			double mrate = 0.001; //originally 1.0
			double crate = 0.5; // originally 0.0
			//You will not need to change the following
			SimpleGeneticAlgorithm ga = new SimpleGeneticAlgorithm(popsize,10,1000,mrate,crate); // 1000 weights
			//Run the GA for 10,000 function calls
			double f = ga.RunSGA(10000,true).GetFitness();// running the genetic algoirthm for 10,000 fitness functions
			System.out.println(f);
		}	
	}
	
	
	public static void ok(double [] weights, int numOfFunctionCalls) {
		
		
		ArrayList<Double> weightsArrayList = new ArrayList<>();
		
		for(int i = 0; i < weights.length; i++) {
			weightsArrayList.add((double) i);
		}
		
		
		ScalesChrome.SetWeights(weightsArrayList);
		
		
			
			ScalesChrome.ClearFC();
			
			int popsize = 100;
			double mrate = 0.001;
			double crate = 0.5;
			
			SimpleGeneticAlgorithm ga2 = new SimpleGeneticAlgorithm(100, 10,1000, mrate, crate);
			
			
			double ok = ga2.RunSGA(numOfFunctionCalls,false).GetFitness();
			
			ScalesChrome ok5 = new ScalesChrome(1000);
			
		
		
		
		
		
		
		
	}
	
	
	
	
	
	
	
}
