import java.util.ArrayList;



public class Lab8 {
	
	public static void main(String args[])
	{
		for(int i=0;i<10;++i)
		{
			int x = CS2004.UI(-1, 1);
			//System.out.println(x);
		}
		
		
		ScalesSolution s = new ScalesSolution("1010x");
		//s.println();
		
		ArrayList<Double> test1 = new ArrayList<Double>();
		
		test1.add(1.00);
		test1.add(2.00);
		test1.add(3.00);
		test1.add(4.00);
		test1.add(10.00);
		
		
		//System.out.println(s.ScalesFitness(test1));
		
		
		
		//Exercise 8.5
		
		ArrayList<Double> test2 = CS2004.ReadNumberFile("C:\\Users\\Raj\\eclipse-workspace\\Lab8\\src\\1000 Primes");
		
		ArrayList<Double> test2TrimmedPrimeNumbers = TrimPrimeNumbers(test2, 8);
		
		
		
		for (int i = 0; i < test2TrimmedPrimeNumbers.size(); i++) {
			System.out.println(test2TrimmedPrimeNumbers.get(i));
		}
		
		ScalesSolution s2 = new ScalesSolution("1010101x");
		s2.println();
		
		System.out.println(s2.ScalesFitness(test2TrimmedPrimeNumbers));
		
	}
	
	//Method made by Raj
	public static ArrayList<Double> TrimPrimeNumbers(ArrayList<Double> arraylist,int numberOfPrimesRequired ){
		
		ArrayList<Double> trimmedArrayList = new ArrayList<Double>(numberOfPrimesRequired);
		
		for(int i = 0; i < numberOfPrimesRequired; i++) {

			trimmedArrayList.add(arraylist.get(i));
		}
		
		
		
		
		return trimmedArrayList;
		
	}

	
	

}
