/* THIS IS FROM LAB 11*/

import java.util.ArrayList;


public class ScalesChrome 
{
	//The number of times the fitness function is evaluated
	private static int fc = 0; /* This relates to the ScalesChrome object*/
	//This is used to store the fitness in
	private double fitness = -1; /* This relates to the ScalesChrome object*/
	//The representation - a binary array
	private ArrayList<Integer> rep = new ArrayList<Integer>(); /* This relates to the ScalesChrome object*/
	//The set of weights we are working with
	private static ArrayList<Double> weights = new ArrayList<Double>();
	//Reset the fitness count back to zero
	public static void ClearFC()/*Start of ClearFC method*/
	{
		fc = 0; /*fc is equal to 0*/
	} /*End of ClearFC method*/
	//Get the number of times the fitness function has been called
	public static int GetFC()
	{
		return(fc);
	}
	//Set the current set of weights
	public static void SetWeights(ArrayList<Double> w)
	{
		weights = w;
	}
	//Create a Chromosome based on string 's'
	public ScalesChrome(String s) /*This is a constructor for the ScalesChrome object, we are passing in String s which will be used to determine the object*/
	{
		boolean ok = true; /* set boolean ok as true*/
		int n = s.length(); /* n is equal to the length of String s */
		for(int i=0;i<n;++i) /* go through a for loop n amount of times (which is the length of string s) */
		{
			char si = s.charAt(i); /* char si is equal to index i of string s  */
			if (si != '0' && si != '1') /* if char si is not equal to '0' and not equal to '1' (keep in mind we are comparing chars not ints) */
			{
				ok = false; /* set boolean ok to false*/
			}/*end if */
			else /* if char si is equal to '0' or '1' (keep in mind we are comparing chars)*/
			{
				rep.add((si == '1')?1:0); /*rep is an arraylist of integers in which */
			}
		}/* end of for loop*/
		if (!ok) /* if ok is false */
		{
			RandomBinary(n); /* call random binary string and pass in n which is the length of String s*/
		}
	}/* end of constructor method*/
	//Create a Chromosome based on array 'r'
	public ScalesChrome(ArrayList<Integer> r) /* This is also a constructor for the ScalesChrome object, in which we are passing an arraylist of integers */
	{
		if (r.size() > weights.size()) /* if size of the arraylist r (which is being passed into the constructor) is greater than the size of the weights arraylist*/
		{
			System.out.println("+++Array size miss match within ScalesChrome"); /*Output error message*/
			System.exit(0); /* Close the program */
		} /* end of if statement */
		for(int i=0;i<r.size();++i) /* go through all the elements of the r arraylist */
		{
			if (r.get(i) < 2) rep.add(r.get(i));/* if r arraylist at index i is less than 2 add from the arraylist r at index i to the rep arraylist  */
		}/* end of for loop */
	}/*end of constructor */
	
	//Create a random binary Chromosome of length 'n' genes/bits
	public ScalesChrome(int n) /* This is also a contrusctor for the ScalesChrome object, we pass in an int into this constructor */
	{
		RandomBinary(n); /* call Randombinary method and pass n */
	}
	//Return the representation (as a pointer) for the Chromosome
	public ArrayList<Integer> GetRep() /* This method used a Chromosome object to be called */
	{
		return(rep); /*The method returns the rep arraylist relating to the object*/
	}/* end of method */
	//Copy and return the representation
	public ArrayList<Integer> CopyRep() /*  This method uses a Chromosome object to be called (I think) */
	{
		ArrayList<Integer> cr = new ArrayList<Integer>(); /* declaring arraylist of integers called cr */
		for(int i=0;i<rep.size();++i) cr.add(rep.get(i)); /* going through all the elements in the rep arraylist and adding the elements to the arraylist called cr */
		return(cr); /* return the cr arraylist*/
	}
	//Create a random binary Chromosome of length 'n' genes/bits
	private void RandomBinary(int n) /* this method takes in an integer */
	{
		rep.clear(); /* remove all the elements in the rep arraylist */
		for(int i=0;i<n;++i) /*for loop to go through */
		{
			Integer y = CS2004.UI(0,1); /* Integer y is equal to the input of CS2004.UI(0,1) (keep in mind CS2004.UI returns a random number of 0's and 1's in any order) (also keep in mind the amount of numbers we will randomly generate depends on the n variable from the for loop )*/
			rep.add(y);/* add Integer y to rep*/
		}/*end of for loop*/
	}/*end of RandomBinary method*/
	//Compute the Scales fitness function
	public double GetFitness() /* This method requires an object to be called, and the method returns a double*/
	{ 
		//We only need to evaluate it once
		if (fitness != -1) return(fitness); /* if fitness is not equal to -1 then return fitness to the value */
		if (weights.size() == 0) /* if weights.size is equal to 0*/
		{
			System.out.println("+++Weights not set within ScalesChrome"); /* output error message*/
			System.exit(1); /* exit the program */
		}/* end of if statement*/
		if (rep.size() > weights.size()) return(-1); /* if rep size is greater than weights size return -1 to the method*/
		double lhs = 0.0,rhs = 0.0; /* lhs is equal to 0.0 and rhs is equal to 0.0 */
		int n = rep.size(); /*n is equal to the size of the rep arraylist*/
		for(int i=0;i<n;++i) /* for loop which happens n amount of times*/
		{
			Integer si = rep.get(i); /* declare Integer si which is equal to index i of arraylist rep */
			double wi = weights.get(i); /* declare double wi which is equal to weights arraylist at index i */
			if (si == 0) lhs += wi; /* if si is equal to 0 then lhs adds the previous value of lhs and the value of wi */
			if (si == 1) rhs += wi; /* if si is equal to 1 then rhs adds the previous value of rhs and the value of wi */
		}/* end of for loop*/
		fitness = Math.abs(lhs-rhs); /* fitness is equal to lhs - rhs*/
		++fc; /*adding one to fc */
		return(fitness); /*return fitness to the method*/
	}/*End of GetFitness method*/
	//Display the fitness and representation
	public void print()
	{
		System.out.print(GetFitness());
		System.out.print(" ");
		System.out.print(rep);
	}
	//Display the fitness and representation followed by a new line
	public void println()
	{
		print();
		System.out.println();
	}
	
	
	
	
	
}