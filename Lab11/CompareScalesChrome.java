/* THIS WAS FROM LAB 11*/
//Compare Scales Chromosomes - used to sort an ArrayList of ScalesChrome
public class CompareScalesChrome implements java.util.Comparator<ScalesChrome>
{
	public int compare(ScalesChrome a, ScalesChrome b) // Method called compare which takes in two ScalesChrome objects and method returns an int
	{
		if (a.GetFitness() < b.GetFitness()) return(-1); // if objects a fitness is less than b's fitness we return -1 (which means a's fitness is better)
		if (a.GetFitness() > b.GetFitness()) return(1); // if object a fitness is greater than a's fitness we return 1 (which means b's fitness is better)
		return(0); // if both a and b scaleschrome object fitnesses are the same we return 0.
	}// end of method
}// end of class