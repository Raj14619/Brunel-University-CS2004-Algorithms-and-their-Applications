import java.util.ArrayList;

public class Lab5 {
	
	
	public static void main(String args[])
	{
		
		ArrayList<Data> array = new ArrayList<Data>();

		
		Data x = new Data("Fred",41);
		Data y = new Data("Jo",43);
		Data z = new Data("Zoe",37);
		Data v = new Data("Harry",78);
		
		array.add(x);
		array.add(y);
		array.add(z);
		
		array.add(2, v);
		

		Data.PrintCollection(array);
		
		
		
		
	}

}
