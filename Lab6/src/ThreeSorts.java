//Java Code for Laboratory 6 – Week 6 Sorting Algorithms
//Brunel University, Computer Science Module CS2004, 2017-2018

import java.lang.reflect.Array;
import java.util.*;
public class ThreeSorts 
{
	private static ArrayList<Integer> CopyArray(ArrayList<Integer> a)
	{
		ArrayList<Integer> resa = new ArrayList<Integer>(a.size());
		for(int i=0;i<a.size();++i) resa.add(a.get(i));
		return(resa);
	}
	public static ArrayList<Integer> SortA(ArrayList<Integer> a)// binary sort
	{
		ArrayList<Integer> array = CopyArray(a);
		int n = a.size(),i;
		boolean noswaps = false;
		
		while (noswaps == false)
		{
			noswaps = true;
			for(i=0;i<n-1;++i)
			{
				if (array.get(i) < array.get(i+1))
				{
					Integer temp = array.get(i);
					array.set(i,array.get(i+1));
					array.set(i+1,temp);
					noswaps = false;
				}
			}
		}
		return(array);
	}
	public static ArrayList<Integer> SortB(ArrayList<Integer> a)
	{
		ArrayList<Integer> array = CopyArray(a);
		Integer[] zero = new Integer[a.size()];
		Integer[] one = new Integer[a.size()];
		int i,b;
		Integer x,p;
		//Change from 8 to 32 for whole integers - will run 4 times slower
		for(b=0;b<8;++b)
		{
			int zc = 0;
			int oc = 0;
			for(i=0;i<array.size();++i)
			{
				x = array.get(i);
				p = 1 << b;
				if ((x & p) == 0)
				{
					zero[zc++] = array.get(i);
				}
				else
				{
					one[oc++] = array.get(i);
				}
			}
			for(i=0;i<oc;++i) array.set(i,one[i]);
			for(i=0;i<zc;++i) array.set(i+oc,zero[i]);
		}
		return(array);
	}
	
	
	public static ArrayList<Integer> SortC(ArrayList<Integer> a)
	{
		ArrayList<Integer> array = CopyArray(a);
		SortC(array,0,array.size()-1);
		return(array);
	}
	public static void SortC(ArrayList<Integer> array,int first,int last)
	{
		if (first < last)
		{
		   int pivot = PivotList(array,first,last);
		   SortC(array,first,pivot-1);
		   SortC(array,pivot+1,last);
		}	
	}
	
	
	private static void Swap(ArrayList<Integer> array,int a,int b)
	{
		Integer temp = array.get(a);
		array.set(a,array.get(b));
		array.set(b,temp);
	}
	private static int PivotList(ArrayList<Integer> array,int first,int last)
	{
		Integer PivotValue = array.get(first);
		int PivotPoint = first;
		for(int index=first+1;index<=last;++index)
		{
			if (array.get(index) > PivotValue)
			{
				PivotPoint = PivotPoint+1;
				Swap(array,PivotPoint,index);
			}
		}
		Swap(array,first,PivotPoint);
		return(PivotPoint);
	}

	
	public static ArrayList<Integer> RandomArray(int n){
		
		
	ArrayList<Integer> Raj = new ArrayList<Integer>();
	
	
	
	
	for(int i = 0; i < n; i++) {
	
		//System.out.println(i);
		
		
	
		Random rand = new Random();

	//creates random integer from 0 to 255
		
	rand.setSeed(System.currentTimeMillis());
	Integer r = Math.abs(rand.nextInt() % 256);
	Raj.add(r);
	
	}
		
		
		return Raj;
	}
	
	public static void ShowArray ( ArrayList<Integer> Array){
		
		for(int i = 0; i < Array.size(); i++) {
			System.out.println(Array.get(i));
		}
		
	}
	
	public static void main(String[] args) {
		
		
		//Get 1000 random integers to Okay Arraylist
		ArrayList<Integer> Okay = RandomArray(1000);
		
		
		
		//ShowArray(Okay);
		
		
		//ArrayList<Integer>Sort1ArrayOutput = SortA(Okay);//SortA meethod sorts from highest to lowest
		//ShowArray(Sort1ArrayOutput);
		
		//ArrayList<Integer>Sort2ArrayOutput = SortB(Okay);//SortB sorts from highest to lowest
		//ShowArray(Sort2ArrayOutput);
		
		//ArrayList<Integer>Sort3ArrayOutput = SortC(Okay);
		//ShowArray(Sort3ArrayOutput);
		
	}
	
	
}