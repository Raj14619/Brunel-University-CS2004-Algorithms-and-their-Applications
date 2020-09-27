import java.io.BufferedReader;
import java.io.FileReader;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

public class TSP 
{
 
 protected static double rmhcFitness;
 protected static double rrhcFitness;
 protected static double shcFitness;
 protected static double saFitness;
 protected static double rmhcMSTEfficiency;
 protected static double rrhcMSTEfficiency;
 protected static double shcMSTEfficiency;
 protected static double saMSTEfficiency;
 protected static double rmhcEfficiency;
 protected static double rrhcEfficiency;
 protected static double shcEfficiency;
 protected static double saEfficiency;
 
 
 /*----------------------------------------------------------------------------------------------------------MAIN-----------------------------------------------------------------------------------------------------------------------------------*/
 public static void main(String [] args) {
  
  Double [][] Matrix = new Double [3][2];
  
  Matrix [0][0] = 5.5;
  Matrix[0][1] = 7.5;
  Matrix[1][0] = 9.5;
  Matrix[1][1] = 10.5;
  Matrix[2][0] = 15.4;
  
  Integer iterations = 1000000000;
  
  //System.out.println(codeRunner(Matrix, iterations));
  
  
  double [][] MatrixPassedIn = ReadArrayFile("H:\\eclipseWorkspace\\Lab15Coderunner\\src\\TSP_48.txt", " ");
  Integer Numofiterations = 1000000;
  
  
  codeRunner(MatrixPassedIn,Numofiterations);
  
 }
 
 /*----------------------------------------------------------------------------------------------------END--OF--MAIN-----------------------------------------------------------------------------------------------------------------------------------*/
 
 
 
 
 
 /*------------------------------------------------=START OF CODERUNNER-----------------------------------------------------*/
 public static ArrayList<Byte> codeRunner(double [][] Matrix, Integer NumberOfIterations) {
  
  //System.out.println("code");
  /*START OF TRY AND CATCH*/
  if(Matrix == null) {
   return null;
  }
  
  
  if(NumberOfIterations == null) {
   return null;
  }
  
  if(NumberOfIterations < 1) {
   return null;
  }
  //Probably wont need these two below if statements
  if(Matrix.length == 0) {
   return null;
  }
  
  if(Matrix[0].length == 0) {
   return null;
  }
  //Probably wont need these two above if statements
  
  System.out.println(Matrix.length);
  System.out.println(Matrix[0].length);
  if(Matrix.length != Matrix[0].length) {
   System.out.println("ok");
   return null;
  }
  
  
  for(int i = 0; i < Matrix.length; i++) {
   for(int j = 0; j < Matrix[i].length; j++) {
    //if(Matrix[i][j] == null) {
    // return null;
    //}
     
   }
   
  }
  /*END OF TRY AND CATCH*/
  
  
  /*START OF CONVERTING DATA TYPES*/
  
  int NewNumberOfIterations = ConvertIterations(NumberOfIterations);
  double [][] NewMatrix = ConvertMatrix(Matrix);
  
  /*END OF CONVERTING DATA TYPES*/
  
  ArrayList<Integer> Tour = new ArrayList<>();
  Tour = PopulateCities(Matrix.length);
  
  
  Collections.shuffle(Tour);
  
  
  //START OF PRESET VALUES - NO TOUCHING
  
  
  double stochasticTemperature = FitnessFunction(Tour.size(),Tour,NewMatrix) / 1000;
  double annealingTemperature = FitnessFunction(Tour.size(), Tour, NewMatrix) / 250;
  double coolingRate = CalculateCoolingRate(stochasticTemperature, NewNumberOfIterations);
  
  
  // END OF PRESET VALUES
  
  
  
  // START OF CALLING ALGORITHMS
  ArrayList<Integer> RRHCTour = RRHC(Tour,NewNumberOfIterations,NewMatrix,false);
  ArrayList<Integer> SATour = SA(Tour,annealingTemperature,NewNumberOfIterations,NewMatrix,coolingRate,false);
  
  ArrayList<Integer>RMHCTour = RMHC(Tour, NewNumberOfIterations, NewMatrix, false);
  ArrayList<Integer>SHCTour = SHC(Tour,NewNumberOfIterations,NewMatrix,stochasticTemperature,false);
  
  
  System.out.println("This is the RMHC fitness  " + rmhcFitness);
  System.out.println("This is the SHC fitness  " + shcFitness);
  System.out.println("This is the RRHC fitness  " + rrhcFitness);
  System.out.println("This is the SA fitness  " + saFitness);
  
  
  // END OF CALLING ALGORITHMS
  
  
  double fitnessOfRMHC = FitnessFunction(RMHCTour.size(), RMHCTour, NewMatrix);
  double fitnessOfSHC = FitnessFunction(SHCTour.size(), SHCTour, NewMatrix);
  double fitnessOfRRHC =FitnessFunction(RRHCTour.size(), RRHCTour, NewMatrix);
  double fitnessOfSATour =FitnessFunction(SATour.size(), SATour, NewMatrix);
  CalculateEfficiencyOfMST(fitnessOfRMHC,NewMatrix);
  
  double [] ArrayOfFitnesses = new double[4];
  ArrayOfFitnesses[0] = fitnessOfRMHC;
  ArrayOfFitnesses[1] = fitnessOfSHC;
  ArrayOfFitnesses[2] = fitnessOfRRHC;
  ArrayOfFitnesses[3] = fitnessOfSATour;
  
  double bestfitness = Double.MAX_VALUE;
  int bestIndex = Integer.MAX_VALUE;
  
  for(int i = 0; i < ArrayOfFitnesses.length; i++) {
   if(ArrayOfFitnesses[i] < bestfitness) {
    bestfitness = ArrayOfFitnesses[i];
    bestIndex = i;
   }
  }
  
  
  
  ArrayList<Integer> BestAlgorithmSolution = new ArrayList<Integer>();
  
  if(bestIndex == 0) {
   BestAlgorithmSolution = (ArrayList<Integer>) RMHCTour.clone();
  }
  if(bestIndex == 1) {
   BestAlgorithmSolution = (ArrayList<Integer>) SHCTour.clone();
  }
  if(bestIndex == 2) {
   BestAlgorithmSolution = (ArrayList<Integer>) RRHCTour.clone();
  }
  if(bestIndex == 3) {
   BestAlgorithmSolution = (ArrayList<Integer>) SATour.clone();
  }
  
  
  ArrayList<Byte> ReturningSolution = new ArrayList<>();
  ReturningSolution = ConvertBestAlgorithmSolutionToByteArrayList(BestAlgorithmSolution);
  
  for(int i = 0; i < ReturningSolution.size(); i++) {
   System.out.println(ReturningSolution.get(i));
  }
  
  
  
  return ReturningSolution;
  
  
 
  
 }
 /*------------------------------------------------=END OF CODERUNNER-----------------------------------------------------*/
 
 
 
 
 
 
 //fix this up later ray
 public static void test(int numOfCities, int numOfIterations) {
  
  double Distance2DArray [][] = ReadArrayFile("C:\\Users\\Raj\\eclipse-workspace\\Lab15Ray\\src\\TSP_48.txt", " ");
  
  ArrayList<Integer> tour = new ArrayList<>();
  
  tour = PopulateCities(numOfCities);
  
  Collections.shuffle(tour);
  
  double stochasticTemperature = FitnessFunction(tour.size(),tour,Distance2DArray) / 1000;
  double annealingTemperature = FitnessFunction(tour.size(), tour, Distance2DArray) / 250;
  double coolingRate = CalculateCoolingRate(stochasticTemperature, numOfIterations);
  
  System.out.println("RMHC");
  RMHC(tour, 100000, Distance2DArray, true);
  System.out.println("This is the RMHC FITNESS " + rmhcFitness);
  
  
  System.out.println();
  
  
  System.out.println("SHC");
  SHC(tour,100000,Distance2DArray,stochasticTemperature,true);
  System.out.println("This is the SHC FITNESS " + shcFitness);
  
  
  System.out.println();
  
  System.out.println("RRHC");
  RRHC(tour,100000,Distance2DArray,true);
  System.out.println("This is the RRHC Fitness " + rrhcFitness);
  
  
  System.out.println();
  
  
  System.out.println("SA");
  SA(tour,annealingTemperature,100000,Distance2DArray,coolingRate,true);
  System.out.println("This is the SA Fitness " + saFitness);
  
  
 }
 
 
 
 
 
 
 //............................................................................................................................................
 //.DDDDDDDDD......OOOOOOO.......... NNN...NNNN....OOOOOOO.....OTTTTTTTTTT.....TTTTTTTTTTT..OOOOOOO....OUUU....UUUU....CCCCCCC...CHHH....HHHH..
 //.DDDDDDDDDD....OOOOOOOOOO........ NNNN..NNNN...OOOOOOOOOO...OTTTTTTTTTT.....TTTTTTTTTTT.OOOOOOOOOO..OUUU....UUUU...CCCCCCCCC..CHHH....HHHH..
 //.DDDDDDDDDDD..OOOOOOOOOOOO....... NNNN..NNNN..OOOOOOOOOOOO..OTTTTTTTTTT.....TTTTTTTTTTTOOOOOOOOOOOO.OUUU....UUUU..CCCCCCCCCCC.CHHH....HHHH..
 //.DDDD...DDDD..OOOOO..OOOOO....... NNNNN.NNNN..OOOOO..OOOOO.....TTTT............TTTT...TOOOO...OOOOO.OUUU....UUUU.UCCCC...CCCC.CHHH....HHHH..
 //.DDDD....DDDDDOOOO....OOOOO...... NNNNN.NNNN.NOOOO....OOOOO....TTTT............TTTT...TOOO.....OOOOOOUUU....UUUU.UCCC....CCC..CHHH....HHHH..
 //.DDDD....DDDDDOOO......OOOO...... NNNNNNNNNN.NOOO......OOOO....TTTT............TTTT...TOOO......OOOOOUUU....UUUU.UCCC.........CHHHHHHHHHHH..
 //.DDDD....DDDDDOOO......OOOO...... NNNNNNNNNN.NOOO......OOOO....TTTT............TTTT...TOOO......OOOOOUUU....UUUU.UCCC.........CHHHHHHHHHHH..
 //.DDDD....DDDDDOOO......OOOO...... NNNNNNNNNN.NOOO......OOOO....TTTT............TTTT...TOOO......OOOOOUUU....UUUU.UCCC.........CHHHHHHHHHHH..
 //.DDDD....DDDDDOOOO....OOOOO...... NNNNNNNNNN.NOOOO....OOOOO....TTTT............TTTT...TOOO.....OOOOOOUUUU...UUUU.UCCC....CCC..CHHH....HHHH..
 //.DDDD...DDDDD.OOOOO..OOOOO....... NNN.NNNNNN..OOOOO..OOOOO.....TTTT............TTTT...TOOOOO..OOOOO..UUUU..UUUUU.UCCCC...CCCCCCHHH....HHHH..
 //.DDDDDDDDDDD..OOOOOOOOOOOO....... NNN..NNNNN..OOOOOOOOOOOO.....TTTT............TTTT....OOOOOOOOOOOO..UUUUUUUUUUU..CCCCCCCCCCC.CHHH....HHHH..
 //.DDDDDDDDDD....OOOOOOOOOO........ NNN..NNNNN...OOOOOOOOOO......TTTT............TTTT.....OOOOOOOOOO...UUUUUUUUUU....CCCCCCCCC..CHHH....HHHH..
 //.DDDDDDDDD.......OOOOOO.......... NNN...NNNN.....OOOOOO........TTTT............TTTT......OOOOOOO.......UUUUUUU......CCCCCCC...CHHH....HHHH..
 //............................................................................................................................................
 
 
 public static ArrayList<Integer> PopulateCities(int numberOfCities) {
  ArrayList<Integer> tour = new ArrayList<>();
  for (int i = 0; i < numberOfCities; i++) {
   tour.add(i);
  }
  return tour;
 }
 
 
 
 public static double FitnessFunction(int NumberOfCitiesToVisit, ArrayList<Integer> Tour, double Distance [][]) {
  
  double s = 0;
  
  for(int i = 0; i < NumberOfCitiesToVisit-1; i++) {// may need to put -1 here
   int a = Tour.get(i);
   int b = Tour.get(i+1);
   s = s + Distance[a][b];
  }
  
  int End_City = Tour.get(NumberOfCitiesToVisit-1); //Ray make sure you know if the -1 should be there
  int Start_City = Tour.get(0);
  s = s + Distance[End_City][Start_City];
  
  return s;
 }
 
 public static ArrayList<Integer> RandomPerm(int n) {// FROM LECTURE SLIDE
  
  ArrayList<Integer> Premutation = new ArrayList<Integer>();
  
  for(int i = 0; i <Premutation.size(); i++) {
   Premutation.add(i);
  }
  
  
  ArrayList<Integer> Tour = new ArrayList<Integer>();
  
  while( Premutation.size() > 0) {
   int i = UI(1,Premutation.size()); //RAY SHOULD THIS BE 0 AND PREMUTATION.SIZE???
   Tour.add(Premutation.get(i));
   Premutation.remove(i);
  }
  
  return Tour;
  
 }
 
 
 
 public static ArrayList<Integer> Swap(ArrayList<Integer> Tour) {
  int i;
  int j;
  i = j = 0;// is this correct ray?
  
  while(i == j) {
   
   i = UI(0,Tour.size()-1);
   j = UI(0,Tour.size()-1);
 
  }
  
  int temp = Tour.get(i);
  int temp2 = Tour.get(j);
  
  Tour.set(i, temp2);
  Tour.set(j, temp);
  
  
  
  return Tour;
  
  
 }
 
 public static double QualityOfATSPSolution(int NumOfCities, ArrayList<Integer> tour, double Matrix [][]){
  
  
  double [][] OutputFromMSTFitnessFunction = MST.PrimsMST(Matrix);
  
  double OutputFromTSPFitnessFunction = FitnessFunction(NumOfCities,tour,Matrix);
  
  
  
  double OverallFitnessFromMST = 0;
  
  for(int i = 0; i < OutputFromMSTFitnessFunction.length; i++) {
   for(int j = 0; j <OutputFromMSTFitnessFunction[i].length; j++) {
    OverallFitnessFromMST += OutputFromMSTFitnessFunction[i][j];
   }
  }
  
  //System.out.println(OverallFitnessFromMST + "r9345=ff9r4ef");
  
  
  
  if(OutputFromTSPFitnessFunction < OverallFitnessFromMST) {
   //THIS IS NOT ALLOWED RETURN FROM METHOD OR AS SUCH
  }
  
  //double MSTDividedByTSP = OverallFitnessFromMST/OutputFromTSPFitnessFunction;
  
  double MSTDividedByTSP = OutputFromTSPFitnessFunction/OverallFitnessFromMST;
  
  double solution = MSTDividedByTSP*100;
  
  //GO TO LECTURE SLIDE NUMBER 16 SLDIE 23 AND CHECK THE LESS THAN OR EQUAL TO 100, YOU MAY NEED TO REDO THIS RAY
  //UPDATE FOR THE ABOVE COMMENT RAY, YOU HAVE COMPLETED THIS
  return solution;
 }
 
 
 public static double PR(double newFitness, double oldFitness, double temperature) {
  double changeInFitness = Math.abs(newFitness - oldFitness);
  changeInFitness = -1 * changeInFitness;
  double prScore = Math.exp(changeInFitness / temperature);
  return prScore;
 }
 
 
 static double CalculateAcceptanceProbability(double newFitness, double oldFitness, double temperature) {
  double deltaFitness = newFitness - oldFitness;
  double p = 1 / (1 + Math.exp(deltaFitness / temperature));
  return p;
 }
 
 static double CalculateEfficiencyOfMST(double fitness, double [][] Distance) {
  double[][] mstArray = MST.PrimsMST(Distance);
  
  double mstFitnessScore = CalculateFitnessOfMST(mstArray);
  double efficiency = (mstFitnessScore / fitness) * 100;
  
  //efficiency = CS2004.RoundDecimal(efficiency);
  
  return efficiency;
 }
 
 private static double CalculateFitnessOfMST(double[][] mstArray) {
  double fitness = 0.0;
  for (int i = 0; i < mstArray.length; i++) {
   for (int j = 0; j < mstArray[i].length; j++) {
    fitness = fitness + mstArray[i][j];
   }
  }
  
  return fitness;
  
//  double fitnessScore = fitness / 2;
//  return fitnessScore;
 }
 
 
 //...........................................................................................................................
 //.....AAAAA.....LLLL..........GGGGGGG......OOOOOOO.....RRRRRRRRRR..IIIII.TTTTTTTTTTHHHHH...HHHH.MMMMMM...MMMMMM..SSSSSSS....
 //.....AAAAA.....LLLL........GGGGGGGGGG....OOOOOOOOOO...RRRRRRRRRRR.IIIII.TTTTTTTTTTHHHHH...HHHH.MMMMMM...MMMMMM.SSSSSSSSS...
 //....AAAAAA.....LLLL.......GGGGGGGGGGGG..OOOOOOOOOOOO..RRRRRRRRRRR.IIIII.TTTTTTTTTTHHHHH...HHHH.MMMMMM...MMMMMM.SSSSSSSSSS..
 //....AAAAAAA....LLLL.......GGGGG..GGGGG..OOOOO..OOOOO..RRRR...RRRRRIIIII....TTTT....HHHH...HHHH.MMMMMMM.MMMMMMMSSSSS..SSSS..
 //...AAAAAAAA....LLLL......GGGGG....GGG..OOOOO....OOOOO.RRRR...RRRRRIIIII....TTTT....HHHH...HHHH.MMMMMMM.MMMMMMMSSSSS........
 //...AAAAAAAA....LLLL......GGGG..........OOOO......OOOO.RRRRRRRRRRR.IIIII....TTTT....HHHHHHHHHHH.MMMMMMM.MMMMMMM.SSSSSSS.....
 //...AAAA.AAAA...LLLL......GGGG..GGGGGGGGOOOO......OOOO.RRRRRRRRRRR.IIIII....TTTT....HHHHHHHHHHH.MMMMMMMMMMMMMMM..SSSSSSSSS..
 //..AAAAAAAAAA...LLLL......GGGG..GGGGGGGGOOOO......OOOO.RRRRRRRR....IIIII....TTTT....HHHHHHHHHHH.MMMMMMMMMMMMMMM....SSSSSSS..
 //..AAAAAAAAAAA..LLLL......GGGGG.GGGGGGGGOOOOO....OOOOO.RRRR.RRRR...IIIII....TTTT....HHHH...HHHH.MMMMMMMMMMMMMMM.......SSSS..
 //..AAAAAAAAAAA..LLLL.......GGGGG....GGGG.OOOOO..OOOOO..RRRR..RRRR..IIIII....TTTT....HHHH...HHHH.MMMM.MMMMM.MMMMSSSS....SSS..
 //.AAAA....AAAA..LLLLLLLLLL.GGGGGGGGGGGG..OOOOOOOOOOOO..RRRR..RRRRR.IIIII....TTTT....HHHH...HHHH.MMMM.MMMMM.MMMMSSSSSSSSSSS..
 //.AAAA.....AAAA.LLLLLLLLLL..GGGGGGGGGG....OOOOOOOOOO...RRRR...RRRRRIIIII....TTTT....HHHH...HHHH.MMMM.MMMMM.MMMM.SSSSSSSSSS..
 //.AAAA.....AAAA.LLLLLLLLLL....GGGGGGG.......OOOOOO.....RRRR....RRRRIIIII....TTTT....HHHH...HHHH.MMMM.MMMMM.MMMM..SSSSSSSS...
 //...........................................................................................................................
 
 
 public static ArrayList<Integer> RRHC(ArrayList<Integer> Tour, int NumberOfIterations, double [][] Distance, boolean PrintReport){
  
  ArrayList<Integer> currentTour;
  ArrayList<Integer> bestTour = new ArrayList<>();
  int numberOfRepeats = 2;
  
  rrhcFitness = FitnessFunction(Tour.size(), Tour, Distance);
  
  double fitness;
  
  for(int i = 0; i < numberOfRepeats; i++) {
   
   Collections.shuffle(Tour);
   currentTour = RMHC(Tour, NumberOfIterations/ numberOfRepeats, Distance, false);
   
   fitness = FitnessFunction(currentTour.size(), currentTour, Distance);
   
   if( fitness < rrhcFitness) {
    
    rrhcFitness = fitness;
    bestTour = (ArrayList<Integer>) currentTour.clone();
    
   }
   
   
   
  }
  
  //rrhcFitness = FitnessFunction(bestTour.size(), bestTour, Distance);
  //rrhcMSTEfficiency = CalculateEfficiencyOfMST(rrhcFitness, Distance);
  //rrhcEfficiency = QualityOfATSPSolution(bestTour.size(), bestTour, Distance);
  
  
  if(PrintReport == true) {
   System.out.println("RRHC TSP Fitness " + rrhcFitness);
   System.out.println("RRHC Optimal Efficiency: " + rrhcEfficiency + "%");
   System.out.println("RRHC MST Efficiency: " + rrhcMSTEfficiency + "%");
  }
  
  
  return bestTour;
 }
 
 
 
 
 //Random mutation hill climbing
 //UNCOMMENT COMMENTS AFTER
 public static ArrayList<Integer> RMHC (ArrayList<Integer> Tour, int NumberOfIterations, double[][] Distance, boolean PrintReport){
  
  double newFitness;
  rmhcFitness = 0;
  
  ArrayList<Integer> oldTour = new ArrayList<>();
  
  for(int i = 0; i < NumberOfIterations; i++) {
   
   oldTour = (ArrayList<Integer>) Tour.clone();
   
   rmhcFitness = FitnessFunction(Tour.size(),Tour,Distance);
   
   Tour = Swap(Tour);//RAY THIS IS CS2004.SmallChange(tour);
   
   newFitness = FitnessFunction(Tour.size(),Tour,Distance); // RAY, THIS IS shcFitness = Calculations.CalculateFitness(tour);
   
   if(newFitness > rmhcFitness) {
    Tour = oldTour;
   }
   
  }
  //RAY COMPLETE THIS AFTER
  //rmhcEfficiency = Calculations.CalculateEfficiency(rmhcFitness);
  //rmhcMSTEfficiency = Calculations.CalculateEfficiencyOfMST(rmhcFitness);
  
  // YOU DID THE BELOW RAY
  
  rmhcFitness = FitnessFunction(Tour.size(), Tour, Distance);
  rmhcMSTEfficiency = CalculateEfficiencyOfMST(rmhcFitness,Distance);
  rmhcEfficiency = QualityOfATSPSolution(Tour.size(), Tour, Distance);
  
  
  if (PrintReport == true) {
  // System.out.println("Fitness: " + rmhcFitness);
  // System.out.println("MST Efficiency: " + rmhcMSTEfficiency + "%");
   
   System.out.println("RMHC TSP Fitness" + rmhcFitness);
   System.out.println("RMHC Optimal Efficiency: " + rmhcEfficiency + "%");
   System.out.println("RMHC MST Efficiency" + rmhcMSTEfficiency + "%");
  }
  return Tour;
  
 }
 
 //Stochastic Hill Climber
 public static ArrayList<Integer> SHC(ArrayList<Integer> Tour, int numberOfIterations, double [][]Distance, double temperature, boolean printReport){
  
  double newFitness;
  shcFitness = 0;
  double p;
  
  ArrayList<Integer> oldTour;
  
   for(int i = 0; i < numberOfIterations; i++) {
    
    oldTour = (ArrayList<Integer>)Tour.clone();
    
    shcFitness = FitnessFunction(Tour.size(),Tour,Distance);
    
    Tour = Swap(Tour);
    
    newFitness = FitnessFunction(Tour.size(),Tour,Distance);
    
   
    
    p = CalculateAcceptanceProbability(newFitness, shcFitness, temperature);
    
    double random = UR(0.0, 1.0);
    
     if(p < random) {
     Tour = oldTour;
     }
    
    
   }
   //RAY FINISH THIS LATER
   shcFitness = FitnessFunction(Tour.size(), Tour, Distance);
   shcMSTEfficiency = CalculateEfficiencyOfMST(shcFitness, Distance);
   shcEfficiency = QualityOfATSPSolution(Tour.size(), Tour, Distance);
   
   if (printReport == true) {
    System.out.println("SHC TSP Fitness: " + shcFitness);
    //System.out.println("MST Efficiency: " + shcMSTEfficiency + "%");
    System.out.println("SHC Optimal Efficiency: " + shcEfficiency + "%");
    System.out.println("SHC MST Effiency" + shcMSTEfficiency + "%");
   }
   
   return Tour;
   
   
 }
 
 
 
 //Simulated annealing algorithm
 
 public static ArrayList<Integer> SA (ArrayList<Integer> Tour, double temperature, int numberOfIterations, double [][]Distance, double coolingRate, boolean printReport) {
  
  ArrayList<Integer> oldTour;
  saFitness = 0;
  double newFitness;
  double p;
  
  for (int i = 0; i < numberOfIterations; i++) {
   
   oldTour = (ArrayList<Integer>) Tour.clone();
   
   saFitness = FitnessFunction(Tour.size(),Tour,Distance);
   
   Tour = Swap(Tour);
   
   
   newFitness = FitnessFunction(Tour.size(),Tour,Distance);
   
   if (newFitness > saFitness) {
    
    p = PR(newFitness, saFitness, temperature);
    
    if(p < UR(0.0,1.0)) {
     Tour = oldTour;
    }
    else {
     //Dont do anything
    }
   }else {
     //Dont do anything
    }
    
    temperature = coolingRate * temperature;
    
  }
  
  //saEfficiency = Calculations.CalculateEfficiency(saFitness);
  //saMSTEfficiency = Calculations.CalculateEfficiencyOfMST(saFitness);
  saEfficiency = QualityOfATSPSolution(Tour.size(), Tour, Distance);
  saFitness = FitnessFunction(Tour.size(), Tour, Distance);
  saMSTEfficiency = CalculateEfficiencyOfMST(saFitness, Distance);
  
  if(printReport == true) {
  // System.out.println("Fitness: " + saFitness);
  // System.out.println("MST Efficiency: " + saMSTEfficiency + "%");
   
   System.out.println("SA TSP Fitness" + saFitness);
   System.out.println("SA Optimal Efficiency: " + saEfficiency + "%");
   System.out.println("SA MST Effiency" + saMSTEfficiency + "%");
  }
  
  return Tour;
 
 }
 
 
 
 
 //...................................................................................................................................................................
 //....CCCCCCC......OOOOOOO.....NNNN...NNNNNVVVV....VVVVVEEEEEEEEEEE.RRRRRRRRRR...TTTTTTTTTTT.... MMMMM...MMMMMM....AAAAA...ATTTTTTTTTTTRRRRRRRRRR..RIII.IXXX...XXXX..
 //...CCCCCCCCC....OOOOOOOOOO...NNNNN..NNNN.VVVV....VVVV.EEEEEEEEEEE.RRRRRRRRRRR..TTTTTTTTTTT.... MMMMM...MMMMMM....AAAAA...ATTTTTTTTTTTRRRRRRRRRRR.RIII..XXXX.XXXXX..
 //..CCCCCCCCCCC..OOOOOOOOOOOO..NNNNN..NNNN.VVVV....VVVV.EEEEEEEEEEE.RRRRRRRRRRR..TTTTTTTTTTT.... MMMMM...MMMMMM...AAAAAA...ATTTTTTTTTTTRRRRRRRRRRR.RIII..XXXXXXXXX...
 //..CCCC...CCCCC.OOOOO..OOOOO..NNNNNN.NNNN.VVVVV..VVVV..EEEE........RRRR...RRRRR....TTTT........ MMMMMM.MMMMMMM...AAAAAAA......TTTT...TRRR....RRRR.RIII...XXXXXXXX...
 //.CCCC.....CCC.OOOOO....OOOOO.NNNNNN.NNNN..VVVV..VVVV..EEEE........RRRR...RRRRR....TTTT........ MMMMMM.MMMMMMM..AAAAAAAA......TTTT...TRRR....RRRR.RIII...XXXXXXX....
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN..VVVV..VVVV..EEEEEEEEEE..RRRRRRRRRRR.....TTTT........ MMMMMM.MMMMMMM..AAAAAAAA......TTTT...TRRRRRRRRRRR.RIII....XXXXX.....
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN..VVVVVVVVV...EEEEEEEEEE..RRRRRRRRRRR.....TTTT........ MMMMMMMMMMMMMM..AAAA.AAAA.....TTTT...TRRRRRRRRRR..RIII....XXXXX.....
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN...VVVVVVVV...EEEEEEEEEE..RRRRRRRR........TTTT........ MMMMMMMMMMMMMM.AAAAAAAAAA.....TTTT...TRRRRRRRR....RIII....XXXXXX....
 //.CCCC.....CCC.OOOOO....OOOOO.NNNNNNNNNNN...VVVVVVVV...EEEE........RRRR.RRRR.......TTTT........ MMMMMMMMMMMMMM.AAAAAAAAAAA....TTTT...TRRR.RRRRR...RIII...XXXXXXXX...
 //..CCCC...CCCCC.OOOOO..OOOOO..NNNN.NNNNNN...VVVVVVV....EEEE........RRRR..RRRR......TTTT........ MMM.MMMMM.MMMM.AAAAAAAAAAA....TTTT...TRRR..RRRRR..RIII..XXXXXXXXX...
 //..CCCCCCCCCCC..OOOOOOOOOOOO..NNNN..NNNNN....VVVVVV....EEEEEEEEEEE.RRRR..RRRRR.....TTTT........ MMM.MMMMM.MMMMAAAA....AAAA....TTTT...TRRR...RRRRR.RIII..XXXX.XXXXX..
 //...CCCCCCCCCC...OOOOOOOOOO...NNNN..NNNNN....VVVVVV....EEEEEEEEEEE.RRRR...RRRRR....TTTT........ MMM.MMMMM.MMMMAAAA.....AAAA...TTTT...TRRR....RRRR.RIII.IXXXX..XXXX..
 //....CCCCCCC.......OOOOOO.....NNNN...NNNN....VVVVV.....EEEEEEEEEEE.RRRR....RRRR....TTTT........ MMM.MMMMM.MMMMAAAA.....AAAA...TTTT...TRRR.....RRRRRIIIIIXXX....XXX..
 //................................................................................................................................................................... 
 
 public static double[][] ConvertMatrix(double [][] Matrix) {
  
  double newMatrix [][] = new double [Matrix.length][Matrix[0].length];
  
  for(int i = 0; i < newMatrix.length; i++) {
   for(int j = 0; j < newMatrix[0].length; j++) {
   
    newMatrix[i][j] = Matrix[i][j];
    
   }
  }
  
  return  newMatrix;
 }
 
 public static double[][] ConvertMatrix(float [][] Matrix) {
  
  double newMatrix [][] = new double [Matrix.length][Matrix[0].length];
  
  for(int i = 0; i < newMatrix.length; i++) {
   for(int j = 0; j < newMatrix[0].length; j++) {
   
    newMatrix[i][j] = Matrix[i][j];
    
   }
  }
  
  return  newMatrix;
 }
 public static double[][] ConvertMatrix(Double [][] Matrix) {
  
  double newMatrix [][] = new double [Matrix.length][Matrix[0].length];
  
  for(int i = 0; i < newMatrix.length; i++) {
   for(int j = 0; j < newMatrix[0].length; j++) {
   
    newMatrix[i][j] = Matrix[i][j];
    
   }
  }
  
  return  newMatrix;
 }
 
 public static double[][] ConvertMatrix(Float [][] Matrix) {
  
  double newMatrix [][] = new double [Matrix.length][Matrix[0].length];
  
  for(int i = 0; i < newMatrix.length; i++) {
   for(int j = 0; j < newMatrix[0].length; j++) {
   
    newMatrix[i][j] = Matrix[i][j];
    
   }
  }
  
  return  newMatrix;
 }
 
 
 //...........................................................................................................................................................................................................
 //....CCCCCCC......OOOOOOO.....NNNN...NNNNNVVVV....VVVVVEEEEEEEEEEE.RRRRRRRRRR...TTTTTTTTTTT.... IIII.TTTTTTTTTTTEEEEEEEEEEE.RRRRRRRRRR......AAAAA...AATTTTTTTTTTTII....OOOOOOO....OONNN...NNNN...SSSSSSS....
 //...CCCCCCCCC....OOOOOOOOOO...NNNNN..NNNN.VVVV....VVVV.EEEEEEEEEEE.RRRRRRRRRRR..TTTTTTTTTTT.... IIII.TTTTTTTTTTTEEEEEEEEEEE.RRRRRRRRRRR.....AAAAA...AATTTTTTTTTTTII...OOOOOOOOOO..OONNN...NNNN..NSSSSSSSS...
 //..CCCCCCCCCCC..OOOOOOOOOOOO..NNNNN..NNNN.VVVV....VVVV.EEEEEEEEEEE.RRRRRRRRRRR..TTTTTTTTTTT.... IIII.TTTTTTTTTTTEEEEEEEEEEE.RRRRRRRRRRR....AAAAAA...AATTTTTTTTTTTII..IOOOOOOOOOOO.OONNNN..NNNN.NNSSSSSSSSS..
 //..CCCC...CCCCC.OOOOO..OOOOO..NNNNNN.NNNN.VVVVV..VVVV..EEEE........RRRR...RRRRR....TTTT........ IIII....TTTT....EEEE........RRRR...RRRRR...AAAAAAA......TTTT...TTII.IIOOO...OOOOO.OONNNNN.NNNN.NNSS...SSSS..
 //.CCCC.....CCC.OOOOO....OOOOO.NNNNNN.NNNN..VVVV..VVVV..EEEE........RRRR...RRRRR....TTTT........ IIII....TTTT....EEEE........RRRR...RRRRR..AAAAAAAA......TTTT...TTII.IIOO.....OOOOOOONNNNN.NNNN.NNSSS........
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN..VVVV..VVVV..EEEEEEEEEE..RRRRRRRRRRR.....TTTT........ IIII....TTTT....EEEEEEEEEE..RRRRRRRRRRR...AAAAAAAA......TTTT...TTII.IIOO......OOOOOONNNNNNNNNN..NSSSSSS.....
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN..VVVVVVVVV...EEEEEEEEEE..RRRRRRRRRRR.....TTTT........ IIII....TTTT....EEEEEEEEEE..RRRRRRRRRRR...AAAA.AAAA.....TTTT...TTII.IIOO......OOOOOONNNNNNNNNN...SSSSSSSS...
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN...VVVVVVVV...EEEEEEEEEE..RRRRRRRR........TTTT........ IIII....TTTT....EEEEEEEEEE..RRRRRRRR.....AAAAAAAAAA.....TTTT...TTII.IIOO......OOOOOONN.NNNNNNN.....SSSSSSS..
 //.CCCC.....CCC.OOOOO....OOOOO.NNNNNNNNNNN...VVVVVVVV...EEEE........RRRR.RRRR.......TTTT........ IIII....TTTT....EEEE........RRRR.RRRR....AAAAAAAAAAA....TTTT...TTII.IIOO.....OOOOOOONN.NNNNNNN........SSSS..
 //..CCCC...CCCCC.OOOOO..OOOOO..NNNN.NNNNNN...VVVVVVV....EEEE........RRRR..RRRR......TTTT........ IIII....TTTT....EEEE........RRRR..RRRR...AAAAAAAAAAA....TTTT...TTII.IIOOOO..OOOOO.OONN..NNNNNN.NNSS...SSSS..
 //..CCCCCCCCCCC..OOOOOOOOOOOO..NNNN..NNNNN....VVVVVV....EEEEEEEEEEE.RRRR..RRRRR.....TTTT........ IIII....TTTT....EEEEEEEEEEE.RRRR..RRRRR.RAAA....AAAA....TTTT...TTII..IOOOOOOOOOOO.OONN..NNNNNN.NNSSSSSSSSS..
 //...CCCCCCCCCC...OOOOOOOOOO...NNNN..NNNNN....VVVVVV....EEEEEEEEEEE.RRRR...RRRRR....TTTT........ IIII....TTTT....EEEEEEEEEEE.RRRR...RRRRRRAAA.....AAAA...TTTT...TTII...OOOOOOOOOO..OONN...NNNNN..NSSSSSSSSS..
 //....CCCCCCC.......OOOOOO.....NNNN...NNNN....VVVVV.....EEEEEEEEEEE.RRRR....RRRR....TTTT........ IIII....TTTT....EEEEEEEEEEE.RRRR....RRRRRAAA.....AAAA...TTTT...TTII....OOOOOOO....OONN....NNNN...SSSSSSS....
 //...........................................................................................................................................................................................................
 
 
 public static int ConvertIterations(byte iterations) {
  
  int newIterations = iterations;
  
  return newIterations;
  
 }
 
 
 public static int ConvertIterations(short iterations) {
  
  int newIterations = iterations;
  
  return newIterations;
  
 }
 
 public static int ConvertIterations(long iterations) {
  
  int newIterations = (int) iterations;
  
  return newIterations;
  
 }
 
 
 public static int ConvertIterations(Byte iterations) {
  
  int newIterations = (int) iterations;
  
  return newIterations;
  
 }
 
 public static int ConvertIterations(Short iterations) {
  
  int newIterations = iterations;
  
  return newIterations;
  
 }
 
 
 public static int ConvertIterations(Integer iterations) {
  
  int newIterations = iterations;
  
  return newIterations;
  
 }
 
 public static int ConvertIterations(Long iterations) {
  
  long iterat = iterations;
  
  int newIterations = (int) iterat;
  return newIterations;
  
 }
 
 public static int ConvertIterations(Float iterations) {
  
  float iterat = iterations;
  
  int newIterations = (int) iterat;
  return newIterations;
  
 }
 
 
 //............................................................................................................................................................................................
 //....CCCCCCC......OOOOOOO.....NNNN...NNNNNVVVV....VVVVVEEEEEEEEEEE.RRRRRRRRRR...TTTTTTTTTTT......SSSSSSS......OOOOOOO.....LLLL.......UUUU...UUUU..TTTTTTTTTTTIIII...OOOOOOO.....NNNN...NNNN..
 //...CCCCCCCCC....OOOOOOOOOO...NNNNN..NNNN.VVVV....VVVV.EEEEEEEEEEE.RRRRRRRRRRR..TTTTTTTTTTT.....SSSSSSSSS....OOOOOOOOOO...LLLL.......UUUU...UUUU..TTTTTTTTTTTIIII..OOOOOOOOOO...NNNNN..NNNN..
 //..CCCCCCCCCCC..OOOOOOOOOOOO..NNNNN..NNNN.VVVV....VVVV.EEEEEEEEEEE.RRRRRRRRRRR..TTTTTTTTTTT.....SSSSSSSSSS..OOOOOOOOOOOO..LLLL.......UUUU...UUUU..TTTTTTTTTTTIIII.OOOOOOOOOOOO..NNNNN..NNNN..
 //..CCCC...CCCCC.OOOOO..OOOOO..NNNNNN.NNNN.VVVVV..VVVV..EEEE........RRRR...RRRRR....TTTT........ SSSS..SSSS..OOOOO..OOOOO..LLLL.......UUUU...UUUU.....TTTT...TIIII.OOOOO..OOOOO..NNNNNN.NNNN..
 //.CCCC.....CCC.OOOOO....OOOOO.NNNNNN.NNNN..VVVV..VVVV..EEEE........RRRR...RRRRR....TTTT........ SSSS.......SOOOO....OOOOO.LLLL.......UUUU...UUUU.....TTTT...TIIIIIOOOO....OOOOO.NNNNNN.NNNN..
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN..VVVV..VVVV..EEEEEEEEEE..RRRRRRRRRRR.....TTTT.........SSSSSSS....SOOO......OOOO.LLLL.......UUUU...UUUU.....TTTT...TIIIIIOOO......OOOO.NNNNNNNNNNN..
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN..VVVVVVVVV...EEEEEEEEEE..RRRRRRRRRRR.....TTTT..........SSSSSSSSS.SOOO......OOOO.LLLL.......UUUU...UUUU.....TTTT...TIIIIIOOO......OOOO.NNNNNNNNNNN..
 //.CCCC.........OOOO......OOOO.NNNNNNNNNNN...VVVVVVVV...EEEEEEEEEE..RRRRRRRR........TTTT............SSSSSSS.SOOO......OOOO.LLLL.......UUUU...UUUU.....TTTT...TIIIIIOOO......OOOO.NNNNNNNNNNN..
 //.CCCC.....CCC.OOOOO....OOOOO.NNNNNNNNNNN...VVVVVVVV...EEEE........RRRR.RRRR.......TTTT...............SSSSSSOOOO....OOOOO.LLLL.......UUUU...UUUU.....TTTT...TIIIIIOOOO....OOOOO.NNNNNNNNNNN..
 //..CCCC...CCCCC.OOOOO..OOOOO..NNNN.NNNNNN...VVVVVVV....EEEE........RRRR..RRRR......TTTT........ SSS....SSSS.OOOOO..OOOOO..LLLL.......UUUU...UUUU.....TTTT...TIIII.OOOOO..OOOOO..NNNN.NNNNNN..
 //..CCCCCCCCCCC..OOOOOOOOOOOO..NNNN..NNNNN....VVVVVV....EEEEEEEEEEE.RRRR..RRRRR.....TTTT........ SSSSSSSSSSS.OOOOOOOOOOOO..LLLLLLLLLL.UUUUUUUUUUU.....TTTT...TIIII.OOOOOOOOOOOO..NNNN..NNNNN..
 //...CCCCCCCCCC...OOOOOOOOOO...NNNN..NNNNN....VVVVVV....EEEEEEEEEEE.RRRR...RRRRR....TTTT.........SSSSSSSSSS...OOOOOOOOOO...LLLLLLLLLL..UUUUUUUUU......TTTT...TIIII..OOOOOOOOOO...NNNN..NNNNN..
 //....CCCCCCC.......OOOOOO.....NNNN...NNNN....VVVVV.....EEEEEEEEEEE.RRRR....RRRR....TTTT..........SSSSSSSS......OOOOOO.....LLLLLLLLLL...UUUUUUU.......TTTT...TIIII....OOOOOO.....NNNN...NNNN..
 //............................................................................................................................................................................................
 
 public static ArrayList<Float> ConvertBestAlgorithmSolutionToFloatArrayList(ArrayList<Integer> BestAlgorithmSolution){
 
  ArrayList<Float> ArrayListDoubleWeights5 = new ArrayList<Float>();
  ArrayListDoubleWeights5 = (ArrayList<Float>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 public static ArrayList<String> ConvertBestAlgorithmSolutionToStringArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<String> ArrayListDoubleWeights5 = new ArrayList<String>();
  ArrayListDoubleWeights5 = (ArrayList<String>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 public static ArrayList<Character> ConvertBestAlgorithmSolutionToCharacterArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<Character> ArrayListDoubleWeights5 = new ArrayList<Character>();
  ArrayListDoubleWeights5 = (ArrayList<Character>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 public static ArrayList<Short> ConvertBestAlgorithmSolutionToShortArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<Short> ArrayListDoubleWeights5 = new ArrayList<Short>();
  ArrayListDoubleWeights5 = (ArrayList<Short>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 public static ArrayList<Integer> ConvertBestAlgorithmSolutionToIntegerArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<Integer> ArrayListDoubleWeights5 = new ArrayList<Integer>();
  ArrayListDoubleWeights5 = (ArrayList<Integer>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 
 public static ArrayList<Long> ConvertBestAlgorithmSolutionToLongArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<Long> ArrayListDoubleWeights5 = new ArrayList<Long>();
  ArrayListDoubleWeights5 = (ArrayList<Long>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 public static ArrayList<Byte> ConvertBestAlgorithmSolutionToByteArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<Byte> ArrayListDoubleWeights5 = new ArrayList<Byte>();
  ArrayListDoubleWeights5 = (ArrayList<Byte>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 
 public static ArrayList<Byte> ConvertBestAlgorithmSolutionToDoubleArrayList(ArrayList<Integer> BestAlgorithmSolution){
  
  ArrayList<Byte> ArrayListDoubleWeights5 = new ArrayList<Byte>();
  ArrayListDoubleWeights5 = (ArrayList<Byte>) BestAlgorithmSolution.clone();
  
  return ArrayListDoubleWeights5;
  
 
 }
 
 
 
 //................................................................................................................................................................................................................................................
 //.IIIII....GGGGGGG....NNNN...NNNN....OOOOOOO.....RRRRRRRRRR...EEEEEEEEEEE......BBBBBBBBBB...EEEEEEEEEEE.LLLL.........OOOOOOO...OOWWW..WWWWW...WWW ....MMMMMM...MMMMMM.EEEEEEEEEEE.TTTTTTTTTTHHHHH...HHHH....OOOOOOO.....DDDDDDDDD.....SSSSSSS....
 //.IIIII..GGGGGGGGGG...NNNNN..NNNN...OOOOOOOOOO...RRRRRRRRRRR..EEEEEEEEEEE......BBBBBBBBBBB..EEEEEEEEEEE.LLLL........OOOOOOOOOO..OWWW..WWWWW..WWWW.....MMMMMM...MMMMMM.EEEEEEEEEEE.TTTTTTTTTTHHHHH...HHHH...OOOOOOOOOO...DDDDDDDDDD...SSSSSSSSS...
 //.IIIII.GGGGGGGGGGGG..NNNNN..NNNN..OOOOOOOOOOOO..RRRRRRRRRRR..EEEEEEEEEEE......BBBBBBBBBBB..EEEEEEEEEEE.LLLL.......OOOOOOOOOOOO.OWWW..WWWWWW.WWWW.....MMMMMM...MMMMMM.EEEEEEEEEEE.TTTTTTTTTTHHHHH...HHHH..OOOOOOOOOOOO..DDDDDDDDDDD..SSSSSSSSSS..
 //.IIIII.GGGGG..GGGGG..NNNNNN.NNNN..OOOOO..OOOOO..RRRR...RRRRR.EEEE.............BBBB...BBBB..EEEE........LLLL.......OOOOO..OOOOO.OWWW.WWWWWWW.WWWW.....MMMMMMM.MMMMMMM.EEEE...........TTTT....HHHH...HHHH..OOOOO..OOOOO..DDDD...DDDD.SSSSS..SSSS..
 //.IIIIIIGGGG....GGG...NNNNNN.NNNN.NOOOO....OOOOO.RRRR...RRRRR.EEEE.............BBBB...BBBB..EEEE........LLLL......LOOOO....OOOOOOWWW.WWWWWWW.WWWW.....MMMMMMM.MMMMMMM.EEEE...........TTTT....HHHH...HHHH.OOOOO....OOOOO.DDDD....DDDDSSSSS........
 //.IIIIIIGGG...........NNNNNNNNNNN.NOOO......OOOO.RRRRRRRRRRR..EEEEEEEEEE.......BBBBBBBBBBB..EEEEEEEEEE..LLLL......LOOO......OOOO.WWWWWWWWWWW.WWW......MMMMMMM.MMMMMMM.EEEEEEEEEE.....TTTT....HHHHHHHHHHH.OOOO......OOOO.DDDD....DDDD.SSSSSSS.....
 //.IIIIIIGGG..GGGGGGGG.NNNNNNNNNNN.NOOO......OOOO.RRRRRRRRRRR..EEEEEEEEEE.......BBBBBBBBBB...EEEEEEEEEE..LLLL......LOOO......OOOO.WWWWWWW.WWWWWWW......MMMMMMMMMMMMMMM.EEEEEEEEEE.....TTTT....HHHHHHHHHHH.OOOO......OOOO.DDDD....DDDD..SSSSSSSSS..
 //.IIIIIIGGG..GGGGGGGG.NNNNNNNNNNN.NOOO......OOOO.RRRRRRRR.....EEEEEEEEEE.......BBBBBBBBBBB..EEEEEEEEEE..LLLL......LOOO......OOOO.WWWWWWW.WWWWWWW......MMMMMMMMMMMMMMM.EEEEEEEEEE.....TTTT....HHHHHHHHHHH.OOOO......OOOO.DDDD....DDDD....SSSSSSS..
 //.IIIIIIGGGG.GGGGGGGG.NNNNNNNNNNN.NOOOO....OOOOO.RRRR.RRRR....EEEE.............BBBB....BBBB.EEEE........LLLL......LOOOO....OOOOO.WWWWWWW.WWWWWWW......MMMMMMMMMMMMMMM.EEEE...........TTTT....HHHH...HHHH.OOOOO....OOOOO.DDDD....DDDD.......SSSS..
 //.IIIII.GGGGG....GGGG.NNNN.NNNNNN..OOOOO..OOOOO..RRRR..RRRR...EEEE.............BBBB....BBBB.EEEE........LLLL.......OOOOO..OOOOO..WWWWWWW.WWWWWWW......MMMM.MMMMM.MMMM.EEEE...........TTTT....HHHH...HHHH..OOOOO..OOOOO..DDDD...DDDDDSSSS....SSS..
 //.IIIII.GGGGGGGGGGGG..NNNN..NNNNN..OOOOOOOOOOOO..RRRR..RRRRR..EEEEEEEEEEE......BBBBBBBBBBBB.EEEEEEEEEEE.LLLLLLLLLL.OOOOOOOOOOOO...WWWWW...WWWWW.......MMMM.MMMMM.MMMM.EEEEEEEEEEE....TTTT....HHHH...HHHH..OOOOOOOOOOOO..DDDDDDDDDDD.SSSSSSSSSSS..
 //.IIIII..GGGGGGGGGG...NNNN..NNNNN...OOOOOOOOOO...RRRR...RRRRR.EEEEEEEEEEE......BBBBBBBBBBB..EEEEEEEEEEE.LLLLLLLLLL..OOOOOOOOOO....WWWWW...WWWWW.......MMMM.MMMMM.MMMM.EEEEEEEEEEE....TTTT....HHHH...HHHH...OOOOOOOOOO...DDDDDDDDDD...SSSSSSSSSS..
 //.IIIII....GGGGGGG....NNNN...NNNN.....OOOOOO.....RRRR....RRRR.EEEEEEEEEEE......BBBBBBBBBB...EEEEEEEEEEE.LLLLLLLLLL....OOOOOO......WWWWW...WWWWW.......MMMM.MMMMM.MMMM.EEEEEEEEEEE....TTTT....HHHH...HHHH.....OOOOOO.....DDDDDDDDD.....SSSSSSSS...
 //................................................................................................................................................................................................................................................
 
 // Not really sure how this works but it does the job 
 static double CalculateCoolingRate(double startingTemperature, int numberOfIterations) {
  double tIter, tValue, coolingRate, powerValue;
  int iter;
  tIter = 0.001;
  iter = numberOfIterations;
  tValue = tIter / startingTemperature;
  powerValue = 1.0 / iter;
  coolingRate = Math.pow(tValue, powerValue);
  return coolingRate;
 }
 
 static public int UI(int aa,int bb)//stolen from one of the labs
 {
   Random rand = new Random();
  int a = Math.min(aa,bb);
  int b = Math.max(aa,bb);
  if (rand == null) 
  {
   rand = new Random();
   rand.setSeed(System.nanoTime());
  }
  int d = b - a + 1;
  int x = rand.nextInt(d) + a;
  return(x);
 }
 
 public static double UR(double a, double b) {
  Random rand = new Random();
  if (rand == null) {
   rand = new Random();
   rand.setSeed(System.nanoTime());
  }
  return ((b - a) * rand.nextDouble() + a);
 }
 
 //Print a 2D double array to the console Window
 static public void PrintArray(double x[][])
 {
  for(int i=0;i<x.length;++i)
  {
   for(int j=0;j<x[i].length;++j)
   {
    System.out.print(x[i][j]);
    System.out.print(" ");
   }
   System.out.println();
  }
 }
 //This method reads in a text file and parses all of the numbers in it
 //This method is for reading in a square 2D numeric array from a text file
 //This code is not very good and can be improved!
 //But it should work!!!
 //'sep' is the separator between columns
 static public double[][] ReadArrayFile(String filename,String sep)
 {
  double res[][] = null;
  try
  {
   BufferedReader input = null;
   input = new BufferedReader(new FileReader(filename));
   String line = null;
   int ncol = 0;
   int nrow = 0;
   
   while ((line = input.readLine()) != null) 
   {
    ++nrow;
    String[] columns = line.split(sep);
    ncol = Math.max(ncol,columns.length);
   }
   res = new double[nrow][ncol];
   input = new BufferedReader(new FileReader(filename));
   int i=0,j=0;
   while ((line = input.readLine()) != null) 
   {
    
    String[] columns = line.split(sep);
    for(j=0;j<columns.length;++j)
    {
     res[i][j] = Double.parseDouble(columns[j]);
    }
    ++i;
   }
  }
  catch(Exception E)
  {
   System.out.println("+++ReadArrayFile: "+E.getMessage());
  }
     return(res);
 }
 //This method reads in a text file and parses all of the numbers in it
 //This code is not very good and can be improved!
 //But it should work!!!
 //It takes in as input a string filename and returns an array list of Integers
 static public ArrayList<Integer> ReadIntegerFile(String filename)
 {
  ArrayList<Integer> res = new ArrayList<Integer>();
  Reader r;
  try
  {
   r = new BufferedReader(new FileReader(filename));
   StreamTokenizer stok = new StreamTokenizer(r);
   stok.parseNumbers();
   stok.nextToken();
   while (stok.ttype != StreamTokenizer.TT_EOF) 
   {
    if (stok.ttype == StreamTokenizer.TT_NUMBER)
    {
     res.add((int)(stok.nval));
    }
    stok.nextToken();
   }
  }
  catch(Exception E)
  {
   System.out.println("+++ReadIntegerFile: "+E.getMessage());
  }
     return(res);
 }
 
 
 
 
 
 
 
 
}//END OF TSP CLASS

//.........................................................................................................................
//.EEEEEEEEEEE.ENNN...NNNN..NDDDDDDDD...........OOOOOOO.....OFFFFFFFFF.FFFFFFFFFF......TTTTTTTTTTT.SSSSSSS....PPPPPPPPP....
//.EEEEEEEEEEE.ENNNN..NNNN..NDDDDDDDDD.........OOOOOOOOOO...OFFFFFFFFF.FFFFFFFFFF......TTTTTTTTTTTSSSSSSSSS...PPPPPPPPPP...
//.EEEEEEEEEEE.ENNNN..NNNN..NDDDDDDDDDD....... OOOOOOOOOOO..OFFFFFFFFF.FFFFFFFFFF......TTTTTTTTTTTSSSSSSSSSS..PPPPPPPPPPP..
//.EEEE........ENNNNN.NNNN..NDDD...DDDD....... OOOO..OOOOO..OFFF.......FFFF...............TTTT...TSSSS..SSSS..PPPP...PPPP..
//.EEEE........ENNNNN.NNNN..NDDD....DDDD..... OOO....OOOOO.OFFF.......FFFF...............TTTT...TSSSS........PPPP...PPPP..
//.EEEEEEEEEE..ENNNNNNNNNN..NDDD....DDDD..... OO......OOOO.OFFFFFFFF..FFFFFFFFF..........TTTT....SSSSSSS.....PPPPPPPPPPP..
//.EEEEEEEEEE..ENNNNNNNNNN..NDDD....DDDD..... OO......OOOO.OFFFFFFFF..FFFFFFFFF..........TTTT.....SSSSSSSSS..PPPPPPPPPP...
//.EEEEEEEEEE..ENNNNNNNNNN..NDDD....DDDD..... OO......OOOO.OFFFFFFFF..FFFFFFFFF..........TTTT.......SSSSSSS..PPPPPPPPP....
//.EEEE........ENNNNNNNNNN..NDDD....DDDD..... OOO....OOOOO.OFFF.......FFFF...............TTTT..........SSSSS.PPPP.........
//.EEEE........ENNN.NNNNNN..NDDD...DDDDD...... OOOO..OOOOO..OFFF.......FFFF...............TTTT...TSSS....SSSS.PPPP.........
//.EEEEEEEEEEE.ENNN..NNNNN..NDDDDDDDDDD....... OOOOOOOOOOO..OFFF.......FFFF...............TTTT...TSSSSSSSSSSS.PPPP.........
//.EEEEEEEEEEE.ENNN..NNNNN..NDDDDDDDDD.........OOOOOOOOOO...OFFF.......FFFF...............TTTT....SSSSSSSSSS..PPPP.........
//.EEEEEEEEEEE.ENNN...NNNN..NDDDDDDDD............OOOOOO.....OFFF.......FFFF...............TTTT.....SSSSSSSS...PPPP.........
//.........................................................................................................................

//....................................................................................................................................................................................................................
//.IIIII....GGGGGGG....NNNN...NNNN....OOOOOOO.....RRRRRRRRRR...EEEEEEEEEEE......BBBBBBBBBB...EEEEEEEEEEE.LLLL.........OOOOOOO...OOWWW..WWWWW...WWW .......CCCCCCC....LLLL..........AAAAA......SSSSSSS.....SSSSSSS.....
//.IIIII..GGGGGGGGGG...NNNNN..NNNN...OOOOOOOOOO...RRRRRRRRRRR..EEEEEEEEEEE......BBBBBBBBBBB..EEEEEEEEEEE.LLLL........OOOOOOOOOO..OWWW..WWWWW..WWWW.......CCCCCCCCC...LLLL..........AAAAA.....SSSSSSSSS...SSSSSSSSS....
//.IIIII.GGGGGGGGGGGG..NNNNN..NNNN..OOOOOOOOOOOO..RRRRRRRRRRR..EEEEEEEEEEE......BBBBBBBBBBB..EEEEEEEEEEE.LLLL.......OOOOOOOOOOOO.OWWW..WWWWWW.WWWW......CCCCCCCCCCC..LLLL.........AAAAAA.....SSSSSSSSSS..SSSSSSSSSS...
//.IIIII.GGGGG..GGGGG..NNNNNN.NNNN..OOOOO..OOOOO..RRRR...RRRRR.EEEE.............BBBB...BBBB..EEEE........LLLL.......OOOOO..OOOOO.OWWW.WWWWWWW.WWWW......CCCC...CCCCC.LLLL.........AAAAAAA...SSSSS..SSSS.SSSSS..SSSS...
//.IIIIIIGGGG....GGG...NNNNNN.NNNN.NOOOO....OOOOO.RRRR...RRRRR.EEEE.............BBBB...BBBB..EEEE........LLLL......LOOOO....OOOOOOWWW.WWWWWWW.WWWW.....CCCC.....CCC..LLLL........AAAAAAAA...SSSSS.......SSSSS.........
//.IIIIIIGGG...........NNNNNNNNNNN.NOOO......OOOO.RRRRRRRRRRR..EEEEEEEEEE.......BBBBBBBBBBB..EEEEEEEEEE..LLLL......LOOO......OOOO.WWWWWWWWWWW.WWW......CCCC..........LLLL........AAAAAAAA....SSSSSSS.....SSSSSSS......
//.IIIIIIGGG..GGGGGGGG.NNNNNNNNNNN.NOOO......OOOO.RRRRRRRRRRR..EEEEEEEEEE.......BBBBBBBBBB...EEEEEEEEEE..LLLL......LOOO......OOOO.WWWWWWW.WWWWWWW......CCCC..........LLLL........AAAA.AAAA....SSSSSSSSS...SSSSSSSSS...
//.IIIIIIGGG..GGGGGGGG.NNNNNNNNNNN.NOOO......OOOO.RRRRRRRR.....EEEEEEEEEE.......BBBBBBBBBBB..EEEEEEEEEE..LLLL......LOOO......OOOO.WWWWWWW.WWWWWWW......CCCC..........LLLL.......AAAAAAAAAA......SSSSSSS.....SSSSSSS...
//.IIIIIIGGGG.GGGGGGGG.NNNNNNNNNNN.NOOOO....OOOOO.RRRR.RRRR....EEEE.............BBBB....BBBB.EEEE........LLLL......LOOOO....OOOOO.WWWWWWW.WWWWWWW......CCCC.....CCC..LLLL.......AAAAAAAAAAA........SSSSS.......SSSSS..
//.IIIII.GGGGG....GGGG.NNNN.NNNNNN..OOOOO..OOOOO..RRRR..RRRR...EEEE.............BBBB....BBBB.EEEE........LLLL.......OOOOO..OOOOO..WWWWWWW.WWWWWWW.......CCCC...CCCCC.LLLL.......AAAAAAAAAAA.SSSS....SSSSSSSS....SSSS..
//.IIIII.GGGGGGGGGGGG..NNNN..NNNNN..OOOOOOOOOOOO..RRRR..RRRRR..EEEEEEEEEEE......BBBBBBBBBBBB.EEEEEEEEEEE.LLLLLLLLLL.OOOOOOOOOOOO...WWWWW...WWWWW........CCCCCCCCCCC..LLLLLLLLLLAAAA....AAAA.SSSSSSSSSSSSSSSSSSSSSSSS..
//.IIIII..GGGGGGGGGG...NNNN..NNNNN...OOOOOOOOOO...RRRR...RRRRR.EEEEEEEEEEE......BBBBBBBBBBB..EEEEEEEEEEE.LLLLLLLLLL..OOOOOOOOOO....WWWWW...WWWWW.........CCCCCCCCCC..LLLLLLLLLLAAAA.....AAAA.SSSSSSSSSS..SSSSSSSSSS...
//.IIIII....GGGGGGG....NNNN...NNNN.....OOOOOO.....RRRR....RRRR.EEEEEEEEEEE......BBBBBBBBBB...EEEEEEEEEEE.LLLLLLLLLL....OOOOOO......WWWWW...WWWWW..........CCCCCCC....LLLLLLLLLLAAAA.....AAAA..SSSSSSSS....SSSSSSSS....
//....................................................................................................................................................................................................................

class CompareEdge implements java.util.Comparator 
{
 public int compare(Object a, Object b) 
 {
  if (((Edge)a).w < ((Edge)b).w) return(-1);
  if (((Edge)a).w > ((Edge)b).w) return(1);
  return(0);
 }
}
class Edge extends Object
{
 public int i,j;
 public double w;
 Edge(int ii,int jj,double ww)
 {
  i = ii;
  j = jj;
  w = ww;
 };
 public void Print()
 {
  System.out.print("(");
  System.out.print(i);
  System.out.print(",");
  System.out.print(j);
  System.out.print(",");
  System.out.print(w);
  System.out.print(")");
 }
};
class MST
{
 //Search for the next applicable edge
 static private Edge LocateEdge(ArrayList<Integer> v,ArrayList<Edge> edges)
 {
  for (Iterator<Edge> it = edges.iterator(); it.hasNext();)
  {
         Edge e = it.next();
   int x = e.i;
   int y = e.j;
   int xv = v.indexOf(x);
   int yv = v.indexOf(y);
   if (xv > -1 && yv == -1)
   {
    return(e);
   }
   if (xv == -1 && yv > -1)
   {
    return(e);
   }
  }
  //Error condition
  return(new Edge(-1,-1,0.0));
 }
 
 //d is a distance matrix, high value edges are more costly
 //Assume that d is symmetric and square
 public static double[][] PrimsMST(double[][] d)
 {
  int i,j,n = d.length;
  double res[][] = new double[n][n];
  //Store edges as an ArrayList
  ArrayList<Edge> edges = new ArrayList<Edge>();
  for(i=0;i<n-1;++i)
  {
   for(j=i+1;j<n;++j)
   {
    //Only non zero edges
    if (d[i][j] != 0.0) edges.add(new Edge(i,j,d[i][j]));
   }
  }
  //Sort the edges by weight
  Collections.sort(edges,new CompareEdge());
  //Don't do anything more if all the edges are zero
  if (edges.size() == 0) return(res);
  //List of variables that have been allocated
  ArrayList<Integer> v = new ArrayList<Integer>();
  //Pick cheapest edge
  v.add(edges.get(0).i);
  //Loop while there are still nodes to connect
  while(v.size() != n)
  {
   Edge e = LocateEdge(v,edges);
   if (v.indexOf(e.i) == -1) v.add(e.i);
   if (v.indexOf(e.j) == -1) v.add(e.j);
   res[e.i][e.j] = e.w;
   res[e.j][e.i] = e.w;
  }
  return(res);
 }
}
