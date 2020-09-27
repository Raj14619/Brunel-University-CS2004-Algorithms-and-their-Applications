import java.util.ArrayList;
import java.util.Random;

public class CannonSolution {
	
	

	public static ArrayList<Double> SmallChangeOperator(double RandomAngle, double RandomMuzzleVelocity, int TargetRange, int NumberOfLoopsToDo ) {
		
		Random Random = new Random();
		
		final double HighesAngleValueMinusLowestAngleValue = 55-25;
		double valueTobeAddedToAngle = 0;
		double percentageValueToAddOnToAngle = 1;
	
		
		
		final double HighestVelocityMinusLowestVelocityValue = 1650-1500;
		double percentageValueToAddOnToVelocity = 1;
		double valueToBeAddedToVelocity = 0;
		
		
		double coinToss;
		
		double NewRandomAngleValue = RandomAngle;
		double NewRandomMuzzleVelocity  = RandomMuzzleVelocity;
		
		double range = Double.MAX_VALUE;
		double fitness = Double.MAX_VALUE;
		
		
		
		
		double bestVelocity = Double.MAX_VALUE;
		double bestRange = Double.MAX_VALUE;
		double bestAngle = Double.MAX_VALUE;
		
		
		for(int i = 0; i < NumberOfLoopsToDo; i++) {
		
			
			
		coinToss = Math.random();
		
		if(coinToss <= 0.495) { 
			
			
			valueTobeAddedToAngle = ((HighesAngleValueMinusLowestAngleValue/100)*(percentageValueToAddOnToAngle)); //
			
			NewRandomAngleValue = RandomAngle+valueTobeAddedToAngle;
			
			
			
			
			percentageValueToAddOnToAngle++;

			if(NewRandomAngleValue > 55) {
				NewRandomAngleValue = 25;
				percentageValueToAddOnToAngle = 1;
			}
			if(NewRandomAngleValue < 25) {
				NewRandomAngleValue = 55;
				percentageValueToAddOnToAngle = 1;
			}
			
			range = Cannon.GetMaxRange(NewRandomAngleValue, NewRandomMuzzleVelocity);
			
		}
			
	
		
		if(coinToss >= 0.496) { // if cointoss is 0 we will change RandomAngle
			
			valueToBeAddedToVelocity = ((HighestVelocityMinusLowestVelocityValue/100)*(percentageValueToAddOnToVelocity));// divide by 100 and multiply by 1 to get 1 percentage to then be added on
			
			NewRandomMuzzleVelocity = RandomMuzzleVelocity+valueToBeAddedToVelocity;
		
			percentageValueToAddOnToVelocity++;
			
			if(NewRandomMuzzleVelocity > 1650) {
				NewRandomMuzzleVelocity = 1500;
				percentageValueToAddOnToVelocity = 1;
			}
			
			if(NewRandomMuzzleVelocity < 1500) {
				NewRandomMuzzleVelocity = 1650;
				percentageValueToAddOnToVelocity = 1;
			}

			range = Cannon.GetMaxRange(NewRandomAngleValue, NewRandomMuzzleVelocity);
			
		}
		
		if (range < bestRange ) {
			bestRange = range;
			bestVelocity = NewRandomMuzzleVelocity;
			bestAngle = NewRandomAngleValue;
		}
		
		
		}
		
		ArrayList<Double> Output =new ArrayList<Double>();
		Output.add(bestRange);
		Output.add(bestVelocity);	
		Output.add(bestAngle);
		
		return Output;
	}
	
	

}
