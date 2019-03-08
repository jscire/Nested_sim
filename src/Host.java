import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Random;


/**
 *
 */
public class Host {

	public static final boolean backwardMutation = true; //are backward mutations allowed? (with same rate mu as forward mutations)



	public static double r0 = 0.6; // wt growth-rate without treatment
	public static double flatness = 15;
	public static double alpha = 0.017; // fitnessCost
	public static double theta1 = 0.3; // a parameter for the shape of the dose-response curve for growth of wild-type strain
	public static double theta2 = 0.6; // a parameter for the shape of the dose-response curve for growth of sensitive strain

	public double r; // effective growth-rate of wt strain
	public double rm; // effective growth-rate of mutant strain

	// parameters for the within-host model
	public static double eta = 0.02; // death rate of the pathogen
	public static double mu = 0.002; // mutation rate towards resistance // was 0.01 in Day and Read // changed from 0.001 to 0.002 to accommodate the change in how mutations occur
	public static double delta = 0.01; // death rate of the immune-system
	public static double kappa = 0.0225; // rate of clearance of pathogen by the immune-system // was 0.0175
	public static double lambda = 0.035; // growth "rate" of the immune-system

	public static double tau = 0.005; // step-size of the tau-leaping algorithm

	public static final int transmittedLoad = 7; // pathogen load transmitted to an infected target
	public static final int basalImmuneSystemLevel = 7;

	public static double infectiosityThreshold = 50; // threshold for the pathogen load above which patients can infect other people
	public static double symptomsThreshold = 100; // threshold for the pathogen load above which patients suffer from symptoms
	public static final double treatmentStartDelay = 1; // how many days after symptoms have set in should a patient wait to receive treatment?
	public static double treatmentLength = Double.POSITIVE_INFINITY; // how many days does the standard treatment last? // was 7 days

	public final int ID;
	public double immuneSystemLevel;
	public double wtPathogenLevel = 0;
	public double mutantPathogenLevel = 0;

	public boolean treated = false;
	public boolean toBeTreated = false;
	public boolean recovered = false;
	public double dose = 0;

	public double treatmentStart = 0;
	public double treatmentEnd = 0;
	public double age = 0;

	public double[] nbEventsInOneStep;

	public Random randomValueGenerator = new Random();// trying that, but may be faster to use the same as the one in the population

	public Host(int ID, double wt, double mut, double immuneLevel){
		this.immuneSystemLevel = immuneLevel;
		this.wtPathogenLevel = wt;
		this.mutantPathogenLevel = mut;
		this.ID = ID;

		randomValueGenerator = new Random();
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];

		setPathogenGrowthRate();
	}

	public Host(int ID) {
		this.immuneSystemLevel = basalImmuneSystemLevel;
		this.ID = ID;
		randomValueGenerator = new Random();
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];

		setPathogenGrowthRate();
	}

	public Host(int ID, Random randGen) {
		this.immuneSystemLevel = basalImmuneSystemLevel;
		this.ID = ID;
		this.randomValueGenerator = randGen;
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];

		setPathogenGrowthRate();
	}

	public Host(int ID, double dose) {
		this.immuneSystemLevel = basalImmuneSystemLevel;
		this.ID = ID;

		this.setDose(dose);
		randomValueGenerator = new Random();
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];
	}

	public Host(int ID, double dose, Random randGen) {
		this.immuneSystemLevel = basalImmuneSystemLevel;
		this.ID = ID;

		this.setDose(dose);
		this.randomValueGenerator = randGen;
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];
	}

	public Host(int ID, double wt, double mut, double immuneLevel, double dose){
		this.ID = ID;
		this.immuneSystemLevel = immuneLevel;
		this.wtPathogenLevel = wt;
		this.mutantPathogenLevel = mut;

		this.dose = dose;
		setPathogenGrowthRate();

		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else
			nbEventsInOneStep = new double[7];
	}

	public void setDose(double d){
		this.dose = d;
		setPathogenGrowthRate();
	}
	public void setToBeTreated(){
		toBeTreated = true;
		treatmentStart = age + treatmentStartDelay;
		treatmentEnd = treatmentStart + treatmentLength;
	}

	public void setTreated(boolean t) {
		this.treated = t;
		if(treated)
			toBeTreated = false;

		setPathogenGrowthRate();
	}

	public void setPathogenGrowthRate(){
		this.r = this.treated? r0/2*(1-Math.tanh(flatness*(dose-theta1))) : r0; // growth-rate for wt-strain with possible treatment
		this.rm = this.treated? r0*(1-alpha)/2*(1-Math.tanh(flatness*(dose-theta2))) : r0*(1-alpha); // growth-rate for mutant-strain with possible treatment
	}

	/**
	 * Simulate one time step of dynamics for the host.
	 * return 0 if no change of state occurred.
	 * return -1 if the host recovered
	 * return +1 if the host went from non-infectious to infectious or from infectious to non-infectious
	 */
	public int simulateOneTimeStep(boolean infectious){
		// update the host's age
		age += tau;

		this.calculateNumberOfEventsInOneStep();

		// to do if a level gets below zero, set it to zero

		wtPathogenLevel += (nbEventsInOneStep[0] - nbEventsInOneStep[2]);
		mutantPathogenLevel += (nbEventsInOneStep[1] + nbEventsInOneStep[6] - nbEventsInOneStep[3]);
		immuneSystemLevel += (nbEventsInOneStep[4] - nbEventsInOneStep[5]);

		if(backwardMutation) {
			wtPathogenLevel += nbEventsInOneStep[7];
		}


		if(!toBeTreated && !treated && (wtPathogenLevel+mutantPathogenLevel)>= Host.symptomsThreshold)
			this.setToBeTreated(); // schedule treatment start
		else if(toBeTreated && age > treatmentStart && age < treatmentEnd)
			this.setTreated(true); // time to start treatment
		else if(treated && age > treatmentEnd && (wtPathogenLevel+mutantPathogenLevel)> Host.symptomsThreshold)	
			treatmentEnd = age + treatmentLength; // renew treatment
		else if(treated && age > treatmentEnd && (wtPathogenLevel+mutantPathogenLevel)<= Host.symptomsThreshold)
			this.setTreated(false); // stop treatment


		if(wtPathogenLevel <1 && mutantPathogenLevel <1)
			return -1; 	
		else if ( !infectious && (wtPathogenLevel+mutantPathogenLevel)>= Host.infectiosityThreshold)
			return 1;
		else if (infectious && (wtPathogenLevel+mutantPathogenLevel)< Host.infectiosityThreshold)
			return 1;
		else
			return 0;

	}

	public int simulateAndGetEmergenceResistance(){
		
		while(this.simulateOneTimeStep(true) >= 0){
			if(mutantPathogenLevel >= Host.infectiosityThreshold)
				return 1;
		}
		return 0;
	}
	
	public double simulateAndGetTimeSpentInfectious(){
		double timeSpentInfectious = 0;
		while(this.simulateOneTimeStep(true) >= 0){
			if((wtPathogenLevel+mutantPathogenLevel)>= Host.infectiosityThreshold)
				timeSpentInfectious += tau;
		}
		return timeSpentInfectious;
	}

	/**
	 * Calculates the within-host changes in one step
	 */
	public void calculateNumberOfEventsInOneStep(){

		nbEventsInOneStep[0] = getPoisson(tau * r * (1-mu)*wtPathogenLevel); // +1 to wtPathogenLevel
		nbEventsInOneStep[1] = backwardMutation? getPoisson(tau * rm * (1-mu) * mutantPathogenLevel): getPoisson(tau * rm * mutantPathogenLevel); // +1 to mutantPathogenLevel
		nbEventsInOneStep[2] = getPoisson(tau * (eta + kappa* immuneSystemLevel)* wtPathogenLevel); // -1 to wtPathogenLevel
		nbEventsInOneStep[3] = getPoisson(tau * (eta + kappa* immuneSystemLevel)* mutantPathogenLevel); // -1 to mutantPathogenLevel
		nbEventsInOneStep[4] = getPoisson(tau * lambda * (mutantPathogenLevel + wtPathogenLevel)); // +1 to immuneSystemLevel
		nbEventsInOneStep[5] = getPoisson(tau * delta * immuneSystemLevel); // -1 to immuneSystemLevel
		nbEventsInOneStep[6] = getPoisson(tau * r * mu * wtPathogenLevel); // +1 to mutantPathogenLevel

		if(backwardMutation)
			nbEventsInOneStep[7] = getPoisson(tau * rm * mu * mutantPathogenLevel); // +1 to wtPathogenLevel
	}

	/**
	 * Based on Donald E. Knuth (1969). Seminumerical Algorithms. The Art of Computer Programming, Volume 2. Addison Wesley.
	 * Adapted from http://stackoverflow.com/questions/9832919/generate-poisson-arrival-in-java#9832977
	 * @param mean
	 * @return
	 */
	public  int getPoisson(double mean) {
		double L = Math.exp(-mean);
		int k = 0;
		double p = 1.0;
		do {
			p = p * randomValueGenerator.nextDouble();
			k++;
		} while (p > L);
		return k - 1;
	}

	public boolean isSymptomatic(){
		return (this.wtPathogenLevel + this.mutantPathogenLevel)>symptomsThreshold ;
	}

	/**
	 * Represents an intial infection with wild-type strain,
	 * it is assumed that the host is not infected before this infection.
	 * @throws Exception 
	 */
	public void getInfectedByWildType() throws Exception {
		if( (wtPathogenLevel+mutantPathogenLevel)> 0)
			throw new Exception("Host is already infected.");

		wtPathogenLevel += transmittedLoad;
		mutantPathogenLevel = 0;
	}

	public boolean getInfectedByFullLoadAndIsInfectedByRes(Host John){
		if ( randomValueGenerator.nextDouble()> John.wtPathogenLevel/(John.wtPathogenLevel + John.mutantPathogenLevel)) {
			mutantPathogenLevel += transmittedLoad; //Infected host gets the full pathogen load from the resistant strain
			return true;
		}else{
			wtPathogenLevel += transmittedLoad;
			return false;
		}
	}

	public boolean equals(Host e){	
		return (this.ID == e.ID);
	}

	public String toString(){
		DecimalFormat df = new DecimalFormat("###0");

		return ("\nHost # " + this.ID
				+ "\nWild-type pathogen level: " + df.format(this.wtPathogenLevel)
				+ "\nMutant pathogen level: " + df.format(this.mutantPathogenLevel)
				+ "\nImmune-system level: " + df.format(this.immuneSystemLevel));

	}

	public static void main(String[] args) {

		System.out.println(Double.MIN_VALUE + "\t" + Double.MAX_VALUE);
		
		Host john  = new Host(0, 7, 0, 7, 0.4);
		int nbOfEmerged = 0;
		for(int i = 0; i < 1000; i++){
			nbOfEmerged += john.simulateAndGetEmergenceResistance();
			john  = new Host(0, 7, 0, 7, 0.4);
			if(i % 1000 ==0)
				System.out.println(i + "\t" + nbOfEmerged);
		}
		System.out.println("Emergence of resistance in " + nbOfEmerged + " simuls.");

	}
}
