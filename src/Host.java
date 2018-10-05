import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Random;


/**
 *
 */
public class Host {

	public static final boolean backwardMutation = true; //are backward mutations allowed? (with same rate mu as forward mutations)

	static double fitnessCost = 0;

	// to do: i may want to add an extra "treatment class" to allow for taking the dose response curve shape as a parameter
	public final double theta1 = 0.3; // a parameter for the shape of the dose-response curve for growth of wild-type strain
	public final double theta2 = 0.6; // a parameter for the shape of the dose-response curve for growth of sensitive strain

	public final double r0 = 0.6; // wt growth-rate without treatment
	public final double rm0 = r0 - fitnessCost; // mutant-strain growth-rate without treatment // was 0.59

	public double r; // effective growth-rate of wt strain
	public double rm; // effective growth-rate of mutant strain

	// parameters for the within-host model
	public final double eta = 0.02; // death rate of the pathogen 
	public static double mu = 0.001; // mutation rate towards resistance // was 0.01
	public final double delta = 0.01; // death rate of the immune-system
	public final double kappa = 0.0225; // rate of clearance of pathogen by the immune-system // was 0.0175 
	public final double lambda = 0.035; // growth "rate" of the immune-system

	/*
	 * c(r = 0.6, rm= 0.59, mu = 0.01, gamma = 0.02, gammam = 0.02, delta = 0.01, kappa = 0.0175,
	 *  theta1 = 0.3, theta2=0.6, c=i, alpha = 0.035,
	 *  infectivityThreshold = 50, symptomaticThreshold = 100,
	 *   treatmentLength=7, treatmentDelay = 1)
	 */

	public final static double tau = 0.005; // step-size of the tau-leaping algorithm

	public final int transmittedLoad = 7; // pathogen load transmitted to an infected target
	public final int basalImmuneSystemLevel = 7;

	public static final int infectiosityThreshold = 50; // threshold for the pathogen load above which patients can infect other people
	public static final int symptomsThreshold = 100; // threshold for the pathogen load above which patients suffer from symptoms
	public static final double treatmentStartDelay = 1; // how many days after symptoms have set in should a patient wait to receive treatment?
	public static final double treatmentLength = 7; // how many days does the standard treatment last?

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

		// TO DO maybe change way this is done when applying treatment
		this.r = r0;
		this.rm = rm0;
		this.ID = ID;

		randomValueGenerator = new Random();
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];
	}

	public Host(int ID) {
		this.immuneSystemLevel = basalImmuneSystemLevel;
		this.ID = ID;

		this.r = r0;
		this.rm = rm0;

		randomValueGenerator = new Random();
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];
	}

	public Host(int ID, Random randGen) {
		this.immuneSystemLevel = basalImmuneSystemLevel;
		this.ID = ID;

		this.r = r0;
		this.rm = rm0;
		this.randomValueGenerator = randGen;
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];
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

		this.r = this.treated? r0/2*(1-Math.tanh(15*(dose-theta1))) : r0; // growth-rate for wt-strain with possible treatment
		this.rm = this.treated? rm0/2*(1-Math.tanh(15*(dose-theta2))) : rm0; // growth-rate for mutant-strain with possible treatment
		if(backwardMutation)
			nbEventsInOneStep = new double[8];
		else 
			nbEventsInOneStep = new double[7];
	}

	public void setDose(double d){
		this.dose = d;

		this.r = this.treated? r0/2*(1-Math.tanh(15*(dose-theta1))) : r0; // growth-rate for wt-strain with possible treatment
		this.rm = this.treated? rm0/2*(1-Math.tanh(15*(dose-theta2))) : rm0; // growth-rate for mutant-strain with possible treatment
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

		this.r = this.treated? r0/2*(1-Math.tanh(15*(dose-theta1))) : r0; // growth-rate for wt-strain with possible treatment
		this.rm = this.treated? rm0/2*(1-Math.tanh(15*(dose-theta2))) : rm0; // growth-rate for mutant-strain with possible treatment
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

		wtPathogenLevel += (nbEventsInOneStep[0] - nbEventsInOneStep[2] - nbEventsInOneStep[6]);
		mutantPathogenLevel += (nbEventsInOneStep[1] + nbEventsInOneStep[6] - nbEventsInOneStep[3]);
		immuneSystemLevel += (nbEventsInOneStep[4] - nbEventsInOneStep[5]);

		if(backwardMutation) {
			mutantPathogenLevel -= nbEventsInOneStep[7];
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


	/**
	 * TODO delete
	 * Alternative to take into account if an individual is symptomatic or not
	 * should replace entirely simulateOneTimeStep (reuse its name) once it works
	 * Simulate one time step of dynamics for the host.
	 * return 0 if no change of state occurred.
	 * return -3 if the host recovered
	 * return -2 if the host went from infectious to non infectious
	 * return -1 if the host went from symptomatic to infectious but not symptomatic
	 * return 0 if no change in state
	 * return 1 if host went from non-infectious to infectious
	 * return 2 if host went from infectious to symptomatic
	 */
	public int simulateOneTimeStepWithSymptoms(boolean infectious, boolean symptomatic){
		// update the host's age
		age += tau;

		this.calculateNumberOfEventsInOneStep();

		// to do if a level gets below zero, set it to zero

		wtPathogenLevel += (nbEventsInOneStep[0] - nbEventsInOneStep[2] - nbEventsInOneStep[6]);
		mutantPathogenLevel += (nbEventsInOneStep[1] + nbEventsInOneStep[6] - nbEventsInOneStep[3]);
		immuneSystemLevel += (nbEventsInOneStep[4] - nbEventsInOneStep[5]);

		if(backwardMutation) {
			mutantPathogenLevel -= nbEventsInOneStep[7];
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
			return -3;
		else if ( !infectious && (wtPathogenLevel+mutantPathogenLevel)>= Host.infectiosityThreshold)
			return 1;
		else if (infectious && (wtPathogenLevel+mutantPathogenLevel)< Host.infectiosityThreshold)
			return 1;
		else
			return 0;

	}


	// to do add treatment
	public int simulateAndGetEndTime (double[][] dyn, double maxTime) {

		int timeStep = 0;
		int nbTimeSteps = (int) Math.floor(maxTime/tau);

		dyn[0] = new double[]{0,wtPathogenLevel, mutantPathogenLevel, immuneSystemLevel}; 
		timeStep ++;

		while(timeStep < nbTimeSteps){
			this.calculateNumberOfEventsInOneStep();

			// update each state variable according to the events that happened during the last step of size tau
			wtPathogenLevel += (nbEventsInOneStep[0] - nbEventsInOneStep[2] - nbEventsInOneStep[6]);
			mutantPathogenLevel += (nbEventsInOneStep[1] + nbEventsInOneStep[6] - nbEventsInOneStep[3]);
			immuneSystemLevel += (nbEventsInOneStep[4] - nbEventsInOneStep[5]);

			if(backwardMutation) {
				mutantPathogenLevel -= nbEventsInOneStep[7];
				wtPathogenLevel += nbEventsInOneStep[7];
			}

			if(wtPathogenLevel <0 || mutantPathogenLevel <0 || (mutantPathogenLevel+wtPathogenLevel)==0 || immuneSystemLevel<0)
				return timeStep;

			dyn[timeStep] = new double[]{timeStep*tau,wtPathogenLevel, mutantPathogenLevel, immuneSystemLevel}; 
			timeStep ++;
		}
		return timeStep;
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
	
	public double simulateAndGetTimeSpentSymptomatic(){
		double timeSpentInfectious = 0;
		while(this.simulateOneTimeStep(true) >= 0){
			if((wtPathogenLevel+mutantPathogenLevel)>= Host.symptomsThreshold)
				timeSpentInfectious += tau;
		}
		return timeSpentInfectious;
	}
	
	public double simulateAndGetFirstTimeSymptomatic(){
		double timeBeforeSymptoms = 0;
		while(this.simulateOneTimeStep(true) >= 0){
			timeBeforeSymptoms += tau;
			if((wtPathogenLevel+mutantPathogenLevel)>= Host.symptomsThreshold) {
				return timeBeforeSymptoms;
			}
		}
		return -1;
	}
	
	public double simulateAndGetFirstTimeInfectious(){
		double timeBeforeInfectious = 0;
		while(this.simulateOneTimeStep(true) >= 0){
			timeBeforeInfectious += tau;
			if((wtPathogenLevel+mutantPathogenLevel)>= Host.infectiosityThreshold) {
				return timeBeforeInfectious;
			}
		}
		return -1;
	}
	
	/**
	 * Calculates the within-host changes in one step
	 */
	public void calculateNumberOfEventsInOneStep(){

		nbEventsInOneStep[0] = getPoisson(tau * r * (1-mu)*wtPathogenLevel); // +1 to wtPathogenLevel
		nbEventsInOneStep[1] = getPoisson(tau * rm * mutantPathogenLevel); // +1 to mutantPathogenLevel
		nbEventsInOneStep[2] = getPoisson(tau * (eta + kappa* immuneSystemLevel)* wtPathogenLevel); // -1 to wtPathogenLevel
		nbEventsInOneStep[3] = getPoisson(tau * (eta + kappa* immuneSystemLevel)* mutantPathogenLevel); // -1 to mutantPathogenLevel
		nbEventsInOneStep[4] = getPoisson(tau * lambda * (mutantPathogenLevel + wtPathogenLevel)); // +1 to immuneSystemLevel
		nbEventsInOneStep[5] = getPoisson(tau * delta * immuneSystemLevel); // -1 to immuneSystemLevel
		nbEventsInOneStep[6] = getPoisson(tau * mu * wtPathogenLevel); // -1 to wtPathogenLevel and +1 to mutantPathogenLevel

		if(backwardMutation)
			nbEventsInOneStep[7] = getPoisson(tau * mu * mutantPathogenLevel); // +1 to wtPathogenLevel and -1 to mutantPathogenLevel
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

	public static double getMaxInCol(double[][] dyn, int n){
		double result = 0;

		for (int i =0; i<dyn.length; i++) {
			if(dyn[i][n] > result)
				result = dyn[i][n];
		}

		return result;
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

	public void getInfectedByMixedLoad(Host John){

		int wildTypeLoad = (int) (Math.round(John.wtPathogenLevel/(John.wtPathogenLevel+John.mutantPathogenLevel) *  transmittedLoad));
		int mutantLoad = transmittedLoad - wildTypeLoad;

		wtPathogenLevel += wildTypeLoad;
		mutantPathogenLevel += mutantLoad;

	}

	public void getInfectedByFullLoad(Host John){
		if ( randomValueGenerator.nextDouble()> John.wtPathogenLevel/(John.wtPathogenLevel + John.mutantPathogenLevel))
			mutantPathogenLevel += transmittedLoad; //Infected host gets the full pathogen load from the resistant strain
		else
			wtPathogenLevel += transmittedLoad;
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

	public double getInfectedByFullLoadAndGetResistantPortion(Host John){
		if ( randomValueGenerator.nextDouble()> John.wtPathogenLevel/(John.wtPathogenLevel + John.mutantPathogenLevel)) {
			mutantPathogenLevel += transmittedLoad; //Infected host gets the full pathogen load from the resistant strain
			return 1;
		}
		else {
			wtPathogenLevel += transmittedLoad;
			return 0;
		}

	}

	public double getInfectedByMixedLoadAndGetResistantPortion(Host John) throws Exception{
		int wildTypeLoad = (int) (Math.round(John.wtPathogenLevel/(John.wtPathogenLevel+John.mutantPathogenLevel) *  transmittedLoad));
		int mutantLoad = transmittedLoad - wildTypeLoad;

		wtPathogenLevel += wildTypeLoad;
		mutantPathogenLevel += mutantLoad;

		return mutantLoad*1.0/transmittedLoad;
	}

	public void recovers(){
		this.recovered = true;
	}

	public void checkRecovery(){
		// TO DO implement
	}

	public boolean equals(Host e){	
		return (this.ID == e.ID);
	}

	public int compareTo(Host e2){
		if(e2.ID > this.ID) return -1;
		if(e2.ID < this.ID) return 1;

		return 0;
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
		
		//		Host michael  = new Host(0, 7,0,7);
		//		Host john  = new Host(1, 90,55,7);
		//		
		//		michael.getInfectedByFullLoad(john);
		//		
		//		System.out.println(michael);

		Host michael  = new Host(1, 0,7,7); // wt pathogen/mutant pathogen/immune-syst

		//		Host michael  = new Host(0, 7,0,7, 0.4);

		//		int count = 0;
		//		while(michael.simulateOneTimeStep(false) > -1){
		//			count++;
		//			if (count % 100 == 0) {
		//				System.out.println("count: " + count + "\n" + michael);
		//			}
		//			
		//		}
		//		
		//		System.out.println(count);


//		double maxTime = 50;
//
//		int nbTimeSteps = (int) Math.floor(maxTime/michael.tau);
//
//		double[][] dynamicsSeries = new double[nbTimeSteps][4]; 
//
//		int nbSteps = michael.simulateAndGetEndTime(dynamicsSeries, maxTime);
//
//		for(int i = 0; i<nbSteps; i++) {
//			System.out.println(Arrays.toString(dynamicsSeries[i]));
//		}
//
//		System.out.println(nbSteps);
//
//		System.out.println(Host.getMaxInCol(dynamicsSeries, 1));
	}
}
