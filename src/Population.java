import java.util.HashSet;
import java.util.Random;


public class Population {

	public static final boolean infectionByMixedPathogenLoad = false;

	public static final int maxTime = 300;

	public int populationSize;

	public int susceptibleCount;
	public int nonInfectiousCount;
	public int infectiousCount;
	public int recoveredCount;
	public int symptomaticCount;

	public static double transmissionRate = 2.5E-5;
	public final int timeStepRatio  = 1; // represents how many within-host time steps of length tau are performed during one population time step

	public HashSet<Host> susceptiblePopulation;
	public HashSet<Host> nonInfectiousPopulation;
	public HashSet<Host> infectiousPopulation;
	public HashSet<Host> recoveredPopulation;

	public HashSet<Host> bufferFromSusceptibleToNonInfectious = new HashSet<Host>(); // the "buffers" are used in simulation step to move hosts from one state to the other while making sure we don't simulate the same host twice in a step
	public HashSet<Host> bufferFromNonInfectiousToInfectious = new HashSet<Host>();
	public HashSet<Host> bufferFromInfectiousToNonInfectious = new HashSet<Host>();
	public HashSet<Host> bufferFromNonInfectiousToRecovered = new HashSet<Host>();


	public Random randomValueGenerator = new Random();

	public Population(int popSize, double dose) {
		populationSize = popSize;

		susceptiblePopulation = new HashSet<Host>(2 * popSize); //values of initialCapacity applied here are pretty much random
		nonInfectiousPopulation = new HashSet<Host>((int) Math.floor(0.2 * popSize));
		infectiousPopulation = new HashSet<Host>((int) Math.floor(0.2 * popSize));
		recoveredPopulation = new HashSet<Host>(popSize);

		for (int i = 0; i < popSize; i++) {
			this.susceptiblePopulation.add(new Host(i, dose, randomValueGenerator));
		}
	}

	public Population(int popSize, double dose, double beta) {
		populationSize = popSize;

		susceptiblePopulation = new HashSet<Host>(2*popSize); //values of initialCapacity applied here are pretty much random
		nonInfectiousPopulation = new HashSet<Host>((int) Math.floor(0.2*popSize));
		infectiousPopulation = new HashSet<Host>((int) Math.floor(0.2*popSize));
		recoveredPopulation = new HashSet<Host>(popSize);

		for (int i =0; i<popSize; i++) {
			this.susceptiblePopulation.add(new Host(i, dose, randomValueGenerator));
		}

		this.transmissionRate = beta;
	}

	public void setTransmissionRate(double beta){
		this.transmissionRate = beta;
	}

	public void addSusceptibleHost(Host John) throws Exception{
		populationSize ++;

		if((John.mutantPathogenLevel + John.wtPathogenLevel) > 0)
			throw new Exception("Host is already carrying the pathogen, it should be susceptible.");

		this.susceptiblePopulation.add(John);
		susceptibleCount = susceptiblePopulation.size();
	}

	/**
	 * Represents an initial infection with the wild-type strain
	 * @throws Exception
	 */
	public void infectWildType() throws Exception{

		for (Host h : susceptiblePopulation) { //the iterator here is just a trick to get one host infected (we don't care which), remove it from susceptibles and add it to non-infectious
			h.getInfectedByWildType();
			nonInfectiousPopulation.add(h);
			susceptiblePopulation.remove(h);

			infectiousCount = infectiousPopulation.size();
			nonInfectiousCount = nonInfectiousPopulation.size();
			recoveredCount = recoveredPopulation.size();
			// to do if needed while debugging, check that susceptibleCount == susceptiblePopulation.size()
			susceptibleCount = populationSize - (infectiousCount + nonInfectiousCount + recoveredCount);

			break; //the loop executes only once
		}
	}

	public void infectWildType(Host michael) throws Exception{

		for (Host h : susceptiblePopulation) {

			if(h.equals(michael)) {
				h.getInfectedByWildType();
				nonInfectiousPopulation.add(h);
				susceptiblePopulation.remove(h);

				infectiousCount = infectiousPopulation.size();
				nonInfectiousCount = nonInfectiousPopulation.size();
				recoveredCount = recoveredPopulation.size();
				// to do if needed while debugging, check that susceptibleCount == susceptiblePopulation.size()
				susceptibleCount = populationSize - (infectiousCount + nonInfectiousCount + recoveredCount);

				break;
			}
		}
	}

	/**
	 * Based on Donald E. Knuth (1969). Seminumerical Algorithms. The Art of Computer Programming, Volume 2. Addison Wesley.
	 * Adapted from http://stackoverflow.com/questions/9832919/generate-poisson-arrival-in-java#9832977
	 * @param mean
	 * @return
	 */
	public int getPoisson(double mean) {
		double L = Math.exp(-mean);
		int k = 0;
		double p = 1.0;
		do {
			p = p * randomValueGenerator.nextDouble();
			k++;
		} while (p > L);
		return k - 1;
	}

	public int[] simulateOnePopulationTimeStepAndTrackInfections() throws Exception{

		for (int i =0 ; i < timeStepRatio ; i++)
			simulateAllWithinHostDynamics();

		int[] res = simulateInfectionEventsInOneStepAndTrackInfections();

		if ((nonInfectiousCount + infectiousCount) < 1)
			return null;

		return res;
	}

	public int[] simulateOnePopulationTimeStepTrackInfectionsAndHost(Host michael) throws Exception{

		for (int i =0 ; i < timeStepRatio ; i++)
			simulateAllWithinHostDynamics();

		int[] res = simulateInfectionEventsInOneStepTrackInfectionsAndHost(michael);

		if ((nonInfectiousCount + infectiousCount) < 1)
			return null;

		return res;
	}

	public void simulateAllWithinHostDynamics() throws Exception{

		// initialize number of symptomatic hosts
		symptomaticCount = 0;

		for(Host h : nonInfectiousPopulation) { // simulate the non-infectious hosts' dynamics for one step

			int simulatedStepResult = h.simulateOneTimeStep(false);
			if (simulatedStepResult == 0) {
				continue; // Host remained non-infectious
			} else {
				if( simulatedStepResult == 1)
					bufferFromNonInfectiousToInfectious.add(h);
				else {
					bufferFromNonInfectiousToRecovered.add(h);
				}
			}
		}

		for (Host h: infectiousPopulation) { // simulate the infectious hosts' dynamics for one step

			int simulatedStepResult = h.simulateOneTimeStep(true);
			if (simulatedStepResult == 0){ // Host remained infectious
				if (h.isSymptomatic())
					symptomaticCount ++;
				continue; }
			else { // Host became non-infectious
				bufferFromInfectiousToNonInfectious.add(h);
			}
		}

		nonInfectiousPopulation.removeAll(bufferFromNonInfectiousToRecovered);
		recoveredPopulation.addAll(bufferFromNonInfectiousToRecovered);
		bufferFromNonInfectiousToRecovered.clear();

		nonInfectiousPopulation.removeAll(bufferFromNonInfectiousToInfectious); // remove the newly infected hosts from the set of non-infectious
		infectiousPopulation.addAll(bufferFromNonInfectiousToInfectious); // add the newly infectious hosts to the set of infectious hosts
		bufferFromNonInfectiousToInfectious.clear(); // clear buffer

		infectiousPopulation.removeAll(bufferFromInfectiousToNonInfectious);
		nonInfectiousPopulation.addAll(bufferFromInfectiousToNonInfectious); // add the newly non-infectious hosts to the set of non-infectious hosts
		bufferFromInfectiousToNonInfectious.clear(); // clear buffer

		infectiousCount = infectiousPopulation.size();
		nonInfectiousCount = nonInfectiousPopulation.size();
		recoveredCount = recoveredPopulation.size();
		susceptibleCount = susceptiblePopulation.size();

		int total = susceptibleCount + infectiousCount + nonInfectiousCount + recoveredCount;

		// to do remove, just for debugging
		if(total != populationSize) {
			throw new Exception("Problem: total amount of individuals is not equal to population size");
		}
	}

	public int[] simulateInfectionEventsInOneStepAndTrackInfections() throws Exception{

		if(infectionByMixedPathogenLoad) {
			throw new Exception("Tracking infections is not compatible with the mixed infections model");
		}

		int nbInfectionEventsInStep = getPoisson(timeStepRatio * Host.tau * transmissionRate * infectiousCount * (populationSize-recoveredCount-1)); // the -1 is here to not allow a host to infect themself

		if(nbInfectionEventsInStep > (populationSize-recoveredCount-1) && nbInfectionEventsInStep>1) {
			throw new Exception("Too many infection events happened at once, maybe reduce timeStepRatio or Host.tau.");
		}

		int[] infections = new int[]{0,0,0,0}; // new infect by res / superinfection by res / new infect by wt / superinfection infect by wt

		for (int i = 0; i<nbInfectionEventsInStep; i++) {

			int randomInfectedHost = randomValueGenerator.nextInt(populationSize-recoveredCount-1); // the -1 is here to not allow a host to infect themself
			int randomInfectingHost = randomValueGenerator.nextInt(infectiousCount);

			Host infectingHost = new Host(0); // find the infecting host among the infectious population
			int countInfecting = 0;
			for (Host h1 : infectiousPopulation) {
				if(countInfecting == randomInfectingHost) {
					infectingHost = h1;
					break; // break the loop once infecting host has been found
				}
				countInfecting ++;
			}

			Host infectedHost = new Host(0); // find the infected host among the non-recovered population 

			int countInfected = 0;
			if(randomInfectedHost < susceptibleCount) { // the newly infected host was previously a susceptible host
				for (Host h2 : susceptiblePopulation) { // iterate n-times (with n randomly drawn) in the susceptible pop to find the newly infected host
					if (countInfected == randomInfectedHost) {
						infectedHost = h2;
						break;
					}
					countInfected ++;
				}

				if(infectedHost.getInfectedByFullLoadAndIsInfectedByRes(infectingHost))
					infections[0] ++;
				else infections[2] ++;

				susceptiblePopulation.remove(infectedHost);
				nonInfectiousPopulation.add(infectedHost); //newly infected host is added to non-infectious population

			} else {
				randomInfectedHost -= susceptibleCount;
				if(randomInfectedHost < nonInfectiousCount) { // the newly infected host was previously a non-infectious host
					for (Host h2 : nonInfectiousPopulation) { // iterate n-times (with n randomly drawn) in the non-infectious pop to find the newly infected host
						if (countInfected == randomInfectedHost) {
							infectedHost = h2;
							break;
						}
						countInfected ++; // note: we don't care to check if the infected host is now infectious, it will be done in the next simulation step 
					}
				} else { // the newly infected host was previously an infectious host
					randomInfectedHost -= nonInfectiousCount;
					for (Host h2 : infectiousPopulation) { // iterate n-times (with n randomly drawn) in the infectious pop to find the newly infected host
						if (h2.equals(infectingHost)) continue; // so that infectedHost and infectingHost cannot be the same

						if (countInfected == randomInfectedHost) {
							infectedHost = h2;
							break;
						}
						countInfected ++;
					}
				}
				if(infectedHost.getInfectedByFullLoadAndIsInfectedByRes(infectingHost))
					infections[1] ++;
				else infections[3] ++;
			}

			//TO DO, think of doing that in a count instead going to the size every time
			infectiousCount = infectiousPopulation.size();
			nonInfectiousCount = nonInfectiousPopulation.size();
			recoveredCount = recoveredPopulation.size();
			susceptibleCount = susceptiblePopulation.size();
		}

		return infections;
	}

	public int[] simulateInfectionEventsInOneStepTrackInfectionsAndHost(Host michael) throws Exception{

		if(infectionByMixedPathogenLoad) {
			throw new Exception("Tracking infections is not compatible with the mixed infections model");
		}

		int nbInfectionEventsInStep = getPoisson(timeStepRatio * Host.tau * transmissionRate * infectiousCount * (populationSize-recoveredCount-1)); // the -1 is here to not allow a host to infect themself

		if(nbInfectionEventsInStep > (populationSize-recoveredCount-1) && nbInfectionEventsInStep>1) {
			throw new Exception("Too many infection events happened at once, maybe reduce timeStepRatio or Host.tau.");
		}

		int[] infections = new int[]{0,0,0,0,0}; // new infect by res / superinfection by res / new infect by wt / superinfection infect by wt / infections by tracked Host

		for (int i = 0; i<nbInfectionEventsInStep; i++) {

			int randomInfectedHost = randomValueGenerator.nextInt(populationSize-recoveredCount-1); // the -1 is here to not allow a host to infect themself
			int randomInfectingHost = randomValueGenerator.nextInt(infectiousCount);

			Host infectingHost = new Host(0); // find the infecting host among the infectious population
			int countInfecting = 0;
			for (Host h1 : infectiousPopulation) {
				if(countInfecting == randomInfectingHost) {
					infectingHost = h1;
					break; // break the loop once infecting host has been found
				}
				countInfecting ++;
			}

			if(infectingHost.equals(michael))
				infections[4] ++;

			Host infectedHost = new Host(0); // find the infected host among the non-recovered population

			int countInfected = 0;
			if(randomInfectedHost < susceptibleCount) { // the newly infected host was previously a susceptible host
				for (Host h2 : susceptiblePopulation) { // iterate n-times (with n randomly drawn) in the susceptible pop to find the newly infected host
					if (countInfected == randomInfectedHost) {
						infectedHost = h2;
						break;
					}
					countInfected ++;
				}

				if(infectedHost.getInfectedByFullLoadAndIsInfectedByRes(infectingHost))
					infections[0] ++;
				else infections[2] ++;

				susceptiblePopulation.remove(infectedHost);
				nonInfectiousPopulation.add(infectedHost); //newly infected host is added to non-infectious population

			} else {
				randomInfectedHost -= susceptibleCount;
				if(randomInfectedHost < nonInfectiousCount) { // the newly infected host was previously a non-infectious host
					for (Host h2 : nonInfectiousPopulation) { // iterate n-times (with n randomly drawn) in the non-infectious pop to find the newly infected host
						if (countInfected == randomInfectedHost) {
							infectedHost = h2;
							break;
						}
						countInfected ++; // note: we don't care to check if the infected host is now infectious, it will be done in the next simulation step
					}
				} else { // the newly infected host was previously an infectious host
					randomInfectedHost -= nonInfectiousCount;
					for (Host h2 : infectiousPopulation) { // iterate n-times (with n randomly drawn) in the infectious pop to find the newly infected host
						if (h2.equals(infectingHost)) continue; // so that infectedHost and infectingHost cannot be the same

						if (countInfected == randomInfectedHost) {
							infectedHost = h2;
							break;
						}
						countInfected ++;
					}
				}
				if(infectedHost.getInfectedByFullLoadAndIsInfectedByRes(infectingHost))
					infections[1] ++;
				else infections[3] ++;
			}

			//TO DO, think of doing that in a count instead going to the size every time
			infectiousCount = infectiousPopulation.size();
			nonInfectiousCount = nonInfectiousPopulation.size();
			recoveredCount = recoveredPopulation.size();
			susceptibleCount = susceptiblePopulation.size();
		}

		return infections;
	}

	public int[] simulatePopulationAndTrackInfectionsAndBurden() throws Exception{

		int[] infections = new int[]{0,0,0,0,0}; // new infect by res / superinfection by res / new infect by wt / superinfection infect by wt // burden (floored to the integer)

		double burden = 0;

		int[] stepResult = simulateOnePopulationTimeStepAndTrackInfections();

		while(stepResult != null){
			for (int i= 0; i < 4; i++){
				infections[i] += stepResult[i];
			}
			burden += (infectiousCount * Host.tau * timeStepRatio);
			stepResult = simulateOnePopulationTimeStepAndTrackInfections();
		}

		infections[4] = (int) burden; // this way, no excessive rounding during each time step.

		return infections;
	}

	public int[] simulatePopulationAndTrackInfectionsAndFirstHostAndTwoBurdenTypes(Host michael) throws Exception{

		int[] infections = new int[]{0,0,0,0,0,0,0}; // new infect by res / superinfection by res / new infect by wt / superinfection infect by wt // burden (floored to the integer)

		double infectiousBurden = 0;
		double symptomaticBurden = 0;

		int[] stepResult = simulateOnePopulationTimeStepTrackInfectionsAndHost(michael);

		while(stepResult != null){
			for (int i= 0; i < 5; i++){
				infections[i] += stepResult[i];
			}
			infectiousBurden += (infectiousCount * Host.tau * timeStepRatio);
			symptomaticBurden += (symptomaticCount * Host.tau * timeStepRatio);
			stepResult = simulateOnePopulationTimeStepTrackInfectionsAndHost(michael);
		}

		infections[5] = (int) infectiousBurden;
		infections[6] = (int) symptomaticBurden;

		return infections;
	}



	public static void main(String[] args) throws Exception {

		int populationSize = 100;
		Population pop = new Population(populationSize -1, 0.4, 2.3E-3); // popSize = 100, R0 = 3
		Host john = new Host(populationSize, 0, 0, 7);
		pop.addSusceptibleHost(john);
		pop.infectWildType(john);

		int[] infectionEvents  = pop.simulatePopulationAndTrackInfectionsAndFirstHostAndTwoBurdenTypes(john);
		for (int i = 0; i < infectionEvents.length; i++) {
			System.out.println(i + " " + infectionEvents[i]);
		}
	}
}
