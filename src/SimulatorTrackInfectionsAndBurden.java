import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Random;


public class SimulatorTrackInfectionsAndBurden {

	static DecimalFormat df = new DecimalFormat("0.00E0");
	static DecimalFormat df2 = new DecimalFormat("0.0##");
	static DecimalFormat df3 = new DecimalFormat("00000");
	static DecimalFormat df4 = new DecimalFormat("000000");

	static int numberOfRunsPerStep = 20;
	static int numberOfSteps = 5000;
	static int populationSize = 10000;
	static double beta = 2.7E-5;
	static double[] dose = new double[50];
	static double mutationRate = 0.001;

	static double doseApplied;
	//public static boolean done = false;

	static Random randomGen = new Random();

	public static void main(String[] args) throws Exception {
		
		int betaChosen = 5;
		double[] betas = new double[]{1.4E-5, 1.5E-5, 1.6E-5, 1.7E-5, 1.8E-5, 1.9E-5, 2E-5, 2.1E-5, 2.3E-5, 2.5E-5, 2.7E-5, 3.2E-5, 3.7E-5, 4.1E-5, 4.6E-5, 5E-5, 5.5E-5, 6E-5};
//		double[] fitnessCosts = new double[]{0, 0.01, 0.02, 0.03};
		double fitnessCost = 0.01;
		double[] mutationRates = new double[]{0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001};
		int mutationRateChosen = 2;
		
		if(args.length == 1){
			betaChosen = Integer.parseInt(args[0]);
			beta = betas[betaChosen];
		} else if(args.length == 2) {
			betaChosen = Integer.parseInt(args[0]);
			mutationRateChosen = Integer.parseInt(args[1]);
			beta = betas[betaChosen];
			mutationRate = mutationRates[mutationRateChosen];
		}
		
		Host.fitnessCost = fitnessCost;
		Host.mu = mutationRate;

		for (int z=0; z<47; z++){
			dose[z] = z*0.7/46;
		}
		dose[47] = 0.8;
		dose[48] = 0.9;
		dose[49] = 1;
		//		System.out.println(Arrays.toString(dose));

		for (int j = 0; j< numberOfSteps; j++) {
			
			doseApplied = dose[randomGen.nextInt(dose.length)]; //draw randomly the dose for this thread

			StringBuilder builder = new StringBuilder(numberOfRunsPerStep);

			for (int i = 0; i < numberOfRunsPerStep; i++) {

				Population pop = new Population(populationSize, 0, beta, doseApplied);
				pop.infectWildType();
				int[] infectionEvents  = pop.simulatePopulationAndTrackInfectionsAndBurden();
//				System.out.println(df.format(burden));
				
				builder.append(df4.format(infectionEvents[4]) + "\t" + df3.format(infectionEvents[0]) + "\t" + df3.format(infectionEvents[1]) + "\t" + 
						df3.format(infectionEvents[2]) + "\t" + df3.format(infectionEvents[3]) + "\t" +df2.format(doseApplied) + "\n");
			}

			File outputFile = new File("Simulations_trackInfectionsAndBurden_R0_" + betaChosen + "_mu_" + mutationRateChosen + ".txt");
			Path path = outputFile.toPath();
			try {
				Files.createFile(path);
			} catch (IOException e) {
				//do nothing
				//					System.out.println("File already exists");
			}

			FileChannel channel;
			channel = FileChannel.open(path,  EnumSet.of(StandardOpenOption.APPEND));

			try {
				// Get an exclusive lock on the whole file
				FileLock lock = channel.lock();
				try{

					String newData = builder.toString();

					if(newData != null) {

						ByteBuffer buf = ByteBuffer.allocate(1024);
						buf.clear();
						buf.put(newData.getBytes());
						buf.flip();

						while(buf.hasRemaining()) {
							channel.write(buf);
						}
						//														System.out.println("Something was added to the buffered writer.");
					} 

				} catch (Exception e) {
					// TODO Auto-generated catch block
					//						Thread.currentThread().interrupt();
				} finally {
					lock.release();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				//					e.printStackTrace();
			} finally {
				channel.close();
			}

		}

	}
}
