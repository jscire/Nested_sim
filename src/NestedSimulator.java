import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.EnumSet;
import java.util.Locale;
import java.util.Random;

public class NestedSimulator {

	static int populationSize = 10000;

	public static void assignParameterValue(String paramName, double paramValue) {

		switch(paramName){
			case "beta":
				Population.transmissionRate = paramValue;
				break;
			case "mu":
				Host.mu = paramValue;
				break;
			case "kappa":
				Host.kappa = paramValue;
				break;
			case "epsilon":
				Host.r0 = paramValue;
				break;
			case "fc":
				Host.alpha = paramValue;
				break;
			case "symptomThreshold":
				Host.symptomsThreshold = paramValue;
				break;
			case "infectiosityThreshold":
				Host.infectiosityThreshold = paramValue;
				break;
			case "flatness":
				Host.flatness = paramValue;
				break;
			case "theta1":
				Host.theta1  =paramValue;
				break;
			case "theta2":
				Host.theta2  =paramValue;
				break;
			case "lambda":
				Host.lambda = paramValue;
				break;
			case "delta":
				Host.delta = paramValue;
				break;
			case "eta":
				Host.eta = paramValue;
				break;
			case "S0":
				populationSize = (int) paramValue;
				break;
			case "tau":
				Host.tau = paramValue;
				break;
			default: break;
		}
	}



	public static void main(String[] args) throws Exception {

		DecimalFormatSymbols custom=new DecimalFormatSymbols();
		custom.setDecimalSeparator('.');

		DecimalFormat df = new DecimalFormat("0.00E00");
		df.setDecimalFormatSymbols(custom);
		DecimalFormat df2 = new DecimalFormat("0.00");
		df2.setDecimalFormatSymbols(custom);
		DecimalFormat df3 = new DecimalFormat("00000");
		df3.setDecimalFormatSymbols(custom);


		int numberOfRunsPerDose = 20;
		int numberOfDosesPerJob = 5000;
		int numberOfDosesExplored=31;
		boolean isInfiniteTreatment=true;

		Random randomGen = new Random();

		String nameOfFirstParameterToVary;
		double firstParameterValue;
		String nameOfSecondParameterToVary="none";
		double secondParameterValue=0;

		if(args.length >= 1){
			nameOfFirstParameterToVary = args[0];
			firstParameterValue = Double.parseDouble(args[1]);
		}
		else {
			nameOfFirstParameterToVary = "beta";
			firstParameterValue =  2.5 * 0.1/populationSize; // keeps a reasonable R0 (around 2.50 with default parameter values
		}

		if(args.length >= 4){
			nameOfSecondParameterToVary = args[2];
			secondParameterValue = Double.parseDouble(args[3]);
		}


		boolean isDoseFixed = false;
		double fixedDose = 0;
		if(args.length >=5) {
			isDoseFixed = true;
			fixedDose = Double.parseDouble(args[4]);
		}

		//Parametrisation of Population and Host
		assignParameterValue(nameOfFirstParameterToVary, firstParameterValue);
		assignParameterValue(nameOfSecondParameterToVary, secondParameterValue);

		//Treatment doses
		double[] doses;
		if (!isDoseFixed) {
			doses = new double[numberOfDosesExplored];
			for (int z=0; z<numberOfDosesExplored; z++){
				doses[z] = z*1.0/(numberOfDosesExplored -1);
			}
		} else {
			doses = new double[]{fixedDose};
		}


		//Output file
		File outputFile;
		if (nameOfSecondParameterToVary.equals("none")) {
			outputFile = new File("NestedSimulation_FullTracking_" + nameOfFirstParameterToVary + ".txt");
		}
		else {
			outputFile = new File("NestedSimulation_FullTracking_" + nameOfFirstParameterToVary + "_" + nameOfSecondParameterToVary + ".txt");
		}

		Path path = outputFile.toPath();
		try {
			Files.createFile(path);
			FileOutputStream fos = new FileOutputStream(outputFile);
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
			if (nameOfSecondParameterToVary.equals("none"))
				bw.write( "NewInfectRes,SuperinfectionRes,NewInfectWT,SuperinfectionWT,FirstInfections,InfectiousBurden,SymptomaticBurden,Dose," + nameOfFirstParameterToVary + ",isInfTreatment");
			else
				bw.write( "NewInfectRes,SuperinfectionRes,NewInfectWT,SuperinfectionWT,FirstInfections,InfectiousBurden,SymptomaticBurden,Dose," + nameOfFirstParameterToVary + "," + nameOfSecondParameterToVary +  ",isInfTreatment");
			bw.newLine();
			bw.close();
		} catch (IOException e) {
			//do nothing
		}


		for (int j = 0; j< numberOfDosesPerJob; j++) {

			double doseApplied = doses[randomGen.nextInt(doses.length)];

			StringBuilder builder = new StringBuilder(numberOfRunsPerDose*10);

			for (int i = 0; i < numberOfRunsPerDose; i++) {

				Population pop = new Population(populationSize -1, doseApplied);
				Host john = new Host(populationSize, 0, 0, 7, doseApplied);
				pop.addSusceptibleHost(john);
				pop.infectWildType(john);
				int[] infectionEvents  = pop.simulatePopulationAndTrackInfectionsAndFirstHostAndTwoBurdenTypes(john);

				if(args.length<=3) {
					builder.append(df3.format(infectionEvents[0]) + "," + df3.format(infectionEvents[1]) + "," + df3.format(infectionEvents[2]) + "," +
							df3.format(infectionEvents[3]) + "," + df.format(infectionEvents[4]) + "," + df.format(infectionEvents[5]) + "," + df.format(infectionEvents[6]) + "," +
							df2.format(doseApplied) + "," + df.format(firstParameterValue) + ","  + isInfiniteTreatment + "\n");
				}
				else {
					builder.append(df3.format(infectionEvents[0]) + "," + df3.format(infectionEvents[1]) + "," + df3.format(infectionEvents[2]) + "," +
							df3.format(infectionEvents[3]) + "," + df3.format(infectionEvents[4]) + "," + df.format(infectionEvents[5]) + "," + df.format(infectionEvents[6]) + "," +
							df2.format(doseApplied) + "," + df.format(firstParameterValue) + "," + df.format(secondParameterValue)  + ","  + isInfiniteTreatment + "\n");
				}
			}

			// Write the latest results to file
			FileChannel channel;
			channel = FileChannel.open(path, EnumSet.of(StandardOpenOption.APPEND));

			try {
				// FileLock lock = channel.lock(); // Get an exclusive lock on the whole file, only works on Windows
				File lockFile = new File("./lockdir_NestedSimulator_" + nameOfFirstParameterToVary);
				boolean hasLock = lockFile.mkdir();
				while(!hasLock) {
					Thread.sleep(100); // retry to acquire the lock every 100ms.
					hasLock = lockFile.mkdir();
				}
				try{

					String newData = builder.toString();

					if (newData != null) {
						ByteBuffer buf = ByteBuffer.allocate(10000);
						buf.clear();
						buf.put(newData.getBytes());
						buf.flip();
						while (buf.hasRemaining()) {
							channel.write(buf);
						}
					}
				}
				catch (Exception e) {
					Thread.currentThread().interrupt();
				}
				finally {
//                    lock.release();
					lockFile.delete();
				}
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			finally {
				channel.close();
			}
		}
		System.exit(0);
	}
}
