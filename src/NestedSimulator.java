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

	public static void assignParameterValue(String paramName, double paramValue) {

		switch(paramName){
			case "beta":
				Population.transmissionRate = paramValue;
			case "mu":
				Host.mu = paramValue;
				break;
			case "kappa":
				Host.kappa = paramValue;
				break;
			case "fc":
				Host.alpha = paramValue;
				break;
			case "symptomThreshold":
				Host.symptomsThreshold = paramValue;
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
		int numberOfDosesExplored=51;
		int populationSize = 1000;

		Random randomGen = new Random();

		String nameOfParameterToVary;
		double parameterValue;

		if(args.length >= 1){
			nameOfParameterToVary = args[0];
			parameterValue = Double.parseDouble(args[1]);
		} else {
			nameOfParameterToVary = "beta";
			parameterValue =  2.5 * 0.1/populationSize; // keeps a reasonable R0 (around 2.50 with default parameter values
		}

		boolean isInfiniteTreatment;
		if(args.length==3) {
			isInfiniteTreatment = Boolean.parseBoolean(args[2]);
		}
		else {
			isInfiniteTreatment = false;
		}

		//Parametrisation of Population and Host
		assignParameterValue(nameOfParameterToVary, parameterValue);

		if(isInfiniteTreatment)
			Host.treatmentLength = Double.POSITIVE_INFINITY;


		//Treatment doses
		double[] doses = new double[numberOfDosesExplored];
		for (int z=0; z<numberOfDosesExplored; z++){
			doses[z] = z*1.0/(numberOfDosesExplored -1);
		}

		//Output file
		File outputFile = new File("NestedSimulation_" + nameOfParameterToVary + ".txt");
		Path path = outputFile.toPath();
		try {
			Files.createFile(path);
			FileOutputStream fos = new FileOutputStream(outputFile);
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
			if (args.length==3)
				bw.write( "NewInfectRes,SuperinfectionRes,NewInfectWT,SuperinfectionWT,InfectiousBurden,Dose," + nameOfParameterToVary + ",NumberOfRuns,isInfTreatment");
			else
				bw.write( "NewInfectRes,SuperinfectionRes,NewInfectWT,SuperinfectionWT,InfectiousBurden,Dose," + nameOfParameterToVary + ",NumberOfRuns");
			bw.newLine();
			bw.close();
		} catch (IOException e) {
			//do nothing
		}


		for (int j = 0; j< numberOfDosesPerJob; j++) {

			double doseApplied = doses[randomGen.nextInt(doses.length)];

			StringBuilder builder = new StringBuilder(numberOfRunsPerDose*10);

			for (int i = 0; i < numberOfRunsPerDose; i++) {

				Population pop = new Population(populationSize, doseApplied);
				pop.infectWildType();
				int[] infectionEvents  = pop.simulatePopulationAndTrackInfectionsAndBurden();

				if(args.length==3) {
					builder.append(df3.format(infectionEvents[0]) + "," + df3.format(infectionEvents[1]) + "," + df3.format(infectionEvents[2]) + "," +
							df3.format(infectionEvents[3]) + "," + df.format(infectionEvents[4]) + "," +df2.format(doseApplied) + "," + df.format(parameterValue) + "," + df3.format(numberOfRunsPerDose) + ","  + isInfiniteTreatment + "\n");
				}
				else {
					builder.append(df3.format(infectionEvents[0]) + "," + df3.format(infectionEvents[1]) + "," + df3.format(infectionEvents[2]) + "," +
							df3.format(infectionEvents[3]) + "," + df.format(infectionEvents[4]) + "," +df2.format(doseApplied) + "," + df.format(parameterValue) + "," + df3.format(numberOfRunsPerDose) + "\n");
				}
			}

			// Write the latest results to file
			FileChannel channel;
			channel = FileChannel.open(path, EnumSet.of(StandardOpenOption.APPEND));

			try {
				// FileLock lock = channel.lock(); // Get an exclusive lock on the whole file, only works on Windows
				File lockFile = new File("./lockdir_NestedSimulator_" + nameOfParameterToVary);
				boolean hasLock = lockFile.mkdir();
				while(!hasLock) {
					Thread.sleep(100); // retry to acquire the lock every 100ms.
					hasLock = lockFile.mkdir();
				}
				try{

					String newData = builder.toString();

					if (newData != null) {
						ByteBuffer buf = ByteBuffer.allocate(100);
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
	}
}
