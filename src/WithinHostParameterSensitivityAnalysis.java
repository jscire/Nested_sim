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

public class WithinHostParameterSensitivityAnalysis {


    public static void assignParameterValue(String paramName, double paramValue) {

        switch(paramName){
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
            case "tau":
                Host.tau = paramValue;
                break;
            default: break;
        }
    }


    public static void main(String[] args) throws Exception {

        DecimalFormatSymbols custom=new DecimalFormatSymbols();
        custom.setDecimalSeparator('.');

        DecimalFormat df = new DecimalFormat("0.00E0");
        df.setDecimalFormatSymbols(custom);
        DecimalFormat df2 = new DecimalFormat("0.00");
        df2.setDecimalFormatSymbols(custom);
        DecimalFormat df3 = new DecimalFormat("00000");
        df3.setDecimalFormatSymbols(custom);
        DecimalFormat df4 = new DecimalFormat("0.0000");
        df4.setDecimalFormatSymbols(custom);

        int numberOfRunsPerDose = 10000;
        int numberOfDosePerJob = 10000;
        int numberOfDosesExplored = 51;

        Random randomGen = new Random();

        String measureToTake;
        String strainOfInterest;

        if(args.length >= 1) {
            measureToTake = args[0];
            if(!measureToTake.equals("RecoveryRateWT") && !measureToTake.equals("RecoveryRateMS") && !measureToTake.contains("ProbabilityEmergence")) {
                throw new IllegalArgumentException("The measure to take is either RecoveryRate or ProbabilityEmergence");
            }
        }
        else {
            measureToTake = "ProbabilityEmergence";
        }

        if(!measureToTake.equals("ProbabilityEmergence"))
            strainOfInterest = measureToTake.substring(measureToTake.length()-2);
        else
            strainOfInterest = "WT";

        String nameOfParameterToVary;
        double parameterValue;

        if(args.length >= 3) {
            nameOfParameterToVary = args[1];
            parameterValue = Double.parseDouble(args[2]);

            if((measureToTake.equals("RecoveryRateWT") || measureToTake.equals("RecoveryRateMS")) && nameOfParameterToVary.equals("mu"))
                throw new IllegalArgumentException("Cannot vary mutation rate and record recovery rate. Mu is set to 0 when recording the recovery rate.");
        }
        else {
            nameOfParameterToVary = "kappa";
            parameterValue = 0.02;
        }

        boolean isInfiniteTreatment;
        if(args.length==4) {
            isInfiniteTreatment = Boolean.parseBoolean(args[3]);
        }
        else {
            isInfiniteTreatment = false;
        }

        //Treatment doses
        double[] doses = new double[numberOfDosesExplored];
        for (int z=0; z<numberOfDosesExplored; z++){
            doses[z] = z*1.0/(numberOfDosesExplored -1);
        }

        //Output file
        File outputFile = new File(measureToTake + "_" + nameOfParameterToVary + ".txt");
        Path path = outputFile.toPath();
        try {
            Files.createFile(path);
            FileOutputStream fos = new FileOutputStream(outputFile);
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
            if (args.length==4)
                bw.write(measureToTake + ",dose," + nameOfParameterToVary + ",NumberOfRuns,isInfTreatment");
            else
                bw.write(measureToTake + ",dose," + nameOfParameterToVary + ",NumberOfRuns");
            bw.newLine();
            bw.close();
        } catch (IOException e) {
            //do nothing
        }


        //Parametrisation of Host
        assignParameterValue(nameOfParameterToVary, parameterValue);

        if(measureToTake.equals("RecoveryRateWT") || measureToTake.equals("RecoveryRateMS"))
            Host.mu = 0;
        if(isInfiniteTreatment)
            Host.treatmentLength = Double.POSITIVE_INFINITY;

        int WTload = strainOfInterest.equals("WT")? Host.transmittedLoad: 0; // WT stands for "wild-type strain"
        int MSload = strainOfInterest.equals("MS")? Host.transmittedLoad: 0; // MS stands for "mutant strain"


        for (int i = 0; i < numberOfDosePerJob; i++) {

            String reportedMeasure;
            double doseApplied = doses[randomGen.nextInt(doses.length)];

            if(measureToTake.contains("ProbabilityEmergence")) {
                int nbOfEmerged = 0;
                for (int j = 0; j < numberOfRunsPerDose; j++) {
                    Host john = new Host(0, WTload, MSload, Host.basalImmuneSystemLevel, doseApplied);
                    nbOfEmerged += john.simulateAndGetEmergenceOtherStrain(strainOfInterest.equals("WT"));
                }
                reportedMeasure = df3.format(nbOfEmerged);
            }
            else { // measure recovery rate
                double totalInfectiousLength = 0;
                for (int j = 0; j < numberOfRunsPerDose; j++) {
                    Host john = new Host(0, WTload, MSload, Host.basalImmuneSystemLevel, doseApplied);
                    totalInfectiousLength += john.simulateAndGetTimeSpentInfectious();
                }
                double meanRecoveryRate = numberOfRunsPerDose * 1.0 / totalInfectiousLength;
                if(totalInfectiousLength == 0)
                    reportedMeasure = "Inf";
                else
                    reportedMeasure = df4.format(meanRecoveryRate);
            }

            // Write the latest results to file
            FileChannel channel;
            channel = FileChannel.open(path, EnumSet.of(StandardOpenOption.APPEND));

            try {
                // FileLock lock = channel.lock(); // Get an exclusive lock on the whole file, only works on Windows
                File lockFile = new File("./lockdir_" + measureToTake + "_" + nameOfParameterToVary);
                boolean hasLock = lockFile.mkdir();
                while(!hasLock) {
                    Thread.sleep(100); // retry to acquire the lock every 100ms.
                    hasLock = lockFile.mkdir();
                }
                try {
                    String newData;
                    if (args.length==4)
                        newData = reportedMeasure + "," + df2.format(doseApplied) + "," + df.format(parameterValue) + "," + df3.format(numberOfRunsPerDose) + "," + isInfiniteTreatment + "\n";
                    else
                        newData = reportedMeasure + "," + df2.format(doseApplied) + "," + df.format(parameterValue) + "," + df3.format(numberOfRunsPerDose) +  "\n";

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
        System.exit(0);
    }
}
