package midesp;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import midesp.methods.MICalculator;
import midesp.objects.Phenotype;
import midesp.objects.SNP;
import midesp.objects.SigFinderResult;
import midesp.objects.LimitedPriorityQueue;
import midesp.objects.MIResult;
import midesp.methods.SignificanceFinder;

public class Main {

	private static Path tpedFile, tfamFile, outFile, snpListFile, discCovariatesFile;
	
	private static boolean isContinuous = false;
	private static boolean isNoAPC = false;
	private static boolean isNoEpi = false;
	private static boolean isPrintAll = false;
	
	private static double keepPercentage = 1;
	
	private static int threadCount = Runtime.getRuntime().availableProcessors() / 2;
	private static int kNext = 30;
	private static int apcAverageNumber = 5000;
	private static double fdr = 0.005;
	
	public static void main(String[] args) {
		if(args.length < 2){
			printHelp();
			return;
		}
		long startTime = System.nanoTime();
		if(!execMIDESP(args)) {
			System.out.println("There was an error while executing MIDESP");
			return;
		}
		else {
			System.out.println("Total runtime of " + (System.nanoTime() - startTime) / 1_000_000_000 / 60 + " minutes");
		}	
	}
	
	public static boolean execMIDESP(String[] args) {
		double start;
		try {
			parseArgs(args);
		}
		catch(IllegalArgumentException e) {
			System.out.println("Error while parsing parameters");
			System.out.println(e.getMessage());
			return false;
		}
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(threadCount));
		List<String> fileSigSNPIDsList = null;
		List<SNP> fileSigSNPList = null;
		List<SNP> snpList;
		Phenotype pheno;
		System.out.println("Reading data from files");
		try {
			snpList = SNP.readTPed(tpedFile);
			pheno = Phenotype.readTFam(tfamFile, isContinuous);
			if(discCovariatesFile != null) {
				pheno.readDiscCovariateFile(discCovariatesFile);
			}
			if(snpListFile != null) {
				fileSigSNPIDsList = Files.lines(snpListFile).collect(Collectors.toList());
			}
		}
		catch(IOException e) {
			System.out.println("Error while reading files");
			System.out.println(e.getMessage());
			e.printStackTrace();
			return false;
		}
		System.out.println("Read phenotypes for " + pheno.getLength() + " samples");
		System.out.println("Read data of " + snpList.size() + " SNPs");
		if(fileSigSNPIDsList != null) {
			List<String> allSNPIDs = snpList.parallelStream().map(SNP::getID).distinct().collect(Collectors.toList());
			long foundCount = fileSigSNPIDsList.parallelStream().filter(snp -> allSNPIDs.contains(snp)).count();
			System.out.println("Using " + foundCount + " SNPs from the given list as important instead of using the SNPs that are significant according to their MI value");
			if(foundCount != fileSigSNPIDsList.size()) {
				System.out.println((fileSigSNPIDsList.size()-foundCount) + " SNPs from the list could not be found in the tped file and will be ignored");
			}
			List<String> tmpList = fileSigSNPIDsList;
			fileSigSNPList = snpList.parallelStream().filter(snp -> tmpList.contains(snp.getID())).collect(Collectors.toList());
		}
		System.out.println("Phenotype = " + (isContinuous ? "Continuous" : "Discrete"));
		if(isContinuous) {
			System.out.println("Number of neighbours for MI = " + kNext);
		}
		System.out.println("FDR = " + fdr);
		System.out.println("Number of samples used for APC = " + apcAverageNumber);
		System.out.println("Number of threads = " + threadCount);
		System.out.println("Calculating single SNP association values");
		start = System.nanoTime();
		List<Double> singleSNPMI = snpList.parallelStream().map(snp ->{
			double mi;
			if(isContinuous) {
				mi = MICalculator.calcMI_ContPheno(pheno, kNext, snp);
			}
			else {
				mi = MICalculator.calcMI_DiscPheno(pheno, snp);
			}
			snp.setMItoPheno(mi);
			return mi;
		}).collect(Collectors.toList());
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		System.out.println("Calculated single SNP association values");
		System.out.println("Calculating pvalues for single SNP association");
		SigFinderResult singleSNPMI_SigFinderResult = SignificanceFinder.findSignificantScores(singleSNPMI, fdr);
		if(singleSNPMI_SigFinderResult == null) {
			System.out.println("Could not calculate pvalues for single SNP association");
			return false;
		}
		for(int i = 0; i < snpList.size(); i++) {
			snpList.get(i).setPValue(singleSNPMI_SigFinderResult.getPValues().get(i));
		}
		System.out.println("Calculated pvalues for single SNP association");
		List<SNP> sigSNPList = singleSNPMI_SigFinderResult.getSignificantIndices().parallelStream().map(idx -> snpList.get(idx)).collect(Collectors.toList());
		System.out.println("Number of significantly associated single SNPs = " + sigSNPList.size());
		System.out.println("Writing significantly associated single SNPs to file");
		PrintWriter sigSNPPW;
		try {
			sigSNPPW = new PrintWriter(Files.newBufferedWriter(Paths.get(outFile.toAbsolutePath().toString() + ".sigSNPs")));
			sigSNPPW.println("SNP Entropy MI PValue");
			if(fileSigSNPList == null) {
				sigSNPList.forEach(snp -> sigSNPPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno() + " " + snp.getPValue()));
			}
			else {
				fileSigSNPList.forEach(snp -> sigSNPPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno() + " " + snp.getPValue()));
			}
			sigSNPPW.close();
		} catch (IOException e) {
			System.out.println("Error while writing SNPs to file");
			e.printStackTrace();
			return false;
		}
		if(isPrintAll) {
			try {
				PrintWriter allSNPsPW = new PrintWriter(Files.newBufferedWriter(Paths.get(outFile.toAbsolutePath().toString() + ".allSNPs")));
				allSNPsPW.println("SNP Entropy MI PValue");
				snpList.forEach(snp -> allSNPsPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno() + " " + snp.getPValue()));
				allSNPsPW.close();
			} catch (IOException e) {
				System.out.println("Error while writing SNPs to file");
				e.printStackTrace();
				return false;
			}
		}
		if(isNoEpi) {
			System.out.println("Stopping without calculating epistatic SNP pairs");
			return true;
		}
		if(isNoAPC) {
			System.out.println("Calculating MI values for SNP pairs");
			long possiblePairCount = (fileSigSNPList == null) ? calcPossiblePairsCount(sigSNPList.size(), snpList.size()) : calcPossiblePairsCount(fileSigSNPList.size(), snpList.size());
			int savePairCount = (int) (possiblePairCount * (keepPercentage / 100.0));
			System.out.println("Saving top " + savePairCount + " pairs (" + keepPercentage + "% from " + possiblePairCount + " calculated pairs)");
			LimitedPriorityQueue<MIResult> topResults = new LimitedPriorityQueue<>(savePairCount);
			List<SNP> tmpList = fileSigSNPList;
			List<SNP> snpWithoutSigList = (fileSigSNPList == null) ? snpList.stream().filter(snp -> !sigSNPList.contains(snp)).collect(Collectors.toList()) : snpList.stream().filter(snp -> !tmpList.contains(snp)).collect(Collectors.toList());
			List<SNP> targetSNPList = (fileSigSNPList == null) ? Stream.concat(sigSNPList.stream(), snpWithoutSigList.stream()).collect(Collectors.toList()) : Stream.concat(tmpList.stream(), snpWithoutSigList.stream()).collect(Collectors.toList());
			if(targetSNPList.size() != snpList.size()) {
				System.out.println("List sizes are not equal");
				return false;
			}
			start = System.nanoTime();
			if(fileSigSNPList == null) {
				for(int i = 0; i < sigSNPList.size(); i++) {
					if(i % 10 == 0) {
						System.out.println("Iteration " + i);
					}
					SNP firstSNP = sigSNPList.get(i);
					List<MIResult> tempResults = targetSNPList.subList(i, targetSNPList.size()).parallelStream().map(secondSNP ->{
						double mi;
						if(isContinuous) {
							mi = MICalculator.calcMI_ContPheno(pheno, kNext, firstSNP, secondSNP);
						}
						else {
							mi = MICalculator.calcMI_DiscPheno(pheno, firstSNP, secondSNP);
						}
						return new MIResult(firstSNP.getID(), secondSNP.getID(), mi);
					}).collect(Collectors.toList());
					topResults.addAll(tempResults);
				}
			}
			else {
				for(int i = 0; i < fileSigSNPList.size(); i++) {
					if(i % 10 == 0) {
						System.out.println("Iteration " + i);
					}
					SNP firstSNP = fileSigSNPList.get(i);
					List<MIResult> tempResults = targetSNPList.subList(i, targetSNPList.size()).parallelStream().map(secondSNP ->{
						double mi;
						if(isContinuous) {
							mi = MICalculator.calcMI_ContPheno(pheno, kNext, firstSNP, secondSNP);
						}
						else {
							mi = MICalculator.calcMI_DiscPheno(pheno, firstSNP, secondSNP);
						}
						return new MIResult(firstSNP.getID(), secondSNP.getID(), mi);
					}).collect(Collectors.toList());
					topResults.addAll(tempResults);
				}
			}
			System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
			System.out.println("Writing results to file");
			PrintWriter outPW;
			try {
				outPW = new PrintWriter(Files.newBufferedWriter(outFile));
				outPW.println("SNP1 + SNP2 MI");
				while(!topResults.isEmpty()) {
					MIResult result = topResults.poll();
					outPW.println(result.toNoAPCString());
				}
				outPW.close();
			}
			catch(IOException e) {
				System.out.println("Error while writing results to file");
				e.printStackTrace();
				return false;
			}
			return true;
		}
		if(singleSNPMI_SigFinderResult.getZeroToLambda1Indices().size() < apcAverageNumber || singleSNPMI_SigFinderResult.getBackgroundIndices().size() < sigSNPList.size()) {
			System.out.println("Not enough SNPs to calculate APC averages. Use a smaller value for -apc or deactivate APC with -noapc.");
			System.out.println("Maximum possible value for current dataset is " + singleSNPMI_SigFinderResult.getZeroToLambda1Indices().size());
			return false;
		}
		System.out.println("Calculating SNP-specific average effects for significantly associated single SNPs");
		start = System.nanoTime();
		sigSNPList.parallelStream().forEach(snp ->{
			List<Integer> randIndices = new ArrayList<>(singleSNPMI_SigFinderResult.getZeroToLambda1Indices());
			Collections.shuffle(randIndices);
			double sum = 0.0;
			for(int i = 0; i < apcAverageNumber; i++) {
				if(isContinuous) {
					sum += MICalculator.calcMI_ContPheno(pheno, kNext, snp, snpList.get(randIndices.get(i)));
				}
				else {
					sum += MICalculator.calcMI_DiscPheno(pheno, snp, snpList.get(randIndices.get(i)));
				}
			}
			snp.setAverageMItoPheno(sum / apcAverageNumber);
		});
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		start = System.nanoTime();
		System.out.println("Calculating single average effect for background SNPs");
		List<SNP> backgroundSNPList = singleSNPMI_SigFinderResult.getBackgroundIndices().parallelStream().map(idx -> snpList.get(idx)).collect(Collectors.toList());
		List<Integer> randSNPIdxList = IntStream.range(0, backgroundSNPList.size()).boxed().collect(Collectors.toList());
		Collections.shuffle(randSNPIdxList);
		double backgroundMeanEffect = randSNPIdxList.subList(0, sigSNPList.size()).parallelStream().mapToDouble(idx ->{
			SNP snp = backgroundSNPList.get(idx);
			List<Integer> randIndices = new ArrayList<>(singleSNPMI_SigFinderResult.getZeroToLambda1Indices());
			Collections.shuffle(randIndices);
			double sum = 0.0;
			for(int i = 0; i < apcAverageNumber; i++) {
				if(isContinuous) {
					sum += MICalculator.calcMI_ContPheno(pheno, kNext, snp, snpList.get(randIndices.get(i)));
				}
				else {
					sum += MICalculator.calcMI_DiscPheno(pheno, snp, snpList.get(randIndices.get(i)));
				}
			}
			return sum / apcAverageNumber;
		}).average().getAsDouble();
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		IntStream.range(0, snpList.size()).parallel().filter(idx -> !singleSNPMI_SigFinderResult.getSignificantIndices().contains(idx)).forEach(idx -> snpList.get(idx).setAverageMItoPheno(backgroundMeanEffect));
		System.out.println("Calculating overall average effect based on all SNPs");
		start = System.nanoTime();
		randSNPIdxList = IntStream.range(0, snpList.size()).boxed().collect(Collectors.toList());
		Collections.shuffle(randSNPIdxList);
		double overallMeanEffect = randSNPIdxList.subList(0, sigSNPList.size()).parallelStream().mapToDouble(idx ->{
			SNP snp = snpList.get(idx);
			List<Integer> randIndices = IntStream.range(0, snpList.size()).boxed().collect(Collectors.toList());
			Collections.shuffle(randIndices);
			double sum = 0.0;
			for(int i = 0; i < apcAverageNumber; i++) {
				if(isContinuous) {
					sum += MICalculator.calcMI_ContPheno(pheno, kNext, snp, snpList.get(randIndices.get(i)));
				}
				else {
					sum += MICalculator.calcMI_DiscPheno(pheno, snp, snpList.get(randIndices.get(i)));
				}
			}
			return sum / apcAverageNumber;
		}).average().getAsDouble();
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		System.out.println("SigSNP-Mean = " + sigSNPList.parallelStream().mapToDouble(snp -> snp.getAverageMItoPheno()).average().getAsDouble());
		System.out.println("Background-Mean = " + backgroundMeanEffect);
		System.out.println("Overall-Mean = " + overallMeanEffect);

		if(fileSigSNPList != null) {
			System.out.println("Calculating SNP-specific average effects for single SNPs given in the list");
			start = System.nanoTime();
			fileSigSNPList.parallelStream().forEach(snp ->{
				List<Integer> randIndices = new ArrayList<>(singleSNPMI_SigFinderResult.getZeroToLambda1Indices());
				Collections.shuffle(randIndices);
				double sum = 0.0;
				for(int i = 0; i < apcAverageNumber; i++) {
					if(isContinuous) {
						sum += MICalculator.calcMI_ContPheno(pheno, kNext, snp, snpList.get(randIndices.get(i)));
					}
					else {
						sum += MICalculator.calcMI_DiscPheno(pheno, snp, snpList.get(randIndices.get(i)));
					}
				}
				snp.setAverageMItoPheno(sum / apcAverageNumber);
			});
			System.out.println("FileSigSNP-Mean = " + fileSigSNPList.parallelStream().mapToDouble(snp -> snp.getAverageMItoPheno()).average().getAsDouble());
			System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		}
		System.out.println("Calculating MI values for SNP pairs");
		long possiblePairCount = (fileSigSNPList == null) ? calcPossiblePairsCount(sigSNPList.size(), snpList.size()) : calcPossiblePairsCount(fileSigSNPList.size(), snpList.size());
		int savePairCount = (int) (possiblePairCount * (keepPercentage / 100.0));
		System.out.println("Saving top " + savePairCount + " pairs (" + keepPercentage + "% from " + possiblePairCount + " calculated pairs)");
		LimitedPriorityQueue<MIResult> topResults = new LimitedPriorityQueue<>(savePairCount);
		List<SNP> tmpList = fileSigSNPList;
		List<SNP> snpWithoutSigList = (fileSigSNPList == null) ? snpList.stream().filter(snp -> !sigSNPList.contains(snp)).collect(Collectors.toList()) : snpList.stream().filter(snp -> !tmpList.contains(snp)).collect(Collectors.toList());
		List<SNP> targetSNPList = (fileSigSNPList == null) ? Stream.concat(sigSNPList.stream(), snpWithoutSigList.stream()).collect(Collectors.toList()) : Stream.concat(tmpList.stream(), snpWithoutSigList.stream()).collect(Collectors.toList());
		if(targetSNPList.size() != snpList.size()) {
			System.out.println("List sizes are not equal");
			return false;
		}
		start = System.nanoTime();
		if(fileSigSNPList == null) {
			for(int i = 0; i < sigSNPList.size(); i++) {
				if(i % 10 == 0) {
					System.out.println("Iteration " + i);
				}
				SNP firstSNP = sigSNPList.get(i);
				List<MIResult> tempResults = targetSNPList.subList(i, targetSNPList.size()).parallelStream().map(secondSNP ->{
					double mi;
					if(isContinuous) {
						mi = MICalculator.calcMI_ContPheno(pheno, kNext, firstSNP, secondSNP);
					}
					else {
						mi = MICalculator.calcMI_DiscPheno(pheno, firstSNP, secondSNP);
					}
					double mi_apc = mi - (firstSNP.getAverageMItoPheno() * secondSNP.getAverageMItoPheno() / overallMeanEffect);
					return new MIResult(firstSNP.getID(), secondSNP.getID(), mi, mi_apc);
				}).collect(Collectors.toList());
				topResults.addAll(tempResults);
			}
		}
		else {
			for(int i = 0; i < fileSigSNPList.size(); i++) {
				if(i % 10 == 0) {
					System.out.println("Iteration " + i);
				}
				SNP firstSNP = fileSigSNPList.get(i);
				List<MIResult> tempResults = targetSNPList.subList(i, targetSNPList.size()).parallelStream().map(secondSNP ->{
					double mi;
					if(isContinuous) {
						mi = MICalculator.calcMI_ContPheno(pheno, kNext, firstSNP, secondSNP);
					}
					else {
						mi = MICalculator.calcMI_DiscPheno(pheno, firstSNP, secondSNP);
					}
					double mi_apc = mi - (firstSNP.getAverageMItoPheno() * secondSNP.getAverageMItoPheno() / overallMeanEffect);
					return new MIResult(firstSNP.getID(), secondSNP.getID(), mi, mi_apc);
				}).collect(Collectors.toList());
				topResults.addAll(tempResults);
			}
		}
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		System.out.println("Writing results to file");
		PrintWriter outPW;
		try {
			outPW = new PrintWriter(Files.newBufferedWriter(outFile));
			outPW.println("SNP1 + SNP2 MI MI_APC");
			while(!topResults.isEmpty()) {
				MIResult result = topResults.poll();
				outPW.println(result);
			}
			outPW.close();
		}
		catch(IOException e) {
			System.out.println("Error while writing results to file");
			e.printStackTrace();
			return false;
		}
		return true;
	}

	private static long calcPossiblePairsCount(int sigSNPCount, int totalSNPCount) {
		long pairCounter = 0;
		for(int i = 1; i <= sigSNPCount; i++) {
			pairCounter += totalSNPCount - (i-1);
		}
		return pairCounter;
	}
	
	private static void parseArgs(String[] args) throws IllegalArgumentException{
		Path tmpOutFile = null;
		for(int i = 0; i < args.length - 2; i++){
			switch(args[i]){
			case "-keep":
				keepPercentage = Double.parseDouble(args[i+1]);
				if(keepPercentage <= 0 || keepPercentage > 100) {
					throw new IllegalArgumentException("Value for -keep needs to be in the intervall (0;100]");
				}
				i++;
				break;
			case "-cont":
				isContinuous = true;
				break;
			case "-apc":
				apcAverageNumber = Integer.parseInt(args[i+1]);
				if(apcAverageNumber < 1) {
					throw new IllegalArgumentException("Value for -apc needs to be greater than 0");
				}
				i++;
				break;
			case "-k":
				kNext = Integer.parseInt(args[i+1]);
				if(kNext < 1) {
					throw new IllegalArgumentException("Value for -k needs to be greater than 0");
				}
				i++;
				break;
			case "-fdr":
				fdr = Double.parseDouble(args[i+1]);
				if(fdr <= 0) {
					throw new IllegalArgumentException("Value for -fdr needs to be positive");
				}
				i++;
				break;
			case "-out":
				tmpOutFile = Paths.get(args[i+1]);
				i++;
				break;
			case "-threads":
				threadCount = Integer.parseInt(args[i+1]);
				i++;
				break;
			case "-list":
				snpListFile = Paths.get(args[i+1]);
				i++;
				break;
			case "-disccovariates":
				discCovariatesFile = Paths.get(args[i+1]);
				i++;
				break;
			case "-noapc":
				isNoAPC = true;
				break;
			case "-noepi":
				isNoEpi = true;
				break;				
			case "-all":
				isPrintAll = true;
				break;
			}	
		}
		tpedFile = Paths.get(args[args.length-2]);
		tfamFile = Paths.get(args[args.length-1]);
		if(tmpOutFile == null) {
			if(isNoAPC) {
				outFile = Paths.get(args[args.length-2] + ".epiNoAPC");
			}
			else {
				outFile = Paths.get(args[args.length-2] + ".epi");
			}
		}
		else {
			outFile = tmpOutFile;
		}
	}
	private static void printHelp(){
		System.out.println("Usage: java -jar MIDESP.jar {Options} tpedFile tfamFile");
		System.out.println("Options:");
		System.out.println("-out\t\tfile\tname of outputfile (default tpedFile.epi)");
		System.out.println("-threads\tnumber\tnumber of threads to use (default = Number_of_Cores / 2)");
		System.out.println("-keep\t\tnumber\tkeep only the top X percentage pairs with highest MI (default = 1)");
		System.out.println("-cont\t\t\tindicate that the phenotype is continuous");
		System.out.println("-k\t\tnumber\tset the value of k for MI estimation for continuous phenotypes (default = 30)");
		System.out.println("-fdr\t\tnumber\tset the value of the false discovery rate for finding significantly associated SNPs (default = 0.005)");
		System.out.println("-apc\t\tnumber\tset the number of samples that should be used to estimate the average effects of the SNPs (default = 5000)");
		System.out.println("-list\t\tfile\tname of file with list of SNP IDs to analyze instead of using the SNPs that are significant according to their MI value");
		System.out.println("-disccovariates\tfile\tname of file that contains discrete covariate variables for the samples as tab-separated list");
		System.out.println("-noapc\t\t\tindicate that the APC should not be applied");
		System.out.println("-noepi\t\t\tindicate that no epistatic SNP pairs should be calculated");
		System.out.println("-all\t\t\twrite an additional file containing the MI values for all SNPs (outputfile.allSNPs)");
	}
	
}