package midesp;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import midesp.methods.MICalculator;
import midesp.methods.MICalculator.SNPBiCalculator;
import midesp.methods.MICalculator.SNPCalculator;
import midesp.objects.Phenotype;
import midesp.objects.SNP;
import midesp.objects.SigFinderResult;
import midesp.objects.TopKHeap_MI;
import midesp.objects.TopKHeap_MIAPC;
import midesp.objects.EntropyCache;
import midesp.objects.GeneralizedBitSet;
import midesp.methods.SignificanceFinder;

public class Main {
	
	private static Path tpedFile, tfamFile, outFile, snpListFile, discCovariatesFile, contCovariatesFile;
	
	private static boolean isContinuous = false;
	private static boolean isNoAPC = false;
	private static boolean isNoEpi = false;
	private static boolean isPrintAll = false;
	
	private static double keepPercentage = 1;
	
	private static int threadCount = Math.max(1, Runtime.getRuntime().availableProcessors() / 2);
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
		Set<String> fileSigSNPIDsSet = null;
		List<SNP> fileSigSNPList = null;
		List<SNP> snpList;
		Phenotype pheno;
		EntropyCache entropyCache;
		if(contCovariatesFile != null) {
			throw new UnsupportedOperationException(
			        "Continuous covariates are not supported in the current implementation."
			);
		}
		System.out.println("Reading data from files");
		try {
			int sampleCount = (int) Files.lines(tfamFile).count();
			entropyCache = new EntropyCache(sampleCount);
			snpList = SNP.readTPed(tpedFile);
			pheno = Phenotype.readTFam(tfamFile, isContinuous);
			if(discCovariatesFile != null) {
				pheno.readDiscCovariateFile(discCovariatesFile);
			}
			if(contCovariatesFile != null) {
				pheno.readContCovariateFile(contCovariatesFile);
			}
			if(snpListFile != null) {
				try(Stream<String> lines = Files.lines(snpListFile)){
					fileSigSNPIDsSet = lines.collect(Collectors.toSet());
				}
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
		if(fileSigSNPIDsSet != null) {
			Set<String> sigIDsSet = fileSigSNPIDsSet;
			long foundCount = snpList.parallelStream().map(SNP::getID).filter(sigIDsSet::contains).count();
			System.out.println("Using " + foundCount + " SNPs from the given list as important instead of using the SNPs that are significant according to their MI value");
			if(foundCount != fileSigSNPIDsSet.size()) {
				System.out.println((fileSigSNPIDsSet.size()-foundCount) + " SNPs from the list could not be found in the tped file and will be ignored");
			}
			fileSigSNPList = snpList.parallelStream().filter(snp -> sigIDsSet.contains(snp.getID())).toList();
		}
		System.out.println("Phenotype = " + (isContinuous ? "Continuous" : "Discrete"));
		if(isContinuous) {
			System.out.println("Number of neighbours for MI = " + kNext);
		}
		System.out.println("FDR = " + fdr);
		System.out.println("Number of samples used for APC = " + apcAverageNumber);
		System.out.println("Number of threads = " + threadCount);
		SNPCalculator singleMICalculator = MICalculator.resolveMICalculator(isContinuous, entropyCache, pheno, kNext);
		SNPBiCalculator pairMICalculator = MICalculator.resolvePairMICalculator(isContinuous, entropyCache, pheno, kNext);
		if(pheno.hasDiscCovariate()) {
			System.out.println("Calculating joint entropy between SNPs and discrete covariate");
			start = System.nanoTime();
			snpList.parallelStream().forEach(snp ->{
				GeneralizedBitSet jointBitSet = GeneralizedBitSet.combineTwo(snp.getBitSet(), pheno.getDiscCovariate());
				snp.setSNPDiscCovariateBitSet(jointBitSet);		
			});
			if(!pheno.isContinuous()) {
				snpList.parallelStream().forEach(snp ->{
					GeneralizedBitSet jointBitSet = GeneralizedBitSet.combineTwo(snp.getBitSet(), pheno.getDiscPhenotype_DiscCovariate());
					snp.setSNPDiscPhenoDiscCovariateBitSet(jointBitSet);			
				});
			}
			System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		}
		System.out.println("Calculating single SNP association values");
		start = System.nanoTime();
		List<Double> singleSNPMI = snpList.parallelStream().map(snp ->{
			double mi = singleMICalculator.compute(snp, snp.getEntropyNats());
			snp.setMItoPheno(mi);
			return mi;
		}).toList();
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		System.out.println("Calculated single SNP association values");
		/**TEMP_TO_DELETE**/
		try(PrintWriter testPW = new PrintWriter(Files.newBufferedWriter(Paths.get(outFile.toAbsolutePath().toString() + ".sigSNPs")))) {
			testPW.println("SNP Entropy MI");
			snpList.forEach(snp -> testPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno()));
		} catch (IOException e) {
			System.out.println("Error while writing SNPs to file");
			e.printStackTrace();
			return false;
		}
		/**TEMP_TO_DELETE**/
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
		List<SNP> sigSNPList = singleSNPMI_SigFinderResult.getSignificantIndices().parallelStream().map(idx -> snpList.get(idx)).toList();
		System.out.println("Number of significantly associated single SNPs = " + sigSNPList.size());
		System.out.println("Writing significantly associated single SNPs to file");
		try(PrintWriter sigSNPPW = new PrintWriter(Files.newBufferedWriter(Paths.get(outFile.toAbsolutePath().toString() + ".sigSNPs")))) {
			sigSNPPW.println("SNP Entropy MI PValue");
			if(fileSigSNPList == null) {
				sigSNPList.forEach(snp -> sigSNPPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno() + " " + snp.getPValue()));
			}
			else {
				fileSigSNPList.forEach(snp -> sigSNPPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno() + " " + snp.getPValue()));
			}
		} catch (IOException e) {
			System.out.println("Error while writing SNPs to file");
			e.printStackTrace();
			return false;
		}
		if(isPrintAll) {
			try(PrintWriter allSNPsPW = new PrintWriter(Files.newBufferedWriter(Paths.get(outFile.toAbsolutePath().toString() + ".allSNPs")))) {
				allSNPsPW.println("SNP Entropy MI PValue");
				snpList.forEach(snp -> allSNPsPW.println(snp.getID() + " " + snp.getEntropyLog2() + " " + snp.getMItoPheno() + " " + snp.getPValue()));
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
			List<SNP> effectiveSigSNPList = (fileSigSNPList != null) ? fileSigSNPList : sigSNPList; 
			long possiblePairCount = calcPossiblePairsCount(effectiveSigSNPList.size(), snpList.size());
			int savePairCount = (int) (possiblePairCount * (keepPercentage / 100.0));
			System.out.println("Saving top " + savePairCount + " pairs (" + keepPercentage + "% from " + possiblePairCount + " calculated pairs)");
			
			Set<SNP> sigSNPSet = new HashSet<>(effectiveSigSNPList);
			List<SNP> snpWithoutSigList = snpList.stream().filter(snp -> !sigSNPSet.contains(snp)).toList();
			List<SNP> targetSNPList = Stream.concat(effectiveSigSNPList.stream(), snpWithoutSigList.stream()).toList();
			if(targetSNPList.size() != snpList.size()) {
				System.out.println("List sizes are not equal");
				return false;
			}
			start = System.nanoTime();
			int chunkSize = (effectiveSigSNPList.size() + threadCount - 1) / threadCount;
			
			List<TopKHeap_MI> partialQueues = IntStream.range(0, threadCount).parallel().mapToObj(worker -> {
				int startIdx = worker * chunkSize;
				int endIdx = Math.min(startIdx + chunkSize, effectiveSigSNPList.size());
				TopKHeap_MI queue = new TopKHeap_MI(savePairCount);
				for (int i = startIdx; i < endIdx; i++) {
					SNP firstSNP = effectiveSigSNPList.get(i);
					for (int j = i; j < targetSNPList.size(); j++) {
						SNP secondSNP = targetSNPList.get(j);
						double mi = pairMICalculator.compute(firstSNP, secondSNP);
						queue.offer(i, j, mi);
					}
				}
				return queue;
			}).collect(Collectors.toCollection(ArrayList::new));
			System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
			System.out.println("Combining results from threads");
			start = System.nanoTime();
			//Fast enough for now but may be improved by parallelization and more advanced merge
			partialQueues.sort(
				    Comparator.comparingDouble(TopKHeap_MI::min).reversed()
				);
			TopKHeap_MI topHeap = partialQueues.get(0);
			for (int i = 1; i < partialQueues.size(); i++) {
				partialQueues.get(i).forEach((snp1, snp2, mi) -> { 
					topHeap.offer(snp1, snp2, mi); 
				});
			}
			System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");

			System.out.println("Writing results to file");
			
			try(PrintWriter outPW = new PrintWriter(Files.newBufferedWriter(outFile))) {
				outPW.println("SNP1 + SNP2 MI");
				topHeap.drain((snp1Index, snp2Index, mi) -> {
				    String snp1 = effectiveSigSNPList.get(snp1Index).getID();
				    String snp2 = targetSNPList.get(snp2Index).getID();
				    outPW.println(snp1 + " + " + snp2 + " " + mi);
				});
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
		List<Integer> candidateIndices = singleSNPMI_SigFinderResult.getZeroToLambda1Indices();
		final int candidateCount = candidateIndices.size();
		SNP[] targetSNPs = new SNP[candidateCount];
		for (int i = 0; i < candidateCount; i++) {
			targetSNPs[i] = snpList.get(candidateIndices.get(i));
		}
		// One permutation buffer per worker thread.
		ThreadLocal<int[]> scratch =
				ThreadLocal.withInitial(() -> {
		            int[] indices = new int[candidateCount];
		            for (int i = 0; i < candidateCount; i++) {
		                indices[i] = i;
		            }
		            return indices;
		        });
		// Stores the swaps so that we can restore the array afterwards.
		ThreadLocal<int[]> swapPositions =
		    ThreadLocal.withInitial(() -> new int[apcAverageNumber]);
		sigSNPList.parallelStream().forEach(snp ->{
			int[] randomIndices = scratch.get();
		    int[] swaps = swapPositions.get();
		    ThreadLocalRandom random = ThreadLocalRandom.current();
		    
		    double sum = 0.0;
		    
			for(int i = 0; i < apcAverageNumber; i++) {
				int j = random.nextInt(i, candidateCount);

		        swaps[i] = j;

		        int tmp = randomIndices[i];
		        randomIndices[i] = randomIndices[j];
		        randomIndices[j] = tmp;

		        sum += pairMICalculator.compute(snp, targetSNPs[randomIndices[i]]);
			}
			 // Restore the permutation so the same scratch array can be reused.
		    for (int i = apcAverageNumber - 1; i >= 0; i--) {
		        int j = swaps[i];

		        int tmp = randomIndices[i];
		        randomIndices[i] = randomIndices[j];
		        randomIndices[j] = tmp;
		    }
			snp.setAverageMItoPheno(sum / apcAverageNumber);
		});
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		start = System.nanoTime();
		System.out.println("Calculating single average effect for background SNPs");
		List<SNP> backgroundSNPList = singleSNPMI_SigFinderResult.getBackgroundIndices().parallelStream().map(idx -> snpList.get(idx)).toList();
		List<Integer> randSNPIdxList = IntStream.range(0, backgroundSNPList.size()).boxed().collect(Collectors.toList());
		Collections.shuffle(randSNPIdxList);
		double backgroundMeanEffect = randSNPIdxList.subList(0, sigSNPList.size()).parallelStream().mapToDouble(idx ->{
			int[] randomIndices = scratch.get();
		    int[] swaps = swapPositions.get();
		    ThreadLocalRandom random = ThreadLocalRandom.current();
		    
		    double sum = 0.0;
		    SNP firstSNP = backgroundSNPList.get(idx);
		    
			for(int i = 0; i < apcAverageNumber; i++) {
				int j = random.nextInt(i, candidateCount);

		        swaps[i] = j;

		        int tmp = randomIndices[i];
		        randomIndices[i] = randomIndices[j];
		        randomIndices[j] = tmp;

		        sum += pairMICalculator.compute(firstSNP, targetSNPs[randomIndices[i]]);
			}
			 // Restore the permutation so the same scratch array can be reused.
		    for (int i = apcAverageNumber - 1; i >= 0; i--) {
		        int j = swaps[i];

		        int tmp = randomIndices[i];
		        randomIndices[i] = randomIndices[j];
		        randomIndices[j] = tmp;
		    }
			return sum / apcAverageNumber;
		}).average().getAsDouble();
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		IntStream.range(0, snpList.size()).parallel().filter(idx -> !singleSNPMI_SigFinderResult.getSignificantIndices().contains(idx)).forEach(idx -> snpList.get(idx).setAverageMItoPheno(backgroundMeanEffect));
		System.out.println("Calculating overall average effect based on all SNPs");
		start = System.nanoTime();
		randSNPIdxList = IntStream.range(0, snpList.size()).boxed().collect(Collectors.toList());
		Collections.shuffle(randSNPIdxList);
		// Recreate the permutation buffers as we now use the total set of all SNPs
		ThreadLocal<int[]> scratch2 =
						ThreadLocal.withInitial(() -> {
				            int[] indices = new int[snpList.size()];
				            for (int i = 0; i < snpList.size(); i++) {
				                indices[i] = i;
				            }
				            return indices;
				        });
		double overallMeanEffect = randSNPIdxList.subList(0, sigSNPList.size()).parallelStream().mapToDouble(idx ->{
			int[] randomIndices = scratch2.get();
		    int[] swaps = swapPositions.get();
		    ThreadLocalRandom random = ThreadLocalRandom.current();
		    
		    double sum = 0.0;
		    SNP firstSNP = snpList.get(idx);
		    
			for(int i = 0; i < apcAverageNumber; i++) {
				int j = random.nextInt(i, snpList.size());

		        swaps[i] = j;

		        int tmp = randomIndices[i];
		        randomIndices[i] = randomIndices[j];
		        randomIndices[j] = tmp;

		        sum += pairMICalculator.compute(firstSNP, snpList.get(randomIndices[i]));
			}
			 // Restore the permutation so the same scratch array can be reused.
		    for (int i = apcAverageNumber - 1; i >= 0; i--) {
		        int j = swaps[i];

		        int tmp = randomIndices[i];
		        randomIndices[i] = randomIndices[j];
		        randomIndices[j] = tmp;
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
				int[] randomIndices = scratch.get();
			    int[] swaps = swapPositions.get();
			    ThreadLocalRandom random = ThreadLocalRandom.current();
			    
			    double sum = 0.0;
			    
				for(int i = 0; i < apcAverageNumber; i++) {
					int j = random.nextInt(i, candidateCount);

			        swaps[i] = j;

			        int tmp = randomIndices[i];
			        randomIndices[i] = randomIndices[j];
			        randomIndices[j] = tmp;

			        sum += pairMICalculator.compute(snp, targetSNPs[randomIndices[i]]);
				}
				 // Restore the permutation so the same scratch array can be reused.
			    for (int i = apcAverageNumber - 1; i >= 0; i--) {
			        int j = swaps[i];

			        int tmp = randomIndices[i];
			        randomIndices[i] = randomIndices[j];
			        randomIndices[j] = tmp;
			    }
				snp.setAverageMItoPheno(sum / apcAverageNumber);
			});
			System.out.println("FileSigSNP-Mean = " + fileSigSNPList.parallelStream().mapToDouble(snp -> snp.getAverageMItoPheno()).average().getAsDouble());
			System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		}
		System.out.println("Calculating MI_APC values for SNP pairs");
		List<SNP> effectiveSigSNPList = (fileSigSNPList != null) ? fileSigSNPList : sigSNPList; 
		long possiblePairCount = calcPossiblePairsCount(effectiveSigSNPList.size(), snpList.size());
		int savePairCount = (int) (possiblePairCount * (keepPercentage / 100.0));
		System.out.println("Saving top " + savePairCount + " pairs (" + keepPercentage + "% from " + possiblePairCount + " calculated pairs)");
		Set<SNP> sigSNPSet = new HashSet<>(effectiveSigSNPList);
		List<SNP> snpWithoutSigList = snpList.stream().filter(snp -> !sigSNPSet.contains(snp)).toList();
		List<SNP> targetSNPList = Stream.concat(effectiveSigSNPList.stream(), snpWithoutSigList.stream()).toList();
		if(targetSNPList.size() != snpList.size()) {
			System.out.println("List sizes are not equal");
			return false;
		}
		start = System.nanoTime();
		int chunkSize = (effectiveSigSNPList.size() + threadCount - 1) / threadCount;
		
		List<TopKHeap_MIAPC> partialQueues = IntStream.range(0, threadCount).parallel().mapToObj(worker -> {
			int startIdx = worker * chunkSize;
			int endIdx = Math.min(startIdx + chunkSize, effectiveSigSNPList.size());
			TopKHeap_MIAPC queue = new TopKHeap_MIAPC(savePairCount);
			for (int i = startIdx; i < endIdx; i++) {
				SNP firstSNP = effectiveSigSNPList.get(i);
				for (int j = i; j < targetSNPList.size(); j++) {
					SNP secondSNP = targetSNPList.get(j);
					double mi = pairMICalculator.compute(firstSNP, secondSNP);
					double mi_apc = mi - (firstSNP.getAverageMItoPheno() * secondSNP.getAverageMItoPheno() / overallMeanEffect);
					queue.offer(i, j, mi, mi_apc);
				}
			}
			return queue;
		}).collect(Collectors.toCollection(ArrayList::new));
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");
		System.out.println("Combining results from threads");
		start = System.nanoTime();
		//Fast enough for now but may be improved by parallelization and more advanced merge
		partialQueues.sort(
			    Comparator.comparingDouble(TopKHeap_MIAPC::min).reversed()
			);
		TopKHeap_MIAPC topHeap = partialQueues.get(0);
		for (int i = 1; i < partialQueues.size(); i++) {
			partialQueues.get(i).forEach((snp1, snp2, mi, mi_apc) -> { 
				topHeap.offer(snp1, snp2, mi, mi_apc); 
			});
		}
		System.out.println("Done in " + (System.nanoTime() - start) / 1_000_000_000 / 60 + " minutes");

		System.out.println("Writing results to file");
		
		try(PrintWriter outPW = new PrintWriter(Files.newBufferedWriter(outFile))) {
			outPW.println("SNP1 + SNP2 MI MI_APC");
			topHeap.drain((snp1Index, snp2Index, mi, mi_apc) -> {
			    String snp1 = effectiveSigSNPList.get(snp1Index).getID();
			    String snp2 = targetSNPList.get(snp2Index).getID();
			    outPW.println(snp1 + " + " + snp2 + " " + mi + " " + mi_apc);
			});
		}
		catch(IOException e) {
			System.out.println("Error while writing results to file");
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	private static long calcPossiblePairsCount(long sigSNPCount, long totalSNPCount) {
		return totalSNPCount * sigSNPCount - (sigSNPCount * (sigSNPCount - 1)) / 2;
	}
	
	private static void parseArgs(String[] args) throws IllegalArgumentException{
		int i = 0;
		while(i < args.length - 2) {
			switch(args[i]) {
			case "-keep" ->{
				keepPercentage = Double.parseDouble(args[++i]);
				if(keepPercentage <= 0 || keepPercentage > 100) {
					throw new IllegalArgumentException("Value for -keep needs to be in the intervall (0, 100]");
				}	
			}
			case "-cont" -> isContinuous = true;
			case "-apc" ->{
				apcAverageNumber = Integer.parseInt(args[++i]);
				if(apcAverageNumber < 1) {
					throw new IllegalArgumentException("Value for -apc needs to be greater than 0");
				}
			}
			case "-k" ->{
				kNext = Integer.parseInt(args[++i]);
				if(kNext < 1) {
					throw new IllegalArgumentException("Value for -k needs to be greater than 0");
				}
			}
			case "-fdr" ->{
				fdr = Double.parseDouble(args[++i]);
				if(fdr <= 0 || fdr > 1) {
					throw new IllegalArgumentException("Value for -fdr needs to be in the interval (0, 1]");
				}
			}
            case "-threads" -> {
            	threadCount = Integer.parseInt(args[++i]);
            	if(threadCount < 1) {
					throw new IllegalArgumentException("Value for -threads needs to be at least 1");
				}
            }
			case "-out" -> outFile = Paths.get(args[++i]);
            case "-list" -> snpListFile = Paths.get(args[++i]);
            case "-disccovariates" -> discCovariatesFile = Paths.get(args[++i]);
            case "-contcovariates" -> contCovariatesFile = Paths.get(args[++i]);
            case "-noapc" -> isNoAPC = true;
            case "-noepi" -> isNoEpi = true;
            case "-all" -> isPrintAll = true;
            default -> throw new IllegalArgumentException("Unknown option: " + args[i]);
			}
			i++;
		}
		tpedFile = Paths.get(args[args.length-2]);
		tfamFile = Paths.get(args[args.length-1]);
		if(outFile == null) {
			outFile = Paths.get(args[args.length-2] + (isNoAPC ? ".epiNoAPC" : ".epi"));
		}
	}
	
	private static void printHelp(){
		System.out.println("""
				Usage: java -jar MIDESP.jar {Options} tpedFile tfamFile
				
				Options:
				-out            file    name of outputfile (default tpedFile.epi)
				-threads        number  number of threads to use (default = Number_of_Cores / 2)
				-keep           number  keep only the top X percentage pairs with highest MI (default = 1)
				-cont                   indicate that the phenotype is continuous
				-k              number  set the value of k for MI estimation for continuous phenotypes (default = 30)
				-fdr            number  set the value of the false discovery rate for finding significantly associated SNPs (default = 0.005)
				-apc            number  set the number of samples that should be used to estimate the average effects of the SNPs (default = 5000)
				-list           file    name of file with list of SNP IDs to analyze instead of using the SNPs that are significant according to their MI value
				-disccovariates	file    name of file that contains discrete covariate variables for the samples as tab-separated list
				-contcovariates	file    name of file that contains continuous covariate variables for the samples as tab-separated list
				-noapc                  indicate that the APC should not be applied
				-noepi                  indicate that no epistatic SNP pairs should be calculated
				-all                    write an additional file containing the MI values for all SNPs (outputfile.allSNPs)
				""");
	}
}