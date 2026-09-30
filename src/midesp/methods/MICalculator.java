package midesp.methods;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

import midesp.objects.EntropyCache;
import midesp.objects.GeneralizedBitSet;
import midesp.objects.Phenotype;
import midesp.objects.Phenotype_Legacy;
import midesp.objects.SNP;
import midesp.objects.SNP_Legacy;

public class MICalculator {

	public static final double logtwo = Math.log(2);
	public static final double singleSNPNormFactor = Math.log(3) / Math.log(2);
	private static final double snpPairNormFactor = Math.log(9) / Math.log(2);
	
	/**
	 * Normalizes the given array and returns it
	 * @param orgVec	array that is normalized
	 * @return
	 */
	public static double[] minMaxNormalization(double[] orgVec) {
		double min = Arrays.stream(orgVec).min().getAsDouble();
		double max = Arrays.stream(orgVec).max().getAsDouble();
		if(min == max) {
			throw new ArithmeticException("All entries have the same value");
		}
		double diff = max - min;
		return Arrays.stream(orgVec).map(val -> (val - min) / diff).toArray();
	}
	
	/**
	 * Calculates entropy of given array and returns it in Nats
	 * @param values	array of values for which the entropy is calculated
	 * @param count	the maximum value inside the array
	 * @return
	 */
	public static double calcEntropyInNats(int[] values, int count) {
		int[] freqVec = new int[count];
		for(int i = 0; i < values.length; i++) {
			freqVec[values[i]]++;
		}
		return calcEntropyInNatsFromFreqs(freqVec, values.length);
	}
	
	/**
	 * Calculates entropy of given frequencies and returns it in Nats
	 * @param freqVec	array containing the frequencies of the values
	 * @param length	total number of values for which the entropy is calculated
	 * @return
	 */
	public static double calcEntropyInNatsFromFreqs(int[] freqVec, int length) {
		double entropy = 0.0;
		for(int x : freqVec){
			if(x != 0){
				entropy += ((double) x/length) * Math.log((double) x/length);		
			}
		}	
		return -entropy;
	}
	
	/**
	 * Calculates entropy of given frequencies and returns it in Nats
	 * @param freqVec	array containing the frequencies of the values
	 * @param cache		EntropyCache containing precalculated results for all possible frequencies
	 * @return
	 */
	public static double calcEntropyInNatsFromFreqs_Cached(int[] freqVec, EntropyCache cache) {
		double entropy = 0.0;
		for(int x : freqVec){
			entropy += cache.get(x);
		}	
		return -entropy;
	}
	
	@FunctionalInterface
	public interface SNPCalculator {
	    double compute(SNP snp, double snpEntropy);
	}
	
	public static SNPCalculator resolveMICalculator(boolean isContinuous, EntropyCache entropyCache, Phenotype pheno) {
	    if (isContinuous) { 
		    if(pheno.hasDiscCovariate()) {
		    	if(pheno.hasContCovariate()) {
		    		throw new UnsupportedOperationException("Continuous phenotypes with both discrete and continuous covariates are currently not supported");
		    	}
		    	throw new UnsupportedOperationException("Continuous phenotypes with discrete covariate are currently not supported");
		    }
		    if(pheno.hasContCovariate()) {
		    	throw new UnsupportedOperationException("Continuous phenotypes with continuous covariate are currently not supported");
		    }
		    throw new UnsupportedOperationException("Continuous phenotypes are currently not supported");

	    } else {
	    	if(pheno.hasDiscCovariate()) {
	    		if(pheno.hasContCovariate()) {
		    		throw new UnsupportedOperationException("Discrete phenotypes with both discrete and continuous covariates are currently not supported");
		    	}
	    		double phenoDiscCovariateJointEntropyNats = pheno.getDiscPhenotype_DiscCovariateJointEntropyNats();
	    		double discCovariateEntropyNats = pheno.getDiscCovariateEntropyNats();
	    		return (snp, snpEntropy) ->
	            calcCMI_OneSNP_DiscPheno_DiscCovariate(
	                snp,
	                phenoDiscCovariateJointEntropyNats,
	                discCovariateEntropyNats
	            );
	    	}
    		if(pheno.hasContCovariate()) {
	    		throw new UnsupportedOperationException("Discrete phenotypes with continuous covariates are currently not supported");
	    	}
    		GeneralizedBitSet discPheno = pheno.getDiscPhenotype();
    		double phenoEntropyNats = pheno.getDiscPhenotypeEntropyNats();
    		return (snp, snpEntropy) ->
    		calcMI_OneSNP_DiscPheno(
    				entropyCache,
    				discPheno,
    				snp.getBitSet(),
    				snpEntropy,
    				phenoEntropyNats
    				);
	    }
	}
	
	@FunctionalInterface
	public interface SNPBiCalculator {
	    double compute(SNP first, SNP second);
	}
	
	public static SNPBiCalculator resolvePairMICalculator(boolean isContinuous, EntropyCache entropyCache, Phenotype pheno) {
	    if (isContinuous) {
	    	if(pheno.hasDiscCovariate()) {
		    	if(pheno.hasContCovariate()) {
		    		throw new UnsupportedOperationException("Continuous phenotypes with both discrete and continuous covariates are currently not supported");
		    	}
		    	throw new UnsupportedOperationException("Continuous phenotypes with discrete covariate are currently not supported");
		    }
		    if(pheno.hasContCovariate()) {
		    	throw new UnsupportedOperationException("Continuous phenotypes with continuous covariate are currently not supported");
		    }
		    throw new UnsupportedOperationException("Continuous phenotypes are currently not supported");
	    } else {
	    	if(pheno.hasDiscCovariate()) {
	    		if(pheno.hasContCovariate()) {
		    		throw new UnsupportedOperationException("Discrete phenotypes with both discrete and continuous covariates are currently not supported");
		    	}
	    		double phenoDiscCovariateJointEntropyNats = pheno.getDiscPhenotype_DiscCovariateJointEntropyNats();
	    		double discCovariateEntropyNats = pheno.getDiscCovariateEntropyNats();
	    		return (first, second) ->
	    			calcCMI_TwoSNPs_DiscPheno_DiscCovariate(
	    					entropyCache, 
	    					first, 
	    					second, 
	    					phenoDiscCovariateJointEntropyNats, 
	    					discCovariateEntropyNats);
	    	}
	    	if(pheno.hasContCovariate()) {
	    		throw new UnsupportedOperationException("Discrete phenotypes with continuous covariates are currently not supported");
	    	}
	    	GeneralizedBitSet discPheno = pheno.getDiscPhenotype();
	    	double entropyNats = pheno.getDiscPhenotypeEntropyNats();
	    	return (first, second) ->
	    	calcMI_TwoSNPs_DiscPheno(
	    			entropyCache,
	    			discPheno,
	    			first.getBitSet(),
	    			second.getBitSet(),
	    			entropyNats
	    			);
	    }
	}
	
	@Deprecated
	public static double calcMI_DiscPheno(EntropyCache cache, Phenotype_Legacy phenotype, int k, SNP_Legacy... snps) {
		int sampleCount;
		int xBitLength, xBitMax;
		int[] xVec;
		int[] xCounts;
		int[] xyCounts;
		int[] yVec;
		int yBitLength;
		int yBitMax;
		double xEntropyInNats;
		double yEntropyInNats;
		double nats;
		double xEntropyInLog2;
		double natsInLog2;
		double mi;
		double normFactor;
		if(snps.length == 1) {
			normFactor = singleSNPNormFactor;
		}
		else if(snps.length == 2) {
			normFactor = snpPairNormFactor;
		}
		else {
			throw new IllegalArgumentException("Invalid number of SNPs for MI calculation");
		}
		//Prepare variables
		sampleCount = snps[0].getLength();
		//x
		if(snps.length == 1) {
			xVec = snps[0].getGenotypes();
			xCounts = snps[0].getGenotypesCounts();
			xBitLength = snps[0].getBitLength();
			xBitMax = snps[0].getBitMax();
			xEntropyInNats = snps[0].getEntropyNats();
		}
		else if(snps.length == 2) {
			int maxValue;
			int snp1BitLength;
			int snp1BitMax, snp2BitMax;
			int[] snp1Genotypes, snp2Genotypes;
			if(snps[0].getBitLength() >= snps[1].getBitLength()) {
				snp1BitLength = snps[0].getBitLength();
				snp1BitMax = snps[0].getBitMax();
				snp1Genotypes = snps[0].getGenotypes();
				snp2BitMax = snps[1].getBitMax();
				snp2Genotypes = snps[1].getGenotypes();
			}
			else {
				snp1BitLength = snps[1].getBitLength();
				snp1BitMax = snps[1].getBitMax();
				snp1Genotypes = snps[1].getGenotypes();
				snp2BitMax = snps[0].getBitMax();
				snp2Genotypes = snps[0].getGenotypes();
			}
			maxValue = snp1BitMax;
			maxValue = maxValue << snp1BitLength;
			maxValue += snp2BitMax;
			xCounts = new int[(int)maxValue+1];
			xVec = new int[sampleCount];
			for(int pos_idx = 0; pos_idx < sampleCount; pos_idx++){
				int byteArray = snp1Genotypes[pos_idx];
				byteArray = byteArray << snp1BitLength;
				byteArray += snp2Genotypes[pos_idx];
				xCounts[byteArray]++;
				xVec[pos_idx] = byteArray;
			}
			xBitLength = (int) Math.ceil(Math.log(maxValue+1) / logtwo);
			xBitMax = maxValue+1;
			xEntropyInNats = calcEntropyInNatsFromFreqs_Cached(xCounts,cache);
		}
		else {
			Arrays.sort(snps, new Comparator<SNP_Legacy>() { 
				@Override
				public int compare(SNP_Legacy arg0, SNP_Legacy arg1) {
					if(arg0.getBitLength() > arg1.getBitLength()) {
						return -1;
					}
					return 1;
				}
			});
			int maxValue = 0;
			for(SNP_Legacy snp : snps){
				maxValue += snp.getBitMax();
				maxValue = maxValue << snp.getBitLength();			
			}
			maxValue = maxValue >> snps[snps.length-1].getBitLength();
			xCounts = new int[(int)maxValue+1];
			xVec = new int[sampleCount];
			for(int pos_idx = 0; pos_idx < sampleCount; pos_idx++){
				int byteArray = 0;
				for(int snp_idx = 0; snp_idx < snps.length-1; snp_idx++){
					byteArray += snps[snp_idx].getGenotypes()[pos_idx];
					byteArray = byteArray << snps[snp_idx].getBitLength();
				}
				byteArray += snps[snps.length-1].getGenotypes()[pos_idx];
				xCounts[byteArray]++;
				xVec[pos_idx] = byteArray;
			}
			xBitLength = (int) Math.ceil(Math.log(maxValue+1) / logtwo);
			xBitMax = maxValue+1;
			xEntropyInNats = calcEntropyInNatsFromFreqs_Cached(xCounts,cache);
		}
		if(phenotype.hasDiscCovariate()) {
			if(phenotype.hasContCovariate()) {
				nats = calcMI_DiscPheno_with_BothCovariates(phenotype, k, xVec, xCounts, xBitLength, xBitMax);
			}
			else {
				nats = -1; //New implementation using GeneralizedBitSets ready
			}
		}
		else {
			if(phenotype.hasContCovariate()) {
				nats = calcMI_DiscPheno_with_ContCovariate(phenotype, k, xVec, xCounts, xBitLength, xBitMax);
			}
			else {
				//y
				yVec = phenotype.getDiscPhenotypeBitValues();
				yBitLength = phenotype.getDiscPhenotypeBitLength();
				yBitMax = phenotype.getDiscPhenotypeBitMax();
				yEntropyInNats = phenotype.getDiscPhenotypeEntropyNats();
				//xy
				int maxValue;
				if(xBitLength > yBitLength) {
					maxValue = xBitMax;
					maxValue = maxValue << xBitLength;
					maxValue += yBitMax;
				}
				else {
					maxValue = yBitMax;
					maxValue = maxValue << yBitLength;
					maxValue += xBitMax;
				}	
				xyCounts = new int[(int)maxValue+1];
				if(xBitLength > yBitLength) {
					for(int i = 0; i < sampleCount; i++) {
						int byteArray = 0;
						byteArray += xVec[i];
						byteArray = byteArray << xBitLength;
						byteArray += yVec[i];
						xyCounts[byteArray]++;
					}
				}
				else {
					for(int i = 0; i < sampleCount; i++) {
						int byteArray = 0;
						byteArray += yVec[i];
						byteArray = byteArray << yBitLength;
						byteArray += xVec[i];	
						xyCounts[byteArray]++;
					}
				}
				//H(X) + H(Y) - H(X,Y)
				nats = xEntropyInNats + yEntropyInNats - calcEntropyInNatsFromFreqs_Cached(xyCounts,cache);
			}
		}
		xEntropyInLog2 = xEntropyInNats / logtwo;
		natsInLog2 = nats / logtwo;
		mi = Math.min(Math.max(natsInLog2, 0.0), xEntropyInLog2);
		return 2 * (mi / (normFactor + xEntropyInLog2));
	}
	
	@Deprecated
	public static double calcMI_ContPheno(EntropyCache cache, Phenotype_Legacy phenotype, int k, SNP_Legacy... snps) {
		int sampleCount;
		int numClasses;
		int xBitLength, xBitMax;
		int[] xVec;
		int[] xCounts;
		double xEntropyInNats;
		double nats;
		double xEntropyInLog2;
		double natsInLog2;
		double mi;
		double normFactor;
		if(snps.length == 1) {
			normFactor = singleSNPNormFactor;
		}
		else if(snps.length == 2) {
			normFactor = snpPairNormFactor;
		}
		else {
			throw new IllegalArgumentException("Invalid number of SNPs for MI calculation");
		}
		//Prepare variables
		sampleCount = snps[0].getLength();
		//x
		if(snps.length == 1) {
			xVec = snps[0].getGenotypes();
			xCounts = snps[0].getGenotypesCounts();
			xBitLength = snps[0].getBitLength();
			xBitMax = snps[0].getBitMax();
			xEntropyInNats = snps[0].getEntropyNats();
		}
		else if(snps.length == 2) {
			int maxValue;
			int snp1BitLength;
			int snp1BitMax, snp2BitMax;
			int[] snp1Genotypes, snp2Genotypes;
			if(snps[0].getBitLength() >= snps[1].getBitLength()) {
				snp1BitLength = snps[0].getBitLength();
				snp1BitMax = snps[0].getBitMax();
				snp1Genotypes = snps[0].getGenotypes();
				snp2BitMax = snps[1].getBitMax();
				snp2Genotypes = snps[1].getGenotypes();
			}
			else {
				snp1BitLength = snps[1].getBitLength();
				snp1BitMax = snps[1].getBitMax();
				snp1Genotypes = snps[1].getGenotypes();
				snp2BitMax = snps[0].getBitMax();
				snp2Genotypes = snps[0].getGenotypes();
			}
			maxValue = snp1BitMax;
			maxValue = maxValue << snp1BitLength;
			maxValue += snp2BitMax;
			xCounts = new int[(int)maxValue+1];
			xVec = new int[sampleCount];
			for(int pos_idx = 0; pos_idx < sampleCount; pos_idx++){
				int byteArray = snp1Genotypes[pos_idx];
				byteArray = byteArray << snp1BitLength;
				byteArray += snp2Genotypes[pos_idx];
				xCounts[byteArray]++;
				xVec[pos_idx] = byteArray;
			}
			xBitLength = (int) Math.ceil(Math.log(maxValue+1) / logtwo);
			xBitMax = maxValue+1;
			xEntropyInNats = calcEntropyInNatsFromFreqs_Cached(xCounts,cache);
		}
		else {
			Arrays.sort(snps, new Comparator<SNP_Legacy>() { 
				@Override
				public int compare(SNP_Legacy arg0, SNP_Legacy arg1) {
					if(arg0.getBitLength() > arg1.getBitLength()) {
						return -1;
					}
					return 1;
				}
			});
			int maxValue = 0;
			for(SNP_Legacy snp : snps){
				maxValue += snp.getBitMax();
				maxValue = maxValue << snp.getBitLength();			
			}
			maxValue = maxValue >> snps[snps.length-1].getBitLength();
			xCounts = new int[(int)maxValue+1];
			xVec = new int[sampleCount];
			for(int pos_idx = 0; pos_idx < sampleCount; pos_idx++){
				int byteArray = 0;
				for(int snp_idx = 0; snp_idx < snps.length-1; snp_idx++){
					byteArray += snps[snp_idx].getGenotypes()[pos_idx];
					byteArray = byteArray << snps[snp_idx].getBitLength();
				}
				byteArray += snps[snps.length-1].getGenotypes()[pos_idx];
				xCounts[byteArray]++;
				xVec[pos_idx] = byteArray;
			}
			xBitLength = (int) Math.ceil(Math.log(maxValue+1) / logtwo);
			xBitMax = maxValue+1;
			xEntropyInNats = calcEntropyInNatsFromFreqs_Cached(xCounts,cache);
		}
		if(phenotype.hasDiscCovariate()) {
			if(phenotype.hasContCovariate()) {
				throw new UnsupportedOperationException("Continuous covariates are not yet supported");
			}
			else {
				nats = calcMI_ContPheno_with_DiscCovariate(phenotype, k, xVec, xCounts, xBitLength, xBitMax);
			}
		}
		else {
			if(phenotype.hasContCovariate()) {
				throw new UnsupportedOperationException("Continuous covariates are not yet supported");
			}
			else {
				numClasses = (int) Arrays.stream(xCounts).filter(i -> i != 0).count();
				//y
				double[] digammaValues = phenotype.getDigammaArray();
				int[][] yClosestNeighbours = phenotype.getClosestNeighborsMat();
				double[][] yClosestNeighboursDist = phenotype.getClosestNeighborsDistMat();
				double y_DigammaSum = 0.0;
				double y_X_DigammaSum = 0.0;
				double n_DigammaAvg = digammaValues[sampleCount];
				double n_X_DigammaSum = 0.0;
				for(int i = 0; i < sampleCount; i++) {
					int currentX = xVec[i];
					int currentK = k;
					if(xCounts[currentX] < k+1) {
						if(xCounts[currentX] == 1) { //Case of no neighbour
							//Correction according to the example code provided by Brian C. Ross	
							y_DigammaSum += digammaValues[numClasses * 2];
							y_X_DigammaSum += digammaValues[1];
							n_X_DigammaSum += digammaValues[1];
							continue;
						}
						else { // Case of less than k neighbours
							currentK = xCounts[currentX]-1; //Set k to max. available neighbour
						}
					}
					//Find distance to k-th neighbour for phenotype
					int tmpCounter = 0;
					int kthNeighbour_Pheno = 0;
					for(int j = 1; j < sampleCount; j++) {
						if(xVec[yClosestNeighbours[i][j]] == currentX) {
							tmpCounter++;
							kthNeighbour_Pheno = j;
							if(tmpCounter == currentK) {
								break;
							}
						}
					}
					double epsilonDist = yClosestNeighboursDist[i][kthNeighbour_Pheno];
					
					//Count samples closer than epsilon(i)
					int nY = kthNeighbour_Pheno;
					int nY_X = tmpCounter;
					for(int j = kthNeighbour_Pheno + 1; j < sampleCount; j++) {
						if(yClosestNeighboursDist[i][j] <= epsilonDist){
							//For both without considering X
							nY++;
							if(xVec[yClosestNeighbours[i][j]] == currentX) {
								//For both considering X
								nY_X++;
							}
						}
						else {
							break;
						}
					}
					y_DigammaSum += digammaValues[nY];
					y_X_DigammaSum += digammaValues[nY_X];
					n_X_DigammaSum += digammaValues[xCounts[currentX]];
				}
				nats = - (y_DigammaSum / sampleCount) + (y_X_DigammaSum / sampleCount) + n_DigammaAvg - (n_X_DigammaSum / sampleCount);	
			}	
		}	
		xEntropyInLog2 = xEntropyInNats / logtwo;
		natsInLog2 = nats / logtwo;
		mi = Math.min(Math.max(natsInLog2, 0.0), xEntropyInLog2);
		return 2 * (mi / (normFactor + xEntropyInLog2));
	}
	public static double calcMI_OneSNP_DiscPheno(EntropyCache cache, GeneralizedBitSet phenotype, GeneralizedBitSet snp, double snpEntropy, double phenoEntropy) {
		int c1Count = snp.getNumClasses();
	    int cYCount = phenotype.getNumClasses();
	    int numWords = snp.getNumWords();
	    double hX1Y = 0.0;

	    for (int c1 = 0; c1 < c1Count; c1++) {
	    	long[] m1 = snp.getMask(c1);
	    	for (int cy = 0; cy < cYCount; cy++) {
	    		long[] mY = phenotype.getMask(cy);

	    		// --- HOT LOOP: PURE IN-REGISTER POPCNT ---
	    		int cellCount = 0;
	    		for (int w = 0; w < numWords; w++) {
	    			cellCount += Long.bitCount(m1[w] & mY[w]);
	    		}
	    		// ----------------------------------------

	    		if (cellCount > 0) {
	    			// Add precalculated p*log(p) for H(X1, Y)
	    			hX1Y -= cache.get(cellCount); 
	    		}
	    	}
	    }

	    double nats = snpEntropy + phenoEntropy - hX1Y;
		double xEntropyInLog2 = snpEntropy / logtwo;
		double natsInLog2 = nats / logtwo;
		double mi = Math.min(Math.max(natsInLog2, 0.0), xEntropyInLog2);
		return 2 * (mi / (singleSNPNormFactor + xEntropyInLog2));
	}
	
	public static double calcCMI_OneSNP_DiscPheno_DiscCovariate(SNP snp, double phenoCovariateJointEntropy, double covariateEntropy) {		
		double nats = snp.getDiscCovariateJointEntropyNats() + phenoCovariateJointEntropy - snp.getDiscPhenoDiscCovariateJointEntropyNats() - covariateEntropy;
		double xGivenZEntropyInLog2 = (snp.getDiscCovariateJointEntropyNats() - covariateEntropy) / logtwo;
		double natsInLog2 = nats / logtwo;
		double cmi = Math.min(Math.max(natsInLog2, 0.0), xGivenZEntropyInLog2);
		return 2 * (cmi / (singleSNPNormFactor + xGivenZEntropyInLog2));
	}
	
	public static double calcMI_TwoSNPs_DiscPheno(EntropyCache cache, GeneralizedBitSet phenotype, GeneralizedBitSet snp1, GeneralizedBitSet snp2, double phenoEntropy) {
		int c1Count = snp1.getNumClasses();
	    int c2Count = snp2.getNumClasses();
	    int cYCount = phenotype.getNumClasses();
	    int numWords = snp1.getNumWords();
	    double hX1X2Y = 0.0;
	    double hX1X2 = 0.0;

	    for (int c1 = 0; c1 < c1Count; c1++) {
	        long[] m1 = snp1.getMask(c1);

	        for (int c2 = 0; c2 < c2Count; c2++) {
	            long[] m2 = snp2.getMask(c2);

	            int x1x2Count = 0;

	            for (int cy = 0; cy < cYCount; cy++) {
	                long[] mY = phenotype.getMask(cy);

	                // --- HOT LOOP: PURE IN-REGISTER POPCNT ---
	                int cellCount = 0;
	                for (int w = 0; w < numWords; w++) {
	                    cellCount += Long.bitCount(m1[w] & m2[w] & mY[w]);
	                }
	                // ----------------------------------------

	                if (cellCount > 0) {
	                    x1x2Count += cellCount;
	                    // Add precalculated p*log(p) for H(X1, X2, Y)
	                    hX1X2Y -= cache.get(cellCount); 
	                }
	            }

	            if (x1x2Count > 0) {
	                // Add precalculated p*log(p) for H(X1, X2)
	                hX1X2 -= cache.get(x1x2Count); 
	            }
	        }
	    }

	    double nats = hX1X2 + phenoEntropy - hX1X2Y;
		double xEntropyInLog2 = hX1X2 / logtwo;
		double natsInLog2 = nats / logtwo;
		double mi = Math.min(Math.max(natsInLog2, 0.0), xEntropyInLog2);
		return 2 * (mi / (snpPairNormFactor + xEntropyInLog2));
	}
	
	public static double calcCMI_TwoSNPs_DiscPheno_DiscCovariate(EntropyCache cache, SNP snp1, SNP snp2, double phenoCovariateJointEntropy, double covariateEntropy) {
		// Calculate H(X1, X2, Z) using the cheaper direction for Z
		boolean dirA_Z = (snp1.getSNPDiscCovariateBitSet().getNumClasses() * snp2.getBitSet().getNumClasses())
				<= (snp2.getSNPDiscCovariateBitSet().getNumClasses() * snp1.getBitSet().getNumClasses());
		GeneralizedBitSet maskZ = dirA_Z ? snp1.getSNPDiscCovariateBitSet() : snp2.getSNPDiscCovariateBitSet();
		GeneralizedBitSet rawForZ = dirA_Z ? snp2.getBitSet() : snp1.getBitSet();

	    int cZCount = maskZ.getNumClasses();
	    int cRawForZCount = rawForZ.getNumClasses();
	    int numWords = rawForZ.getNumWords();
	    double hX1X2Z = 0.0;

	    for (int cZ = 0; cZ < cZCount; cZ++) {
	    	long[] mZ = maskZ.getMask(cZ);
	    	for (int cR = 0; cR < cRawForZCount; cR++) {
	    		long[] mR = rawForZ.getMask(cR);

	    		// --- HOT LOOP: PURE IN-REGISTER POPCNT ---
	    		int cellCount = 0;
	    		for (int w = 0; w < numWords; w++) {
	    			cellCount += Long.bitCount(mZ[w] & mR[w]);
	    		}
	    		// ----------------------------------------

	    		if (cellCount > 0) {
	    			// Add precalculated p*log(p) for H(X1,X2,Z)
	    			hX1X2Z -= cache.get(cellCount); 
	    		}
	    	}
	    }
	    
	    // Calculate H(X1, X2, Z) using the cheaper direction for Z
	    boolean dirA_YZ = (snp1.getSNPDiscPhenoDiscCovariateBitSet().getNumClasses() * snp2.getBitSet().getNumClasses())
	    		<= (snp2.getSNPDiscPhenoDiscCovariateBitSet().getNumClasses() * snp1.getBitSet().getNumClasses());
	    GeneralizedBitSet maskYZ = dirA_YZ ? snp1.getSNPDiscPhenoDiscCovariateBitSet() : snp2.getSNPDiscPhenoDiscCovariateBitSet();
	    GeneralizedBitSet rawForYZ = dirA_YZ ? snp2.getBitSet() : snp1.getBitSet();

	    int cYZCount = maskYZ.getNumClasses();
	    int cRawForYZCount = rawForYZ.getNumClasses();
	    double hX1X2YZ = 0.0;

	    for (int cYZ = 0; cYZ < cYZCount; cYZ++) {
	    	long[] mYZ = maskYZ.getMask(cYZ);
	    	for (int cR = 0; cR < cRawForYZCount; cR++) {
	    		long[] mR = rawForYZ.getMask(cR);

	    		// --- HOT LOOP: PURE IN-REGISTER POPCNT ---
	    		int cellCount = 0;
	    		for (int w = 0; w < numWords; w++) {
	    			cellCount += Long.bitCount(mYZ[w] & mR[w]);
	    		}
	    		// ----------------------------------------

	    		if (cellCount > 0) {
	    			// Add precalculated p*log(p) for H(X1,X2,Y,Z)
	    			hX1X2YZ -= cache.get(cellCount); 
	    		}
	    	}
	    }

	    double nats = hX1X2Z + phenoCovariateJointEntropy - hX1X2YZ - covariateEntropy;
		double xGivenZEntropyInLog2 = Math.max(0.0, (hX1X2Z - covariateEntropy) / logtwo);
		double natsInLog2 = nats / logtwo;
		double cmi = Math.min(Math.max(natsInLog2, 0.0), xGivenZEntropyInLog2);
		return 2 * (cmi / (snpPairNormFactor + xGivenZEntropyInLog2));
	}
	
	/**
	 * Calculates MI(X;Y|W)
	 * @param phenotype	as discrete Y
	 * @param snps	as discrete X 
	 * @param contCovariate as continuous W
	 * @return
	 */
	public static double calcMI_DiscPheno_with_ContCovariate(Phenotype_Legacy phenotype, int k, int[] xVec, int[] xCounts, int xBitLength, int xBitMax) {
		int sampleCount;
		int[] xyCounts;
		int[] xyVec;
		int yBitLength = phenotype.getDiscPhenotypeBitLength();
		int yBitMax = phenotype.getDiscPhenotypeBitMax();
		int[] yVec = phenotype.getDiscPhenotypeBitValues();
		int[] yCounts = phenotype.getDiscPhenotypeBitCounts();
		double yEntropy = phenotype.getDiscPhenotypeEntropyNats();
		int[][] wClosestNeighbors = phenotype.getContCovariate_ClosestNeighborsMat();
		double[][] wClosestNeighborsDist = phenotype.getContCovariate_ClosestNeighborsDistMat();
		//Prepare variables
		sampleCount = xVec.length;
		double xEntropy = calcEntropyInNatsFromFreqs(xCounts, sampleCount);
		//xy
		int maxValue;
		if(xBitLength > yBitLength) {
			maxValue = xBitMax;
			maxValue = maxValue << xBitLength;
			maxValue += yBitMax;
		}
		else {
			maxValue = yBitMax;
			maxValue = maxValue << yBitLength;
			maxValue += xBitMax;
		}
		xyCounts = new int[(int)maxValue+1];
		xyVec = new int[sampleCount];
		if(xBitLength > yBitLength) {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += xVec[i];
				byteArray = byteArray << xBitLength;
				byteArray += yVec[i];
				xyCounts[byteArray]++;
				xyVec[i] = byteArray;
			}
		}
		else {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += yVec[i];
				byteArray = byteArray << yBitLength;
				byteArray += xVec[i];	
				xyCounts[byteArray]++;
				xyVec[i] = byteArray;
			}
		}
		double xyEntropy = calcEntropyInNatsFromFreqs(xyCounts, sampleCount);
		double[] digammaValues = phenotype.getDigammaArray();
		double w_DigammaSum = 0;
		double w_X_DigammaSum = 0;
		double w_Y_DigammaSum = 0;
		double w_XY_DigammaSum = 0;
		double n_DigammaSum = 0;
		double n_X_DigammaSum = 0;
		double n_Y_DigammaSum = 0;
		double n_XY_DigammaSum = 0;
		
		for(int i = 0; i < sampleCount; i++) {
			int currentX = xVec[i];
			int currentY = yVec[i];
			int currentXY = xyVec[i];
			int currentK = k;
			if(xyCounts[currentXY] < k+1) {
				if(xyCounts[currentXY] == 1) { //Case of no neighbour
					//Correction according to the example code provided by Brian C. Ross	
					//Slight adjustment to limit the number of classes for w_X and w_Y to the actual samples with classes currentX and currentY
					int numClassesWithX = (int) IntStream.range(0, sampleCount).filter(idx -> xVec[idx] == currentX).map(idx -> xyVec[idx]).distinct().count();
					int numClassesWithY = (int) IntStream.range(0, sampleCount).filter(idx -> yVec[idx] == currentY).map(idx -> xyVec[idx]).distinct().count();
					int numClasses = (int) Arrays.stream(xyCounts).filter(xy -> xy != 0).count();
					w_DigammaSum += digammaValues[numClasses * 2];
					w_X_DigammaSum += digammaValues[numClassesWithX * 2 > xCounts[currentX] ? xCounts[currentX] : numClassesWithX * 2];
					w_Y_DigammaSum += digammaValues[numClassesWithY * 2 > yCounts[currentY] ? yCounts[currentY] : numClassesWithY * 2];
					w_XY_DigammaSum += digammaValues[1];
					n_DigammaSum += digammaValues[sampleCount];
					n_X_DigammaSum += digammaValues[xCounts[currentX]];
					n_Y_DigammaSum += digammaValues[yCounts[currentY]];
					n_XY_DigammaSum += digammaValues[1];
					continue;					
				}
				else { // Case of less than k neighbours
					currentK = xyCounts[currentXY]-1; //Set k to max. available neighbour
				}
			}
			//Find distance to k-th neighbour for covariate
			int nW_X = 0;
			int nW_Y = 0;
			int tmpCounter = 0;
			int kthNeighbour_Covariate = 0;
			for(int j = 1; j < sampleCount; j++) {
				if(xVec[wClosestNeighbors[i][j]] == currentX) {
					nW_X++;
				}
				if(yVec[wClosestNeighbors[i][j]] == currentY) {
					nW_Y++;
				}
				if(xyVec[wClosestNeighbors[i][j]] == currentXY) {
					tmpCounter++;
					kthNeighbour_Covariate = j;
					if(tmpCounter == currentK) {
						break;
					}
				}
			}

			double epsilonDist = wClosestNeighborsDist[i][kthNeighbour_Covariate];
			//Count samples closer than epsilon(i)
			int nW = kthNeighbour_Covariate;
			int nW_XY = tmpCounter;
			for(int j = kthNeighbour_Covariate + 1; j < sampleCount; j++) {
				if(wClosestNeighborsDist[i][j] <= epsilonDist) {
					nW++;
					if(xVec[wClosestNeighbors[i][j]] == currentX) {
						nW_X++;
					}
					if(yVec[wClosestNeighbors[i][j]] == currentY) {
						nW_Y++;
					}
					if(xyVec[wClosestNeighbors[i][j]] == currentXY) {
						//For both considering X
						nW_XY++;
					}
				}
				else {
					break;
				}
			}
			w_DigammaSum += digammaValues[nW];
			w_X_DigammaSum += digammaValues[nW_X];
			w_Y_DigammaSum += digammaValues[nW_Y];
			w_XY_DigammaSum += digammaValues[nW_XY];
			n_DigammaSum += digammaValues[sampleCount];
			n_X_DigammaSum += digammaValues[xCounts[currentX]];
			n_Y_DigammaSum += digammaValues[yCounts[currentY]];
			n_XY_DigammaSum += digammaValues[xyCounts[currentXY]];
		}
		return (w_DigammaSum / sampleCount) - (n_DigammaSum / sampleCount) + (w_XY_DigammaSum / sampleCount) - (n_XY_DigammaSum / sampleCount) - (w_Y_DigammaSum / sampleCount) + (n_Y_DigammaSum / sampleCount) - (w_X_DigammaSum / sampleCount) + (n_X_DigammaSum / sampleCount) + xEntropy + yEntropy - xyEntropy; 
	}
	
	/**
	 * Calculates MI(X;Y|V,W)
	 * @param phenotype	as discrete Y
	 * @param snps	as discrete X 
	 * @param discCovariate as discrete V
	 * @param contCovariate as continuous W
	 * @return
	 */
	public static double calcMI_DiscPheno_with_BothCovariates(Phenotype_Legacy phenotype, int k, int[] xVec, int[] xCounts, int xBitLength, int xBitMax) {
		int sampleCount;
		int[] xvCounts;
		int[] xvVec;
		int[] xyvCounts;
		int[] xyvVec;
		int vBitLength = phenotype.getDiscCovariateBitLength();
		int vBitMax = phenotype.getDiscCovariateBitMax();
		int[] vVec = phenotype.getDiscCovariateBitValues();
		int[] vCounts = phenotype.getDiscCovariateBitCounts();
		int yvBitLength = phenotype.getDiscPhenotype_DiscCovariate_BitLength();
		int yvBitMax = phenotype.getDiscPhenotype_DiscCovariate_BitMax();
		int[] yvVec = phenotype.getDiscPhenotype_DiscCovariate_BitValues();
		int[] yvCounts = phenotype.getDiscPhenotype_DiscCovariate_BitCounts();
		double vEntropy = phenotype.getDiscCovariateEntropyNats();
		double yvEntropy = phenotype.getDiscPhenotype_DiscCovariateJointEntropyNats();
		int[][] wClosestNeighbors = phenotype.getContCovariate_ClosestNeighborsMat();
		double[][] wClosestNeighborsDist = phenotype.getContCovariate_ClosestNeighborsDistMat();
		//Prepare variables
		sampleCount = xVec.length;
		//xv
		int maxValue;
		if(xBitLength > vBitLength) {
			maxValue = xBitMax;
			maxValue = maxValue << xBitLength;
			maxValue += vBitMax;
		}
		else {
			maxValue = vBitMax;
			maxValue = maxValue << vBitLength;
			maxValue += xBitMax;
		}
		xvCounts = new int[(int)maxValue+1];
		xvVec = new int[sampleCount];
		if(xBitLength > vBitLength) {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += xVec[i];
				byteArray = byteArray << xBitLength;
				byteArray += vVec[i];
				xvCounts[byteArray]++;
				xvVec[i] = byteArray;
			}
		}
		else {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += vVec[i];
				byteArray = byteArray << vBitLength;
				byteArray += xVec[i];	
				xvCounts[byteArray]++;
				xvVec[i] = byteArray;
			}
		}
		double xvEntropy = calcEntropyInNatsFromFreqs(xvCounts, sampleCount);
		//xyv
		if(xBitLength > yvBitLength) {
			maxValue = xBitMax;
			maxValue = maxValue << xBitLength;
			maxValue += yvBitMax;
		}
		else {
			maxValue = yvBitMax;
			maxValue = maxValue << yvBitLength;
			maxValue += xBitMax;
		}
		xyvCounts = new int[(int)maxValue+1];
		xyvVec = new int[sampleCount];
		if(xBitLength > yvBitLength) {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += xVec[i];
				byteArray = byteArray << xBitLength;
				byteArray += yvVec[i];
				xyvCounts[byteArray]++;
				xyvVec[i] = byteArray;
			}
		}
		else {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += yvVec[i];
				byteArray = byteArray << yvBitLength;
				byteArray += xVec[i];	
				xyvCounts[byteArray]++;
				xyvVec[i] = byteArray;
			}
		}
		double xyvEntropy = calcEntropyInNatsFromFreqs(xyvCounts, sampleCount);

		double[] digammaValues = phenotype.getDigammaArray();
		double w_YV_DigammaSum = 0;
		double w_XYV_DigammaSum = 0;
		double w_V_DigammaSum = 0;
		double w_XV_DigammaSum = 0;
		double n_YV_DigammaSum = 0;
		double n_XYV_DigammaSum = 0;
		double n_V_DigammaSum = 0;
		double n_XV_DigammaSum = 0;
		
		for(int i = 0; i < sampleCount; i++) {
			int currentV = vVec[i];
			int currentXV = xvVec[i];
			int currentYV = yvVec[i];
			int currentXYV = xyvVec[i];
			int currentK = k;
			if(xyvCounts[currentXYV] < k+1) {
				if(xyvCounts[currentXYV] == 1) { //Case of no neighbour
					//Correction according to the example code provided by Brian C. Ross	
					//Slight adjustment to limit the number of classes for w_V, w_YV and w_XV to the actual samples with classes currentV, currentYV and currentXV
					int numClassesWithV = (int) IntStream.range(0, sampleCount).filter(idx -> vVec[idx] == currentV).map(idx -> xyvVec[idx]).distinct().count();
					int numClassesWithYV = (int) IntStream.range(0, sampleCount).filter(idx -> yvVec[idx] == currentYV).map(idx -> xyvVec[idx]).distinct().count();
					int numClassesWithXV = (int) IntStream.range(0, sampleCount).filter(idx -> xvVec[idx] == currentXV).map(idx -> xyvVec[idx]).distinct().count();
					w_V_DigammaSum += digammaValues[numClassesWithV * 2 > vCounts[currentV] ? vCounts[currentV] : numClassesWithV * 2];
					w_XV_DigammaSum += digammaValues[numClassesWithXV * 2 > xvCounts[currentXV] ? xvCounts[currentXV] : numClassesWithXV * 2];
					w_YV_DigammaSum += digammaValues[numClassesWithYV * 2 > yvCounts[currentYV] ? yvCounts[currentYV] : numClassesWithYV * 2];
					w_XYV_DigammaSum += digammaValues[1];
					n_V_DigammaSum += digammaValues[vCounts[currentV]];
					n_XV_DigammaSum += digammaValues[xvCounts[currentXV]];
					n_YV_DigammaSum += digammaValues[yvCounts[currentYV]];
					n_XYV_DigammaSum += digammaValues[1];
					continue;					
				}
				else { // Case of less than k neighbours
					currentK = xyvCounts[currentXYV]-1; //Set k to max. available neighbour
				}
			}
			//Find distance to k-th neighbour for covariate
			int nW_V = 0;
			int nW_XV = 0;
			int nW_YV = 0;
			int tmpCounter = 0;
			int kthNeighbour_Covariate = 0;
			for(int j = 1; j < sampleCount; j++) {
				if(vVec[wClosestNeighbors[i][j]] == currentV) {
					nW_V++;
				}
				if(xvVec[wClosestNeighbors[i][j]] == currentXV) {
					nW_XV++;
				}
				if(yvVec[wClosestNeighbors[i][j]] == currentYV) {
					nW_YV++;
				}
				if(xyvVec[wClosestNeighbors[i][j]] == currentXYV) {
					tmpCounter++;
					kthNeighbour_Covariate = j;
					if(tmpCounter == currentK) {
						break;
					}
				}
			}

			double epsilonDist = wClosestNeighborsDist[i][kthNeighbour_Covariate];
			//Count samples closer than epsilon(i)
			int nW_XYV = tmpCounter;
			for(int j = kthNeighbour_Covariate + 1; j < sampleCount; j++) {
				if(wClosestNeighborsDist[i][j] <= epsilonDist) {
					if(vVec[wClosestNeighbors[i][j]] == currentV) {
						nW_V++;
					}
					if(xvVec[wClosestNeighbors[i][j]] == currentXV) {
						nW_XV++;
					}
					if(yvVec[wClosestNeighbors[i][j]] == currentYV) {
						nW_YV++;
					}
					if(xyvVec[wClosestNeighbors[i][j]] == currentXYV) {
						nW_XYV++;
					}
				}
				else {
					break;
				}
			}
			w_V_DigammaSum += digammaValues[nW_V];
			w_XV_DigammaSum += digammaValues[nW_XV];
			w_YV_DigammaSum += digammaValues[nW_YV];
			w_XYV_DigammaSum += digammaValues[nW_XYV];
			n_V_DigammaSum += digammaValues[vCounts[currentV]];
			n_XV_DigammaSum += digammaValues[xvCounts[currentXV]];
			n_YV_DigammaSum += digammaValues[yvCounts[currentYV]];
			n_XYV_DigammaSum += digammaValues[xyvCounts[currentXYV]];
		}
		
		return - (w_YV_DigammaSum / sampleCount) + (n_YV_DigammaSum / sampleCount) + (w_XYV_DigammaSum / sampleCount) - (n_XYV_DigammaSum / sampleCount) + (w_V_DigammaSum / sampleCount) - (n_V_DigammaSum / sampleCount) - (w_XV_DigammaSum / sampleCount) + (n_XV_DigammaSum / sampleCount) + xvEntropy + yvEntropy - xyvEntropy - vEntropy; 
	}
	
	/**
	 * Calculates MI(X;Y|V)
	 * @param phenotype	as continuous Y
	 * @param snps	as discrete X 
	 * @param discCovariate as discrete V
	 * @return
	 */
	public static double calcMI_ContPheno_with_DiscCovariate(Phenotype_Legacy phenotype, int k, int[] xVec, int[] xCounts, int xBitLength, int xBitMax) {
		int sampleCount;
		int[] xvCounts;
		int[] xvVec;
		int[] vVec = phenotype.getDiscCovariateBitValues();
		int[] vCounts = phenotype.getDiscCovariateBitCounts();
		int vBitLength = phenotype.getDiscCovariateBitLength();
		int vBitMax = phenotype.getDiscCovariateBitMax();
		//Prepare variables
		sampleCount = xVec.length;
		//xv
		int maxValue;
		if(xBitLength > vBitLength) {
			maxValue = xBitMax;
			maxValue = maxValue << xBitLength;
			maxValue += vBitMax;
		}
		else {
			maxValue = vBitMax;
			maxValue = maxValue << vBitLength;
			maxValue += xBitMax;
		}
		xvCounts = new int[(int)maxValue+1];
		xvVec = new int[sampleCount];
		if(xBitLength > vBitLength) {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += xVec[i];
				byteArray = byteArray << xBitLength;
				byteArray += vVec[i];
				xvCounts[byteArray]++;
				xvVec[i] = byteArray;
			}
		}
		else {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += vVec[i];
				byteArray = byteArray << vBitLength;
				byteArray += xVec[i];	
				xvCounts[byteArray]++;
				xvVec[i] = byteArray;
			}
		}
		//y
		double[] digammaValues = phenotype.getDigammaArray();
		int[][] yClosestNeighbours = phenotype.getClosestNeighborsMat();
		double[][] yClosestNeighboursDist = phenotype.getClosestNeighborsDistMat();
		double y_V_DigammaSum = 0;
		double y_XV_DigammaSum = 0;
		double n_V_DigammaSum = 0;
		double n_XV_DigammaSum = 0;
		for(int i = 0; i < sampleCount; i++) {
			int currentV = vVec[i];
			int currentXV = xvVec[i];
			int currentK = k;
			if(xvCounts[currentXV] < k+1) {
				if(xvCounts[currentXV] == 1) { //Case of no neighbour
					//Correction according to the example code provided by Brian C. Ross	
					//Slight adjustment to limit the number of classes for y_V to the actual samples with class currentV
					int numClassesWithV = (int) IntStream.range(0, sampleCount).filter(idx -> vVec[idx] == currentV).map(idx -> xvVec[idx]).distinct().count();
					y_V_DigammaSum += digammaValues[numClassesWithV * 2 > vCounts[currentV] ? vCounts[currentV] : numClassesWithV * 2];
					y_XV_DigammaSum += digammaValues[1];
					n_V_DigammaSum += digammaValues[vCounts[currentV]];
					n_XV_DigammaSum += digammaValues[1];
					continue;					
				}
				else { // Case of less than k neighbours
					currentK = xvCounts[currentXV]-1; //Set k to max. available neighbour
				}
			}
			//Find distance to k-th neighbour for phenotype
			int nY_V = 0;
			int tmpCounter = 0;
			int kthNeighbour_Pheno = 0;
			for(int j = 1; j < sampleCount; j++) {
				if(vVec[yClosestNeighbours[i][j]] == currentV) {
					nY_V++;
				}
				if(xvVec[yClosestNeighbours[i][j]] == currentXV) {
					tmpCounter++;
					kthNeighbour_Pheno = j;
					if(tmpCounter == currentK) {
						break;
					}
				}
			}
			double epsilonDist = yClosestNeighboursDist[i][kthNeighbour_Pheno];
			//Count samples closer than epsilon(i)
			int nY_XV = tmpCounter;
			for(int j = kthNeighbour_Pheno + 1; j < sampleCount; j++) {
				if(yClosestNeighboursDist[i][j] <= epsilonDist){
					if(vVec[yClosestNeighbours[i][j]] == currentV) {
						//For both considering V
						nY_V++;
					}
					if(xvVec[yClosestNeighbours[i][j]] == currentXV) {
						//For both considering X and V
						nY_XV++;
					}
				}
				else {
					break;
				}
			}
			y_V_DigammaSum += digammaValues[nY_V];
			y_XV_DigammaSum += digammaValues[nY_XV];
			n_V_DigammaSum += digammaValues[vCounts[currentV]];
			n_XV_DigammaSum += digammaValues[xvCounts[currentXV]];
		}
		return - (y_V_DigammaSum / sampleCount) + (n_V_DigammaSum / sampleCount) + (y_XV_DigammaSum / sampleCount) - (n_XV_DigammaSum / sampleCount);
	}
}
