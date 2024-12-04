package midesp.methods;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

import org.apache.commons.math3.special.Gamma;

import midesp.objects.Phenotype;
import midesp.objects.SNP;

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
	
	public static double calcDigamma(double value) {
		return Gamma.digamma(value);
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
	 * Calculates NMI(X;Y)
	 * @param phenotype	as discrete Y
	 * @param snps	as discrete X 
	 * @return
	 */
	public static double calcMI_DiscPheno(Phenotype phenotype, SNP... snps) {
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
			xEntropyInNats = calcEntropyInNatsFromFreqs(xCounts,sampleCount);
		}
		else {
			Arrays.sort(snps, new Comparator<SNP>() { 
				@Override
				public int compare(SNP arg0, SNP arg1) {
					if(arg0.getBitLength() > arg1.getBitLength()) {
						return -1;
					}
					return 1;
				}
			});
			int maxValue = 0;
			for(SNP snp : snps){
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
			xEntropyInNats = calcEntropyInNatsFromFreqs(xCounts,sampleCount);
		}
		if(phenotype.hasDiscCovariate()) {
			if(phenotype.hasContCovariate()) {
				throw new UnsupportedOperationException("Continuous covariates are not yet supported");
			}
			else {
				nats = calcMI_DiscPheno_with_DiscCovariate(phenotype, xVec, xCounts, xBitLength, xBitMax);
			}
		}
		else {
			if(phenotype.hasContCovariate()) {
				throw new UnsupportedOperationException("Continuous covariates are not yet supported");
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
				nats = xEntropyInNats + yEntropyInNats - calcEntropyInNatsFromFreqs(xyCounts,sampleCount);
			}
		}
		xEntropyInLog2 = xEntropyInNats / logtwo;
		natsInLog2 = nats / logtwo;
		mi = Math.min(Math.max(natsInLog2, 0.0), xEntropyInLog2);
		return 2 * (mi / (normFactor + xEntropyInLog2));
	}
	
	/**
	 * Calculates MI(X;Y|V)
	 * @param phenotype	as discrete Y
	 * @param snps	as discrete X 
	 * @param discCovariate as discrete V
	 * @return
	 */
	public static double calcMI_DiscPheno_with_DiscCovariate(Phenotype phenotype, int[] xVec, int[] xCounts, int xBitLength, int xBitMax) {
		int sampleCount;
		int[] xvCounts;
		int[] xyvCounts;
		int[] vVec = phenotype.getDiscCovariateBitValues();
		int[] yvVec = phenotype.getDiscPhenotype_DiscCovariate_BitValues();
		int vBitLength = phenotype.getDiscCovariateBitLength();
		int vBitMax = phenotype.getDiscCovariateBitMax();
		int yvBitLength = phenotype.getDiscPhenotype_DiscCovariate_BitLength();
		int yvBitMax = phenotype.getDiscPhenotype_DiscCovariate_BitMax();
		double vEntropy = phenotype.getDiscCovariateEntropyNats();
		double yvEntropy = phenotype.getDiscPhenotype_DiscCovariateJointEntropyNats();
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
		if(xBitLength > vBitLength) {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += xVec[i];
				byteArray = byteArray << xBitLength;
				byteArray += vVec[i];
				xvCounts[byteArray]++;
			}
		}
		else {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += vVec[i];
				byteArray = byteArray << vBitLength;
				byteArray += xVec[i];	
				xvCounts[byteArray]++;
			}
		}
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
		if(xBitLength > yvBitLength) {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += xVec[i];
				byteArray = byteArray << xBitLength;
				byteArray += yvVec[i];
				xyvCounts[byteArray]++;
			}
		}
		else {
			for(int i = 0; i < sampleCount; i++) {
				int byteArray = 0;
				byteArray += yvVec[i];
				byteArray = byteArray << yvBitLength;
				byteArray += xVec[i];	
				xyvCounts[byteArray]++;
			}
		}
		//H(X,V) + H(Y,V) - H(X,Y,V) - H(V)
		return calcEntropyInNatsFromFreqs(xvCounts,sampleCount) + yvEntropy - calcEntropyInNatsFromFreqs(xyvCounts,sampleCount) - vEntropy;
	}
	
	/**
	 * Calculates NMI(X;Y)
	 * @param phenotype	as continuous Y
	 * @param snps	as discrete X 
	 * @return
	 */
	public static double calcMI_ContPheno(Phenotype phenotype, int k, SNP... snps) {
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
			xEntropyInNats = calcEntropyInNatsFromFreqs(xCounts,sampleCount);
		}
		else {
			Arrays.sort(snps, new Comparator<SNP>() { 
				@Override
				public int compare(SNP arg0, SNP arg1) {
					if(arg0.getBitLength() > arg1.getBitLength()) {
						return -1;
					}
					return 1;
				}
			});
			int maxValue = 0;
			for(SNP snp : snps){
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
			xEntropyInNats = calcEntropyInNatsFromFreqs(xCounts,sampleCount);
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
	
	/**
	 * Calculates MI(X;Y|V)
	 * @param phenotype	as continuous Y
	 * @param snps	as discrete X 
	 * @param discCovariate as discrete V
	 * @return
	 */
	public static double calcMI_ContPheno_with_DiscCovariate(Phenotype phenotype, int k, int[] xVec, int[] xCounts, int xBitLength, int xBitMax) {
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
