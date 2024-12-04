package midesp.objects;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import midesp.methods.MICalculator;


public class Phenotype{

	private String id;
	private int length;
	private boolean isContinuous;
	private boolean hasDiscCovariate;
	private boolean hasContCovariate;
	private int[] discPhenotypeVec;
	private int[] discCovariate_bitValues;
	private int[] discCovariate_bitCounts;
	private int[] discPhenotype_discCovariate_bitValues;
	private int[] bitValues;
	private int bitLength;
	private int bitMax;
	private int discCovariate_bitLength;
	private int discCovariate_bitMax;
	private int discPhenotype_discCovariate_bitLength;
	private int discPhenotype_discCovariate_bitMax;
	private double discPhenotypeEntropyNats;
	private double discCovariateEntropyNats;
	private double discPhenotype_discCovariate_JointEntropyNats;
	private double[] contPhenotypeVec;
	private double[] digammaValuesArray;
	private int[][] closestNeighborsMat;
	private double[][] closestNeighborsDistMat;
	
	public Phenotype(String id, int length, boolean continuous) {
		this.id = id;
		this.length = length;
		isContinuous = continuous;
		hasDiscCovariate = false;
		hasContCovariate = false;
		if(isContinuous) {
			contPhenotypeVec = new double[length];
		}
		else {
			discPhenotypeVec = new int[length];
		}
	}
	
	public int getLength() {
		return length;
	}
	
	public boolean isContinuous() {
		return isContinuous;
	}
	
	public boolean hasDiscCovariate() {
		return hasDiscCovariate;
	}
	
	public boolean hasContCovariate() {
		return hasContCovariate;
	}
	
	public double[] getContPhenotype() {
		return contPhenotypeVec;
	}
	
	public int[] getDiscPhenotype() {
		return discPhenotypeVec;
	}
	
	public double getDiscPhenotypeEntropyNats() {
		return discPhenotypeEntropyNats;
	}
	
	public int getDiscPhenotypeBitMax() {
		return bitMax;
	}
	
	public int getDiscPhenotypeBitLength() {
		return bitLength;
	}
	
	public int[] getDiscPhenotypeBitValues() {
		return bitValues;
	}
	
	public int[][] getClosestNeighborsMat(){
		return closestNeighborsMat;
	}
	
	public double[][] getClosestNeighborsDistMat(){
		return closestNeighborsDistMat;
	}
	
	public double[] getDigammaArray() {
		return digammaValuesArray;
	}
	
	public int[] getDiscCovariateBitValues() {
		return discCovariate_bitValues;
	}

	public int[] getDiscCovariateBitCounts() {
		return discCovariate_bitCounts;
	}
	public int getDiscCovariateBitMax() {
		return discCovariate_bitMax;
	}
	
	public int getDiscCovariateBitLength() {
		return discCovariate_bitLength;
	}
	
	public double getDiscCovariateEntropyNats() {
		return discCovariateEntropyNats;
	}
	
	public double getDiscPhenotype_DiscCovariateJointEntropyNats() {
		return discPhenotype_discCovariate_JointEntropyNats;
	}
	
	public int getDiscPhenotype_DiscCovariate_BitLength() {
		return discPhenotype_discCovariate_bitLength;
	}

	public int getDiscPhenotype_DiscCovariate_BitMax() {
		return discPhenotype_discCovariate_bitMax;
	}

	public int[] getDiscPhenotype_DiscCovariate_BitValues() {
		return discPhenotype_discCovariate_bitValues;
	}

	public void setValueAt(int idx, String value) {
		if(isContinuous) {
			contPhenotypeVec[idx] = Double.parseDouble(value);
		}
		else {
			discPhenotypeVec[idx] = Integer.parseInt(value);
		}
	}	
	
	public void parseValues() {	
		if(isContinuous) {
			//Create sorted lists of neighbors for each entry
			closestNeighborsMat = IntStream.range(0, length).mapToObj(index ->{
				double currentValue = contPhenotypeVec[index];
				List<Pair<Integer,Double>> pairList = new ArrayList<>();
				for(int i = 0; i < length; i++) {
					if(i != index) {
						pairList.add(new Pair<>(i, Math.abs(currentValue - contPhenotypeVec[i])));
					}
				}
				pairList.sort(Comparator.comparing(o -> o.getSecond()));
				return Stream.concat(Stream.of(index), pairList.stream().mapToInt(pair -> pair.getFirst()).boxed()).mapToInt(i->i).toArray();
			}).toArray(int[][]::new);
			//Calculate distances between the entry and its neighbors
			closestNeighborsDistMat = IntStream.range(0, length).mapToObj(i ->{
				return Arrays.stream(closestNeighborsMat[i]).mapToDouble(j -> {
					return Math.abs(contPhenotypeVec[i] - contPhenotypeVec[j]);
				}).toArray();
			}).toArray(double[][]::new);
			//Precalculate all possible digamma values needed for this phenotype 
			digammaValuesArray = IntStream.range(0, length+1).mapToDouble(i -> MICalculator.calcDigamma(i)).toArray();
		}
		else {
			Map<Integer,Byte> bitMap = new HashMap<>();
			byte counter = 0;
			bitValues = new int[length];
			for(int i = 0; i < length; i++) {
				if(!bitMap.containsKey(discPhenotypeVec[i])) {
					bitMap.put(discPhenotypeVec[i], counter);
					counter++;
				}
				bitValues[i] = bitMap.get(discPhenotypeVec[i]);
			}
			bitLength = (int) Math.ceil(Math.log(counter) / MICalculator.logtwo);
			bitMax = counter-1;
			discPhenotypeEntropyNats = MICalculator.calcEntropyInNats(bitValues,counter);
		}
	}
	
	public static Phenotype readTFam(Path tfamFile, boolean isContinuous) throws IOException {
		List<String> values = Files.lines(tfamFile).map(line -> line.split(" ")[5]).collect(Collectors.toList());
		Phenotype pheno = new Phenotype("Phenotype", values.size(), isContinuous);
		for(int i = 0; i < values.size(); i++) {
			pheno.setValueAt(i, values.get(i));
		}
		pheno.parseValues();
		return pheno;
	}
	
	public void readDiscCovariateFile(Path covariateFile) throws IOException{
		List<String[]> covariateList = Files.lines(covariateFile).map(str -> str.split("\t")).toList();
		if(covariateList.size() != this.length) {
			throw new IOException("Number of values for covariate (" + covariateList.size() + ") is different from number of samples (" + this.length + ")");
		}
		int covariateCount = covariateList.get(0).length;
		for(int i = 1; i < covariateList.size(); i++) {
			if(covariateCount != covariateList.get(i).length) {
				throw new IOException("Number of covariates in line " + (i+1) +" (" + covariateList.get(i).length + ") is different from number of covariates in line 1 (" + covariateCount + ")");
			}
		}
		
		int[][] covariateMat = new int[this.length][covariateCount];
		for(int i = 0; i < covariateCount; i++) {
			Map<String, Integer> valueToNumber = new HashMap<>();
			for(int j = 0; j < this.length; j++) {
				String value = covariateList.get(j)[i];
				if(!valueToNumber.containsKey(value)) {
					valueToNumber.put(value, valueToNumber.size());
				}
				covariateMat[j][i] = valueToNumber.get(value);
			}
		}
		this.discCovariate_bitValues = new int[this.length];
		//Combine the different discrete covariates to a single covariate
		Map<String, Byte> valueToNumber = new HashMap<>();
		byte counter = 0;
		for(int i = 0; i < this.length; i++) {
			String combinedValue = Arrays.toString(covariateMat[i]);
			if(!valueToNumber.containsKey(combinedValue)) {
				valueToNumber.put(combinedValue, counter);
				counter++;
			}
			this.discCovariate_bitValues[i] = valueToNumber.get(combinedValue);
		}
		discCovariate_bitLength = (int) Math.ceil(Math.log(counter) / MICalculator.logtwo);
		discCovariate_bitMax = counter-1;
		discCovariateEntropyNats = MICalculator.calcEntropyInNats(discCovariate_bitValues,counter);
		discCovariate_bitCounts = new int[counter];
		for(int i = 0; i < length; i++) {
			discCovariate_bitCounts[discCovariate_bitValues[i]]++;
		}
		if(!isContinuous) {
			discPhenotype_discCovariate_bitValues = new int[this.length];
			valueToNumber = new HashMap<>();
			counter = 0;
			for(int i = 0; i < this.length; i++) {
				String combinedValue = discCovariate_bitValues[i] + "_" + bitValues[i];
				if(!valueToNumber.containsKey(combinedValue)) {
					valueToNumber.put(combinedValue, counter);
					counter++;
				}
				discPhenotype_discCovariate_bitValues[i] = valueToNumber.get(combinedValue); 
			}
			discPhenotype_discCovariate_bitLength = (int) Math.ceil(Math.log(counter) / MICalculator.logtwo);
			discPhenotype_discCovariate_bitMax = counter-1;
			discPhenotype_discCovariate_JointEntropyNats = MICalculator.calcEntropyInNats(discPhenotype_discCovariate_bitValues, counter);
		}
		hasDiscCovariate = true;
	}
	
	@Override
	public String toString() {
		return "Phenotype [id=" + id + "]";
	}
}