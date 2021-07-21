package midesp.objects;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import midesp.methods.MICalculator;


public class Phenotype{

	private String id;
	private int length;
	private boolean isContinuous;
	private int[] discPhenotypeVec;
	private int[] bitValues;
	private int bitLength;
	private int bitMax;
	private double discPhenotypeEntropyNats;
	private double[] contPhenotypeVec;
	private double[] digammaValuesArray;
	private int[][] closestNeighborsMat;
	private double[][] closestNeighborsDistMat;
	
	public Phenotype(String id, int length, boolean continuous) {
		this.id = id;
		this.length = length;
		isContinuous = continuous;
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
	
	public void setValueAt(int idx, String value) {
		if(isContinuous) {
			contPhenotypeVec[idx] = Double.parseDouble(value);
		}
		else {
			discPhenotypeVec[idx] = Integer.parseInt(value);
		}
	}	
	
	public void parseValues() {	
		Random rng = new Random();
		if(isContinuous) {
			rng.setSeed(42);
			//Check for duplicate values and add noise to guarantee uniqueness of values
			while(Arrays.stream(contPhenotypeVec).distinct().count() != length) {
				//Adding small noise as suggested by Kraskov et al in "Estimating Mutual Information"
				for(int i = 0; i < length; i++) {
					contPhenotypeVec[i] += rng.nextGaussian() * 1e-10;
				}
			}
			//Normalize the phenotype to the interval [0;1]
			contPhenotypeVec = MICalculator.minMaxNormalization(contPhenotypeVec);
			//Create sorted lists of neighbors for each entry
			closestNeighborsMat = Arrays.stream(contPhenotypeVec).mapToObj(value ->{
				return IntStream.range(0, length).boxed().sorted((o1, o2) -> {
					if(Math.abs(value - contPhenotypeVec[o1]) < Math.abs(value - contPhenotypeVec[o2])) {
						return -1;
					}	
					else {
						return 1;
					}
				}).mapToInt(i -> i).toArray();
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
	
	@Override
	public String toString() {
		return "Phenotype [id=" + id + "]";
	}
}