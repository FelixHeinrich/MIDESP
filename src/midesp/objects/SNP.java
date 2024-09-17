package midesp.objects;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import midesp.methods.MICalculator;

public class SNP {

	private String id;
	private int length;
	private int[] genotypesCounts;
	private int[] genotypesArray;
	private int bitLength;
	private int bitMax;
	private double entropyNats;
	private double miToPheno;
	private double averageMiToPheno;
	private double pvalue;
	
	public SNP(String id, int length) {
		this.id = id;
		this.length = length;
		genotypesArray = new int[length];
	}

	public String getID() {
		return id;
	}
	
	public int getLength() {
		return length;
	}
	
	public int[] getGenotypes() {
		return genotypesArray;
	}
	
	public void setGenotypeAt(int idx, int genotype) {
		genotypesArray[idx] = genotype;
	}
	
	public int[] getGenotypesCounts() {
		return genotypesCounts;
	}
	
	public int getBitLength() {
		return bitLength;
	}
	
	public int getBitMax() {
		return bitMax;
	}
	
	public double getEntropyNats() {
		return entropyNats;
	}
	
	public double getEntropyLog2() {
		return entropyNats / MICalculator.logtwo;
	}
	
	public double getMItoPheno() {
		return miToPheno;
	}
	
	public void setMItoPheno(double mi) {
		miToPheno = mi;
	}
	
	public double getAverageMItoPheno() {
		return averageMiToPheno;
	}
	
	public void setAverageMItoPheno(double mi) {
		averageMiToPheno = mi;
	}
	
	public double getPValue() {
		return pvalue;
	}
	
	public void setPValue(double p) {
		pvalue = p;
	}
	
	public void parseValues(byte counter) {
		bitLength = (int) Math.ceil(Math.log(counter) / MICalculator.logtwo);
		bitMax = counter-1;
		genotypesCounts = new int[counter];
		for(int i = 0; i < length; i++) {
			genotypesCounts[genotypesArray[i]]++;
		}
		entropyNats = MICalculator.calcEntropyInNatsFromFreqs(genotypesCounts, length);
	}
	
	public static List<SNP> readTPed(Path tpedFile) throws IOException{
		return Files.lines(tpedFile).parallel().map(line ->{
			String[] tmpArr = line.split(" ");
			SNP tmpSNP = new SNP(tmpArr[1], (tmpArr.length-4)/2);
			Map<String,Byte> gtMap = new HashMap<>();
			byte counter = 0;
			for(int i = 0; i < tmpArr.length-4; i+=2) {
				String value = tmpArr[i+4]+tmpArr[i+5];
				if(!gtMap.containsKey(value)) {
					gtMap.put(value,counter);
					counter++;
				}
				tmpSNP.setGenotypeAt(i / 2, gtMap.get(value));
			}
			tmpSNP.parseValues(counter);
			return tmpSNP;
		}).collect(Collectors.toList());
	}
	
	@Override
	public String toString() {
		return "SNP [id=" + id + "]";
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(genotypesArray);
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SNP other = (SNP) obj;
		if (!Arrays.equals(genotypesArray, other.genotypesArray))
			return false;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}
}