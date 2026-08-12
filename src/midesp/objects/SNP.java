package midesp.objects;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import midesp.methods.MICalculator;

public class SNP {

	private final String id;
	private final int length;
	private final int[] genotypesArray;
	private int[] genotypesCounts;
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
		try(Stream<String> lines = Files.lines(tpedFile)) {
			return lines.parallel().map(line ->{
				String[] tmpArr = line.split(" ");
				SNP tmpSNP = new SNP(tmpArr[1], (tmpArr.length-4)/2);
				Map<String,Byte> gtMap = new HashMap<>();
				byte counter = 0;
				for(int i = 0; i < tmpArr.length-4; i+=2) {
					String value = tmpArr[i+4]+tmpArr[i+5];
					Byte mappedValue = gtMap.get(value);
					if(mappedValue == null) {
						mappedValue = counter++;
						gtMap.put(value,mappedValue);
					}
					tmpSNP.setGenotypeAt(i / 2, mappedValue);
				}
				tmpSNP.parseValues(counter);
				return tmpSNP;
			}).collect(Collectors.toMap(
					SNP::getID,
					snp -> snp,
					(existing, duplicate) ->{
						throw new UncheckedIOException(new IOException("Duplicate SNP ID found in tped file: " + existing.getID()));
					},
					LinkedHashMap::new
			)).values().stream().toList();	
		} catch (UncheckedIOException e) {
	        throw e.getCause();
	    }
	}
	
	@Override
	public String toString() {
		return "SNP [id=" + id + "]";
	}
	
	@Override
	public int hashCode() {
		return id.hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!(obj instanceof SNP other))
			return false;
		return Objects.equals(id, other.id);
	}
}