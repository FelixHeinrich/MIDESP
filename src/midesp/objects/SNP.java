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
	
	private GeneralizedBitSet bitSet;
	private GeneralizedBitSet snpDiscCovariateBitSet;
	private GeneralizedBitSet snpDiscPhenoDiscCovariateBitSet;

	private int[] genotypesCounts;
	private int[] genotypesValues;
	private double entropyNats;
	private double discCovariate_JointEntropyNats;
	private double discPhenoDiscCovariate_JointEntropyNats;
	private double miToPheno;
	private double averageMiToPheno;
	private double pvalue;
	
	public SNP(String id, int length) {
		this.id = id;
		this.length = length;
	}

	public String getID() {
		return id;
	}
	
	public int getLength() {
		return length;
	}
	
	public GeneralizedBitSet getBitSet() {
        return bitSet;
    }

	public GeneralizedBitSet getSNPDiscCovariateBitSet() {
		return snpDiscCovariateBitSet;
	}

	public void setSNPDiscCovariateBitSet(GeneralizedBitSet bitSet) {
		snpDiscCovariateBitSet = bitSet;
		discCovariate_JointEntropyNats = MICalculator.calcEntropyInNatsFromFreqs(bitSet.getClassCounts(), length);
	}
	
	public GeneralizedBitSet getSNPDiscPhenoDiscCovariateBitSet() {
		return snpDiscPhenoDiscCovariateBitSet;
	}

	public void setSNPDiscPhenoDiscCovariateBitSet(GeneralizedBitSet bitSet) {
		snpDiscPhenoDiscCovariateBitSet = bitSet;
		discPhenoDiscCovariate_JointEntropyNats = MICalculator.calcEntropyInNatsFromFreqs(bitSet.getClassCounts(), length);
	}
	
	public int[] getGenotypesCounts() {
		return genotypesCounts;
	}
	
	public int[] getGenotypesValues() {
		return genotypesValues;
	}
	
	public double getEntropyNats() {
		return entropyNats;
	}
	
	public double getDiscCovariateJointEntropyNats() {
		return discCovariate_JointEntropyNats;
	}
	
	public double getDiscPhenoDiscCovariateJointEntropyNats() {
		return discPhenoDiscCovariate_JointEntropyNats;
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
	
	public void initBitSet(int[] rawGenotypes, int numClasses) {
		this.bitSet = new GeneralizedBitSet(rawGenotypes, numClasses);
        this.genotypesCounts = this.bitSet.getClassCounts();
        this.genotypesValues = this.bitSet.getClassValues(length);
		entropyNats = MICalculator.calcEntropyInNatsFromFreqs(genotypesCounts, length);
	}
	
	public static List<SNP> readTPed(Path tpedFile) throws IOException{
		try(Stream<String> lines = Files.lines(tpedFile)) {
			return lines.parallel().map(line ->{
				String[] tmpArr = line.split(" ");
				int sampleCount = (tmpArr.length - 4) / 2;
				SNP tmpSNP = new SNP(tmpArr[1], sampleCount);
				int[] rawGenotypes = new int[sampleCount];
				Map<String,Byte> gtMap = new HashMap<>();
				byte counter = 0;
				for(int i = 0; i < tmpArr.length-4; i+=2) {
					String value = tmpArr[i+4]+tmpArr[i+5];
					Byte mappedValue = gtMap.get(value);
					if(mappedValue == null) {
						mappedValue = counter++;
						gtMap.put(value,mappedValue);
					}
					rawGenotypes[i / 2] = mappedValue;
				}
				tmpSNP.initBitSet(rawGenotypes, counter);
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