package midesp.objects;

public class EntropyCache {

	private final double[] terms;
	
	public EntropyCache(int sampleCount) {
		terms = new double[sampleCount + 1];
		
		for(int x = 1; x <= sampleCount; x++) {
			double p = (double) x / sampleCount;
			terms[x] = p * Math.log(p);
		}
	}
	
	public double get(int frequency) {
		return terms[frequency];
	}
}
