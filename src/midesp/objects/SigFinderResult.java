package midesp.objects;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SigFinderResult {

	private double lambda_1;
	private double lambda_2;
	private double tau;
	
	private List<Double> pvaluesList;
	
	private List<Integer> sigIdxList, backgroundIdxList, noiseIdxList;
	
	private List<Integer> zeroToLambda1IdxList, tauToLambda1IdxList;
	
	public double getLambda_1() {
		return lambda_1;
	}
	
	public double getLambda_2() {
		return lambda_2;
	}
	
	public double getTau() {
		return tau;
	}
	
	public List<Double> getPValues(){
		return pvaluesList;
	}
	
	public List<Integer> getSignificantIndices(){
		return sigIdxList;
	}
	
	public List<Integer> getBackgroundIndices(){
		return backgroundIdxList;
	}
	
	public List<Integer> getNoiseIndices(){
		return noiseIdxList;
	}
	
	public List<Integer> getZeroToLambda1Indices(){
		return zeroToLambda1IdxList;
	}
	
	public List<Integer> getTauToLambda1Indices(){
		return tauToLambda1IdxList;
	}
	
	public SigFinderResult(List<Double> pvalues, double lambda_1, double lambda_2, double tau) {
		this.pvaluesList = pvalues.stream().collect(Collectors.toList());
		this.lambda_1 = lambda_1;
		this.lambda_2 = lambda_2;
		this.tau = tau;
		
		sigIdxList = IntStream.range(0, pvaluesList.size()).filter(idx -> pvaluesList.get(idx) <= tau).boxed().collect(Collectors.toList());
		backgroundIdxList = IntStream.range(0, pvaluesList.size()).filter(idx -> pvaluesList.get(idx) >= this.lambda_1 && pvaluesList.get(idx) <= this.lambda_2).boxed().collect(Collectors.toList());
		noiseIdxList = IntStream.range(0, pvaluesList.size()).filter(idx -> pvaluesList.get(idx) > this.lambda_2).boxed().collect(Collectors.toList());
		
		zeroToLambda1IdxList = IntStream.range(0, pvaluesList.size()).filter(idx -> pvaluesList.get(idx) <= this.lambda_1).boxed().collect(Collectors.toList());
		tauToLambda1IdxList = IntStream.range(0, pvaluesList.size()).filter(idx -> pvaluesList.get(idx) <= this.lambda_1 && pvaluesList.get(idx) >= tau).boxed().collect(Collectors.toList());
	}
}