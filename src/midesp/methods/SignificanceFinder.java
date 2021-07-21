package midesp.methods;

import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.BetaDistribution;

import midesp.objects.Pair;
import midesp.objects.SigFinderResult;


public class SignificanceFinder {

	
	public static SigFinderResult findSignificantScores(List<Double> valuesList, double fdr) {
		double mean = calcMean(valuesList);
		double variance = calcVariance(valuesList, mean);
		double alpha = estAlpha(mean, variance);
        double beta = estBeta(mean, variance);
        List<Double> pValuesList = calcPValues(valuesList, alpha, beta);
        if(performChiSquareTest(pValuesList)) {
        	Pair<Pair<Double,Double>,Double> gammaValues = getGamma(pValuesList);
        	Pair<Double,Double> boundaries = gammaValues.getFirst();
        	double gamma = gammaValues.getSecond();
        	double tau = getTau(pValuesList, fdr, boundaries.getFirst(), gamma);
        	if(tau <= 0) {
        		Pair<List<Double>,Pair<Double,Pair<Double,Double>>> result = varianceDecrease(mean, variance, fdr, valuesList);
        		pValuesList = result.getFirst();
        		tau = result.getSecond().getFirst();
        		boundaries = result.getSecond().getSecond();
        	}
        	else {
        		Pair<List<Double>,Pair<Double,Pair<Double,Double>>> result = varianceIncrease(mean, variance, fdr, valuesList, boundaries);
        		pValuesList = result.getFirst();
        		tau = result.getSecond().getFirst();
        		boundaries = result.getSecond().getSecond();
        	}
        	if(tau < 0) {
        		System.out.println("Failed to compute lambda_1 and lambda_2 with sufficient distance between each other.");
        		return null;
        	}
        	return new SigFinderResult(pValuesList, boundaries.getFirst(), boundaries.getSecond(), tau);
        }
        else {
        	System.out.println("ChiSquareTest reported a uniform distribution of pvalues");
        	return null;
        }
	}
	
	private static double calcMean(List<Double> valuesList) {
		return valuesList.parallelStream().mapToDouble(Double::doubleValue).average().getAsDouble();
	}
	
	private static double calcVariance(List<Double> valuesList, double mean) {
		double tempSum = valuesList.parallelStream().mapToDouble(value -> Math.pow(value - mean, 2)).sum();
		return tempSum / (valuesList.size() - 1); 
	}
	
	private static double estAlpha(double mean, double variance) {
		return mean * (((mean * (1 - mean)) / variance) - 1);
	}
	
	private static double estBeta(double mean, double variance) {
		return (1 - mean) * (((mean * (1 - mean)) / variance) - 1);
	}
	
	private static List<Double> calcPValues(List<Double> valuesList, double alpha, double beta) {
		BetaDistribution betaDist = new BetaDistribution(alpha, beta);
		return valuesList.parallelStream().map(value -> 1.0 - betaDist.cumulativeProbability(value)).collect(Collectors.toList());
	}
	
	private static Pair<Pair<Double,Double>,Double> getGamma(List<Double> valuesList) {
		Pair<Double,Double> boundaries = estBoundaries(valuesList);
		if(boundaries.getSecond() - boundaries.getFirst() < 0.18) {
			return null;
		}
		long counter = valuesList.parallelStream().filter(value -> value >= boundaries.getFirst() && value <= boundaries.getSecond()).count();
		double gamma = counter / (valuesList.size() * (boundaries.getSecond() - boundaries.getFirst()));
		return new Pair<>(boundaries, gamma);
	}
	
	private static Pair<Double,Double> estBoundaries(List<Double> valuesList ){
		double lambda_1 = 0.2;
		double lambda_2 = 0.8;
		double gamma = 0.0;
		double last_gamma = 0.0;
		boolean accept = false;
		while(!accept && lambda_1 < lambda_2) {
			double temp_lambda_1 = lambda_1;
			double temp_lambda_2 = lambda_2;
			long counter = valuesList.parallelStream().filter(value -> value >= temp_lambda_1 && value <= temp_lambda_2).count();
			gamma = counter / (valuesList.size() * (lambda_2 - lambda_1));
			double variance = last_gamma - gamma;
			if(variance < 0.01 && variance > -0.01) {
				lambda_1 -= 0.01;
				lambda_2 += 0.01;
				accept = true;
			}
			else {
				last_gamma = gamma;
				lambda_1 += 0.01;
				lambda_2 -= 0.01;
			}
		}
		
		return new Pair<>(lambda_1, lambda_2);
	}
	
	private static double getTau(List<Double> valuesList, double fdr, double lambda_1, double gamma) {
		if(gamma < 0) {
			return -1;
		}
		double tau = 0.0;
		double tests = 10000;
		double diff = 0.0001;
		for(int i = 0; i < tests; i++) {
			tau += 1 / tests;
			double currentFDR = 0; 
			if(tau < lambda_1) {
				double tempTau = tau;
				long counter = valuesList.parallelStream().filter(value -> value < tempTau).count();
				currentFDR = gamma * valuesList.size() * (tau / counter);
				if(Math.abs(currentFDR - fdr) < diff) {
					return tau;
				}
			}
		}
		return -1;
	}
	
	private static boolean performChiSquareTest(List<Double> valuesList) {
		int[] interval = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (int i = 0; i < valuesList.size(); ++i) {
            double p = valuesList.get(i);
            for (int j = 1; j < 21; ++j) {
                if (p == 1) {
                    interval[19]++;
                    break;
                } else if (((j - 1) * 0.05) <= p && p < (j * 0.05)) {
                    interval[j - 1]++;
                    break;
                }
            }
        }
        double n_0 = (1 / (double) interval.length) * valuesList.size();
        double chiSquare = 0.0;
        for (int i = 0; i < interval.length; ++i) {
            chiSquare += (Math.pow((interval[i] - n_0), 2) / n_0);
        }
        if (chiSquare < 45.0) {//38.58) {
            return false;
        }
        return true;
	}
	
	private static Pair<List<Double>,Pair<Double,Pair<Double,Double>>> varianceDecrease(double mean, double variance, double fdr, List<Double> valuesList) {
		while(true) {
			variance -= variance * 0.1;
			double alpha = estAlpha(mean, variance);
			double beta = estBeta(mean, variance);
			List<Double> pValuesList = calcPValues(valuesList, alpha, beta);
			if(performChiSquareTest(pValuesList)) {
				Pair<Pair<Double,Double>,Double> gammaValues = getGamma(pValuesList);
	        	Pair<Double,Double> boundaries = gammaValues.getFirst();
	        	double gamma = gammaValues.getSecond();
	        	double tau = getTau(pValuesList, fdr, boundaries.getFirst(), gamma);
	        	if(tau > 0) {
	        		return new Pair<>(pValuesList,new Pair<>(tau, boundaries));
	        	}
			}
		}
	}
	
	private static Pair<List<Double>,Pair<Double,Pair<Double,Double>>> varianceIncrease(double mean, double variance, double fdr, List<Double> valuesList, Pair<Double, Double> boundariesOrg) {
		boolean alphaCheck = true;
		double tau = 0;
		double lastTau = 0;
		while(true) {
			variance += variance * 0.01;
			double alpha = estAlpha(mean, variance);
			if(alpha < 40) {
				variance -= variance * 0.01;
				alpha = estAlpha(mean, variance);
				alphaCheck = false;
			}
			double beta = estBeta(mean, variance);
			List<Double> pValuesList = calcPValues(valuesList, alpha, beta);
			if(performChiSquareTest(pValuesList)) {
				Pair<Pair<Double,Double>,Double> gammaValues = getGamma(pValuesList);
	        	Pair<Double,Double> boundaries = gammaValues.getFirst();
	        	double gamma = gammaValues.getSecond();
	        	if(!checkGamma(boundaries, boundariesOrg)) {
	        		variance -= variance * 0.01;
	        		alpha = estAlpha(mean, variance);
	        		beta = estBeta(mean, variance);
	        		pValuesList = calcPValues(valuesList, alpha, beta);
	        		gammaValues = getGamma(pValuesList);
	        		boundaries = gammaValues.getFirst();
		        	gamma = gammaValues.getSecond();
	        		tau = getTau(pValuesList, fdr, boundaries.getFirst(), gamma);
	        	}
	        	tau = getTau(pValuesList, fdr, boundaries.getFirst(), gamma);
	        	if(tau <= 0) {
	        		return new Pair<>(pValuesList,new Pair<>(lastTau, boundaries));
	        	}
			}
			if(!alphaCheck) {
        	    Pair<Pair<Double,Double>,Double> gammaValues = getGamma(pValuesList);
        	    Pair<Double,Double> boundaries = gammaValues.getFirst();
				return new Pair<>(pValuesList,new Pair<>(tau, boundaries));
			}
			lastTau = tau;
		}
	}
	
	private static boolean checkGamma(Pair<Double,Double> boundariesOrg, Pair<Double,Double> boundariesNew) {
		double epsilon = 0.05;
		if(Math.abs(boundariesOrg.getFirst() - boundariesNew.getFirst()) < epsilon && Math.abs(boundariesOrg.getSecond() - boundariesNew.getSecond()) < epsilon) {
			return true;
		}			
		return false;
	}
}