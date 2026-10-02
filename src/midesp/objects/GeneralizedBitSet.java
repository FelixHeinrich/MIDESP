package midesp.objects;

import java.util.Arrays;

public class GeneralizedBitSet {
	
	private long[][] masks;
	private int[] classCounts;

	public GeneralizedBitSet(int[] values, int numClasses) {
		int sampleCount = values.length;
        int numWords = (sampleCount + 63) / 64;
        masks = new long[numClasses][numWords];
        classCounts = new int[numClasses];
        for (int i = 0; i < sampleCount; i++) {
            int val = values[i];
            masks[val][i >>> 6] |= (1L << (i & 63));
            classCounts[val]++;
        }
    }
    
	public GeneralizedBitSet(long[][] masks, int[] classCounts) {
	    this.masks = masks;
	    this.classCounts = classCounts;
	}
	
    public int getNumClasses() {
        return masks.length;
    }

    public int getNumWords() {
        return masks.length > 0 ? masks[0].length : 0;
    }

    public int[] getClassCounts() {
        return classCounts;
    }
    
    public long[] getMask(int classIdx) {
        return masks[classIdx];
    }
    
    public static GeneralizedBitSet combine(SNP[] snps) {
        GeneralizedBitSet combined = snps[0].getBitSet();
        for (int i = 1; i < snps.length; i++) {
            combined = combineTwo(combined, snps[i].getBitSet());
        }
        return combined;
    }
    
    public static GeneralizedBitSet combineTwo(GeneralizedBitSet b1, GeneralizedBitSet b2) {
        int c1Count = b1.getNumClasses();
        int c2Count = b2.getNumClasses();
        int maxClasses = c1Count * c2Count;
        int numWords = b1.getNumWords();

        int[] b1Counts = b1.getClassCounts();
        int[] b2Counts = b2.getClassCounts();

        // Pre-allocate upper bound
        long[][] combinedMasks = new long[maxClasses][numWords];
        int[] classCounts = new int[maxClasses];

        int rowIdx = 0;

        for (int c1 = 0; c1 < c1Count; c1++) {
            if (b1Counts[c1] == 0) continue; // Skip zero marginals
            long[] mask1 = b1.getMask(c1);

            for (int c2 = 0; c2 < c2Count; c2++) {
                if (b2Counts[c2] == 0) continue; // Skip zero marginals
                long[] mask2 = b2.getMask(c2);

                // Write directly into current slot
                long[] targetMask = combinedMasks[rowIdx];
                int count = 0;

                for (int w = 0; w < numWords; w++) {
                    long combinedBitmask = mask1[w] & mask2[w];
                    targetMask[w] = combinedBitmask;
                    count += Long.bitCount(combinedBitmask);
                }
                // Only advance rowIdx if this class actually exists in the data!
                // If count == 0, the next iteration overwrites combinedMasks[rowIdx].
                if (count > 0) {
                    classCounts[rowIdx] = count;
                    rowIdx++;
                }
            }
        }

        long[][] activeMasks = Arrays.copyOf(combinedMasks, rowIdx);
        int[] activeCounts = Arrays.copyOf(classCounts, rowIdx);

        return new GeneralizedBitSet(activeMasks, activeCounts);
    }
    
    public int[] getClassValues(int sampleCount) {
    	int[] values  = new int[sampleCount];

    	for (int c = 0; c < masks.length; c++) {
    		long[] mask = masks[c];

    		for (int w = 0; w < mask.length; w++) {
    			long word = mask[w];

    			while (word != 0L) {
    				int bit = Long.numberOfTrailingZeros(word);
    				int sampleIdx = (w << 6) + bit;

    				if (sampleIdx < sampleCount) {
    					values[sampleIdx] = c;
    				}

    				word &= word - 1;
    			}
    		}
    	}
    	return values;
    }
}
