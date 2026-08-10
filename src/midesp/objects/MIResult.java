package midesp.objects;

public record MIResult(String snp1, String snp2, double mi, double mi_apc) implements Comparable<MIResult>{
	
	public MIResult(String snp1, String snp2, double mi) {
		this(snp1, snp2, mi, 0.0);
	}

	@Override
	public int compareTo(MIResult o) {
		int cmp = Double.compare(this.mi_apc, o.mi_apc);
		if(cmp != 0) {
			return cmp;
		}
		return Double.compare(this.mi, o.mi);
	}
	
	public String toNoAPCString() {
		return snp1 + " + " + snp2 + " " + mi;
	}
	
	@Override
	public String toString() {
		return snp1 + " + " + snp2 + " " + mi + " " + mi_apc;
	}
}
