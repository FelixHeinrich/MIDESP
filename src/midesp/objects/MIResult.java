package midesp.objects;

public class MIResult implements Comparable<MIResult>{

	private String snp1;
	private String snp2;
	private double mi;
	
	private double mi_apc;
	
	public MIResult(String snp1, String snp2, double mi, double mi_apc) {
		this.snp1 = snp1;
		this.snp2 = snp2;
		this.mi = mi;
		this.mi_apc = mi_apc;
	}
	
	public MIResult(String snp1, String snp2, double mi) {
		this(snp1, snp2, mi, 0.0);
	}
	
	public String getSNP1() {
		return snp1;
	}
	
	public String getSNP2() {
		return snp2;
	}
	
	public double getMI() {
		return mi;
	}
	
	public double getMI_apc() {
		return mi_apc;
	}
	
	@Override
	public String toString() {
		return snp1 + " + " + snp2 + " " + mi + " " + mi_apc;
	}
	
	public String toNoAPCString() {
		return snp1 + " + " + snp2 + " " + mi;
	}
	
	@Override
	public int compareTo(MIResult o) {
		if(this.mi_apc < o.mi_apc) {
			return -1;
		}
		else if(this.mi_apc == o.mi_apc) {
			if(this.mi < o.mi) {
				return -1;
			} else if(this.mi == o.mi) {
				return 0;
			}
			return 1;
		}
		return 1;
	}

}