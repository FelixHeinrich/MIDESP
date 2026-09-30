package midesp.objects;

public class TopKHeap_MIAPC {
	private final int[] snp1;
    private final int[] snp2;
    private final double[] mi;
    private final double[] miApc;

    private int size;

    public TopKHeap_MIAPC(int capacity) {
        this.snp1 = new int[capacity];
        this.snp2 = new int[capacity];
        this.mi = new double[capacity];
        this.miApc = new double[capacity];
    }

    public int size() {
        return size;
    }

    public boolean isFull() {
        return size == mi.length;
    }

    public double min() {
        return size == 0 ? Double.NEGATIVE_INFINITY : miApc[0];
    }

    public void offer(int snp1Index, int snp2Index, double mi_score, double miApc_score) {
    	if (size < mi.length) {
            snp1[size] = snp1Index;
            snp2[size] = snp2Index;
            mi[size] = mi_score;
            miApc[size] = miApc_score;

            siftUp(size);
            size++;
            return;
        }

    	if (isBetter(mi_score, miApc_score, mi[0], miApc[0])) {
    		snp1[0] = snp1Index;
            snp2[0] = snp2Index;
            mi[0] = mi_score;
            miApc[0] = miApc_score;

            siftDown(0);
        } 
    }
    
    private boolean isBetter(double miA, double miApcA, double miB, double miApcB) {
        if (miApcA > miApcB) {
            return true;
        }
        if (miApcA < miApcB) {
            return false;
        }
        return miA > miB;
    }

    
    private boolean isBetter(int a, int b) {        
        if (miApc[a] > miApc[b]) return true;
        if (miApc[a] < miApc[b]) return false;
        
        return mi[a] > mi[b];
    }

    private void siftUp(int child) {
        while (child > 0) {
            int parent = (child - 1) >>> 1;

            if (!isBetter(parent, child)) {
                break;
            }

            swap(parent, child);
            child = parent;
        }
    }

    private void siftDown(int parent) {
        int half = size >>> 1;

        while (parent < half) {
            int left = (parent << 1) + 1;
            int right = left + 1;

            int smallest = left;

            if (right < size && isBetter(left, right)) {
                smallest = right;
            }

            if (!isBetter(parent, smallest)) {
                break;
            }

            swap(parent, smallest);
            parent = smallest;
        }
    }

    private void swap(int a, int b) {
        int tmpSnp = snp1[a];
        snp1[a] = snp1[b];
        snp1[b] = tmpSnp;

        tmpSnp = snp2[a];
        snp2[a] = snp2[b];
        snp2[b] = tmpSnp;

        double doubleTemp = mi[a];
        mi[a] = mi[b];
        mi[b] = doubleTemp;
        
        doubleTemp = miApc[a];
        miApc[a] = miApc[b];
        miApc[b] = doubleTemp;
    }
    
    public interface ResultConsumer {
        void accept(int snp1Index, int snp2Index, double mi, double mi_apc);
    }

    public void forEach(ResultConsumer consumer) {
        for (int i = 0; i < size; i++) {
            consumer.accept(snp1[i], snp2[i], mi[i], miApc[i]);
        }
    }
    
    public void drain(ResultConsumer consumer) {

        while (size > 0) {

            int currentSnp1 = snp1[0];
            int currentSnp2 = snp2[0];
            double currentMi = mi[0];
            double currentMiApc = miApc[0];

            size--;

            if (size > 0) {
                snp1[0] = snp1[size];
                snp2[0] = snp2[size];
                mi[0] = mi[size];
                miApc[0] = miApc[size];

                siftDown(0);
            }

            consumer.accept(currentSnp1, currentSnp2, currentMi, currentMiApc);
        }
    }
}
