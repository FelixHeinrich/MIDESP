package midesp.objects;

public class TopKHeap_MI {
	private final int[] snp1;
    private final int[] snp2;
    private final double[] mi;

    private int size;

    public TopKHeap_MI(int capacity) {
        this.snp1 = new int[capacity];
        this.snp2 = new int[capacity];
        this.mi = new double[capacity];    
    }

    public int size() {
        return size;
    }

    public boolean isFull() {
        return size == mi.length;
    }

    public double min() {
        return size == 0 ? Double.NEGATIVE_INFINITY : mi[0];
    }

    public void offer(int snp1Index, int snp2Index, double score) {
    	if (size < mi.length) {
            snp1[size] = snp1Index;
            snp2[size] = snp2Index;
            mi[size] = score;

            siftUp(size);
            size++;
            return;
        }

        if (score <= mi[0]) {
        	return;
        }
        snp1[0] = snp1Index;
        snp2[0] = snp2Index;
        mi[0] = score;

        siftDown(0);
    }

    private void siftUp(int child) {
        while (child > 0) {
            int parent = (child - 1) >>> 1;

            if (mi[parent] <= mi[child]) {
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

            if (right < size && mi[right] < mi[left]) {
                smallest = right;
            }

            if (mi[parent] <= mi[smallest]) {
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

        double tmpMi = mi[a];
        mi[a] = mi[b];
        mi[b] = tmpMi;
    }
    
    public interface ResultConsumer {
        void accept(int snp1Index, int snp2Index, double mi);
    }

    public void forEach(ResultConsumer consumer) {
        for (int i = 0; i < size; i++) {
            consumer.accept(snp1[i], snp2[i], mi[i]);
        }
    }
    
    public void drain(ResultConsumer consumer) {

        while (size > 0) {

            int currentSnp1 = snp1[0];
            int currentSnp2 = snp2[0];
            double currentMi = mi[0];

            size--;

            if (size > 0) {
                snp1[0] = snp1[size];
                snp2[0] = snp2[size];
                mi[0] = mi[size];

                siftDown(0);
            }

            consumer.accept(currentSnp1, currentSnp2, currentMi);
        }
    }
}
