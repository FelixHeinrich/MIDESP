package midesp.objects;

import java.util.Collection;
import java.util.PriorityQueue;

public class LimitedPriorityQueue<E extends Comparable<E>> extends PriorityQueue<E> {

	private static final long serialVersionUID = 1L;
	
	private final int limit;
	
	public LimitedPriorityQueue(int limit) {
		super();
		this.limit = limit;
	}
	
	public int getLimit() {
		return limit;
	}
	
	@Override
	public boolean add(E element) {
		if(this.size() == limit) {
			E head = this.peek();
			if(head.compareTo(element) < 0) {
				this.remove();
			}
			else {
				return false;
			}
		}
		return super.add(element);
	}
	
	@Override
	public boolean addAll(Collection<? extends E> col) {
		for(E element : col) {
			this.add(element);
		}
		return true;
	}
}