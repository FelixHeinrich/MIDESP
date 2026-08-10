package midesp.objects;

import java.util.Collection;
import java.util.Objects;
import java.util.PriorityQueue;

public class LimitedPriorityQueue<E extends Comparable<E>> extends PriorityQueue<E> {

	private static final long serialVersionUID = 1L;
	
	private final int limit;
	
	public LimitedPriorityQueue(int limit) {
		super(Math.max(1, limit));
		if(limit <= 0) {
			throw new IllegalArgumentException("Limit must be greater than 0");
		}
		this.limit = limit;
	}
	
	public int getLimit() {
		return limit;
	}
	
	@Override
	public boolean offer(E element) {
		Objects.requireNonNull(element, "Null elements are not permitted");
		
		if(size() >= limit) {
			E head = peek();
			if(head.compareTo(element) < 0) {
				poll();
			}
			else {
				return false;
			}
		}
		
		return super.offer(element);
	}
	
	@Override
	public boolean add(E element) {
		return offer(element);
	}
	
	@Override
	public boolean addAll(Collection<? extends E> col) {
		boolean modified = false;
		for(E element : col) {
			if(offer(element)) {
				modified = true;
			}
		}
		return modified;
	}
}