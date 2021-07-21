package midesp.objects;

import java.util.Objects;

public class Pair<E,F> {
	private E first;
	private F second;
	
	public Pair(E first, F second){
		this.first = first;
		this.second = second;
	}
	
	public Pair() {
		this.first = null;
		this.second = null;
	}
	
	public E getFirst() {
		return first;
	}
	
	public void setFirst(E first) {
		this.first = first;
	}
	
	public F getSecond() {
		return second;
	}
	
	public void setSecond(F second) {
		this.second = second;
	}
	
	public String toString(){
		return first.toString() + "\t" + second.toString();
	}
	
    @Override
    public boolean equals(Object o){
    	if(o instanceof Pair<?,?>){
    		Pair<?,?> pair2 = (Pair<?,?>) o;
    		if(this.getFirst().equals(pair2.getFirst()) && this.getSecond().equals(pair2.getSecond())){
    			return true;
    		}
    		else{
    			return false;
    		} 
    	}
    	else{
    		return false;
    	}
    }
    
    public boolean equals(Pair<E,F> pair2){
    	if(this.getFirst().equals(pair2.getFirst()) && this.getSecond().equals(pair2.getSecond())){
    		return true;
    	}
    	else{
    		return false;
    	} 
    }
    
    public int hashCode(){
    	return Objects.hash(this.first,this.second);
    }
}