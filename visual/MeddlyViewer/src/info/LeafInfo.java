package info;

import java.util.PriorityQueue;

/**
 * 
 * @author Coleman
 *
 */
public class LeafInfo {
	private PriorityQueue<Long> level;
	private PriorityQueue<Long> anc;

	/**
	 * 
	 */
	public LeafInfo() {
		this.level = new PriorityQueue<Long>();
		this.anc = new PriorityQueue<Long>();
	}

	/**
	 * 
	 * @param level
	 */
	public void addLevel(long level) {
		this.level.add(level);
	}

	/**
	 * 
	 * @param anc
	 */
	public void addAnc(long anc) {
		this.anc.add(anc);
	}

	/**
	 * 
	 * @return
	 */
	public Long getLevel() {
		return this.level.poll();
	}

	/**
	 * 
	 * @return
	 */
	public Long getAnc() {
		return this.anc.poll();
	}

	/**
	 * 
	 * @return
	 */
	public boolean hasNext() {
		if (level.peek() == null || anc.peek() == null)
			return false;

		return true;
	}

}
