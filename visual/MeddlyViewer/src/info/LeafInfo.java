package info;

import java.util.PriorityQueue;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * 
 * @author Coleman
 *
 */
public class LeafInfo {
	private LinkedBlockingQueue<Integer> level;
	private LinkedBlockingQueue<Integer> anc;

	/**
	 * 
	 */
	public LeafInfo() {
		this.level = new LinkedBlockingQueue<Integer>();
		this.anc = new LinkedBlockingQueue<Integer>();
	}

	/**
	 * 
	 * @param level
	 */
	public void addLevel(int level) {
		this.level.add(level);
	}

	/**
	 * 
	 * @param anc
	 */
	public void addAnc(int anc) {
		this.anc.add(anc);
	}

	/**
	 * 
	 * @return
	 */
	public int getLevel() {
		return this.level.poll();
	}

	/**
	 * 
	 * @return
	 */
	public int getAnc() {
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
