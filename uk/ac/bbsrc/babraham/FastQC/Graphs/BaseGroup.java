/**
 * Copyright Copyright 2010-11 Simon Andrews
 *
 *    This file is part of FastQC.
 *
 *    FastQC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    FastQC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with FastQC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
package uk.ac.bbsrc.babraham.FastQC.Graphs;

import java.util.Vector;

/**
 * The base group class is a simple way to create a set of bins into
 * which positions within a read can be put such that early positions
 * get a group to themselves and later positions get averaged so that
 * general trends can still be observed.
 *
 */
public class BaseGroup {

	private int lowerCount;
	private int upperCount;
	
	public static BaseGroup [] makeBaseGroups (int maxLength) {
		
		int startingBase = 1;
		int interval = 1;
		
		// The user can set a parameter which forces us to not group
		// bases together, so check if this is present before starting
		boolean noIncrement = false;
		
		if (System.getProperty("fastqc.nogroup") != null && System.getProperty("fastqc.nogroup").equals("true")) {
			noIncrement = true;
		}

		
		
		
		Vector<BaseGroup> groups = new Vector<BaseGroup>();
		
		while (startingBase <= maxLength) {
			
			int endBase = startingBase+(interval-1);
			if (endBase > maxLength) endBase = maxLength;
			
			BaseGroup bg = new BaseGroup(startingBase, endBase);
			groups.add(bg);
			
			startingBase += interval;
						
			// See if we need to increase the interval
			if (!noIncrement) {
				if (startingBase == 10 && maxLength > 75) {
					interval = 5;
				}
				if (startingBase == 50 && maxLength > 200) {
					interval = 10;
				}
				if (startingBase == 100 && maxLength > 300) {
					interval = 50;
				}
				if (startingBase == 500 && maxLength > 1000) {
					interval = 100;
				}
				if (startingBase == 1000 && maxLength > 2000) {
					interval = 500;
				}
			}
			
		}
		
		return groups.toArray(new BaseGroup[0]);
		
	}
	
	
	/**
	 * 
	 * @param lowerCount
	 * @param upperCount
	 */
	
	private BaseGroup (int lowerCount, int upperCount) {
		this.lowerCount = lowerCount;
		this.upperCount = upperCount;
	}
	
	public int lowerCount () {
		return lowerCount;
	}
	
	public int upperCount () {
		return upperCount;
	}
	
	public boolean containsValue (int value) {
		return value>=lowerCount && value<=upperCount;
	}
	
	public String toString () {
		if (lowerCount == upperCount) {
			return ""+lowerCount;
		}
		else {
			return ""+lowerCount+"-"+upperCount;
		}
	}
	
}
