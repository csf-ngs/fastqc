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
package uk.ac.bbsrc.babraham.FastQC.Analysis;

import java.util.LinkedList;

import uk.ac.bbsrc.babraham.FastQC.Modules.QCModule;
import uk.ac.bbsrc.babraham.FastQC.Sequence.SequenceFile;

public class AnalysisQueue implements Runnable, AnalysisListener{

	private static AnalysisQueue instance = new AnalysisQueue();
	
	private LinkedList<AnalysisRunner>queue = new LinkedList<AnalysisRunner>();
	
	private int availableSlots = 1;
	private int usedSlots = 0;
	
	public static AnalysisQueue getInstance () {
		return instance;
	}
	
	private AnalysisQueue () {
		
		if (System.getProperty("fastqc.threads") != null) {
			try {
				availableSlots = Integer.parseInt(System.getProperty("fastqc.threads"));
			}
			catch (NumberFormatException nfe) {
				System.err.println("Thead limit '"+System.getProperty("fastqc.threads")+" wasn't a number");
				availableSlots = 1;
			}
			
			if (availableSlots < 1) {
				System.err.println("Thead limit '"+System.getProperty("fastqc.threads")+" must be > 1");
				availableSlots = 1;
			}
			
		}
		
		Thread t = new Thread(this);
		t.start();
	}
	
	public void addToQueue (AnalysisRunner runner) {
		queue.add(runner);
	}

	public void run() {

		while (true) {
//			System.err.println("Status available="+availableSlots+" used="+usedSlots+" queue="+queue.size());
			if (availableSlots > usedSlots && queue.size() > 0) {
				++usedSlots;
				AnalysisRunner currentRun = queue.getFirst();
				queue.removeFirst();
				currentRun.addAnalysisListener(this);
				Thread t = new Thread(currentRun);
				t.start();
			}
			
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {}
		}
	}

	public void analysisComplete(SequenceFile file, QCModule[] results) {
		--usedSlots;
	}

	public void analysisUpdated(SequenceFile file, int sequencesProcessed, int percentComplete) {}

	public void analysisExceptionReceived(SequenceFile file, Exception e) {
		--usedSlots;
	}

	public void analysisStarted(SequenceFile file) {}
	
	
}
