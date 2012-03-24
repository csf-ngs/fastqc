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
package uk.ac.bbsrc.babraham.FastQC.Modules;

import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import uk.ac.bbsrc.babraham.FastQC.Graphs.LineGraph;
import uk.ac.bbsrc.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.bbsrc.babraham.FastQC.Sequence.Sequence;

public class DuplicationLevel implements QCModule {

	private OverRepresentedSeqs overrepresentedModule;
	private double [] counts = null;
	private double maxCount = 0;
	private double duplicatePercentTotal = 0;
	private String [] labels;
	private static final DecimalFormat df = new DecimalFormat("#.##");
	
	protected DuplicationLevel (OverRepresentedSeqs overrepresentedModule) {
		this.overrepresentedModule = overrepresentedModule;
	}
	
	public String description() {
		return "Plots the number of sequences which are duplicated to different levels";
	}

	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	protected synchronized void calculateLevels () {
		
		if (counts != null) return;
		
		counts = new double[10];
		int [] tempCounts = new int [10];
		labels = new String [10];
		maxCount = 0;
				
		Iterator<String> it = overrepresentedModule.sequences.keySet().iterator();
				
		while (it.hasNext()) {
			int thisCount = overrepresentedModule.sequences.get(it.next());
			if (thisCount > 10) thisCount = 10;
			thisCount -= 1; // Convert count to array index
			++tempCounts[thisCount];
		}
		
		// To prevent infinite values later we need to ensure we have a count
		// of at least 1 in the unique sequences
		if (tempCounts[0] == 0) tempCounts[0] = 1;
	
		for (int i=0;i<tempCounts.length;i++) {
			counts[i] = (tempCounts[i]*100d)/tempCounts[0];
			if (counts[i] > maxCount) maxCount = counts[i];
			labels[i] = ""+(i+1);
			if (i == 9) {
				labels[i] = labels[i]+"+";
			}
		}

		// We work out the duplicate total percentage by
		// extrapolating the unique count forward to the
		// total number of sequences in the file.
				
		// First calculate the proportion of sequences we actually measured
		duplicatePercentTotal = ((double)overrepresentedModule.count) / Math.min(overrepresentedModule.countAtUniqueLimit,overrepresentedModule.count);
		
		// Now we extrapolate our unique count to cover the sequences
		// we didn't see (assuming they come at the same rate
		duplicatePercentTotal *= tempCounts[0];
		
//		System.err.println("Extrapolated unique count from "+tempCounts[0]+" to "+duplicatePercentTotal);
		
		// Now we subtract that from the total count to get the total non-unique count
		duplicatePercentTotal = overrepresentedModule.count - duplicatePercentTotal;
		
//		System.out.println("Non unique count is "+duplicatePercentTotal);
		
		// Now we express this as a percentage
		duplicatePercentTotal *= 100;
		duplicatePercentTotal /= overrepresentedModule.count;
		
//		System.out.println("Non unique percentage is "+duplicatePercentTotal);
		
		// In corner cases we could get silly answers
		if (duplicatePercentTotal>100) duplicatePercentTotal = 100;
		if (duplicatePercentTotal<0) duplicatePercentTotal=0;

	}

	public JPanel getResultsPanel() {
		if (counts == null) calculateLevels();

		return new LineGraph(new double [][] {counts}, 0d, maxCount, "Sequence Duplication Level",new String [] {"%Duplicate relative to unique"}, labels, "Sequence Duplication Level >= "+df.format(duplicatePercentTotal)+"%");
	}
	
	public void makeReport(HTMLReportArchive report) throws IOException {
		if (counts == null) calculateLevels();

		ZipOutputStream zip = report.zipFile();
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/duplication_levels.png"));
		
		BufferedImage b = new BufferedImage(800,600,BufferedImage.TYPE_INT_RGB);
		Graphics g = b.getGraphics();
		
		LineGraph lg = new LineGraph(new double [][] {counts}, 0d, maxCount,"Sequence Duplication Level", new String [] {"%Duplicate relative to unique"}, labels, "Sequence Duplication Level >= "+df.format(duplicatePercentTotal)+"%");
		lg.paint(g,800,600);

		ImageIO.write((BufferedImage)(b),"PNG",zip);

		StringBuffer sb = report.htmlDocument();
		
		sb.append("<p><img class=\"indented\" src=\"Images/duplication_levels.png\" alt=\"Duplication level graph\"></p>\n");
		
		sb = report.dataDocument();
		
		sb.append("#Total Duplicate Percentage\t");
		sb.append(duplicatePercentTotal);
		sb.append("\n");
		
		sb.append("#Duplication Level\tRelative count\n");
		for (int i=0;i<labels.length;i++) {
			sb.append(labels[i]);
			if (i == labels.length-1) {
				sb.append("+");
			}
			sb.append("\t");
			sb.append(counts[i]);
			sb.append("\n");
		}
		
	}

	public String name() {
		return "Sequence Duplication Levels";
	}

	public void processSequence(Sequence sequence) {
		// We don't need to do anything since we use 
		// the data structure from the overrepresented sequences
		// module.
	}

	public boolean raisesError() {
		if (counts == null) calculateLevels();
		
		// Anything over 50% duplicate gets us a error
		if (duplicatePercentTotal > 50) {
			return true;
		}
		
		return false;
	}

	public boolean raisesWarning() {
		if (counts == null) calculateLevels();

		// Anything over 20% duplicate gets us a warning
		if (duplicatePercentTotal > 20) {
			return true;
		}
		
		return false;
	}

	public void reset() {
		counts = null;
	}

}
