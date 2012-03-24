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
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import uk.ac.bbsrc.babraham.FastQC.Graphs.BaseGroup;
import uk.ac.bbsrc.babraham.FastQC.Graphs.LineGraph;
import uk.ac.bbsrc.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.bbsrc.babraham.FastQC.Sequence.Sequence;

public class PerBaseGCContent implements QCModule {

	public long [] gcCounts = new long [0];
	public long [] atCounts = new long [0];
	private double [] percentages = null;
	private String [] xCategories = new String[0];
	private boolean calculated = false;
	private double maxDeviation;
	
	
	public JPanel getResultsPanel() {
		
		if (!calculated) getPercentages();

		return new LineGraph(new double [][] {percentages}, 0d, 100d, "Position in read (bp)", new String [] {"%GC"}, xCategories, "GC content across all bases");
	}

	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	private synchronized void getPercentages () {

		BaseGroup [] groups = BaseGroup.makeBaseGroups(gcCounts.length);
		
		xCategories = new String[groups.length];

		percentages = new double [groups.length];

		
		// For error calling we also want to work out what the
		// maximum deviation from the mean GC is to see if we
		// have any great level of variation
		double mean = 0;

		long total;
		long gcCount;

		for (int i=0;i<groups.length;i++) {
						
			xCategories[i] = groups[i].toString();

			gcCount = 0;
			total = 0;
			
			for (int bp=groups[i].lowerCount()-1;bp<groups[i].upperCount();bp++) {		
				gcCount += gcCounts[bp];
				total += gcCounts[bp];
				total += atCounts[bp];
			}
			
			percentages[i] = 100*(gcCount/(double)total);
			mean += percentages[i];
		}
		
		mean /= percentages.length;
		
		// We now calculate the maximumn percentage deviation
		maxDeviation = Math.abs(mean - percentages[0]);
		for (int i=1;i<percentages.length;i++) {
			double diff = Math.abs(mean-percentages[i]);
			if (diff > maxDeviation) {
				maxDeviation = diff;
			}
		}
		
		calculated = true;
	}
	
	public void processSequence(Sequence sequence) {
		calculated = false;
		char [] seq = sequence.getSequence().toCharArray();

		if (gcCounts.length < seq.length) {
			// We need to expand the size of the data structures
			
			long [] gcCountsNew = new long [seq.length];
			long [] atCountsNew = new long [seq.length];

			for (int i=0;i<gcCounts.length;i++) {
				gcCountsNew[i] = gcCounts[i];
				atCountsNew[i] = atCounts[i];
			}
			
			gcCounts = gcCountsNew;
			atCounts = atCountsNew;
		}
		
		for (int i=0;i<seq.length;i++) {
			if (seq[i] == 'G' || seq[i] == 'C') {
				++gcCounts[i];
			}
			if (seq[i] == 'A' || seq[i] == 'T' || seq[i] == 'U') {
				++atCounts[i];
			}
		}
		
	}
	
	public void reset () {
		gcCounts = new long[0];
		atCounts = new long[0];
	}

	public String description() {
		return "Shows the GC content of all bases at a given position in a sequencing run";
	}

	public String name() {
		return "Per base GC content";
	}

	public boolean raisesError() {
		if (!calculated) getPercentages();
		
		return maxDeviation > 10;
	}

	public boolean raisesWarning() {
		if (!calculated) getPercentages();

		return maxDeviation > 5;
	}

	public void makeReport(HTMLReportArchive report) throws IOException {
		if (!calculated) getPercentages();
		ZipOutputStream zip = report.zipFile();
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/per_base_gc_content.png"));

		BufferedImage b = new BufferedImage(Math.max(800, percentages.length*15),600,BufferedImage.TYPE_INT_RGB);
		Graphics g = b.getGraphics();

		LineGraph lg = new LineGraph(new double [][] {percentages}, 0d, 100d, "Position in read (bp)", new String [] {"%GC"}, xCategories, "GC content across all bases");
		lg.paint(g,b.getWidth(),b.getHeight());

		ImageIO.write((BufferedImage)(b),"PNG",zip);

		StringBuffer sb = report.htmlDocument();

		sb.append("<p><img class=\"indented\" src=\"Images/per_base_gc_content.png\" alt=\"Per base GC content graph\"></p>\n");	
		
		sb = report.dataDocument();
		sb.append("#Base\t%GC\n");
		for (int i=0;i<xCategories.length;i++) {
			sb.append(xCategories[i]);
			sb.append("\t");
			sb.append(percentages[i]);
			sb.append("\n");
		}
	} 

}
