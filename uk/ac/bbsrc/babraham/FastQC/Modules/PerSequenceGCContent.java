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
import java.util.Hashtable;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

import uk.ac.bbsrc.babraham.FastQC.Graphs.LineGraph;
import uk.ac.bbsrc.babraham.FastQC.Modules.GCModel.GCModel;
import uk.ac.bbsrc.babraham.FastQC.Modules.GCModel.GCModelValue;
import uk.ac.bbsrc.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.bbsrc.babraham.FastQC.Sequence.Sequence;
import uk.ac.bbsrc.babraham.FastQC.Statistics.NormalDistribution;

public class PerSequenceGCContent implements QCModule {

	private double [] gcDistribution = new double[101];
	private double [] theoreticalDistribution  = new double[101];
	private int [] xCategories = new int[0];
	private double max = 0;
	private double deviationPercent;
	private boolean calculated = false;
	
	private Hashtable<Integer, GCModel> cachedModels = new Hashtable<Integer, GCModel>();
	
	public JPanel getResultsPanel() {
	
		if (!calculated) calculateDistribution();
				
		return new LineGraph(new double [][] {gcDistribution,theoreticalDistribution}, 0d, max, "Mean GC content (%)", new String [] {"GC count per read","Theoretical Distribution"}, xCategories, "GC distribution over all sequences");
	}
	
	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	
	private synchronized void calculateDistribution () {
		max = 0;
		xCategories = new int[gcDistribution.length];
		double totalCount = 0;
		
		
		// We use the mode to calculate the theoretical distribution
		// so that we cope better with skewed distributions.
		int firstMode = 0;
		double modeCount = 0;
		
		for (int i=0;i<gcDistribution.length;i++) {
			xCategories[i] = i;
			totalCount += gcDistribution[i];
			
			if (gcDistribution[i] > modeCount) {
				modeCount = gcDistribution[i];
				firstMode = i;
			}
			if (gcDistribution[i] > max) max = gcDistribution[i];
		}
				
		// The mode might not be a very good measure of the centre
		// of the distribution either due to duplicated vales or
		// several very similar values next to each other.  We therefore
		// average over adjacent points which stay above 95% of the modal
		// value

		double mode = 0;
		int modeDuplicates = 0;
		
		boolean fellOffTop = true;

		for (int i=firstMode;i<gcDistribution.length;i++) {
			if (gcDistribution[i] > gcDistribution[firstMode] - (gcDistribution[firstMode]/10)) {
				mode += i;
				modeDuplicates++;
			}
			else {
				fellOffTop = false;
				break;
			}
		}

		boolean fellOffBottom = true;
		
		for (int i=firstMode-1;i>=0;i--) {
			if (gcDistribution[i] > gcDistribution[firstMode] - (gcDistribution[firstMode]/10)) {
				mode += i;
				modeDuplicates++;
			}
			else {
				fellOffBottom = false;
				break;
			}
		}

		if (fellOffBottom || fellOffTop) {
			// If the distribution is so skewed that 95% of the mode
			// is off the 0-100% scale then we keep the mode as the 
			// centre of the model
			mode = firstMode;
		}
		else {
			mode /= modeDuplicates;
		}
		
		
		
		// We can now work out a theoretical distribution
		double stdev = 0;
		
		for (int i=0;i<gcDistribution.length;i++) {
			stdev += Math.pow((i-mode),2) * gcDistribution[i];
		}
		
		stdev /= totalCount-1;
		
		stdev = Math.sqrt(stdev);
		
		NormalDistribution nd = new NormalDistribution(mode, stdev);
		
		deviationPercent = 0;
		
		for (int i=0;i<theoreticalDistribution.length;i++) {
			double probability = nd.getZScoreForValue(i);
			theoreticalDistribution[i] = probability*totalCount;
			
			if (theoreticalDistribution[i] > max) {
				max = theoreticalDistribution[i];
			}
			
			deviationPercent += Math.abs(theoreticalDistribution[i]-gcDistribution[i]);
		}
		
		deviationPercent /= totalCount;
		deviationPercent *= 100;
		
//		System.out.println("Percentage deviation from normality is "+deviationPercent);
		
		
		calculated = true;
	}

	public void processSequence(Sequence sequence) {
		
		// Because we keep a model around for every possible sequence length we
		// encounter we need to reduce the number of models.  We can do this by
		// rounding off the sequence once we get above a certain size
		
		
		char [] seq = truncateSequence(sequence);
		
		if (seq.length == 0) return; // Ignore empty sequences
		
		
		int thisSeqGCCount = 0;
		for (int i=0;i<seq.length;i++) {
			if (seq[i] == 'G' || seq[i] == 'C') {
				++thisSeqGCCount;
			}
		}

		if (! cachedModels.containsKey(seq.length)) {
			cachedModels.put(seq.length,new GCModel(seq.length));
		}

		GCModelValue [] values = cachedModels.get(seq.length).getModelValues(thisSeqGCCount);

		for (int i=0;i<values.length;i++) {
			gcDistribution[values[i].percentage()] += values[i].increment();
		}
		
	}
	
	private char [] truncateSequence (Sequence sequence) {
		
		String seq = sequence.getSequence();
		
		// TODO: We should return a random chunk of sequence, rather
		// than the start.
		
		if (seq.length() > 1000) {
			int length = (seq.length()/1000)*1000;
			return seq.substring(0, length).toCharArray();
		}
		if (seq.length() > 100) {
			int length = (seq.length()/100)*100;
			return seq.substring(0, length).toCharArray();
		}

		return seq.toCharArray();		
		
	}
	
	public void reset () {
		gcDistribution = new double[101];
	}

	public String description() {
		return "Shows the distribution of GC contents for whole sequences";
	}

	public String name() {
		return "Per sequence GC content";
	}

	public boolean raisesError() {
		if (!calculated) calculateDistribution();

		return deviationPercent > 30;
	}

	public boolean raisesWarning() {
		if (!calculated) calculateDistribution();

		return deviationPercent > 15;
	}

	public void makeReport(HTMLReportArchive report) throws IOException {
		ZipOutputStream zip = report.zipFile();
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/per_sequence_gc_content.png"));

		BufferedImage b = new BufferedImage(800,600,BufferedImage.TYPE_INT_RGB);
		Graphics g = b.getGraphics();

		LineGraph lg = new LineGraph(new double [][] {gcDistribution,theoreticalDistribution}, 0d, max, "Mean GC content (%)", new String [] {"GC count per read","Theoretical Distribution"}, xCategories, "GC distribution over all sequences");

		lg.paint(g,800,600);

		ImageIO.write((BufferedImage)(b),"PNG",zip);

		StringBuffer sb = report.htmlDocument();

		sb.append("<p><img class=\"indented\" src=\"Images/per_sequence_gc_content.png\" alt=\"Per sequence GC content graph\"></p>\n");

		sb = report.dataDocument();
		sb.append("#GC Content\tCount\n");
		for (int i=0;i<xCategories.length;i++) {
			sb.append(xCategories[i]);
			sb.append("\t");
			sb.append(gcDistribution[i]);
			sb.append("\n");
		}
	}
}
