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

import java.awt.BorderLayout;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import uk.ac.bbsrc.babraham.FastQC.Graphs.BaseGroup;
import uk.ac.bbsrc.babraham.FastQC.Graphs.LineGraph;
import uk.ac.bbsrc.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.bbsrc.babraham.FastQC.Sequence.Sequence;

public class KmerContent implements QCModule {

	private Hashtable<String, Kmer> kmers = new Hashtable<String, Kmer>((int)Math.pow(4, MAX_KMER_SIZE));
	private long gCount = 0;
	private long aCount = 0;
	private long tCount = 0;
	private long cCount = 0;
	
	private int longestSequence = 0;
	private long [][] totalKmerCounts = new long [0][0];
	private long skipCount = 0;
	
	private static final int MIN_KMER_SIZE = 5;
	private static final int MAX_KMER_SIZE = 5;
	
	public boolean calculated = false;
	
	private Kmer [] enrichedKmers = null;
	private double [][] enrichments = null;
	private String [] xCategories = new String[0];
	private String [] xLabels = new String[0];
	
	BaseGroup [] groups;

	public boolean ignoreFilteredSequences() {
		return true;
	}
	
	public JPanel getResultsPanel() {
		
		if (!calculated) calculateEnrichment();
		JPanel returnPanel = new JPanel();
		returnPanel.setLayout(new BorderLayout());
		returnPanel.add(new JLabel("Overrepresented Kmers",JLabel.CENTER),BorderLayout.NORTH);
		
		JSplitPane splitPanel = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
		
		if (enrichedKmers.length > 0) {
			TableModel model = new ResultsTable(enrichedKmers);
			splitPanel.setBottomComponent(new JScrollPane(new JTable(model)));
			splitPanel.setTopComponent(new LineGraph(enrichments, 0d, 100d, "Position in read (bp)", xLabels, xCategories, "Relative enrichment over read length"));
			returnPanel.add(splitPanel,BorderLayout.CENTER);
		}
		else {
			returnPanel.add(new JLabel("There are no overrepresented Kmers",JLabel.CENTER),BorderLayout.CENTER);
		}
		
		return returnPanel;
	}
	
	/**
	 * This method simply keeps a count of the number of Kmers of a given size
	 * seen at each position within the run.  We can use this later on to calculate
	 * the enrichment of the Kmers we actually count.
	 * 
	 * We take in the Kmer sequence even though this isn't used in the total counts
	 * we do this because we don't want to count Kmers with Ns in them, but we do
	 * need to ensure that the data structure is expanded to the right size, and if
	 * we have libraries where later positions are Ns in all sequences then our
	 * data structure ends up too short and we crash. 
	 * 
	 * @param position Position within the read.  0 indexed
	 * @param kmerLength Actual length of the Kmer analysed
	 */
	private void addKmerCount (int position,int kmerLength, String kmer) {
	
		
		if (position >= totalKmerCounts.length) {
			// We need to expand the array
			long [][] newCounts = new long[position+1][];
			for (int i=0;i<totalKmerCounts.length;i++) {
				newCounts[i] = totalKmerCounts[i];
			}
			for (int i=totalKmerCounts.length;i<newCounts.length;i++) {
				newCounts[i] = new long[MAX_KMER_SIZE];
			}
			
			totalKmerCounts = newCounts;
		}
		
		if (kmer.indexOf("N") >=0) return;

		++totalKmerCounts[position][kmerLength-1];
		
	}

	private synchronized void calculateEnrichment () {
		
		// For each kmer we work out the number of times we should have
		// seen this by chance, and then the obs/exp for each of the
		// kmers.  We can then filter for a specific level of enrichment
		// and show only those Kmers
		
		float totalBases = gCount+aCount+tCount+cCount;
		
		float gProb = ((float)gCount)/totalBases;
		float aProb = ((float)aCount)/totalBases;
		float tProb = ((float)tCount)/totalBases;
		float cProb = ((float)cCount)/totalBases;
		
		// We'll be grouping together positions later so make up the groups now
		groups = BaseGroup.makeBaseGroups((longestSequence-MIN_KMER_SIZE)+1);

		Vector<Kmer>enrichedKmers = new Vector<Kmer>();
				
		Iterator<Kmer> rawKmers = kmers.values().iterator();
		
		KMER: while (rawKmers.hasNext()) {
			Kmer k = rawKmers.next();
			
			long totalKmerCount = 0;

			for (int i=0;i<totalKmerCounts.length;i++) {
				totalKmerCount += totalKmerCounts[i][k.sequence().length()-1];
			}
			
			
			float prob = 1;
			char [] chars = k.sequence().toCharArray();
			for (int c=0;c<chars.length;c++) {
				switch (chars[c]) {
				case 'G': prob *= gProb;break;
				case 'A': prob *= aProb;break;
				case 'T': prob *= tProb;break;
				case 'C': prob *= cProb;break;
				default: continue KMER; // Ignore Kmers containing non-GATC chars
				}
			}
			// Now work out how many of these kmers we should have seen
			float predicted = prob * totalKmerCount;
						
			k.setObsExp(k.count()/predicted);
				
			// We shall also calculate the positional variation in obs/exp
				
			float [] obsExpPositions = new float[groups.length];
				
			long [] positionCounts = k.getPositions();
				
			for (int g=0;g<groups.length;g++) {
				// This is a summation of the number of Kmers of this length which
				// fall into this base group
				long totalGroupCount = 0;
				
				// This is a summation of the number of hit Kmers which fall within
				// this base group.
				long totalGroupHits = 0;
				for (int p=groups[g].lowerCount()-1;p<groups[g].upperCount() && p < positionCounts.length ;p++) {
					totalGroupCount += totalKmerCounts[p][chars.length-1];
					totalGroupHits += positionCounts[p];
				}
			
				// We used to have a filter here which provided a default value where there
				// were fewer than 1000 observations at a given position.  This caused breakage
				// in small files where every position had fewer than 1000 observations and since
				// I can't see why this exclusion was there in the first place, I've now removed
				// it.  If the reason for its original inclusion reappears then we can put it back
				// (with a more sensible default value), but for the moment this appears to be the
				// better fix.
				
				predicted = prob * totalGroupCount;
				obsExpPositions[g] = totalGroupHits/predicted;
			}
			k.setObsExpPositions(obsExpPositions);
			
			if (k.obsExp() > 3 || k.maxObsExp() > 5) {
				enrichedKmers.add(k);
			}			
			
		}
		
		Kmer [] finalKMers = enrichedKmers.toArray(new Kmer[0]);
		Arrays.sort(finalKMers);				
		
		// Now we take the enrichment positions for the top 6 hits and
		// record these so we can plot them on a line graph
		enrichments = new double [Math.min(6, finalKMers.length)][];
		xLabels = new String[enrichments.length];
		
		xCategories = new String [groups.length];
		
		for (int i=0;i<xCategories.length;i++) {
			xCategories[i] = groups[i].toString();
		}
		
		for (int k=0;k<enrichments.length;k++) {
			enrichments[k] = new double[groups.length];
			
			float [] obsExpPos = finalKMers[k].getObsExpPositions();
			
			// Find the max enrichment for this sequence
			float max = 0;
			for (int i=0;i<obsExpPos.length;i++) {								
				if (obsExpPos[i] > max) max = obsExpPos[i]; 
			}
			
			for (int g=0;g<groups.length;g++) {				
				enrichments[k][g] = obsExpPos[g]/max * 100; 
			}
			
			xLabels[k] = finalKMers[k].sequence();
			
		}
		
		this.enrichedKmers = finalKMers;		
		
		calculated = true;
	}
	
		
	public void processSequence(Sequence sequence) {
		calculated = false;
		
		++skipCount;
		if (skipCount % 5 != 0) return;
		
		char [] seq = sequence.getSequence().toCharArray();

		if (seq.length > longestSequence) {
			longestSequence = seq.length;
		}
				
		for (int i=0;i<seq.length;i++) {
			switch (seq[i]) {
			case 'G': ++gCount;break;
			case 'A': ++aCount;break;
			case 'T': ++tCount;break;
			case 'C': ++cCount;break;
			}
		}
		
		// Now we go through all of the Kmers to count these
		for (int kmerSize=MIN_KMER_SIZE;kmerSize<=MAX_KMER_SIZE;kmerSize++) {
			for (int i=0;i<=seq.length-kmerSize;i++) {
				
				String kmer = sequence.getSequence().substring(i, i+kmerSize);
				
				// Add to the counts before skipping Kmers containing Ns (see
				// explanation in addKmerCount for the reasoning).
				addKmerCount(i, kmerSize, kmer);
				
				// Skip Kmers containing N
				if (kmer.indexOf("N") >=0) continue;

				if (kmers.containsKey(kmer)) {
					kmers.get(kmer).incrementCount(i);
				}
				else {
					kmers.put(kmer, new Kmer(kmer,i,(seq.length-kmerSize)+1));
				}

			}
		}
	}
	
	public void reset () {
		calculated = false;
		gCount = 0;
		aCount = 0;
		tCount = 0;
		cCount = 0;
		totalKmerCounts = new long[0][0];
		longestSequence = 0;
		skipCount = 0;
		enrichedKmers = null;
		kmers.clear();
	}

	public String description() {
		return "Identifies short sequences which are overrepresented";
	}

	public String name() {
		return "Kmer Content";
	}

	public boolean raisesError() {
		if (!calculated) calculateEnrichment();
		
		// We raise an error if the most enriched kmer is seen more than 100 times
		// more frequently than we expect.
		
		if (enrichedKmers.length > 0 && enrichedKmers[0].maxObsExp() > 10) return true;
		return false;
	}

	public boolean raisesWarning() {
		if (!calculated) calculateEnrichment();
		
		// We raise a warning if there are any enriched kmers
		if (enrichedKmers.length > 0) return true;
		return false;
	}

	public void makeReport(HTMLReportArchive report) throws IOException {
		if (!calculated) calculateEnrichment();
		
		if (enrichedKmers.length > 0) {
			ZipOutputStream zip = report.zipFile();
			zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/kmer_profiles.png"));

			BufferedImage b = new BufferedImage(800,600,BufferedImage.TYPE_INT_RGB);
			Graphics g = b.getGraphics();

			LineGraph lg = new LineGraph(enrichments, 0d, 100d, "Position in read (bp)", xLabels, xCategories, "Relative enrichment over read length");
			lg.paint(g,800,600);

			ImageIO.write((BufferedImage)(b),"PNG",zip);

			StringBuffer sb = report.htmlDocument();

			sb.append("<p><img class=\"indented\" src=\"Images/kmer_profiles.png\" alt=\"Kmer graph\"></p>\n");

		}
		
		
		ResultsTable table = new ResultsTable(enrichedKmers);
		
		StringBuffer b = report.htmlDocument();
		StringBuffer d = report.dataDocument();
		
		if (enrichedKmers.length == 0) {
			b.append("<p>No overrepresented Kmers</p>\n");
		}
		
		else {
			b.append("<table>\n");
			// Do the headers
			b.append("<tr>\n");
			d.append("#");
			for (int c=0;c<table.getColumnCount();c++) {
				b.append("<th>");
				b.append(table.getColumnName(c));
				d.append(table.getColumnName(c));
				b.append("</th>\n");
				if (c<table.getColumnCount()-1) {
					d.append("\t");
				}
			}
			b.append("</tr>\n");
			d.append("\n");
			
			// Do the rows
			for (int r=0;r<table.getRowCount();r++) {
				b.append("<tr>\n");
				for (int c=0;c<table.getColumnCount();c++) {
					b.append("<td>");
					b.append(table.getValueAt(r, c));
					d.append(table.getValueAt(r, c));
					b.append("</td>\n");
					if (c<table.getColumnCount()-1) {
						d.append("\t");
					}
				}
				b.append("</tr>\n");
				d.append("\n");
			}
			
			b.append("</table>\n");
		}	
	}
	
	private class Kmer implements Comparable<Kmer>{
		
		private String sequence;
		private long count = 0;
		private float obsExp = 0;
		private float [] obsExpPositions = null;
		private long [] positions = new long[0];
		
		public Kmer (String sequence, int position, int seqLength) {

			// Do this slightly convoluted dance to try to avoid
			// keeping the whole original sequence in memory
			char [] chars = sequence.toCharArray();
			this.sequence = new String(chars);
			count = 1;
			positions = new long[seqLength];
			++positions[position];
		}
		
		public void incrementCount (int position) {
			++count;
			
			if (position >= positions.length) {
				long [] newPositions = new long[position+1];
				for (int i=0;i<positions.length;i++) {
					newPositions[i] = positions[i];
				}
				positions = newPositions;
			}
			
			++positions[position];
			
		}
		
		public long [] getPositions () {
			return positions;
		}
		
		public String sequence () {
			return sequence;
		}
		
		public long count () {
			return count;
		}
		
		public void setObsExp (float oe) {
			this.obsExp = oe;
		}
		
		public void setObsExpPositions (float [] oePositions) {
			this.obsExpPositions = oePositions;
		}
		
		public float [] getObsExpPositions () {
			return obsExpPositions;
		}
		
		public float obsExp () {
			return obsExp;
		}
		
		public float maxObsExp () {
			float max = 0;
			for (int i=0;i<obsExpPositions.length;i++) {
				if (obsExpPositions[i]>max) max = obsExpPositions[i];
			}
			return max;
		}

		public int maxPosition () {
			float max = 0;
			int position = 0;
			for (int i=0;i<obsExpPositions.length;i++) {
				if (obsExpPositions[i]>max) {
					max = obsExpPositions[i];
					position = i+1;
				}
			}
			
			if (position == 0) {
				System.err.println("No value > 0 for "+sequence);
				position = 1;
			}
			
			return position;
		}

		public int compareTo(Kmer o) {
			return Float.compare(o.obsExp(), obsExp());
		}
	}
	
	private class ResultsTable extends AbstractTableModel {
		
		private Kmer [] kmers;
		
		public ResultsTable (Kmer [] kmers) {
			this.kmers = kmers;
		}
		
		
		// Sequence - Count - Obs/Exp
		public int getColumnCount() {
			return 5;
		}

		public int getRowCount() {
			return kmers.length;
		}

		public Object getValueAt(int rowIndex, int columnIndex) {
			switch (columnIndex) {
				case 0: return kmers[rowIndex].sequence();
				case 1: return kmers[rowIndex].count()*5;
				case 2: return kmers[rowIndex].obsExp();
				case 3: return kmers[rowIndex].maxObsExp();
				case 4: return groups[kmers[rowIndex].maxPosition()-1].toString();
						
			}
			return null;
		}
		
		public String getColumnName (int columnIndex) {
			switch (columnIndex) {
				case 0: return "Sequence";
				case 1: return "Count";
				case 2: return "Obs/Exp Overall";
				case 3: return "Obs/Exp Max";
				case 4: return "Max Obs/Exp Position";
			}
			return null;
		}
		
		public Class<?> getColumnClass (int columnIndex) {
			switch (columnIndex) {
			case 0: return String.class;
			case 1: return Integer.class;
			case 2: return Float.class;
			case 3: return Float.class;
			case 4: return String.class;
		}
		return null;
			
		}
	}

}
