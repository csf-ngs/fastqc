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
package uk.ac.bbsrc.babraham.FastQC;

import java.awt.BorderLayout;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import javax.swing.ImageIcon;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;

import uk.ac.bbsrc.babraham.FastQC.Analysis.AnalysisRunner;
import uk.ac.bbsrc.babraham.FastQC.Analysis.OfflineRunner;
import uk.ac.bbsrc.babraham.FastQC.Dialogs.WelcomePanel;
import uk.ac.bbsrc.babraham.FastQC.FileFilters.BAMFileFilter;
import uk.ac.bbsrc.babraham.FastQC.FileFilters.CasavaFastQFileFilter;
import uk.ac.bbsrc.babraham.FastQC.FileFilters.FastQFileFilter;
import uk.ac.bbsrc.babraham.FastQC.FileFilters.MappedBAMFileFilter;
import uk.ac.bbsrc.babraham.FastQC.FileFilters.SequenceFileFilter;
import uk.ac.bbsrc.babraham.FastQC.Modules.BasicStats;
import uk.ac.bbsrc.babraham.FastQC.Modules.KmerContent;
import uk.ac.bbsrc.babraham.FastQC.Modules.OverRepresentedSeqs;
import uk.ac.bbsrc.babraham.FastQC.Modules.PerBaseGCContent;
import uk.ac.bbsrc.babraham.FastQC.Modules.NContent;
import uk.ac.bbsrc.babraham.FastQC.Modules.PerBaseQualityScores;
import uk.ac.bbsrc.babraham.FastQC.Modules.PerBaseSequenceContent;
import uk.ac.bbsrc.babraham.FastQC.Modules.PerSequenceGCContent;
import uk.ac.bbsrc.babraham.FastQC.Modules.PerSequenceQualityScores;
import uk.ac.bbsrc.babraham.FastQC.Modules.QCModule;
import uk.ac.bbsrc.babraham.FastQC.Modules.SequenceLengthDistribution;
import uk.ac.bbsrc.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.bbsrc.babraham.FastQC.Results.ResultsPanel;
import uk.ac.bbsrc.babraham.FastQC.Sequence.SequenceFactory;
import uk.ac.bbsrc.babraham.FastQC.Sequence.SequenceFile;
import uk.ac.bbsrc.babraham.FastQC.Sequence.SequenceFormatException;
import uk.ac.bbsrc.babraham.FastQC.Utilities.CasavaBasename;

public class FastQCApplication extends JFrame {	
	
	public static final String VERSION = "0.10.0";
	
	private JTabbedPane fileTabs;
	private WelcomePanel welcomePanel;
	private File lastUsedDir = null;
	
	public FastQCApplication () {
			setTitle("FastQC");
			setIconImage(new ImageIcon(ClassLoader.getSystemResource("uk/ac/bbsrc/babraham/FastQC/Resources/fastqc_icon.png")).getImage());
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	//		setSize(1280, 720);
			setSize(800,600);
			setLocationRelativeTo(null);
			
			welcomePanel = new WelcomePanel();
			
			fileTabs = new JTabbedPane(JTabbedPane.TOP);
			setContentPane(welcomePanel);
			
			setJMenuBar(new FastQCMenuBar(this));
			
		}

	public void close () {
		if (fileTabs.getSelectedIndex() >=0) {
			fileTabs.remove(fileTabs.getSelectedIndex());
		}
		if (fileTabs.getTabCount() == 0) {
			setContentPane(welcomePanel);
			validate();
			repaint();
		}
	}
	
	public void closeAll () {
		fileTabs.removeAll();
		setContentPane(welcomePanel);
		validate();
		repaint();
	}
	
	public void openFile () {
		JFileChooser chooser;
		
		if (lastUsedDir == null) {
			chooser = new JFileChooser();
		}
		else {
			chooser = new JFileChooser(lastUsedDir);
		}
		chooser.setMultiSelectionEnabled(true);
		SequenceFileFilter sff = new SequenceFileFilter();
		chooser.addChoosableFileFilter(sff);
		chooser.addChoosableFileFilter(new FastQFileFilter());
		chooser.addChoosableFileFilter(new CasavaFastQFileFilter());
		chooser.addChoosableFileFilter(new BAMFileFilter());
		chooser.addChoosableFileFilter(new MappedBAMFileFilter());
		chooser.setFileFilter(sff);
		int result = chooser.showOpenDialog(this);
		if (result == JFileChooser.CANCEL_OPTION) return;
	
		// See if they forced a file format
		FileFilter chosenFilter = chooser.getFileFilter();
		if (chosenFilter instanceof FastQFileFilter) {
			System.setProperty("fastqc.sequence_format", "fastq");
		}
		if (chosenFilter instanceof CasavaFastQFileFilter) {
			System.setProperty("fastqc.sequence_format", "fastq");
			System.setProperty("fastqc.casava","true");
		}
		else if (chosenFilter instanceof BAMFileFilter) {
			System.setProperty("fastqc.sequence_format", "bam");
		}
		else if (chosenFilter instanceof MappedBAMFileFilter) {
			System.setProperty("fastqc.sequence_format", "bam_mapped");
		}
		
		// If we're still showing the welcome panel switch this out for
		// the file tabs panel
		if (fileTabs.getTabCount() == 0) {
			setContentPane(fileTabs);
			validate();
			repaint();
		}
		
		File [] files = chooser.getSelectedFiles();
		
		File [][] fileGroups;
		
		// See if we need to group together files from a casava group
		if (System.getProperty("fastqc.casava") != null && System.getProperty("fastqc.casava").equals("true")) {
			fileGroups = CasavaBasename.getCasavaGroups(files);
		}
		else {
			fileGroups = new File [files.length][1];
			for (int f=0;f<files.length;f++) {
				fileGroups[f][0] = files[f];
			}
		}

	
		for (int i=0;i<fileGroups.length;i++) {
			File [] filesToProcess = fileGroups[i];
			lastUsedDir = filesToProcess[0].getParentFile();
			SequenceFile sequenceFile;
			
			
			try {
				if (filesToProcess.length > 1) {
					sequenceFile = SequenceFactory.getSequenceFile(filesToProcess, 0);
				}
				else {
					sequenceFile = SequenceFactory.getSequenceFile(filesToProcess[0], 0);
				}
			}
			catch (SequenceFormatException e) {
				JPanel errorPanel = new JPanel();
				errorPanel.setLayout(new BorderLayout());
				errorPanel.add(new JLabel("File format error: "+e.getLocalizedMessage(), JLabel.CENTER),BorderLayout.CENTER);
				fileTabs.addTab(filesToProcess[0].getName(), errorPanel);
				e.printStackTrace();
				continue;
			}
			catch (IOException e) {
				System.err.println("File broken");
				e.printStackTrace();
				JOptionPane.showMessageDialog(this, "Couldn't read file:"+e.getLocalizedMessage(), "Error reading file", JOptionPane.ERROR_MESSAGE);
				continue;
			}
					
			AnalysisRunner runner = new AnalysisRunner(sequenceFile);
			ResultsPanel rp = new ResultsPanel(sequenceFile);
			runner.addAnalysisListener(rp);
			fileTabs.addTab(sequenceFile.name(), rp);
			
			OverRepresentedSeqs os = new OverRepresentedSeqs();
			
			QCModule [] module_list = new QCModule [] {
				new BasicStats(),
				new PerBaseQualityScores(),
				new PerSequenceQualityScores(),
				new PerBaseSequenceContent(),
				new PerBaseGCContent(), 
				new PerSequenceGCContent(),
				new NContent(),
				new SequenceLengthDistribution(),
				os.duplicationLevelModule(),
				os,
				new KmerContent()
			};
	
			runner.startAnalysis(module_list);
		}
	}

	public void saveReport () {
		JFileChooser chooser;
		
		if (lastUsedDir == null) {
			chooser = new JFileChooser();
		}
		else {
			chooser = new JFileChooser(lastUsedDir);
		}
		
		if (fileTabs.getSelectedComponent() == null) {
			JOptionPane.showMessageDialog(this, "No FastQ files are open yet", "Can't save report", JOptionPane.ERROR_MESSAGE);
			return;
		}
		chooser.setSelectedFile(new File(((ResultsPanel)fileTabs.getSelectedComponent()).sequenceFile().getFile().getName().replaceAll(".gz$","").replaceAll(".bz2$","").replaceAll(".txt$","").replaceAll(".fastq$", "").replaceAll(".sam$", "").replaceAll(".bam$", "")+"_fastqc.zip"));
		chooser.setMultiSelectionEnabled(false);
		chooser.setFileFilter(new FileFilter() {
		
			public String getDescription() {
				return "Zip files";
			}
		
			public boolean accept(File f) {
				if (f.isDirectory() || f.getName().toLowerCase().endsWith(".zip")) {
					return true;
				}
				else {
					return false;
				}
			}
		
		});
	
		File reportFile;
		while (true) {
			int result = chooser.showSaveDialog(this);
			if (result == JFileChooser.CANCEL_OPTION) return;
			
			reportFile = chooser.getSelectedFile();
			if (! reportFile.getName().toLowerCase().endsWith(".zip")) {
				reportFile = new File(reportFile.getAbsoluteFile()+".zip");
			}
			
			// Check if we're overwriting something
			if (reportFile.exists()) {
				int reply = JOptionPane.showConfirmDialog(this, reportFile.getName()+" already exists.  Overwrite?", "Overwrite existing file?", JOptionPane.YES_NO_OPTION);
				if (reply == JOptionPane.NO_OPTION) {
					continue;
				}
				else {
					break;
				}
			}
			else {
				break;
			}
		}
		
		
		ResultsPanel selectedPanel = (ResultsPanel)fileTabs.getSelectedComponent();
		
		try {
			new HTMLReportArchive(selectedPanel.sequenceFile(), selectedPanel.modules(), reportFile, 0);
		} 
		catch (IOException e) {
			JOptionPane.showMessageDialog(this, "Failed to create archive: "+e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
        System.out.println("arguments: " + Arrays.toString(args));

		// See if we just have to print out the version
		if (System.getProperty("fastqc.show_version") != null && System.getProperty("fastqc.show_version").equals("true")) {
			System.out.println("FastQC v"+VERSION);
			System.exit(0);
		}
		
		if (args.length > 0) {
			// Set headless to true so we don't get problems
			// with people working without an X display.
			System.setProperty("java.awt.headless", "true");
			
			// The non-interactive default is to uncompress the
			// reports after they have been generated
			if (System.getProperty("fastqc.unzip") == null || ! System.getProperty("fastqc.unzip").equals("false")) {
				System.setProperty("fastqc.unzip", "true");
			}
			if(args.length == 2){
                int read = Integer.parseInt(args[0]);
                String[] files = new String[]{args[1]};
                new OfflineRunner(files, read);
            }
			System.exit(0);
		}
		
		else {
			try {
				UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			} catch (Exception e) {}
			
	
			// The interactive default is to not uncompress the
			// reports after they have been generated
			if (System.getProperty("fastqc.unzip") == null || !System.getProperty("fastqc.unzip").equals("true")) {
				System.setProperty("fastqc.unzip", "false");
			}
	
			FastQCApplication app = new FastQCApplication();
	
			app.setVisible(true);
		}
	}	

}
