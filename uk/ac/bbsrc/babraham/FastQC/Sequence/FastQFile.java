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
package uk.ac.bbsrc.babraham.FastQC.Sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.itadaki.bzip2.BZip2InputStream;

public class FastQFile implements SequenceFile {

	private Sequence nextSequence = null;
	private File file;
	private long fileSize = 0;
	private int count = 0;
	
	private boolean casavaMode = false;
	
	// We actually read our final data from this buffered reader
	private BufferedReader br;
	
	// We keep the file stream around just so we can see how far through
	// the file we've got.  We don't read from this directly, but it's the
	// only way to access the file pointer.
	private FileInputStream fis;
	
	private String name;
	private boolean isColorspace = false;

	protected FastQFile(File file) throws SequenceFormatException, IOException {
		this.file = file;
		fileSize = file.length();
		name = file.getName();
		
		if (System.getProperty("fastqc.casava") != null && System.getProperty("fastqc.casava").equals("true")) {
			casavaMode = true;
		}
		
		fis = new FileInputStream(file);
		if (file.getName().toLowerCase().endsWith(".gz")) {
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis)));
		} 
		else if (file.getName().toLowerCase().endsWith(".bz2")) {
			br = new BufferedReader(new InputStreamReader(new BZip2InputStream(fis,false)));
		} 
		
		else {
			br = new BufferedReader(new InputStreamReader(fis));
		}
		readNext();
	}

	public String name() {
		return name;
	}

	public int getPercentComplete() {
		try {
			int percent = (int) (((double)fis.getChannel().position()/ fileSize)*100);
			return percent;
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		return 0;
	}

	public boolean isColorspace() {
		return isColorspace;
	}

	public void setIsColorspace(boolean isColorspace) {
		this.isColorspace = isColorspace;
	}

	public boolean hasNext() {
		return nextSequence != null;
	}

	public Sequence next() throws SequenceFormatException {
		Sequence seq = nextSequence;
		readNext();
		return seq;
	}

	private void readNext() throws SequenceFormatException {
		try {
			// First line should be the id
			
			// We might have blank lines between entries or at the end
			// so allow for this
			String id;
			
			while (true) {
				id = br.readLine();

				if (id == null) {
					nextSequence = null;
					return;
				}
				if (id.length()==0) {
					continue;
				}				
				
				break;
			}
			

			if (!id.startsWith("@")) {
				nextSequence = null;
				throw new SequenceFormatException("ID line didn't start with '@'");
			}

			// Then the sequence
			String seq = br.readLine();
			// Then another id which we don't need
			String midLine = br.readLine();
			if (!midLine.startsWith("+")) {
				throw new SequenceFormatException("Midline '"+midLine+"' didn't start with '+'");
			}
			// Then the quality string
			String quality = br.readLine();

			if (nextSequence == null && seq != null) {
				checkColorspace(seq);
			}

			if (quality == null) {
				nextSequence = null;
				br.close();
			} 
			else {
				++count;
				if (isColorspace()) {
					nextSequence = new Sequence(this,
							convertColorspaceToBases(seq.toUpperCase()), seq.toUpperCase(), quality, id);
				} 
				else {
					nextSequence = new Sequence(this, seq.toUpperCase(),quality, id);
				}
				
				// If we're running in --casava mode then we will flag any sequences which
				// are marked as being filtered.
				if (casavaMode) {
					
					// This is the test illumina suggest, but it's a bit flakey, and I'm not
					// sure it's not going to catch things it shouldn't.
					if (id.indexOf(":Y:") > 0) {
						nextSequence.setIsFiltered(true);
					}
				}

			}


		} 
		catch (IOException ioe) {
			nextSequence = null;
			ioe.printStackTrace();
		}
	}

	private void checkColorspace(String seq) {
		// Some basecalled files can be all dots, which leads to them
		// being identified as colorspace data. This check should find
		// only true colorspace files.
		String regex = "^[GATCNgatcn][\\.0123456]+$";
		Pattern pattern = Pattern.compile(regex);
		Matcher matcher = pattern.matcher(seq);
		if (matcher.find()) {
			isColorspace = true;
		} else {
			isColorspace = false;
		}
	}

	private String convertColorspaceToBases(String s) {

		char[] cs = s.toUpperCase().toCharArray();
		
		// We've had a crash report where a file contained a zero length
		// colorspace entry.  This is completely invalid, but we should
		// handle it anyway.
		if (cs.length == 0) {
			return "";
		}
		
		char[] bp = new char[cs.length - 1];

		char refBase;

		for (int i = 1; i < cs.length; i++) {
			if (i == 1) {
				refBase = cs[i - 1];
			} else {
				refBase = bp[i - 2];
			}
			if (!(refBase == 'G' || refBase == 'A' || refBase == 'T' || refBase == 'C')) {
				throw new IllegalArgumentException("Refbase was " + refBase
						+ " at position " + i);
			}
			switch (cs[i]) {
			case ('0'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'G';
					break;
				case ('A'):
					bp[i - 1] = 'A';
					break;
				case ('T'):
					bp[i - 1] = 'T';
					break;
				case ('C'):
					bp[i - 1] = 'C';
					break;
				}
				break;
			case ('1'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'T';
					break;
				case ('A'):
					bp[i - 1] = 'C';
					break;
				case ('T'):
					bp[i - 1] = 'G';
					break;
				case ('C'):
					bp[i - 1] = 'A';
					break;
				}
				break;

			case ('2'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'A';
					break;
				case ('A'):
					bp[i - 1] = 'G';
					break;
				case ('T'):
					bp[i - 1] = 'C';
					break;
				case ('C'):
					bp[i - 1] = 'T';
					break;
				}
				break;

			case ('3'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'C';
					break;
				case ('A'):
					bp[i - 1] = 'T';
					break;
				case ('T'):
					bp[i - 1] = 'A';
					break;
				case ('C'):
					bp[i - 1] = 'G';
					break;
				}
				break;

			case ('.'):
			case ('4'):
			case ('5'):
			case ('6'):
				for (; i < cs.length; i++) {
					bp[i - 1] = 'N';
				}
				break;
			default:
				throw new IllegalArgumentException("Unexpected cs char "
						+ cs[i]);
			}
		}

		return new String(bp);
	}

	public void remove() {
		// No action here
	}

	public File getFile() {
		return file;
	}

}
