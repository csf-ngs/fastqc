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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;

public class BAMFile implements SequenceFile {

	private File file;
	private boolean onlyMapped;
	private long fileSize = 0;
	private long recordSize = 0;
	private int count = 0;
	private int rawCount = 0;
	
	// We keep the file stream around just so we can see how far through
	// the file we've got.  We don't read from this directly, but it's the
	// only way to access the file pointer.
	private FileInputStream fis;

	private SAMFileReader br;
	private String name;
	private Sequence nextSequence = null;
	Iterator<SAMRecord> it;
	
	
	protected BAMFile (File file, boolean onlyMapped) throws SequenceFormatException, IOException {
		this.file = file;
		fileSize = file.length();
		name = file.getName();
		this.onlyMapped = onlyMapped;

		SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);

		fis = new FileInputStream(file);
		
		br = new SAMFileReader(fis);
		
		it = br.iterator();
		readNext();
	}
	
	public String name () {
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

	public boolean isColorspace () {
		return false;
	}
		
	public boolean hasNext() {
		return nextSequence != null;
	}

	public Sequence next () throws SequenceFormatException {
		Sequence returnSeq = nextSequence;
		readNext();
		return returnSeq;
	}
	
	private void readNext() throws SequenceFormatException {
		
		SAMRecord record;
		
		while (true) {
			
			if (!it.hasNext()) {
				nextSequence = null;
				return;
			}
		
			try {
				record = it.next();
			}
			catch (SAMFormatException sfe) {
				throw new SequenceFormatException(sfe.getMessage());
			}
			++rawCount;
		
			// We skip over entries with no mapping if that's what the user asked for
			if (onlyMapped && record.getReadUnmappedFlag()) {
				continue;
			}
			else {
				break;
			}
		}
		
		++count;

		if (recordSize == 0) {
			recordSize = (record.getReadLength()*2)+150;
			if (br.isBinary()) {
				recordSize /= 4;
			}
		}

		String sequence = record.getReadString();
		String qualities = record.getBaseQualityString();

		// BAM/SAM files always show sequence relative to the top strand of
		// the mapped reference so if this sequence maps to the reverse strand
		// we need to reverse complement the sequence and reverse the qualities
		// to get the original orientation of the read.
		if (record.getReadNegativeStrandFlag()) {
			sequence = reverseComplement(sequence);
			qualities = reverse(qualities);
		}

		nextSequence = new Sequence(this, sequence, qualities, record.getReadName());

	}

	
	private String reverseComplement (String sequence) {
		
		char [] letters = reverse(sequence).toUpperCase().toCharArray();
		char [] rc = new char[letters.length];
		
		for (int i=0;i<letters.length;i++) {
			switch(letters[i]) {
			case 'G': rc[i] = 'C';break;
			case 'A': rc[i] = 'T';break;
			case 'T': rc[i] = 'A';break;
			case 'C': rc[i] = 'G';break;
			default: rc[i] = letters[i];
			}
		}
	
		return new String(rc);

	}
	
	private String reverse (String sequence) {
		char [] starting = sequence.toCharArray();
		char [] reversed = new char[starting.length];
		
		for (int i=0;i<starting.length;i++) {
			reversed[reversed.length-(1+i)] = starting[i];
		}
		
		return new String(reversed);
	}

	public File getFile() {
		return file;
	}
	
}
