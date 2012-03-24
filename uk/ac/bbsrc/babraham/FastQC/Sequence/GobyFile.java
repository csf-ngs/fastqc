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

//import edu.cornell.med.icb.goby.reads.ReadsReader;
//import edu.cornell.med.icb.goby.reads.Reads.ReadEntry;

public class GobyFile implements SequenceFile {

	private File file;
	private long fileSize = 0;
	private int rawCount = 0;
	
	// We keep the file stream around just so we can see how far through
	// the file we've got.  We don't read from this directly, but it's the
	// only way to access the file pointer.
	private FileInputStream fis;

//	private ReadsReader br;
	private String name;	
	
	protected GobyFile (File file) throws SequenceFormatException, IOException {
		this.file = file;
		fileSize = file.length();
		name = file.getName();

		fis = new FileInputStream(file);
		
//		br = new ReadsReader(fis);
		
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
//		if (br.hasNext()) return true;
//		else {
//			try {
//				br.close();
//			} 
//			catch (IOException e) {
//				e.printStackTrace();
//			}
			return false;
//		}
	}

	public Sequence next () throws SequenceFormatException {
//		ReadEntry entry = br.next();
		
		++rawCount;
		
//		String sequence = entry.getSequence().toStringUtf8();
//		
//		if (entry.getQualityScores().isEmpty()) {
//			throw new SequenceFormatException("No quality scores associates with these data");
//		}
//		
//		String qualities = entry.getQualityScores().toStringUtf8();
//		
//		// The bytes here are encoded without an offset so we'll apply a sanger offset of 33
//		// so the rest of the pipeline understands them.
//		byte [] qualBytes = qualities.getBytes();
//
//		for (int i=0;i<qualBytes.length;i++) {
//			qualBytes[i] += 33;
//		}
//		
//		qualities = new String(qualBytes);
//		
//		Sequence returnSeq = new Sequence(this, sequence, qualities, "");
//		return returnSeq;
		
		return null;
	}

	public File getFile() {
		return file;
	}
	
}
