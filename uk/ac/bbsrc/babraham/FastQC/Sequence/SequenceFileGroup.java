package uk.ac.bbsrc.babraham.FastQC.Sequence;

import java.io.File;

import uk.ac.bbsrc.babraham.FastQC.Utilities.CasavaBasename;
import uk.ac.bbsrc.babraham.FastQC.Utilities.NameFormatException;

public class SequenceFileGroup implements SequenceFile {

	private SequenceFile [] files;
	private File groupFile;
	private int currentIndex = 0;
	
	public SequenceFileGroup (SequenceFile [] files) {
		this.files = files;
		try {
			String baseName = CasavaBasename.getCasavaBasename(files[0].name());
			if (files[0].getFile().getParent() == null) {
				groupFile = new File(baseName);
			}
			else {
				groupFile = new File(files[0].getFile().getParent()+"/"+baseName);				
			}
		}
		catch (NameFormatException nfe) {
			groupFile = files[0].getFile();
		}
	}
	
	public File getFile() {
		return groupFile;
	}

	public int getPercentComplete() {
		return ((100*currentIndex)/files.length)+(files[currentIndex].getPercentComplete()/files.length);
	}

	public boolean hasNext() {
		if (files[currentIndex].hasNext()) {
			return true;
		}
		else {
			if (currentIndex < files.length-1) {
				++currentIndex;
			}
			return files[currentIndex].hasNext();
		}
	}

	public boolean isColorspace() {
		return files[currentIndex].isColorspace();
	}

	public String name() {
		return groupFile.getName();
	}

	public Sequence next() throws SequenceFormatException {
		return files[currentIndex].next();
	}

}
