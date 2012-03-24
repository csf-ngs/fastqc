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
package uk.ac.bbsrc.babraham.FastQC.Report;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import uk.ac.bbsrc.babraham.FastQC.FastQCApplication;
import uk.ac.bbsrc.babraham.FastQC.Modules.QCModule;
import uk.ac.bbsrc.babraham.FastQC.Sequence.SequenceFile;

public class HTMLReportArchive {

	private StringBuffer html = new StringBuffer();
	private StringBuffer data = new StringBuffer();
	private QCModule [] modules;
	private ZipOutputStream zip;
	private SequenceFile sequenceFile;
	private byte [] buffer = new byte[1024];
	private File file;
	
	public HTMLReportArchive (SequenceFile sequenceFile, QCModule [] modules, File file) throws IOException {
		this.sequenceFile = sequenceFile;
		this.modules = modules;
		this.file = file;
		zip = new ZipOutputStream(new FileOutputStream(file));
		zip.putNextEntry(new ZipEntry(folderName()+"/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Icons/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Images/"));
		startDocument();
		for (int m=0;m<modules.length;m++) {
			html.append("<div class=\"module\"><h2 id=\"M");
			html.append(m);
			html.append("\">");
			
			// Add an icon before the module name
			if (modules[m].raisesError()) {
				html.append("<img src=\"Icons/error.png\" alt=\"[FAIL]\"> ");
			}
			else if (modules[m].raisesWarning()) {
				html.append("<img src=\"Icons/warning.png\" alt=\"[WARN]\"> ");
			}
			else {
				html.append("<img src=\"Icons/tick.png\" alt=\"[OK]\"> ");
			}

			
			
			html.append(modules[m].name());
			data.append(">>");
			data.append(modules[m].name());
			data.append("\t");
			if (modules[m].raisesError()) {
				data.append("fail");
			}
			else if (modules[m].raisesWarning()) {
				data.append("warn");
			}
			else {
				data.append("pass");
			}
			data.append("\n");
			html.append("</h2>\n");
			modules[m].makeReport(this);
			data.append(">>END_MODULE\n");
			
			html.append("</div>\n");
		}
		closeDocument();
		
		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_report.html"));
		zip.write(html.toString().getBytes());
		zip.closeEntry();
		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_data.txt"));
		zip.write(data.toString().getBytes());
		zip.closeEntry();
		zip.close();
		
		if (System.getProperty("fastqc.unzip").equals("true")) {
			unzipZipFile(file);
		}
		
		
	}
	
	private void unzipZipFile (File file) throws IOException {
		ZipFile zipFile = new ZipFile(file);
		Enumeration<? extends ZipEntry> entries = zipFile.entries();
		int size;
		byte [] buffer = new byte[1024];
		
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			
//			System.out.println("Going to extract '"+entry.getName()+"'");
			
			if (entry.isDirectory()) {
				File dir = new File(file.getParent()+"/"+entry.getName());
				if (dir.exists() && dir.isDirectory()) continue; // Don't need to do anything
				if (dir.exists() && ! dir.isDirectory()) throw new IOException ("File exists with dir name "+dir.getName());
				if (!dir.mkdir()) throw new IOException("Failed to make dir for "+dir.getName());
				continue;
			}

			BufferedInputStream bis = new BufferedInputStream(zipFile.getInputStream(entry));
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(file.getParent()+"/"+entry.getName()),buffer.length);
			while ((size = bis.read(buffer,0,buffer.length)) != -1) {
				bos.write(buffer,0,size);
			}
			bos.flush();
			bos.close();
			bis.close();
		}
	}

	public StringBuffer htmlDocument () {
		return html;
	}
	
	public StringBuffer dataDocument() {
		return data;
	}
	
	public String folderName () {
		return file.getName().replaceAll(".zip$", "");
	}
	
	public ZipOutputStream zipFile () {
		return zip;
	}
	
	private void startDocument () throws IOException {
		
		// Just put the fastQC version at the start of the text report
		data.append("##FastQC\t");
		data.append(FastQCApplication.VERSION);
		data.append("\n");
		
		// Add in the icon files for pass/fail/warn
		File templatesDir = new File(ClassLoader.getSystemResource("Templates/Icons").getFile().toString().replaceAll("%20", " "));
		String [] names = templatesDir.list();
		for (int n=0;n<names.length;n++) {
			if (names[n].toLowerCase().endsWith(".png") || names[n].toLowerCase().endsWith(".jpg") || names[n].toLowerCase().endsWith(".jpeg")) {
				zip.putNextEntry(new ZipEntry(folderName()+"/Icons/"+names[n]));
				FileInputStream fileIn = new FileInputStream(new File(ClassLoader.getSystemResource("Templates/Icons/"+names[n]).getFile().toString().replaceAll("%20", " ")));
				int len;
				while ((len = fileIn.read(buffer)) > 0) { 
					zip.write(buffer, 0, len); 
				} 
				fileIn.close();
			}
		}		
		

		SimpleDateFormat df = new SimpleDateFormat("EEE d MMM yyyy");
		addTemplate(sequenceFile.name(),df.format(new Date()));

		html.append("<h2>Summary</h2>\n<ul>\n");
		
		StringBuffer summaryText = new StringBuffer();
		
		for (int m=0;m<modules.length;m++) {
			html.append("<li>");
			html.append("<img src=\"");
			if (modules[m].raisesError()) {
				html.append("Icons/error.png\" alt=\"[FAIL]\"> ");
				summaryText.append("FAIL");
			}
			else if (modules[m].raisesWarning()) {
				html.append("Icons/warning.png\" alt=\"[WARNING]\"> ");				
				summaryText.append("WARN");
			}
			else {
				html.append("Icons/tick.png\" alt=\"[PASS]\"> ");
				summaryText.append("PASS");
			}
			summaryText.append("\t");
			summaryText.append(modules[m].name());
			summaryText.append("\t");
			summaryText.append(sequenceFile.name());
			summaryText.append(System.getProperty("line.separator"));
			
			html.append("<a href=\"#M");
			html.append(m);
			html.append("\">");
			html.append(modules[m].name());
			html.append("</a></li>\n");
			
		}
		html.append("</ul>\n</div>\n<div class=\"main\">\n");

		zip.putNextEntry(new ZipEntry(folderName()+"/summary.txt"));
		zip.write(summaryText.toString().getBytes());

	}
	
	private void closeDocument () {
		
		html.append("</div><div class=\"footer\">Produced by <a href=\"http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/\">FastQC</a> (version ");
		html.append(FastQCApplication.VERSION);
		html.append(")</div>\n");
		
		html.append("</body></html>");
	}
	
	private void addTemplate (String filename, String date) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(new File(ClassLoader.getSystemResource("Templates/header_template.html").getFile().toString().replaceAll("%20", " "))));
		String line;
		while ((line = br.readLine())!=null) {
			
			line = line.replaceAll("@@FILENAME@@", filename);
			line = line.replaceAll("@@DATE@@", date);
			
			html.append(line);
			html.append("\n");
		}
		
		
		br.close();
	
	}
	
	
	
}
