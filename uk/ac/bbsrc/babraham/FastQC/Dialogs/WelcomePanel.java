package uk.ac.bbsrc.babraham.FastQC.Dialogs;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.JLabel;
import javax.swing.JPanel;

public class WelcomePanel extends JPanel {

	public WelcomePanel () {
		setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		
		gbc.gridx=1;
		gbc.gridy=1;
		gbc.weightx=0.5;
		gbc.weighty=0.99;
		
		add(new JPanel(),gbc);
		gbc.gridy++;
		gbc.weighty=0.01;
		
		gbc.insets = new Insets(10, 10, 10, 10);
		gbc.fill = GridBagConstraints.NONE;
		
		add(new FastQCTitlePanel(),gbc);
		
		gbc.gridy++;
		gbc.weighty=0.5;
		
		add(new JLabel("Use File > Open to select the sequence file you want to check"),gbc);
		
		gbc.gridy++;
		gbc.weighty=0.99;
		add(new JPanel(),gbc);
		
	}
}
