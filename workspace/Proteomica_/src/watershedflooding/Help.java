package watershedflooding;

import ij.IJ;
import ij.gui.GUI;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;

public class Help extends JFrame
  implements ActionListener, WindowListener
{
  private JButton bnClose = new JButton("Close");
  private GridBagLayout layout = new GridBagLayout();
  private GridBagConstraints constraint = new GridBagConstraints();

  public Help()
  {
    super("Help on Watershed");

    JTextArea help = new JTextArea(40, 50);
    help.append("Watershed on graylevel images\n\nhttp://bigwww.epfl.ch/sage/soft/watershed\n\nDaniel Sage\nBiomedical Imaging Group (BIG)\nEcole Polytechnique Fédérale de Lausane, Lausanne, Switzerland.\n\n20 February 2007\n\nGAUSSIAN BLUR\nThe watershed applies on graylevel images is known to oversegment the image\ndue to the noise of the natural images.\nTo avoid the oversegmentation, the noise should be reduce. One method to reduce the.\nnoise is to blur the image. We provide, here, an efficient gaussian blurring tuning with.\nthe parameter 'Radius'.\n\nWATERSHED\nSelect an image and click on 'Start Watershed' to begin the watershed process.\nThe program start to flooding from the 'min' value to the 'max' value.\nBy default, the program start to flood the dark area (in principle low values) and.\ncontinue to the brightest values (in principle high values) considering 4 neigborhood.\nThe default settings could be change by the user'.\nThe user can always stop the process by clicking on 'Stop Watershed'.\nthe parameter 'Radius'.\n\nDISPLAY\nAfter the watershed process, several views of the segmented image could be obtained.\nChoose a display and click on 'Show'.\nThe composite image is a mixture of graylevel and labels, in the basins the pixel value is.\nthe graylevel (integer part) and the label (decimal part) and over the dams the pixel value is.\n0 (integer part) and the graylevel (decimal part).\n");

    help.setForeground(new Color(0, 32, 128));
    help.setBackground(getBackground());

    JPanel panelHelp = new JPanel();
    panelHelp.setLayout(this.layout);
    addComponent(panelHelp, 0, 0, 8, 1, 9, help);
    panelHelp.setBorder(BorderFactory.createEtchedBorder());

    this.bnClose.addActionListener(this);
    addWindowListener(this);

    JPanel panel = new JPanel();
    panel.setLayout(this.layout);
    addComponent(panel, 0, 0, 8, 1, 9, panelHelp);
    addComponent(panel, 1, 0, 1, 1, 9, new JLabel("  "));
    addComponent(panel, 1, 4, 1, 1, 9, this.bnClose);
    addComponent(panel, 1, 5, 1, 1, 9, new JLabel("  "));

    setLayout(this.layout);
    JPanel pnMain = new JPanel();
    pnMain.setLayout(this.layout);
    addComponent(pnMain, 0, 0, 1, 1, 9, panel);
    add(pnMain);
    pack();
    setResizable(true);
    GUI.center(this);
    setVisible(true);
  }

  private void addComponent(JPanel pn, int row, int col, int width, int height, int space, JComponent comp)
  {
    this.constraint.gridx = col;
    this.constraint.gridy = row;
    this.constraint.gridwidth = width;
    this.constraint.gridheight = height;
    this.constraint.anchor = 18;
    this.constraint.insets = new Insets(space, space, space, space);
    this.constraint.weightx = (IJ.isMacintosh() ? 90.0D : 100.0D);
    this.constraint.fill = 2;
    this.layout.setConstraints(comp, this.constraint);
    pn.add(comp);
  }

  public synchronized void actionPerformed(ActionEvent e)
  {
    if (e.getSource() == this.bnClose) {
      dispose();
    }
    notify();
  }

  public void windowActivated(WindowEvent e)
  {
  }

  public void windowClosing(WindowEvent e)
  {
    dispose();
  }

  public void windowClosed(WindowEvent e)
  {
  }

  public void windowDeactivated(WindowEvent e)
  {
  }

  public void windowDeiconified(WindowEvent e)
  {
  }

  public void windowIconified(WindowEvent e)
  {
  }

  public void windowOpened(WindowEvent e)
  {
  }
}