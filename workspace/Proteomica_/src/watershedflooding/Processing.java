package watershedflooding;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.WindowManager;
import imageware.Builder;
import imageware.ImageWare;
import java.awt.Cursor;
import javax.swing.JButton;

public class Processing
  implements Runnable
{
  private int MAXBASINS = 200000;
  private Thread thread = null;
  private final int JOB_SMOOTH = 1;
  private final int JOB_WATERSHED = 2;
  private final int JOB_DISPLAY = 3;
  private String job = "";

  private double sigmaXY = 0.0D;
  private double sigmaZ = 0.0D;
  private boolean animation = false;
  private boolean conn4 = false;
  private boolean progression = false;
  private int minLevel = 0;
  private int maxLevel = 255;
  private int[] display;
  private JButton bnWatershed;
  private ImageWare output = null;
  private ImageWare image = null;
  private ImagePlus imp;

  public Processing(ImagePlus imp)
  {
    this.imp = imp;
  }

  public void paramsSmooth(double sigmaXY, double sigmaZ)
  {
    this.sigmaXY = sigmaXY;
    this.sigmaZ = sigmaZ;
  }

  public void paramsWatershed(ImageWare image, JButton bnWatershed, boolean animation, boolean conn4, boolean progression, int minLevel, int maxLevel)
  {
    this.image = image;
    this.bnWatershed = bnWatershed;
    this.animation = animation;
    this.conn4 = conn4;
    this.progression = progression;
    this.minLevel = minLevel;
    this.maxLevel = maxLevel;
  }

  public void paramsDisplay(int[] display)
  {
    this.display = display;
  }

  public void start(String job)
  {
    if (this.thread != null) {
      return;
    }
    this.job = job.toLowerCase();
    this.thread = new Thread(this);
    this.thread.setPriority(1);
    this.thread.start();
  }

  public void run()
  {
    Cursor cursor = IJ.getInstance().getCursor();
    IJ.getInstance().setCursor(new Cursor(3));

    if (this.job.equals("smooth")) {
      smooth();
    }
    else if (this.job.equals("start watershed")) {
      watershed();
    }
    else if (this.job.equals("display")) {
      display();
    }
    this.thread = null;
    IJ.getInstance().setCursor(cursor);
  }

  private void smooth()
  {
    ImagePlus imp = WindowManager.getCurrentImage();

    if (imp == null) {
      IJ.showMessage("Image required.");
      return;
    }

    if (imp.getType() != 0) {
      IJ.showMessage("8-bit image required.");
      return;
    }

    double step = 1.0D / imp.getStackSize();
    ImageWare original = Builder.create(imp);
    ImageWare smoothImage = original.duplicate();
    smoothImage.smoothGaussian(this.sigmaXY, this.sigmaXY, this.sigmaZ);
    smoothImage.show("Smooth");
  }

  private void watershed()
  {
    double t = System.currentTimeMillis();

    int nx = this.image.getSizeX();
    int ny = this.image.getSizeY();
    int nz = this.image.getSizeZ();
    ImageWare inputSlice = Builder.create(nx, ny, 1, this.image.getType());
    this.output = Builder.create(nx, ny, nz, 3);

    if (this.bnWatershed != null)
      this.bnWatershed.setText("Stop Watershed");
    Watershed ws = new Watershed(this.progression);

    if (this.animation) {
      ws.enableAnimation();
    }

    for (int z = 0; z < nz; z++) {
      if ((this.progression) && (nz > 1))
        IJ.write(">>>> Start slice: " + z);
      this.image.getXY(0, 0, z, inputSlice);
      ImageWare outputSlice = ws.doWatershed(inputSlice, this.conn4, this.minLevel, this.maxLevel);
      this.output.putXY(0, 0, z, outputSlice);
      if ((this.progression) && (nz > 1))
        IJ.write(">>>> End slice: " + z);
    }
    if (this.bnWatershed != null) {
      this.bnWatershed.setText("Start Watershed");
    }
    if (this.progression) {
      IJ.write("Watershed time: " + IJ.d2s(System.currentTimeMillis() - t) + " ms");
    }
    display();
  }

  private void display()
  {
    if (this.output == null)
      return;
    for (int k = 0; k < this.display.length; k++)
      switch (this.display[k]) {
      case 0:
        WatershedDisplay.showBinary(this.output);
        break;
      case 1:
        WatershedDisplay.showDams(this.output);
        break;
      case 2:
        WatershedDisplay.showRedDams(this.output, this.imp);
        break;
      case 3:
        this.output.show("Catchment Basins");
        break;
      case 4:
        WatershedDisplay.showBasins(this.output, WatershedDisplay.createLUTColor(this.MAXBASINS));
        break;
      case 5:
        WatershedDisplay.showComposite(this.output, this.imp);
        break;
      case 6:
        this.image.show("Watershed input");
      }
  }
}