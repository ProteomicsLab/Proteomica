package imageware;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.process.ColorProcessor;

public class Display
{
  public static void show(String title, ImageWare ds)
  {
    new ImagePlus(title, ds.buildImageStack()).show();
  }

  public static void showColor(String title, ImageWare red, ImageWare green, ImageWare blue)
  {
    new ImagePlus(title, buildColor(red, green, blue)).show();
  }

  public static void show(String title, ImageWare ds, double magnification)
  {
    ImagePlus imp = new ImagePlus(title, ds.buildImageStack());
    imp.show();
    ImageWindow win = imp.getWindow();
    ImageCanvas canvas = win.getCanvas();
    canvas.setMagnification(magnification);
    canvas.setDrawingSize((int)Math.ceil(ds.getWidth() * magnification), (int)Math.ceil(ds.getHeight() * magnification));

    win.pack();
    imp.updateAndRepaintWindow();
  }

  public static void showColor(String title, ImageWare red, ImageWare green, ImageWare blue, double magnification)
  {
    ImagePlus imp = new ImagePlus(title, buildColor(red, green, blue));
    imp.show();
    ImageWindow win = imp.getWindow();
    ImageCanvas canvas = win.getCanvas();
    canvas.setMagnification(magnification);
    canvas.setDrawingSize((int)Math.ceil(red.getWidth() * magnification), (int)Math.ceil(red.getHeight() * magnification));

    win.pack();
    imp.updateAndRepaintWindow();
  }

  private static ImageStack buildColor(ImageWare red, ImageWare green, ImageWare blue)
  {
    if (!red.isSameSize(green)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to create a ImageStack the channel are not the same size.\n[" + red.getSizeX() + "," + red.getSizeY() + "," + red.getSizeZ() + "] != " + "[" + green.getSizeX() + "," + green.getSizeY() + "," + green.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    if (!red.isSameSize(blue)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to create a ImageStack the channel are not the same size.\n[" + red.getSizeX() + "," + red.getSizeY() + "," + red.getSizeZ() + "] != " + "[" + blue.getSizeX() + "," + blue.getSizeY() + "," + blue.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    int nx = red.getSizeX();
    int ny = red.getSizeY();
    int nz = red.getSizeZ();
    int nxy = nx * ny;
    ImageStack imagestack = new ImageStack(nx, ny);

    byte[] r = new byte[nxy];
    byte[] g = new byte[nxy];
    byte[] b = new byte[nxy];
    for (int z = 0; z < nz; z++) {
      ColorProcessor cp = new ColorProcessor(nx, ny);
      switch (red.getType()) {
      case 4:
        double[] dpixred = red.getSliceDouble(z);
        for (int k = 0; k < nxy; k++)
          r[k] = ((byte)(int)dpixred[k]);
        break;
      case 3:
        float[] fpixred = red.getSliceFloat(z);
        for (int k = 0; k < nxy; k++)
          r[k] = ((byte)(int)fpixred[k]);
        break;
      case 2:
        short[] spixred = red.getSliceShort(z);
        for (int k = 0; k < nxy; k++)
          r[k] = ((byte)spixred[k]);
        break;
      case 1:
        r = red.getSliceByte(z);
      }

      switch (green.getType()) {
      case 4:
        double[] dpixgreen = green.getSliceDouble(z);
        for (int k = 0; k < nxy; k++)
          g[k] = ((byte)(int)dpixgreen[k]);
        break;
      case 3:
        float[] fpixgreen = green.getSliceFloat(z);
        for (int k = 0; k < nxy; k++)
          g[k] = ((byte)(int)fpixgreen[k]);
        break;
      case 2:
        short[] spixgreen = green.getSliceShort(z);
        for (int k = 0; k < nxy; k++)
          g[k] = ((byte)spixgreen[k]);
        break;
      case 1:
        g = green.getSliceByte(z);
      }

      switch (blue.getType()) {
      case 4:
        double[] dpixblue = blue.getSliceDouble(z);
        for (int k = 0; k < nxy; k++)
          b[k] = ((byte)(int)dpixblue[k]);
        break;
      case 3:
        float[] fpixblue = blue.getSliceFloat(z);
        for (int k = 0; k < nxy; k++)
          b[k] = ((byte)(int)fpixblue[k]);
        break;
      case 2:
        short[] spixblue = blue.getSliceShort(z);
        for (int k = 0; k < nxy; k++)
          b[k] = ((byte)spixblue[k]);
        break;
      case 1:
        b = blue.getSliceByte(z);
      }

      cp.setRGB(r, g, b);
      imagestack.addSlice("" + z, cp);
    }
    return imagestack;
  }
}