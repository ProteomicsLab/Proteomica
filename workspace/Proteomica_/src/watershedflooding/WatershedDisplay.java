package watershedflooding;

import ij.ImagePlus;
import imageware.Builder;
import imageware.Display;
import imageware.ImageWare;

public class WatershedDisplay
{
  public static byte[][] createLUTColor(int MAXBASINS)
  {
    byte[][] LUTColor = new byte[3][MAXBASINS];
    for (int c = 0; c < 3; c++) {
      for (int n = 0; n < 10; n++)
        LUTColor[c][n] = 0;
      for (int n = 10; n < MAXBASINS; n++)
        LUTColor[c][n] = ((byte)(int)(10.0D + 235.0D * Math.random()));
    }
    return LUTColor;
  }

  public static void showBasins(ImageWare output, byte[][] LUTColor)
  {
    int nx = output.getSizeX();
    int ny = output.getSizeY();
    int nz = output.getSizeZ();
    int size = nx * ny;
    ImageWare r = Builder.create(nx, ny, nz, 1);
    ImageWare g = Builder.create(nx, ny, nz, 1);
    ImageWare b = Builder.create(nx, ny, nz, 1);

    for (int z = 0; z < nz; z++) {
      byte[] pixr = r.getSliceByte(z);
      byte[] pixg = g.getSliceByte(z);
      byte[] pixb = b.getSliceByte(z);
      float[] pixout = output.getSliceFloat(z);
      for (int index = 0; index < size; index++) {
        int v = (int)pixout[index];
        pixr[index] = LUTColor[0][v];
        pixg[index] = LUTColor[1][v];
        pixb[index] = LUTColor[2][v];
      }
    }
    Display.showColor("Colorized Catchment Basins", r, g, b);
  }

  public static void showBinary(ImageWare output)
  {
    int nx = output.getSizeX();
    int ny = output.getSizeY();
    int nz = output.getSizeZ();
    ImageWare out = Builder.create(nx, ny, nz, 1);

    for (int z = 0; z < nz; z++)
      for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
          if (output.getPixel(x, y, z) < 2.0D)
            out.putPixel(x, y, z, 255.0D);
          if (output.getPixel(x - 1, y, z) < 2.0D)
            out.putPixel(x, y - 1, z, 255.0D);
          if (output.getPixel(x, y, z) < 2.0D)
            out.putPixel(x + 1, y, z, 255.0D);
          if (output.getPixel(x, y + 1, z) < 2.0D)
            out.putPixel(x, y, z, 255.0D);
        }
    out.show("Binary watershed lines");
  }

  public static void showDams(ImageWare output)
  {
    int nx = output.getSizeX();
    int ny = output.getSizeY();
    int nz = output.getSizeZ();
    ImageWare out = Builder.create(nx, ny, nz, 1);
    out.fillConstant(255.0D);
    for (int z = 0; z < nz; z++)
      for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++)
          if (output.getPixel(x, y, z) == 0.0D)
            out.putPixel(x, y, z, 0.0D);
    out.show("Watershed Lines");
  }

  public static void showRedDams(ImageWare output, ImagePlus first)
  {
    ImageWare firstData = Builder.create(first);
    int size = output.getSizeX() * output.getSizeY();
    int nz = output.getSizeZ();
    ImageWare original = firstData.convert(1);
    ImageWare over = original.duplicate();
    byte b255 = -1;
    byte b0 = 0;
    for (int z = 0; z < nz; z++) {
      float[] pixOutput = output.getSliceFloat(z);
      byte[] pixOver = over.getSliceByte(z);
      byte[] pixOriginal = original.getSliceByte(z);
      for (int index = 0; index < size; index++)
        if (pixOutput[index] == 0.0D) {
          pixOver[index] = b255;
          pixOriginal[index] = b0;
        }
    }
    Display.showColor("Dams", over, original, original);
  }

  public static void showComposite(ImageWare output, ImagePlus first)
  {
    ImageWare firstData = Builder.wrap(first);
    int nx = output.getSizeX();
    int ny = output.getSizeY();
    int nz = output.getSizeZ();
    int size = nx * ny;
    ImageWare combine = Builder.create(nx, ny, nz, 3);
    byte b255 = -1;
    byte b0 = 0;

    for (int z = 0; z < nz; z++) {
      float[] pixOutput = output.getSliceFloat(z);
      byte[] pixFirst = firstData.getSliceByte(z);
      float[] pixCombine = combine.getSliceFloat(z);
      for (int index = 0; index < size; index++)
        if (pixOutput[index] == 0.0D)
          pixCombine[index] = ((pixFirst[index] & 0xFF) / 10000.0F);
        else
          pixCombine[index] = ((pixFirst[index] & 0xFF) + pixOutput[index] / 10000.0F);
    }
    combine.show("Composite");
  }
}