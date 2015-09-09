package watershedflooding;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.text.TextPanel;
import ij.text.TextWindow;
import imageware.Builder;
import imageware.ImageWare;

public class WatershedMeasurements
{
  public static void measure()
  {
    ImagePlus imp = WindowManager.getCurrentImage();
    measure(imp);
  }

  public static void measure(ImagePlus imp) {
    TextWindow tw = new TextWindow("Watershed Measurements", "", 400, 400);
    String msg = "Slice\tBasins\tXG\tYG\t";
    msg = msg + "Minimum X\tMinimun Y\tMinimum Value\t";
    msg = msg + "Volume\tArea\tPerimeter\tMax\t";
    tw.getTextPanel().setColumnHeadings(msg);

    int nx = imp.getWidth();
    int ny = imp.getHeight();
    int nz = imp.getStack().getSize();
    ImageWare volume = Builder.wrap(imp.getStack());
    ImageWare combine = Builder.create(nx, ny, 1, volume.getType());

    for (int k = 0; k < nz; k++) {
      volume.getXY(0, 0, k, combine);
      ImageWare basins = Builder.create(nx, ny, 1, 3);
      ImageWare image = Builder.create(nx, ny, 1, 3);

      for (int x = 1; x < nx - 1; x++)
        for (int y = 1; y < ny - 1; y++) {
          double p = combine.getPixel(x, y, 0);
          if ((int)p == 0) {
            basins.putPixel(x, y, 0, 0.0D);
            image.putPixel(x, y, 0, Math.round(p * 10000.0D));
          }
          else {
            basins.putPixel(x, y, 0, Math.round((p - (int)p) * 10000.0D));
            image.putPixel(x, y, 0, (int)p);
          }
        }
      measureSlice(basins, image, tw, k + 1);
    }
  }

  private static void measureSlice(ImageWare basins, ImageWare image, TextWindow tw, int slice)
  {
    int nx = basins.getWidth();
    int ny = basins.getHeight();

    int maxClass = (int)basins.getMaximum();
    int minClass = (int)basins.getMinimum();
    int nb = maxClass - minClass + 1;
    double[] xg = new double[nb];
    double[] yg = new double[nb];
    int[] area = new int[nb];
    float[] volume = new float[nb];
    int[] perimeter = new int[nb];
    float[] min = new float[nb];
    int[] minX = new int[nb];
    int[] minY = new int[nb];
    float[] max = new float[nb];
    float[] barycenterX = new float[nb];
    float[] barycenterY = new float[nb];

    int[] status = new int[nb];

    for (int l = 0; l < nb; l++) {
      xg[l] = 0.0D;
      yg[l] = 0.0D;
      area[l] = 0;
      volume[l] = 0.0F;
      perimeter[l] = 0;
      min[l] = 3.4028235E+38F;
      minX[l] = 0;
      minY[l] = 0;
      max[l] = 1.4E-45F;
      barycenterX[l] = 0.0F;
      barycenterY[l] = 0.0F;
    }

    double[][] arr = new double[3][3];

    for (int x = 1; x < nx - 1; x++) {
      for (int y = 1; y < ny - 1; y++) {
        int v = (int)image.getPixel(x, y, 0);
        basins.getNeighborhoodXY(x, y, 0, arr, (byte)2);
        int l = (int)arr[1][1];
        xg[l] += x;
        yg[l] += y;
        area[l] += 1;
        volume[l] += 255 - v;
        barycenterX[l] += (255 - v) * x;
        barycenterY[l] += (255 - v) * y;
        if (v < min[l]) {
          min[l] = v;
          minX[l] = x;
          minY[l] = y;
        }
        if (v > max[l]) {
          max[l] = v;
        }
        for (int k = 0; k < 3; k++)
          for (int kk = 0; kk < 3; kk++)
            if (arr[k][kk] == 0.0D) {
              perimeter[l] += 1;
              break;
            }
      }
    }
    double[][] arr5 = new double[5][5];

    for (int l = 0; l < nb; l++) {
      if (area[l] > 0) {
        xg[l] /= area[l];
        yg[l] /= area[l];
      }

    }

    for (int l = 1; l < nb; l++)
      if (area[l] > 0) {
        String msg = "" + slice + "\t" + l + "\t" + IJ.d2s(xg[l]) + "\t" + IJ.d2s(yg[l]) + "\t";
        msg = msg + IJ.d2s(minX[l], 0) + "\t" + IJ.d2s(minY[l], 0) + "\t" + IJ.d2s(min[l], 0) + "\t";
        msg = msg + IJ.d2s(volume[l], 0) + "\t" + IJ.d2s(area[l], 0) + "\t" + IJ.d2s(perimeter[l], 0) + "\t" + IJ.d2s(max[l], 0) + "\t";
        tw.append(msg);
      }
  }
}