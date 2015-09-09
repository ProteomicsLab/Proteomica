package imageware;

import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import java.awt.Point;

public class ImageAccess
{
  public static final int PATTERN_SQUARE_3x3 = 0;
  public static final int PATTERN_CROSS_3x3 = 1;
  private ImageWare imageware = null;
  private int nx = 0;
  private int ny = 0;

  public ImageAccess(double[][] array)
  {
    if (array == null) {
      throw new ArrayStoreException("Constructor: array == null.");
    }
    this.imageware = Builder.create(array);
    this.nx = this.imageware.getSizeX();
    this.ny = this.imageware.getSizeY();
  }

  public ImageAccess(ImageProcessor ip)
  {
    if (ip == null) {
      throw new ArrayStoreException("Constructor: ImageProcessor == null.");
    }
    ImagePlus imp = new ImagePlus("", ip);
    if ((ip instanceof ByteProcessor))
      this.imageware = Builder.create(imp, 4);
    else if ((ip instanceof ShortProcessor))
      this.imageware = Builder.create(imp, 4);
    else if ((ip instanceof FloatProcessor))
      this.imageware = Builder.create(imp, 4);
    this.nx = this.imageware.getSizeX();
    this.ny = this.imageware.getSizeY();
  }

  public ImageAccess(ColorProcessor cp, int colorPlane)
  {
    if (cp == null) {
      throw new ArrayStoreException("Constructor: ColorProcessor == null.");
    }
    if (colorPlane < 0) {
      throw new ArrayStoreException("Constructor: colorPlane < 0.");
    }
    if (colorPlane > 2) {
      throw new ArrayStoreException("Constructor: colorPlane > 2.");
    }
    this.nx = cp.getWidth();
    this.ny = cp.getHeight();
    int size = this.nx * this.ny;
    ImagePlus imp = new ImagePlus("", cp);
    this.imageware = new DoubleSet(imp.getStack(), (byte)colorPlane);
  }

  public ImageAccess(int nx, int ny)
  {
    this.imageware = new DoubleSet(nx, ny, 1);
    this.nx = nx;
    this.ny = ny;
  }

  public ImageWare getDataset()
  {
    return this.imageware;
  }

  public int getWidth()
  {
    return this.nx;
  }

  public int getHeight()
  {
    return this.ny;
  }

  public double getMaximum()
  {
    return this.imageware.getMaximum();
  }

  public double getMinimum()
  {
    return this.imageware.getMinimum();
  }

  public double getMean()
  {
    return this.imageware.getMean();
  }

  public double[][] getArrayPixels()
  {
    double[][] array = new double[this.nx][this.ny];
    this.imageware.getXY(0, 0, 0, array);
    return array;
  }

  public double[] getPixels()
  {
    return this.imageware.getSliceDouble(0);
  }

  public FloatProcessor createFloatProcessor()
  {
    FloatProcessor fp = new FloatProcessor(this.nx, this.ny);
    double[] pixels = getPixels();
    int size = pixels.length;
    float[] fsrc = new float[size];
    for (int k = 0; k < size; k++)
      fsrc[k] = ((float)pixels[k]);
    fp.setPixels(fsrc);
    return fp;
  }

  public ByteProcessor createByteProcessor()
  {
    ByteProcessor bp = new ByteProcessor(this.nx, this.ny);
    double[] pixels = getPixels();
    int size = pixels.length;
    byte[] bsrc = new byte[size];

    for (int k = 0; k < size; k++) {
      double p = pixels[k];
      if (p < 0.0D)
        p = 0.0D;
      if (p > 255.0D)
        p = 255.0D;
      bsrc[k] = ((byte)(int)p);
    }
    bp.setPixels(bsrc);
    return bp;
  }

  public ImageAccess duplicate()
  {
    double[][] array = new double[this.nx][this.ny];
    this.imageware.getXY(0, 0, 0, array);
    ImageAccess ia = new ImageAccess(array);
    return ia;
  }

  public double getPixel(int x, int y)
  {
    return this.imageware.getPixel(x, y, 0, (byte)2);
  }

  public double getInterpolatedPixel(double x, double y)
  {
    return this.imageware.getInterpolatedPixel(x, y, 0.0D, (byte)2);
  }

  public void getColumn(int x, double[] column)
  {
    if (x < 0) {
      throw new IndexOutOfBoundsException("getColumn: x < 0.");
    }
    if (x >= this.nx) {
      throw new IndexOutOfBoundsException("getColumn: x >= nx.");
    }
    if (column == null) {
      throw new ArrayStoreException("getColumn: column == null.");
    }

    if (column.length != this.ny) {
      throw new ArrayStoreException("getColumn: column.length != ny.");
    }
    this.imageware.getBlockY(x, 0, 0, column, (byte)2);
  }

  public void getColumn(int x, int y, double[] column)
  {
    if (x < 0) {
      throw new IndexOutOfBoundsException("getColumn: x < 0.");
    }
    if (x >= this.nx) {
      throw new IndexOutOfBoundsException("getColumn: x >= nx.");
    }
    if (column == null) {
      throw new ArrayStoreException("getColumn: column == null.");
    }
    this.imageware.getBlockY(x, y, 0, column, (byte)2);
  }

  public void getRow(int y, double[] row)
  {
    if (y < 0) {
      throw new IndexOutOfBoundsException("getRow: y < 0.");
    }
    if (y >= this.ny) {
      throw new IndexOutOfBoundsException("getRow: y >= ny.");
    }
    if (row == null) {
      throw new ArrayStoreException("getColumn: row == null.");
    }
    if (row.length != this.nx) {
      throw new ArrayStoreException("getColumn: row.length != nx.");
    }
    this.imageware.getBlockX(0, y, 0, row, (byte)2);
  }

  public void getRow(int x, int y, double[] row)
  {
    if (y < 0) {
      throw new IndexOutOfBoundsException("getRow: y < 0.");
    }
    if (y >= this.ny) {
      throw new IndexOutOfBoundsException("getRow: y >= ny.");
    }
    if (row == null) {
      throw new ArrayStoreException("getRow: row == null.");
    }
    this.imageware.getBlockX(x, y, 0, row, (byte)2);
  }

  public void getNeighborhood(int x, int y, double[][] neigh)
  {
    this.imageware.getNeighborhoodXY(x, y, 0, neigh, (byte)2);
  }

  public void getPattern(int x, int y, double[] neigh, int pattern)
  {
    if (neigh == null) {
      throw new ArrayStoreException("getPattern: neigh == null.");
    }

    double[][] block = new double[3][3];
    this.imageware.getNeighborhoodXY(x, y, 0, block, (byte)2);

    switch (pattern) {
    case 0:
      if (neigh.length != 9) {
        throw new ArrayStoreException("getPattern: neigh.length != 9.");
      }
      neigh[0] = block[0][0];
      neigh[1] = block[1][0];
      neigh[2] = block[2][0];
      neigh[3] = block[0][1];
      neigh[4] = block[1][1];
      neigh[5] = block[2][1];
      neigh[6] = block[0][2];
      neigh[7] = block[1][2];
      neigh[8] = block[2][2];
      break;
    case 1:
      if (neigh.length != 5) {
        throw new ArrayStoreException("getPattern: neigh.length != 5");
      }
      neigh[0] = block[1][0];
      neigh[1] = block[0][1];
      neigh[2] = block[1][1];
      neigh[3] = block[2][1];
      neigh[4] = block[0][1];
      break;
    default:
      throw new ArrayStoreException("getPattern: unexpected pattern.");
    }
  }

  public void getSubImage(int x, int y, ImageAccess output)
  {
    if (output == null) {
      throw new ArrayStoreException("getSubImage: output == null.");
    }
    if (x < 0) {
      throw new ArrayStoreException("getSubImage: Incompatible image size");
    }
    if (y < 0) {
      throw new ArrayStoreException("getSubImage: Incompatible image size");
    }
    if (x >= this.nx) {
      throw new ArrayStoreException("getSubImage: Incompatible image size");
    }
    if (y >= this.ny) {
      throw new ArrayStoreException("getSubImage: Incompatible image size");
    }
    int nxcopy = output.getWidth();
    int nycopy = output.getHeight();
    double[][] neigh = new double[nxcopy][nycopy];
    this.imageware.getBlockXY(x, y, 0, neigh, (byte)2);
    output.putArrayPixels(neigh);
  }

  public void putPixel(int x, int y, double value)
  {
    if (x < 0) {
      throw new IndexOutOfBoundsException("putPixel: x < 0");
    }
    if (x >= this.nx) {
      throw new IndexOutOfBoundsException("putPixel: x >= nx");
    }
    if (y < 0) {
      throw new IndexOutOfBoundsException("putPixel:  y < 0");
    }
    if (y >= this.ny) {
      throw new IndexOutOfBoundsException("putPixel:  y >= ny");
    }
    this.imageware.putPixel(x, y, 0, value);
  }

  public void putColumn(int x, double[] column)
  {
    if (x < 0) {
      throw new IndexOutOfBoundsException("putColumn: x < 0.");
    }
    if (x >= this.nx) {
      throw new IndexOutOfBoundsException("putColumn: x >= nx.");
    }
    if (column == null) {
      throw new ArrayStoreException("putColumn: column == null.");
    }
    if (column.length != this.ny) {
      throw new ArrayStoreException("putColumn: column.length != ny.");
    }
    this.imageware.putBoundedY(x, 0, 0, column);
  }

  public void putColumn(int x, int y, double[] column)
  {
    if (x < 0) {
      throw new IndexOutOfBoundsException("putColumn: x < 0.");
    }
    if (x >= this.nx) {
      throw new IndexOutOfBoundsException("putColumn: x >= nx.");
    }
    if (column == null) {
      throw new ArrayStoreException("putColumn: column == null.");
    }
    this.imageware.putBoundedY(x, y, 0, column);
  }

  public void putRow(int y, double[] row)
  {
    if (y < 0) {
      throw new IndexOutOfBoundsException("putRow: y < 0.");
    }
    if (y >= this.ny) {
      throw new IndexOutOfBoundsException("putRow: y >= ny.");
    }
    if (row == null) {
      throw new ArrayStoreException("putRow: row == null.");
    }
    if (row.length != this.nx) {
      throw new ArrayStoreException("putRow: row.length != nx.");
    }
    this.imageware.putBoundedX(0, y, 0, row);
  }

  public void putRow(int x, int y, double[] row)
  {
    if (y < 0) {
      throw new IndexOutOfBoundsException("putRow: y < 0.");
    }
    if (y >= this.ny) {
      throw new IndexOutOfBoundsException("putRow: y >= ny.");
    }
    if (row == null) {
      throw new ArrayStoreException("putRow: row == null.");
    }
    this.imageware.putBoundedX(x, y, 0, row);
  }

  public void putArrayPixels(double[][] array)
  {
    if (array == null) {
      throw new IndexOutOfBoundsException("putArrayPixels: array == null.");
    }
    this.imageware.putBoundedXY(0, 0, 0, array);
  }

  public void putSubImage(int x, int y, ImageAccess input)
  {
    if (input == null) {
      throw new ArrayStoreException("putSubImage: input == null.");
    }
    if (x < 0) {
      throw new IndexOutOfBoundsException("putSubImage: x < 0.");
    }
    if (y < 0) {
      throw new IndexOutOfBoundsException("putSubImage: y < 0.");
    }
    if (x >= this.nx) {
      throw new IndexOutOfBoundsException("putSubImage: x >= nx.");
    }
    if (y >= this.ny) {
      throw new IndexOutOfBoundsException("putSubImage: y >= ny.");
    }

    int nxcopy = input.getWidth();
    int nycopy = input.getHeight();
    double[][] sub = input.getArrayPixels();
    this.imageware.putBoundedXY(x, y, 0, sub);
  }

  public void setConstant(double constant)
  {
    this.imageware.fillConstant(constant);
  }

  public void normalizeContrast()
  {
    this.imageware.rescale();
  }

  public void show(String title, Point loc)
  {
    FloatProcessor fp = createFloatProcessor();
    fp.resetMinAndMax();
    ImagePlus impResult = new ImagePlus(title, fp);
    ImageWindow window = impResult.getWindow();
    window.setLocation(loc.x, loc.y);
    impResult.show();
  }

  public void show(String title)
  {
    this.imageware.show(title);
  }

  public void abs()
  {
    this.imageware.abs();
  }

  public void sqrt()
  {
    this.imageware.sqrt();
  }

  public void pow(double a)
  {
    this.imageware.pow(a);
  }

  public void add(double constant)
  {
    this.imageware.add(constant);
  }

  public void multiply(double constant)
  {
    this.imageware.multiply(constant);
  }

  public void subtract(double constant)
  {
    this.imageware.add(-constant);
  }

  public void divide(double constant)
  {
    if (constant == 0.0D) {
      throw new ArrayStoreException("divide: Divide by 0");
    }
    this.imageware.multiply(1.0D / constant);
  }

  public void add(ImageAccess im1, ImageAccess im2)
  {
    if (im1.getWidth() != this.nx) {
      throw new ArrayStoreException("add: incompatible size.");
    }

    if (im1.getHeight() != this.ny) {
      throw new ArrayStoreException("add: incompatible size.");
    }

    if (im2.getWidth() != this.nx) {
      throw new ArrayStoreException("add: incompatible size.");
    }

    if (im2.getHeight() != this.ny) {
      throw new ArrayStoreException("add: incompatible size.");
    }
    this.imageware.copy(im1.getDataset());
    this.imageware.add(im2.getDataset());
  }

  public void multiply(ImageAccess im1, ImageAccess im2)
  {
    if (im1.getWidth() != this.nx) {
      throw new ArrayStoreException("multiply: incompatible size.");
    }

    if (im1.getHeight() != this.ny) {
      throw new ArrayStoreException("multiply: incompatible size.");
    }

    if (im2.getWidth() != this.nx) {
      throw new ArrayStoreException("multiply: incompatible size.");
    }

    if (im2.getHeight() != this.ny) {
      throw new ArrayStoreException("multiply: incompatible size.");
    }

    this.imageware.copy(im1.getDataset());
    this.imageware.multiply(im2.getDataset());
  }

  public void subtract(ImageAccess im1, ImageAccess im2)
  {
    if (im1.getWidth() != this.nx) {
      throw new ArrayStoreException("subtract: incompatible size.");
    }

    if (im1.getHeight() != this.ny) {
      throw new ArrayStoreException("subtract: incompatible size.");
    }

    if (im2.getWidth() != this.nx) {
      throw new ArrayStoreException("subtract: incompatible size.");
    }

    if (im2.getHeight() != this.ny) {
      throw new ArrayStoreException("subtract: incompatible size.");
    }

    this.imageware.copy(im1.getDataset());
    this.imageware.subtract(im2.getDataset());
  }

  public void divide(ImageAccess im1, ImageAccess im2)
  {
    if (im1.getWidth() != this.nx) {
      throw new ArrayStoreException("divide: incompatible size.");
    }

    if (im1.getHeight() != this.ny) {
      throw new ArrayStoreException("divide: incompatible size.");
    }

    if (im2.getWidth() != this.nx) {
      throw new ArrayStoreException("divide: incompatible size.");
    }

    if (im2.getHeight() != this.ny) {
      throw new ArrayStoreException("divide: incompatible size.");
    }

    this.imageware.copy(im1.getDataset());
    this.imageware.divide(im2.getDataset());
  }
}