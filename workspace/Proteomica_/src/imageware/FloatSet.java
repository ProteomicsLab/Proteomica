package imageware;

import ij.ImagePlus;
import ij.ImageStack;
import java.awt.Image;
import java.io.PrintStream;

public class FloatSet extends FloatProcess
  implements ImageWare
{
  protected FloatSet(int nx, int ny, int nz)
  {
    super(nx, ny, nz); } 
  protected FloatSet(Image image, int mode) { super(image, mode); } 
  protected FloatSet(ImageStack stack, int mode) {
    super(stack, mode); } 
  protected FloatSet(ImageStack stack, byte chan) { super(stack, chan); } 
  protected FloatSet(byte[] array, int mode) {
    super(array, mode); } 
  protected FloatSet(byte[][] array, int mode) { super(array, mode); } 
  protected FloatSet(byte[][][] array, int mode) { super(array, mode); } 
  protected FloatSet(short[] array, int mode) { super(array, mode); } 
  protected FloatSet(short[][] array, int mode) { super(array, mode); } 
  protected FloatSet(short[][][] array, int mode) { super(array, mode); } 
  protected FloatSet(float[] array, int mode) { super(array, mode); } 
  protected FloatSet(float[][] array, int mode) { super(array, mode); } 
  protected FloatSet(float[][][] array, int mode) { super(array, mode); } 
  protected FloatSet(double[] array, int mode) { super(array, mode); } 
  protected FloatSet(double[][] array, int mode) { super(array, mode); } 
  protected FloatSet(double[][][] array, int mode) { super(array, mode); }


  public ImageWare duplicate()
  {
    ImageWare out = new FloatSet(this.nx, this.ny, this.nz);

    for (int z = 0; z < this.nz; z++) {
      float[] outdata = (float[])((FloatSet)out).data[z];
      System.arraycopy(this.data[z], 0, outdata, 0, this.nxy);
    }
    return out;
  }

  public ImageWare replicate()
  {
    return new FloatSet(this.nx, this.ny, this.nz);
  }

  public ImageWare replicate(int type)
  {
    switch (type) {
    case 1:
      return new ByteSet(this.nx, this.ny, this.nz);
    case 2:
      return new ShortSet(this.nx, this.ny, this.nz);
    case 3:
      return new FloatSet(this.nx, this.ny, this.nz);
    case 4:
      return new DoubleSet(this.nx, this.ny, this.nz);
    }
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + "].\n" + "-------------------------------------------------------\n");
  }

  public void copy(ImageWare source)
  {
    if (this.nx != source.getSizeX()) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to copy because it is not the same size (" + this.nx + " != " + source.getSizeX() + ").\n" + "-------------------------------------------------------\n");
    }

    if (this.ny != source.getSizeY()) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to copy because it is not the same size (" + this.ny + " != " + source.getSizeY() + ").\n" + "-------------------------------------------------------\n");
    }

    if (this.nz != source.getSizeZ()) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to copy because it is not the same size (" + this.nz + " != " + source.getSizeZ() + ").\n" + "-------------------------------------------------------\n");
    }

    if (getType() != source.getType()) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to copy because it is not the same type (" + getType() + " != " + source.getType() + ").\n" + "-------------------------------------------------------\n");
    }

    for (int z = 0; z < this.nz; z++) {
      float[] src = (float[])((FloatSet)source).data[z];
      System.arraycopy(src, 0, this.data[z], 0, this.nxy);
    }
  }

  public ImageWare convert(int type)
  {
    if (type == 3)
      return duplicate();
    ImageWare out = null;
    switch (type)
    {
    case 1:
      out = new ByteSet(this.nx, this.ny, this.nz);

      for (int z = 0; z < this.nz; z++) {
        float[] slice = (float[])this.data[z];
        byte[] outslice = (byte[])((ByteSet)out).data[z];
        for (int k = 0; k < this.nxy; k++) {
          outslice[k] = ((byte)(int)slice[k]);
        }
      }

      break;
    case 2:
      out = new ShortSet(this.nx, this.ny, this.nz);

      for (int z = 0; z < this.nz; z++) {
        float[] slice = (float[])this.data[z];
        short[] outslice = (short[])((ShortSet)out).data[z];
        for (int k = 0; k < this.nxy; k++) {
          outslice[k] = ((short)(int)slice[k]);
        }
      }

      break;
    case 3:
      out = new FloatSet(this.nx, this.ny, this.nz);

      for (int z = 0; z < this.nz; z++) {
        float[] slice = (float[])this.data[z];
        float[] outslice = (float[])((FloatSet)out).data[z];
        for (int k = 0; k < this.nxy; k++) {
          outslice[k] = slice[k];
        }
      }

      break;
    case 4:
      out = new DoubleSet(this.nx, this.ny, this.nz);

      for (int z = 0; z < this.nz; z++) {
        float[] slice = (float[])this.data[z];
        double[] outslice = (double[])((DoubleSet)out).data[z];
        for (int k = 0; k < this.nxy; k++) {
          outslice[k] = slice[k];
        }
      }

      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + "].\n" + "-------------------------------------------------------\n");
    }

    return out;
  }

  public void printInfo()
  {
    System.out.println("ImageWare object information");
    System.out.println("Dimension: " + getDimension());
    System.out.println("Size: [" + this.nx + ", " + this.ny + ", " + this.nz + "]");
    System.out.println("TotalSize: " + getTotalSize());
    System.out.println("Type: " + getTypeToString());
    System.out.println("Maximun: " + getMaximum());
    System.out.println("Minimun: " + getMinimum());
    System.out.println("Mean: " + getMean());
    System.out.println("Norm1: " + getNorm1());
    System.out.println("Norm2: " + getNorm2());
    System.out.println("Total: " + getTotal());
    System.out.println("");
  }

  public void show()
  {
    String title = getTypeToString();
    switch (getDimension()) {
    case 1:
      title = title + " line";
      break;
    case 2:
      title = title + " image";
      break;
    case 3:
      title = title + " volume";
    }

    ImagePlus imp = new ImagePlus(title, buildImageStack());
    imp.show();
  }

  public void show(String title)
  {
    ImagePlus imp = new ImagePlus(title, buildImageStack());
    imp.show();
  }

  public double getMinimum()
  {
    double min = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        if (slice[k] < min)
          min = slice[k];
    }
    return min;
  }

  public double getMaximum()
  {
    double max = -1.79769313486231E+308D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        if (slice[k] > max)
          max = slice[k];
    }
    return max;
  }

  public double getMean()
  {
    return getTotal() / (this.nz * this.nxy);
  }

  public double getNorm1()
  {
    double norm = 0.0D;
    double value = 0.0D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        value = slice[k];
        norm += (value > 0.0D ? value : -value);
      }
    }
    return norm;
  }

  public double getNorm2()
  {
    double norm = 0.0D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        norm += slice[k] * slice[k];
    }
    return norm;
  }

  public double getTotal()
  {
    double total = 0.0D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        total += slice[k];
    }
    return total;
  }

  public double[] getMinMax()
  {
    double max = -1.79769313486231E+308D;
    double min = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        if (slice[k] > max)
          max = slice[k];
        if (slice[k] < min)
          min = slice[k];
      }
    }
    double[] minmax = { min, max };
    return minmax;
  }
}