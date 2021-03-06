package imageware;

import ij.ImageStack;
import ij.process.FloatProcessor;
import java.awt.Image;
import java.util.Random;

public class DoublePointwise extends DoubleAccess
  implements Pointwise
{
  protected DoublePointwise(int nx, int ny, int nz)
  {
    super(nx, ny, nz); } 
  protected DoublePointwise(Image image, int mode) { super(image, mode); } 
  protected DoublePointwise(ImageStack stack, int mode) {
    super(stack, mode); } 
  protected DoublePointwise(ImageStack stack, byte chan) { super(stack, chan); } 
  protected DoublePointwise(byte[] array, int mode) {
    super(array, mode); } 
  protected DoublePointwise(byte[][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(byte[][][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(short[] array, int mode) { super(array, mode); } 
  protected DoublePointwise(short[][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(short[][][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(float[] array, int mode) { super(array, mode); } 
  protected DoublePointwise(float[][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(float[][][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(double[] array, int mode) { super(array, mode); } 
  protected DoublePointwise(double[][] array, int mode) { super(array, mode); } 
  protected DoublePointwise(double[][][] array, int mode) { super(array, mode); }


  public void fillConstant(double value)
  {
    double typedValue = value;
    double[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = typedValue;
    }
  }

  public void fillRamp()
  {
    int off = 0;
    double[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (off + k);
      off += this.nxy;
    }
  }

  public void fillGaussianNoise(double amplitude)
  {
    Random rnd = new Random();
    double[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (rnd.nextGaussian() * amplitude);
    }
  }

  public void fillUniformNoise(double amplitude)
  {
    Random rnd = new Random();
    double[] slice = null;
    amplitude *= 2.0D;
    for (int z = 0; z < this.nz; z++) {
      slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((rnd.nextDouble() - 0.5D) * amplitude);
    }
  }

  public void fillSaltPepper(double amplitudeSalt, double amplitudePepper, double percentageSalt, double percentagePepper)
  {
    Random rnd = new Random();
    double[] slice = null;

    if (percentageSalt > 0.0D) {
      double nbSalt = this.nxy * this.nz / percentageSalt;
      for (int k = 0; k < nbSalt; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((double[])this.data[z])[index] = (rnd.nextDouble() * amplitudeSalt);
      }
    }
    if (percentagePepper > 0.0D) {
      double nbPepper = this.nxy * this.nz / percentagePepper;
      for (int k = 0; k < nbPepper; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((double[])this.data[z])[index] = (-rnd.nextDouble() * amplitudeSalt);
      }
    }
  }

  public ImageStack buildImageStack()
  {
    ImageStack imagestack = new ImageStack(this.nx, this.ny);
    for (int z = 0; z < this.nz; z++)
    {
      FloatProcessor ip = new FloatProcessor(this.nx, this.ny);
      float[] pix = (float[])ip.getPixels();
      for (int k = 0; k < this.nxy; k++)
        pix[k] = ((float)((double[])(double[])this.data[z])[k]);
      imagestack.addSlice("" + z, ip);
    }
    return imagestack;
  }

  public void invert()
  {
    double max = -1.79769313486231E+308D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        if (slice[k] > max)
          max = slice[k];
      }
    }
    double cst = max;
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (max - slice[k]);
    }
  }

  public void negate()
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (-slice[k]);
    }
  }

  public void clip()
  {
    clip(0.0D, 255.0D);
  }

  public void clip(double minLevel, double maxLevel)
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];

      double min = minLevel;
      double max = maxLevel;
      for (int k = 0; k < this.nxy; k++) {
        double value = slice[k];
        if (value < min)
          slice[k] = min;
        if (value > max)
          slice[k] = max;
      }
    }
  }

  public void rescale()
  {
    double maxImage = -1.79769313486231E+308D;
    double minImage = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        if (slice[k] > maxImage)
          maxImage = slice[k];
        if (slice[k] < minImage)
          minImage = slice[k];
      }
    }
    double a;
    if (minImage - maxImage == 0.0D) {
      a = 1.0D;
      minImage = 128.0D;
    }
    else {
      a = 255.0D / (maxImage - minImage);
    }
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (a * (slice[k] - minImage));
    }
  }

  public void rescale(double minLevel, double maxLevel)
  {
    double maxImage = -1.79769313486231E+308D;
    double minImage = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        if (slice[k] > maxImage)
          maxImage = slice[k];
        if (slice[k] < minImage)
          minImage = slice[k];
      }
    }
    double a;
    if (minImage - maxImage == 0.0D) {
      a = 1.0D;
      minImage = (maxLevel - minLevel) / 2.0D;
    }
    else {
      a = (maxLevel - minLevel) / (maxImage - minImage);
    }
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (a * (slice[k] - minImage) + minLevel);
    }
  }

  public void rescaleCenter(double minLevel, double maxLevel)
  {
    double maxImage = -1.79769313486231E+308D;
    double minImage = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        if (slice[k] > maxImage)
          maxImage = slice[k];
        if (slice[k] < minImage)
          minImage = slice[k];
      }
    }
    double center = (maxLevel + minLevel) / 2.0D;
    double a;
    if (minImage - maxImage == 0.0D) {
      a = 1.0D;
      minImage = (maxLevel - minLevel) / 2.0D;
    }
    else
    {      
      if (Math.abs(maxImage) > Math.abs(minImage))
        a = (maxLevel - center) / Math.abs(maxImage);
      else
        a = (center - minLevel) / Math.abs(minImage);
    }
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (a * (slice[k] - minImage) + center);
    }
  }

  public void abs()
  {
    double zero = 0.0D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        if (slice[k] < zero)
          slice[k] = (-slice[k]);
    }
  }

  public void log()
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = Math.log(slice[k]);
    }
  }

  public void exp()
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = Math.exp(slice[k]);
    }
  }

  public void sqrt()
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = Math.sqrt(slice[k]);
    }
  }

  public void sqr()
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] *= slice[k];
    }
  }

  public void pow(double a)
  {
    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = Math.pow(slice[k], a);
    }
  }

  public void add(double constant)
  {
    double cst = constant;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] += cst;
    }
  }

  public void multiply(double constant)
  {
    double cst = constant;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] *= cst;
    }
  }

  public void subtract(double constant)
  {
    double cst = constant;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] -= cst;
    }
  }

  public void divide(double constant)
  {
    if (constant == 0.0D) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to divide because the constant is 0.\n-------------------------------------------------------\n");
    }

    double cst = constant;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] /= cst;
    }
  }

  public void threshold(double thresholdValue)
  {
    threshold(thresholdValue, 0.0D, 255.0D);
  }

  public void threshold(double thresholdValue, double minLevel, double maxLevel)
  {
    double low = minLevel;
    double high = maxLevel;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (slice[k] > thresholdValue ? high : low);
    }
  }

  public void thresholdSoft(double thresholdValue)
  {
    double zero = 0.0D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        double pixel = slice[k];
        slice[k] = (pixel > thresholdValue ? pixel - thresholdValue : pixel <= -thresholdValue ? pixel + thresholdValue : zero);
      }
    }
  }

  public void thresholdHard(double thresholdValue)
  {
    double zero = 0.0D;

    for (int z = 0; z < this.nz; z++) {
      double[] slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        double pixel = slice[k];
        if ((pixel > -thresholdValue) && (pixel < thresholdValue))
          slice[k] = zero;
      }
    }
  }

  public void addGaussianNoise(double amplitude)
  {
    Random rnd = new Random();
    double[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] += rnd.nextGaussian() * amplitude;
    }
  }

  public void addUniformNoise(double amplitude)
  {
    Random rnd = new Random();
    double[] slice = null;
    amplitude *= 2.0D;
    for (int z = 0; z < this.nz; z++) {
      slice = (double[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] += (rnd.nextDouble() - 0.5D) * amplitude;
    }
  }

  public void addSaltPepper(double amplitudeSalt, double amplitudePepper, double percentageSalt, double percentagePepper)
  {
    Random rnd = new Random();
    double[] slice = null;

    if (percentageSalt > 0.0D) {
      double nbSalt = this.nxy * this.nz / percentageSalt;
      for (int k = 0; k < nbSalt; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((double[])this.data[z])[index] += rnd.nextDouble() * amplitudeSalt;
      }
    }
    if (percentagePepper > 0.0D) {
      double nbPepper = this.nxy * this.nz / percentagePepper;
      for (int k = 0; k < nbPepper; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((double[])this.data[z])[index] -= rnd.nextDouble() * amplitudeSalt;
      }
    }
  }
}