package imageware;

import ij.ImageStack;
import ij.process.FloatProcessor;
import java.awt.Image;
import java.util.Random;

public class FloatPointwise extends FloatAccess
  implements Pointwise
{
  protected FloatPointwise(int nx, int ny, int nz)
  {
    super(nx, ny, nz); } 
  protected FloatPointwise(Image image, int mode) { super(image, mode); } 
  protected FloatPointwise(ImageStack stack, int mode) {
    super(stack, mode); } 
  protected FloatPointwise(ImageStack stack, byte chan) { super(stack, chan); } 
  protected FloatPointwise(byte[] array, int mode) {
    super(array, mode); } 
  protected FloatPointwise(byte[][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(byte[][][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(short[] array, int mode) { super(array, mode); } 
  protected FloatPointwise(short[][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(short[][][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(float[] array, int mode) { super(array, mode); } 
  protected FloatPointwise(float[][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(float[][][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(double[] array, int mode) { super(array, mode); } 
  protected FloatPointwise(double[][] array, int mode) { super(array, mode); } 
  protected FloatPointwise(double[][][] array, int mode) { super(array, mode); }


  public void fillConstant(double value)
  {
    float typedValue = (float)value;
    float[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = typedValue;
    }
  }

  public void fillRamp()
  {
    int off = 0;
    float[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (off + k);
      off += this.nxy;
    }
  }

  public void fillGaussianNoise(double amplitude)
  {
    Random rnd = new Random();
    float[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)(rnd.nextGaussian() * amplitude));
    }
  }

  public void fillUniformNoise(double amplitude)
  {
    Random rnd = new Random();
    float[] slice = null;
    amplitude *= 2.0D;
    for (int z = 0; z < this.nz; z++) {
      slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)((rnd.nextDouble() - 0.5D) * amplitude));
    }
  }

  public void fillSaltPepper(double amplitudeSalt, double amplitudePepper, double percentageSalt, double percentagePepper)
  {
    Random rnd = new Random();
    float[] slice = null;

    if (percentageSalt > 0.0D) {
      double nbSalt = this.nxy * this.nz / percentageSalt;
      for (int k = 0; k < nbSalt; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((float[])this.data[z])[index] = ((float)(rnd.nextDouble() * amplitudeSalt));
      }
    }
    if (percentagePepper > 0.0D) {
      double nbPepper = this.nxy * this.nz / percentagePepper;
      for (int k = 0; k < nbPepper; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((float[])this.data[z])[index] = ((float)(-rnd.nextDouble() * amplitudeSalt));
      }
    }
  }

  public ImageStack buildImageStack()
  {
    ImageStack imagestack = new ImageStack(this.nx, this.ny);
    for (int z = 0; z < this.nz; z++) {
      FloatProcessor ip = new FloatProcessor(this.nx, this.ny);
      float[] pix = (float[])ip.getPixels();
      for (int k = 0; k < this.nxy; k++)
        pix[k] = ((float[])(float[])this.data[z])[k];
      imagestack.addSlice("" + z, ip);
    }
    return imagestack;
  }

  public void invert()
  {
    double max = -1.79769313486231E+308D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        if (slice[k] > max)
          max = slice[k];
      }
    }
    double cst = max;
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)(max - slice[k]));
    }
  }

  public void negate()
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)-slice[k]);
    }
  }

  public void clip()
  {
    clip(0.0D, 255.0D);
  }

  public void clip(double minLevel, double maxLevel)
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];

      float min = (float)minLevel;
      float max = (float)maxLevel;
      for (int k = 0; k < this.nxy; k++) {
        float value = slice[k];
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
      float[] slice = (float[])this.data[z];
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
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)(a * (slice[k] - minImage)));
    }
  }

  public void rescale(double minLevel, double maxLevel)
  {
    double maxImage = -1.79769313486231E+308D;
    double minImage = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
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
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)(a * (slice[k] - minImage) + minLevel));
    }
  }

  public void rescaleCenter(double minLevel, double maxLevel)
  {
    double maxImage = -1.79769313486231E+308D;
    double minImage = 1.7976931348623157E+308D;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
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
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)(a * (slice[k] - minImage) + center));
    }
  }

  public void abs()
  {
    float zero = 0.0F;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        if (slice[k] < zero)
          slice[k] = (-slice[k]);
    }
  }

  public void log()
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)Math.log(slice[k]));
    }
  }

  public void exp()
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)Math.exp(slice[k]));
    }
  }

  public void sqrt()
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)Math.sqrt(slice[k]));
    }
  }

  public void sqr()
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] *= slice[k];
    }
  }

  public void pow(double a)
  {
    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = ((float)Math.pow(slice[k], a));
    }
  }

  public void add(double constant)
  {
    float cst = (float)constant;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] += cst;
    }
  }

  public void multiply(double constant)
  {
    float cst = (float)constant;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] *= cst;
    }
  }

  public void subtract(double constant)
  {
    float cst = (float)constant;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] -= cst;
    }
  }

  public void divide(double constant)
  {
    if (constant == 0.0D) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to divide because the constant is 0.\n-------------------------------------------------------\n");
    }

    float cst = (float)constant;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
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
    float low = (float)minLevel;
    float high = (float)maxLevel;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] = (slice[k] > thresholdValue ? high : low);
    }
  }

  public void thresholdSoft(double thresholdValue)
  {
    float zero = 0.0F;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++) {
        double pixel = slice[k];
        slice[k] = (pixel > thresholdValue ? (float)(pixel - thresholdValue) : pixel <= -thresholdValue ? (float)(pixel + thresholdValue) : zero);
      }
    }
  }

  public void thresholdHard(double thresholdValue)
  {
    float zero = 0.0F;

    for (int z = 0; z < this.nz; z++) {
      float[] slice = (float[])this.data[z];
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
    float[] slice = null;
    for (int z = 0; z < this.nz; z++) {
      slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] += (float)(rnd.nextGaussian() * amplitude);
    }
  }

  public void addUniformNoise(double amplitude)
  {
    Random rnd = new Random();
    float[] slice = null;
    amplitude *= 2.0D;
    for (int z = 0; z < this.nz; z++) {
      slice = (float[])this.data[z];
      for (int k = 0; k < this.nxy; k++)
        slice[k] += (float)((rnd.nextDouble() - 0.5D) * amplitude);
    }
  }

  public void addSaltPepper(double amplitudeSalt, double amplitudePepper, double percentageSalt, double percentagePepper)
  {
    Random rnd = new Random();
    float[] slice = null;

    if (percentageSalt > 0.0D) {
      double nbSalt = this.nxy * this.nz / percentageSalt;
      for (int k = 0; k < nbSalt; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((float[])this.data[z])[index] += (float)(rnd.nextDouble() * amplitudeSalt);
      }
    }
    if (percentagePepper > 0.0D) {
      double nbPepper = this.nxy * this.nz / percentagePepper;
      for (int k = 0; k < nbPepper; k++) {
        int index = (int)(rnd.nextDouble() * this.nxy);
        int z = (int)(rnd.nextDouble() * this.nz);
        ((float[])this.data[z])[index] -= (float)(rnd.nextDouble() * amplitudeSalt);
      }
    }
  }
}