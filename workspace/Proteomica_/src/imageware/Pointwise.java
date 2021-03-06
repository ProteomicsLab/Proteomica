package imageware;

import ij.ImageStack;

public abstract interface Pointwise extends Access
{
  public abstract void fillConstant(double paramDouble);

  public abstract void fillRamp();

  public abstract void fillGaussianNoise(double paramDouble);

  public abstract void fillUniformNoise(double paramDouble);

  public abstract void fillSaltPepper(double paramDouble1, double paramDouble2, double paramDouble3, double paramDouble4);

  public abstract ImageStack buildImageStack();

  public abstract void invert();

  public abstract void negate();

  public abstract void rescale();

  public abstract void clip();

  public abstract void clip(double paramDouble1, double paramDouble2);

  public abstract void rescale(double paramDouble1, double paramDouble2);

  public abstract void rescaleCenter(double paramDouble1, double paramDouble2);

  public abstract void abs();

  public abstract void log();

  public abstract void exp();

  public abstract void sqrt();

  public abstract void sqr();

  public abstract void pow(double paramDouble);

  public abstract void add(double paramDouble);

  public abstract void multiply(double paramDouble);

  public abstract void subtract(double paramDouble);

  public abstract void divide(double paramDouble);

  public abstract void threshold(double paramDouble);

  public abstract void threshold(double paramDouble1, double paramDouble2, double paramDouble3);

  public abstract void thresholdHard(double paramDouble);

  public abstract void thresholdSoft(double paramDouble);

  public abstract void addGaussianNoise(double paramDouble);

  public abstract void addUniformNoise(double paramDouble);

  public abstract void addSaltPepper(double paramDouble1, double paramDouble2, double paramDouble3, double paramDouble4);
}