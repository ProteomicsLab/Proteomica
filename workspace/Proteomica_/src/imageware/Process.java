package imageware;

public abstract interface Process extends Pointwise
{
  public abstract void smoothGaussian(double paramDouble);

  public abstract void smoothGaussian(double paramDouble1, double paramDouble2, double paramDouble3);

  public abstract void max(ImageWare paramImageWare);

  public abstract void min(ImageWare paramImageWare);

  public abstract void add(ImageWare paramImageWare);

  public abstract void multiply(ImageWare paramImageWare);

  public abstract void subtract(ImageWare paramImageWare);

  public abstract void divide(ImageWare paramImageWare);
}