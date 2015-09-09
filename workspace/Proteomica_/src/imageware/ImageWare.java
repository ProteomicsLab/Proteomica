package imageware;

public abstract interface ImageWare extends Process
{
  public static final int UNDEFINED_TYPE = 0;
  public static final int BYTE = 1;
  public static final int SHORT = 2;
  public static final int FLOAT = 3;
  public static final int DOUBLE = 4;
  public static final byte UNDEFINED_BOUNDARY = 0;
  public static final byte NONE = 1;
  public static final byte MIRROR = 2;
  public static final byte PERIODIC = 3;
  public static final int UNDEFINED = 0;
  public static final int CREATE = 1;
  public static final int WRAP = 2;
  public static final byte RED = 0;
  public static final byte GREEN = 1;
  public static final byte BLUE = 2;

  public abstract ImageWare duplicate();

  public abstract ImageWare replicate();

  public abstract ImageWare replicate(int paramInt);

  public abstract void copy(ImageWare paramImageWare);

  public abstract ImageWare convert(int paramInt);

  public abstract void printInfo();

  public abstract void show();

  public abstract void show(String paramString);

  public abstract double getMinimum();

  public abstract double getMaximum();

  public abstract double getMean();

  public abstract double getNorm1();

  public abstract double getNorm2();

  public abstract double getTotal();

  public abstract double[] getMinMax();
}