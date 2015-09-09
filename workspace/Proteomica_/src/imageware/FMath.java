package imageware;

public final class FMath
{
  public static int floor(float f)
  {
    if (f >= 0.0F) return (int)f;

    int iAdd = (int)f - 1;
    return (int)(f - iAdd) + iAdd;
  }

  public static int floor(double d)
  {
    if (d >= 0.0D) return (int)d;

    int iAdd = (int)d - 1;
    return (int)(d - iAdd) + iAdd;
  }

  public static int ceil(float f)
  {
    float mf = -f;
    if (mf >= 0.0F) return -(int)mf;

    int iAdd = (int)mf - 1;
    return -((int)(mf - iAdd) + iAdd);
  }

  public static int ceil(double d)
  {
    double md = -d;
    if (md >= 0.0D) return -(int)md;

    int iAdd = (int)md - 1;
    return -((int)(md - iAdd) + iAdd);
  }

  public static int round(float f)
  {
    float f05 = f + 0.5F;
    if (f05 >= 0.0D) return (int)f05;

    int iAdd = (int)f05 - 1;
    return (int)(f05 - iAdd) + iAdd;
  }

  public static int round(double d)
  {
    double d05 = d + 0.5D;
    if (d05 >= 0.0D) return (int)d05;

    int iAdd = (int)d05 - 1;
    return (int)(d05 - iAdd) + iAdd;
  }

  public static float min(float f1, float f2)
  {
    return f1 < f2 ? f1 : f2;
  }

  public static double min(double d1, double d2)
  {
    return d1 < d2 ? d1 : d2;
  }

  public static float max(float f1, float f2)
  {
    return f1 > f2 ? f1 : f2;
  }

  public static double max(double d1, double d2)
  {
    return d1 > d2 ? d1 : d2;
  }
}