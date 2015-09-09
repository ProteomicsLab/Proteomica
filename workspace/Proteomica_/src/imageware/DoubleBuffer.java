package imageware;

import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import java.awt.Image;
import java.awt.image.ImageObserver;
import java.awt.image.PixelGrabber;

public class DoubleBuffer
  implements Buffer
{
  protected Object[] data = null;
  protected int nx = 0;
  protected int ny = 0;
  protected int nz = 0;
  protected int nxy = 0;

  protected DoubleBuffer(int nx, int ny, int nz)
  {
    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
      throw_constructor(nx, ny, nz);
    allocate();
  }

  protected DoubleBuffer(Image image, int mode)
  {
    if (image == null) {
      throw_constructor();
    }
    ImageObserver observer = null;
    this.nx = image.getWidth(observer);
    this.ny = image.getHeight(observer);
    this.nz = 1;
    this.nxy = (this.nx * this.ny);
    byte[] pixels = new byte[this.nxy];
    PixelGrabber pg = new PixelGrabber(image, 0, 0, this.nx, this.ny, false);
    try {
      pg.grabPixels();
      pixels = (byte[])pg.getPixels();
    }
    catch (Exception e) {
      throw_constructor();
    }
    allocate();
    for (int k = 0; k < this.nxy; k++)
      ((double[])this.data[0])[k] = (pixels[k] & 0xFF);
  }

  protected DoubleBuffer(ImageStack stack, int mode)
  {
    if (stack == null) {
      throw_constructor();
    }
    this.nx = stack.getWidth();
    this.ny = stack.getHeight();
    this.nz = stack.getSize();
    this.nxy = (this.nx * this.ny);
    switch (mode) {
    case 2:
      throw_constructor();
      break;
    case 1:
      allocate();
      ImageProcessor ip = stack.getProcessor(1);
      if ((ip instanceof ByteProcessor)) {
        Object[] vol = stack.getImageArray();
        for (int z = 0; z < this.nz; z++) {
          byte[] slice = (byte[])vol[z];
          for (int k = 0; k < this.nxy; k++) {
            ((double[])this.data[z])[k] = (slice[k] & 0xFF);
          }
        }
      }
      else if ((ip instanceof ShortProcessor)) {
        Object[] vol = stack.getImageArray();
        for (int z = 0; z < this.nz; z++) {
          short[] slice = (short[])vol[z];
          for (int k = 0; k < this.nxy; k++) {
            ((double[])this.data[z])[k] = (slice[k] & 0xFFFF);
          }
        }
      }
      else if ((ip instanceof FloatProcessor)) {
        Object[] vol = stack.getImageArray();
        for (int z = 0; z < this.nz; z++) {
          float[] slice = (float[])vol[z];
          for (int k = 0; k < this.nxy; k++) {
            ((double[])this.data[z])[k] = slice[k];
          }
        }
      }
      else if ((ip instanceof ColorProcessor))
      {
        for (int z = 0; z < this.nz; z++) {
          ColorProcessor cp = (ColorProcessor)stack.getProcessor(z + 1);
          int[] pixels = (int[])cp.getPixels();
          for (int k = 0; k < this.nxy; k++) {
            int c = pixels[k];
            double r = (c & 0xFF0000) >> 16;
            double g = (c & 0xFF00) >> 8;
            double b = c & 0xFF;
            ((double[])this.data[z])[k] = ((r + g + b) / 3.0D);
          }
        }
      }
      else {
        throw_constructor();
      }
      break;
    default:
      throw_constructor();
    }
  }

  protected DoubleBuffer(ImageStack stack, byte channel)
  {
    if (stack == null) {
      throw_constructor();
    }
    this.nx = stack.getWidth();
    this.ny = stack.getHeight();
    this.nz = stack.getSize();
    this.nxy = (this.nx * this.ny);
    allocate();
    ImageProcessor ip = stack.getProcessor(1);
    if ((ip instanceof ByteProcessor)) {
      Object[] vol = stack.getImageArray();
      for (int z = 0; z < this.nz; z++) {
        byte[] slice = (byte[])vol[z];
        for (int k = 0; k < this.nxy; k++) {
          ((double[])this.data[z])[k] = (slice[k] & 0xFF);
        }
      }
    }
    else if ((ip instanceof ShortProcessor)) {
      Object[] vol = stack.getImageArray();
      for (int z = 0; z < this.nz; z++) {
        short[] slice = (short[])vol[z];
        for (int k = 0; k < this.nxy; k++) {
          ((double[])this.data[z])[k] = (slice[k] & 0xFFFF);
        }
      }
    }
    else if ((ip instanceof FloatProcessor)) {
      Object[] vol = stack.getImageArray();
      for (int z = 0; z < this.nz; z++) {
        float[] slice = (float[])vol[z];
        for (int k = 0; k < this.nxy; k++) {
          ((double[])this.data[z])[k] = slice[k];
        }
      }
    }
    else if ((ip instanceof ColorProcessor))
    {
      for (int z = 0; z < this.nz; z++) {
        ColorProcessor cp = (ColorProcessor)stack.getProcessor(z + 1);
        int[] pixels = (int[])cp.getPixels();
        switch (channel) {
        case 0:
          for (int k = 0; k < this.nxy; k++) {
            ((double[])this.data[z])[k] = ((pixels[k] & 0xFF0000) >> 16);
          }
          break;
        case 1:
          for (int k = 0; k < this.nxy; k++) {
            ((double[])this.data[z])[k] = ((pixels[k] & 0xFF00) >> 8);
          }
          break;
        case 2:
          for (int k = 0; k < this.nxy; k++) {
            ((double[])this.data[z])[k] = (pixels[k] & 0xFF);
          }
          break;
        default:
          throw_constructor();
        }
      }
    }
    else {
      throw_constructor();
    }
  }

  protected DoubleBuffer(byte[] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = 1;
    this.nz = 1;
    allocate();
    putX(0, 0, 0, array);
  }

  protected DoubleBuffer(byte[][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = 1;
    allocate();
    putXY(0, 0, 0, array);
  }

  protected DoubleBuffer(byte[][][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = array[0][0].length;
    allocate();
    putXYZ(0, 0, 0, array);
  }

  protected DoubleBuffer(short[] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = 1;
    this.nz = 1;
    allocate();
    putX(0, 0, 0, array);
  }

  protected DoubleBuffer(short[][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = 1;
    allocate();
    putXY(0, 0, 0, array);
  }

  protected DoubleBuffer(short[][][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = array[0][0].length;
    allocate();
    putXYZ(0, 0, 0, array);
  }

  protected DoubleBuffer(float[] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = 1;
    this.nz = 1;
    allocate();
    putX(0, 0, 0, array);
  }

  protected DoubleBuffer(float[][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = 1;
    allocate();
    putXY(0, 0, 0, array);
  }

  protected DoubleBuffer(float[][][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = array[0][0].length;
    allocate();
    putXYZ(0, 0, 0, array);
  }

  protected DoubleBuffer(double[] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = 1;
    this.nz = 1;
    allocate();
    putX(0, 0, 0, array);
  }

  protected DoubleBuffer(double[][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = 1;
    allocate();
    putXY(0, 0, 0, array);
  }

  protected DoubleBuffer(double[][][] array, int mode)
  {
    if (array == null) {
      throw_constructor();
    }
    this.nx = array.length;
    this.ny = array[0].length;
    this.nz = array[0][0].length;
    allocate();
    putXYZ(0, 0, 0, array);
  }

  public int getType()
  {
    return 4;
  }

  public String getTypeToString()
  {
    return "Double";
  }

  public int getDimension()
  {
    int dims = 0;
    dims += (this.nx > 1 ? 1 : 0);
    dims += (this.ny > 1 ? 1 : 0);
    dims += (this.nz > 1 ? 1 : 0);
    return dims;
  }

  public int[] getSize()
  {
    int[] size = { this.nx, this.ny, this.nz };
    return size;
  }

  public int getSizeX()
  {
    return this.nx;
  }

  public int getSizeY()
  {
    return this.ny;
  }

  public int getSizeZ()
  {
    return this.nz;
  }

  public int getWidth()
  {
    return this.nx;
  }

  public int getHeight()
  {
    return this.ny;
  }

  public int getDepth()
  {
    return this.nz;
  }

  public int getTotalSize()
  {
    return this.nxy * this.nz;
  }

  public boolean isSameSize(ImageWare imageware)
  {
    if (this.nx != imageware.getSizeX())
      return false;
    if (this.ny != imageware.getSizeY())
      return false;
    if (this.nz != imageware.getSizeZ())
      return false;
    return true;
  }

  public void putX(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    double[] buf = new double[bnx];
    buffer.getX(0, 0, 0, buf);
    putX(x, y, z, buf);
  }

  public void putY(int x, int y, int z, ImageWare buffer)
  {
    int bny = buffer.getSizeY();
    double[] buf = new double[bny];
    buffer.getY(0, 0, 0, buf);
    putY(x, y, z, buf);
  }

  public void putZ(int x, int y, int z, ImageWare buffer)
  {
    int bnz = buffer.getSizeZ();
    double[] buf = new double[bnz];
    buffer.getZ(0, 0, 0, buf);
    putZ(x, y, z, buf);
  }

  public void putXY(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    int bny = buffer.getSizeY();
    double[][] buf = new double[bnx][bny];
    buffer.getXY(0, 0, 0, buf);
    putXY(x, y, z, buf);
  }

  public void putXZ(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    int bnz = buffer.getSizeZ();
    double[][] buf = new double[bnx][bnz];
    buffer.getXZ(0, 0, 0, buf);
    putXZ(x, y, z, buf);
  }

  public void putYZ(int x, int y, int z, ImageWare buffer)
  {
    int bny = buffer.getSizeY();
    int bnz = buffer.getSizeZ();
    double[][] buf = new double[bny][bnz];
    buffer.getYZ(0, 0, 0, buf);
    putYZ(x, y, z, buf);
  }

  public void putXYZ(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    int bny = buffer.getSizeY();
    int bnz = buffer.getSizeZ();
    double[][][] buf = new double[bnx][bny][bnz];
    buffer.getXYZ(0, 0, 0, buf);
    putXYZ(x, y, z, buf);
  }

  public void putX(int x, int y, int z, byte[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      for (int i = 0; i < leni; i++) {
        tmp[offset] = (buffer[i] & 0xFF);
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "No check", buffer, x, y, z);
    }
  }

  public void putX(int x, int y, int z, short[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      for (int i = 0; i < leni; i++) {
        tmp[offset] = (buffer[i] & 0xFFFF);
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "No check", buffer, x, y, z);
    }
  }

  public void putX(int x, int y, int z, float[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      for (int i = 0; i < leni; i++) {
        tmp[offset] = buffer[i];
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "No check", buffer, x, y, z);
    }
  }

  public void putX(int x, int y, int z, double[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      System.arraycopy(buffer, 0, tmp, offset, leni);
    }
    catch (Exception e) {
      throw_put("X", "No check", buffer, x, y, z);
    }
  }

  public void putY(int x, int y, int z, byte[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        tmp[offset] = (buffer[i] & 0xFF);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "No check", buffer, x, y, z);
    }
  }

  public void putY(int x, int y, int z, short[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        tmp[offset] = (buffer[i] & 0xFFFF);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "No check", buffer, x, y, z);
    }
  }

  public void putY(int x, int y, int z, float[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        tmp[offset] = buffer[i];
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "No check", buffer, x, y, z);
    }
  }

  public void putY(int x, int y, int z, double[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        tmp[offset] = buffer[i];
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "No check", buffer, x, y, z);
    }
  }

  public void putZ(int x, int y, int z, byte[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        ((double[])this.data[z])[offset] = (buffer[i] & 0xFF);
        z++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "No check", buffer, x, y, z);
    }
  }

  public void putZ(int x, int y, int z, short[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        ((double[])this.data[z])[offset] = (buffer[i] & 0xFFFF);
        z++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "No check", buffer, x, y, z);
    }
  }

  public void putZ(int x, int y, int z, float[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        ((double[])this.data[z])[offset] = buffer[i];
        z++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "No check", buffer, x, y, z);
    }
  }

  public void putZ(int x, int y, int z, double[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        ((double[])this.data[z])[offset] = buffer[i];
        z++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "No check", buffer, x, y, z);
    }
  }

  public void putXY(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          tmp[offset] = (buffer[i][j] & 0xFF);

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_put("XY", "No check", buffer, x, y, z);
    }
  }

  public void putXY(int x, int y, int z, short[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          tmp[offset] = (buffer[i][j] & 0xFFFF);

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_put("XY", "No check", buffer, x, y, z);
    }
  }

  public void putXY(int x, int y, int z, float[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          tmp[offset] = buffer[i][j];

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_put("XY", "No check", buffer, x, y, z);
    }
  }

  public void putXY(int x, int y, int z, double[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          tmp[offset] = buffer[i][j];

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_put("XY", "No check", buffer, x, y, z);
    }
  }

  public void putXZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + j * this.nx;
        for (int i = 0; i < leni; offset++) {
          ((double[])this.data[z])[offset] = (buffer[i][j] & 0xFF);

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_put("YZ", "No check", buffer, x, y, z);
    }
  }

  public void putXZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + j * this.nx;
        for (int i = 0; i < leni; offset++) {
          ((double[])this.data[z])[offset] = (buffer[i][j] & 0xFFFF);

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_put("YZ", "No check", buffer, x, y, z);
    }
  }

  public void putXZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + j * this.nx;
        for (int i = 0; i < leni; offset++) {
          ((double[])this.data[z])[offset] = buffer[i][j];

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_put("YZ", "No check", buffer, x, y, z);
    }
  }

  public void putXZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + j * this.nx;
        for (int i = 0; i < leni; offset++) {
          ((double[])this.data[z])[offset] = buffer[i][j];

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_put("YZ", "No check", buffer, x, y, z);
    }
  }

  public void putYZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          ((double[])this.data[z])[offset] = (buffer[i][j] & 0xFF);

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_put("XZ", "No check", buffer, x, y, z);
    }
  }

  public void putYZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          ((double[])this.data[z])[offset] = (buffer[i][j] & 0xFFFF);

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_put("XZ", "No check", buffer, x, y, z);
    }
  }

  public void putYZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          ((double[])this.data[z])[offset] = buffer[i][j];

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_put("XZ", "No check", buffer, x, y, z);
    }
  }

  public void putYZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          ((double[])this.data[z])[offset] = buffer[i][j];

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_put("XZ", "No check", buffer, x, y, z);
    }
  }

  public void putXYZ(int x, int y, int z, byte[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            tmp[offset] = (buffer[i][j][k] & 0xFF);

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_put("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void putXYZ(int x, int y, int z, short[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            tmp[offset] = (buffer[i][j][k] & 0xFFFF);

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_put("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void putXYZ(int x, int y, int z, float[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            tmp[offset] = buffer[i][j][k];

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_put("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void putXYZ(int x, int y, int z, double[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            tmp[offset] = buffer[i][j][k];

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_put("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void getX(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    double[] buf = new double[bnx];
    getX(x, y, z, buf);
    buffer.putX(0, 0, 0, buf);
  }

  public void getY(int x, int y, int z, ImageWare buffer)
  {
    int bny = buffer.getSizeY();
    double[] buf = new double[bny];
    getY(x, y, z, buf);
    buffer.putY(0, 0, 0, buf);
  }

  public void getZ(int x, int y, int z, ImageWare buffer)
  {
    int bnz = buffer.getSizeZ();
    double[] buf = new double[bnz];
    getZ(x, y, z, buf);
    buffer.putZ(0, 0, 0, buf);
  }

  public void getXY(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    int bny = buffer.getSizeY();
    double[][] buf = new double[bnx][bny];
    getXY(x, y, z, buf);
    buffer.putXY(0, 0, 0, buf);
  }

  public void getXZ(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    int bnz = buffer.getSizeZ();
    double[][] buf = new double[bnx][bnz];
    getXZ(x, y, z, buf);
    buffer.putXZ(0, 0, 0, buf);
  }

  public void getYZ(int x, int y, int z, ImageWare buffer)
  {
    int bny = buffer.getSizeY();
    int bnz = buffer.getSizeZ();
    double[][] buf = new double[bny][bnz];
    getYZ(x, y, z, buf);
    buffer.putYZ(0, 0, 0, buf);
  }

  public void getXYZ(int x, int y, int z, ImageWare buffer)
  {
    int bnx = buffer.getSizeX();
    int bny = buffer.getSizeY();
    int bnz = buffer.getSizeZ();
    double[][][] buf = new double[bnx][bny][bnz];
    getXYZ(x, y, z, buf);
    buffer.putXYZ(0, 0, 0, buf);
  }

  public void getX(int x, int y, int z, byte[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      for (int i = 0; i < leni; i++) {
        buffer[i] = ((byte)(int)tmp[offset]);
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getX(int x, int y, int z, short[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      for (int i = 0; i < leni; i++) {
        buffer[i] = ((short)(int)tmp[offset]);
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getX(int x, int y, int z, float[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      for (int i = 0; i < leni; i++) {
        buffer[i] = ((float)tmp[offset]);
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getX(int x, int y, int z, double[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];

      System.arraycopy(tmp, offset, buffer, 0, leni);
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getY(int x, int y, int z, byte[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((byte)(int)tmp[offset]);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getY(int x, int y, int z, short[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((short)(int)tmp[offset]);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getY(int x, int y, int z, float[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((float)tmp[offset]);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getY(int x, int y, int z, double[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      double[] tmp = (double[])this.data[z];
      for (int i = 0; i < leni; i++) {
        buffer[i] = tmp[offset];
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("X", "No check", buffer, x, y, z);
    }
  }

  public void getZ(int x, int y, int z, byte[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((byte)(int)((double[])(double[])this.data[z])[offset]);
        z++;
      }
    }
    catch (Exception e) {
      throw_get("Y", "No check", buffer, x, y, z);
    }
  }

  public void getZ(int x, int y, int z, short[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((short)(int)((double[])(double[])this.data[z])[offset]);
        z++;
      }
    }
    catch (Exception e) {
      throw_get("Y", "No check", buffer, x, y, z);
    }
  }

  public void getZ(int x, int y, int z, float[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((float)((double[])(double[])this.data[z])[offset]);
        z++;
      }
    }
    catch (Exception e) {
      throw_get("Y", "No check", buffer, x, y, z);
    }
  }

  public void getZ(int x, int y, int z, double[] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      for (int i = 0; i < leni; i++) {
        buffer[i] = ((double[])(double[])this.data[z])[offset];
        z++;
      }
    }
    catch (Exception e) {
      throw_get("Y", "No check", buffer, x, y, z);
    }
  }

  public void getXY(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((byte)(int)tmp[offset]);

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_get("XY", "No check", buffer, x, y, z);
    }
  }

  public void getXY(int x, int y, int z, short[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((short)(int)tmp[offset]);

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_get("XY", "No check", buffer, x, y, z);
    }
  }

  public void getXY(int x, int y, int z, float[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((float)tmp[offset]);

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_get("XY", "No check", buffer, x, y, z);
    }
  }

  public void getXY(int x, int y, int z, double[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      double[] tmp = (double[])this.data[z];
      for (int j = 0; j < lenj; j++) {
        offset = x + (y + j) * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = tmp[offset];

          i++;
        }
      }
    }
    catch (Exception e)
    {
      throw_get("XY", "No check", buffer, x, y, z);
    }
  }

  public void getXZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + y * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[z])[offset]);

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_get("XZ", "No check", buffer, x, y, z);
    }
  }

  public void getXZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + y * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[z])[offset]);

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_get("XZ", "No check", buffer, x, y, z);
    }
  }

  public void getXZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + y * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((float)((double[])(double[])this.data[z])[offset]);

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_get("XZ", "No check", buffer, x, y, z);
    }
  }

  public void getXZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; z++) {
        offset = x + y * this.nx;
        for (int i = 0; i < leni; offset++) {
          buffer[i][j] = ((double[])(double[])this.data[z])[offset];

          i++;
        }
        j++;
      }

    }
    catch (Exception e)
    {
      throw_get("XZ", "No check", buffer, x, y, z);
    }
  }

  public void getYZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[z])[offset]);

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_get("YZ", "No check", buffer, x, y, z);
    }
  }

  public void getYZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[z])[offset]);

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_get("YZ", "No check", buffer, x, y, z);
    }
  }

  public void getYZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          buffer[i][j] = ((float)((double[])(double[])this.data[z])[offset]);

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_get("YZ", "No check", buffer, x, y, z);
    }
  }

  public void getYZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      for (int j = 0; j < lenj; offset = x + this.nx * y) {
        for (int i = 0; i < leni; offset += this.nx) {
          buffer[i][j] = ((double[])(double[])this.data[z])[offset];

          i++;
        }
        j++; z++;
      }

    }
    catch (Exception e)
    {
      throw_get("YZ", "No check", buffer, x, y, z);
    }
  }

  public void getXYZ(int x, int y, int z, byte[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            buffer[i][j][k] = ((byte)(int)tmp[offset]);

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_get("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void getXYZ(int x, int y, int z, short[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            buffer[i][j][k] = ((short)(int)tmp[offset]);

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_get("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void getXYZ(int x, int y, int z, float[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            buffer[i][j][k] = ((float)tmp[offset]);

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_get("XYZ", "No check", buffer, x, y, z);
    }
  }

  public void getXYZ(int x, int y, int z, double[][][] buffer)
  {
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; z++) {
        double[] tmp = (double[])this.data[z];
        for (int j = 0; j < lenj; j++) {
          offset = x + (j + y) * this.nx;
          for (int i = 0; i < leni; offset++) {
            buffer[i][j][k] = tmp[offset];

            i++;
          }
        }
        k++;
      }

    }
    catch (Exception e)
    {
      throw_get("XYZ", "No check", buffer, x, y, z);
    }
  }

  protected void throw_constructor()
  {
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to create a double imageware.\n-------------------------------------------------------\n");
  }

  protected void throw_constructor(int nx, int ny, int nz)
  {
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to create a double imageware " + nx + "," + ny + "," + nz + "].\n" + "-------------------------------------------------------\n");
  }

  protected void throw_get(String direction, String border, Object buffer, int x, int y, int z)
  {
    int leni = 0;
    int lenj = 0;
    int lenk = 0;
    String type = " unknown type";
    if ((buffer instanceof byte[])) {
      leni = ((byte[])buffer).length;
      type = " 1D byte";
    }
    else if ((buffer instanceof short[])) {
      leni = ((short[])buffer).length;
      type = " 1D short";
    }
    else if ((buffer instanceof float[])) {
      leni = ((float[])buffer).length;
      type = " 1D float";
    }
    else if ((buffer instanceof double[])) {
      leni = ((double[])buffer).length;
      type = " 1D double";
    }
    else if ((buffer instanceof byte[][])) {
      leni = ((byte[][])buffer).length;
      lenj = ((byte[][])(byte[][])buffer)[0].length;
      type = " 2D byte";
    }
    else if ((buffer instanceof short[][])) {
      leni = ((short[][])buffer).length;
      lenj = ((short[][])(short[][])buffer)[0].length;
      type = " 2D short";
    }
    else if ((buffer instanceof float[][])) {
      leni = ((float[][])buffer).length;
      lenj = ((float[][])(float[][])buffer)[0].length;
      type = " 2D float";
    }
    else if ((buffer instanceof double[][])) {
      leni = ((double[][])buffer).length;
      lenj = ((double[][])(double[][])buffer)[0].length;
      type = " 2D double";
    }
    else if ((buffer instanceof byte[][][])) {
      leni = ((byte[][][])buffer).length;
      lenj = ((byte[][][])(byte[][][])buffer)[0].length;
      lenk = ((byte[][][])(byte[][][])buffer)[0][0].length;
      type = " 3D byte";
    }
    else if ((buffer instanceof short[][][])) {
      leni = ((short[][][])buffer).length;
      lenj = ((short[][][])(short[][][])buffer)[0].length;
      lenk = ((short[][][])(short[][][])buffer)[0][0].length;
      type = " 3D short";
    }
    else if ((buffer instanceof float[][][])) {
      leni = ((float[][][])buffer).length;
      lenj = ((float[][][])(float[][][])buffer)[0].length;
      lenk = ((float[][][])(float[][][])buffer)[0][0].length;
      type = " 3D float";
    }
    else if ((buffer instanceof double[][][])) {
      leni = ((double[][][])buffer).length;
      lenj = ((double[][][])(double[][][])buffer)[0].length;
      lenk = ((double[][][])(double[][][])buffer)[0][0].length;
      type = " 3D double";
    }
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to get a" + type + " buffer [" + (leni == 0 ? "" : new StringBuilder().append("").append(leni).toString()) + (lenj == 0 ? "" : new StringBuilder().append(",").append(lenj).toString()) + (lenk == 0 ? "" : new StringBuilder().append(",").append(lenk).toString()) + "] \n" + "from the double imageware [" + this.nx + "," + this.ny + "," + this.nz + "]\n" + "at the position (" + x + "," + y + "," + z + ") in direction " + direction + "\n" + "using " + border + ".\n" + "-------------------------------------------------------\n");
  }

  protected void throw_put(String direction, String border, Object buffer, int x, int y, int z)
  {
    int leni = 0;
    int lenj = 0;
    int lenk = 0;
    String type = " unknown type";
    if ((buffer instanceof byte[])) {
      leni = ((byte[])buffer).length;
      type = " 1D byte";
    }
    else if ((buffer instanceof short[])) {
      leni = ((short[])buffer).length;
      type = " 1D short";
    }
    else if ((buffer instanceof float[])) {
      leni = ((float[])buffer).length;
      type = " 1D float";
    }
    else if ((buffer instanceof double[])) {
      leni = ((double[])buffer).length;
      type = " 1D double";
    }
    else if ((buffer instanceof byte[][])) {
      leni = ((byte[][])buffer).length;
      lenj = ((byte[][])(byte[][])buffer)[0].length;
      type = " 2D byte";
    }
    else if ((buffer instanceof short[][])) {
      leni = ((short[][])buffer).length;
      lenj = ((short[][])(short[][])buffer)[0].length;
      type = " 2D short";
    }
    else if ((buffer instanceof float[][])) {
      leni = ((float[][])buffer).length;
      lenj = ((float[][])(float[][])buffer)[0].length;
      type = " 2D float";
    }
    else if ((buffer instanceof double[][])) {
      leni = ((double[][])buffer).length;
      lenj = ((double[][])(double[][])buffer)[0].length;
      type = " 2D double";
    }
    else if ((buffer instanceof byte[][][])) {
      leni = ((byte[][][])buffer).length;
      lenj = ((byte[][][])(byte[][][])buffer)[0].length;
      lenk = ((byte[][][])(byte[][][])buffer)[0][0].length;
      type = " 3D byte";
    }
    else if ((buffer instanceof short[][][])) {
      leni = ((short[][][])buffer).length;
      lenj = ((short[][][])(short[][][])buffer)[0].length;
      lenk = ((short[][][])(short[][][])buffer)[0][0].length;
      type = " 3D short";
    }
    else if ((buffer instanceof float[][][])) {
      leni = ((float[][][])buffer).length;
      lenj = ((float[][][])(float[][][])buffer)[0].length;
      lenk = ((float[][][])(float[][][])buffer)[0][0].length;
      type = " 3D float";
    }
    else if ((buffer instanceof double[][][])) {
      leni = ((double[][][])buffer).length;
      lenj = ((double[][][])(double[][][])buffer)[0].length;
      lenk = ((double[][][])(double[][][])buffer)[0][0].length;
      type = " 3D double";
    }
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to put a" + type + " buffer [" + (leni == 0 ? "" : new StringBuilder().append("").append(leni).toString()) + (lenj == 0 ? "" : new StringBuilder().append(",").append(lenj).toString()) + (lenk == 0 ? "" : new StringBuilder().append(",").append(lenk).toString()) + "] \n" + "into the double imageware [" + this.nx + "," + this.ny + "," + this.nz + "]\n" + "at the position (" + x + "," + y + "," + z + ") in direction " + direction + "\n" + "using " + border + ".\n" + "-------------------------------------------------------\n");
  }

  public Object[] getVolume()
  {
    return this.data;
  }

  public byte[] getSliceByte(int z)
  {
    return null;
  }

  public short[] getSliceShort(int z)
  {
    return null;
  }

  public float[] getSliceFloat(int z)
  {
    return null;
  }

  public double[] getSliceDouble(int z)
  {
    return (double[])this.data[z];
  }

  private void allocate()
  {
    try
    {
      this.data = new Object[this.nz];
      this.nxy = (this.nx * this.ny);
      for (int z = 0; z < this.nz; z++)
        this.data[z] = new double[this.nxy];
    }
    catch (Exception e) {
      throw_constructor(this.nx, this.ny, this.nz);
    }
  }
}