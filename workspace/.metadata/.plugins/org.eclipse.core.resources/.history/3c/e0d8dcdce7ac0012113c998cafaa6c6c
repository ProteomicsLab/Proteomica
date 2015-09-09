package imageware;

import ij.ImageStack;
import java.awt.Image;

public class DoubleAccess extends DoubleBuffer
  implements Access
{
  protected DoubleAccess(int nx, int ny, int nz)
  {
    super(nx, ny, nz); } 
  protected DoubleAccess(Image image, int mode) { super(image, mode); } 
  protected DoubleAccess(ImageStack stack, int mode) {
    super(stack, mode); } 
  protected DoubleAccess(ImageStack stack, byte chan) { super(stack, chan); } 
  protected DoubleAccess(byte[] array, int mode) {
    super(array, mode); } 
  protected DoubleAccess(byte[][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(byte[][][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(short[] array, int mode) { super(array, mode); } 
  protected DoubleAccess(short[][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(short[][][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(float[] array, int mode) { super(array, mode); } 
  protected DoubleAccess(float[][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(float[][][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(double[] array, int mode) { super(array, mode); } 
  protected DoubleAccess(double[][] array, int mode) { super(array, mode); } 
  protected DoubleAccess(double[][][] array, int mode) { super(array, mode); }


  public double getPixel(int x, int y, int z)
  {
    if (x >= this.nx) return 0.0D;
    if (y >= this.ny) return 0.0D;
    if (z >= this.nz) return 0.0D;
    if (x < 0) return 0.0D;
    if (y < 0) return 0.0D;
    if (z < 0) return 0.0D;
    return ((double[])(double[])this.data[z])[(x + y * this.nx)];
  }

  public double getPixel(int x, int y, int z, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to put a pixel \nat the position (" + x + "," + y + "," + z + ".\n" + "-------------------------------------------------------\n");
    }

    int xp = x;
    while (xp < 0)
      xp += xperiod;
    while (xp >= this.nx) {
      xp = xperiod - xp;
      xp = xp < 0 ? -xp : xp;
    }
    int yp = y;
    while (yp < 0)
      yp += yperiod;
    while (yp >= this.ny) {
      yp = yperiod - yp;
      yp = yp < 0 ? -yp : yp;
    }
    int zp = z;
    while (zp < 0)
      zp += zperiod;
    while (zp >= this.nz) {
      zp = zperiod - zp;
      zp = zp < 0 ? -zp : zp;
    }
    return ((double[])(double[])this.data[zp])[(xp + yp * this.nx)];
  }

  public double getInterpolatedPixel(double x, double y, double z)
  {
    if (x > this.nx - 1) return 0.0D;
    if (y > this.ny - 1) return 0.0D;
    if (z > this.nz - 1) return 0.0D;
    if (x < 0.0D) return 0.0D;
    if (y < 0.0D) return 0.0D;
    if (z < 0.0D) return 0.0D;
    double output = 0.0D;

    int i = x >= 0.0D ? (int)x : (int)x - 1;
    int j = y >= 0.0D ? (int)y : (int)y - 1;
    int k = z >= 0.0D ? (int)z : (int)z - 1;
    boolean fi = i == this.nx - 1;
    boolean fj = j == this.ny - 1;
    boolean fk = k == this.nz - 1;
    int index = i + j * this.nx;
    switch (getDimension()) {
    case 1:
      double v1_0 = ((double[])(double[])this.data[k])[index];
      double v1_1 = fi ? v1_0 : ((double[])(double[])this.data[k])[(index + 1)];
      double dx1 = x - i;
      return v1_1 * dx1 - v1_0 * (dx1 - 1.0D);
    case 2:
      double v2_00 = ((double[])(double[])this.data[k])[index];
      double v2_10 = fi ? v2_00 : ((double[])(double[])this.data[k])[(index + 1)];
      double v2_01 = fj ? v2_00 : ((double[])(double[])this.data[k])[(index + this.nx)];
      double v2_11 = fi ? v2_01 : fj ? v2_00 : ((double[])(double[])this.data[k])[(index + 1 + this.nx)];
      double dx2 = x - i;
      double dy2 = y - j;
      return dx2 * (v2_11 * dy2 - v2_10 * (dy2 - 1.0D)) - (dx2 - 1.0D) * (v2_01 * dy2 - v2_00 * (dy2 - 1.0D));
    case 3:
      double v3_000 = ((double[])(double[])this.data[k])[index];
      double v3_100 = fi ? v3_000 : ((double[])(double[])this.data[k])[(index + 1)];
      double v3_010 = fj ? v3_000 : ((double[])(double[])this.data[k])[(index + this.nx)];
      double v3_110 = fi ? v3_010 : fj ? v3_000 : ((double[])(double[])this.data[k])[(index + 1 + this.nx)];
      double v3_001 = fk ? v3_000 : ((double[])(double[])this.data[(k + 1)])[index];
      double v3_011 = fk ? v3_010 : fj ? v3_000 : ((double[])(double[])this.data[(k + 1)])[(index + 1)];
      double v3_101 = fk ? v3_100 : fi ? v3_000 : ((double[])(double[])this.data[(k + 1)])[(index + this.nx)];
      double v3_111 = fk ? v3_110 : fj ? v3_100 : fi ? v3_000 : ((double[])(double[])this.data[(k + 1)])[(index + 1 + this.nx)];
      double dx3 = x - i;
      double dy3 = y - j;
      double dz3 = z - k;
      double z1 = dx3 * (v3_110 * dy3 - v3_100 * (dy3 - 1.0D)) - (dx3 - 1.0D) * (v3_010 * dy3 - v3_000 * (dy3 - 1.0D));
      double z2 = dx3 * (v3_111 * dy3 - v3_101 * (dy3 - 1.0D)) - (dx3 - 1.0D) * (v3_011 * dy3 - v3_001 * (dy3 - 1.0D));
      return z2 * dz3 - z1 * (dz3 - 1.0D);
    }
    return output;
  }

  public double getInterpolatedPixel(double x, double y, double z, byte boundaryConditions)
  {
    double output = 0.0D;
    int i = x >= 0.0D ? (int)x : (int)x - 1;
    int j = y >= 0.0D ? (int)y : (int)y - 1;
    int k = z >= 0.0D ? (int)z : (int)z - 1;
    double dx1;
    switch (getDimension()) {
    case 1:
      double v1_0 = getPixel(i, j, k, boundaryConditions);
      double v1_1 = getPixel(i + 1, j, k, boundaryConditions);
      dx1 = x - i;
    case 2:
      double v2_00 = getPixel(i, j, k, boundaryConditions);
      double v2_10 = getPixel(i + 1, j, k, boundaryConditions);
      double v2_01 = getPixel(i, j + 1, k, boundaryConditions);
      double v2_11 = getPixel(i + 1, j + 1, k, boundaryConditions);
      double dx2 = x - i;
      double dy2 = y - j;
      return dx2 * (v2_11 * dy2 - v2_10 * (dy2 - 1.0D)) - (dx2 - 1.0D) * (v2_01 * dy2 - v2_00 * (dy2 - 1.0D));
    case 3:
      double v3_000 = getPixel(i, j, k, boundaryConditions);
      double v3_100 = getPixel(i + 1, j, k, boundaryConditions);
      double v3_010 = getPixel(i, j + 1, k, boundaryConditions);
      double v3_110 = getPixel(i + 1, j + 1, k, boundaryConditions);
      double v3_001 = getPixel(i, j, k + 1, boundaryConditions);
      double v3_011 = getPixel(i + 1, j, k + 1, boundaryConditions);
      double v3_101 = getPixel(i, j + 1, k + 1, boundaryConditions);
      double v3_111 = getPixel(i + 1, j + 1, k + 1, boundaryConditions);
      double dx3 = x - i;
      double dy3 = y - j;
      double dz3 = z - k;
      double z1 = dx3 * (v3_110 * dy3 - v3_100 * (dy3 - 1.0D)) - (dx3 - 1.0D) * (v3_010 * dy3 - v3_000 * (dy3 - 1.0D));
      double z2 = dx3 * (v3_111 * dy3 - v3_101 * (dy3 - 1.0D)) - (dx3 - 1.0D) * (v3_011 * dy3 - v3_001 * (dy3 - 1.0D));
      return z2 * dz3 - z1 * (dz3 - 1.0D);
    }
    return output;
  }

  public void putPixel(int x, int y, int z, double value)
  {
    if (x >= this.nx) return;
    if (y >= this.ny) return;
    if (z >= this.nz) return;
    if (x < 0) return;
    if (y < 0) return;
    if (z < 0) return;
    ((double[])this.data[z])[(x + y * this.nx)] = value;
  }

  public void getBoundedX(int x, int y, int z, byte[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((byte)(int)tmp[offset]);
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedX(int x, int y, int z, short[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((short)(int)tmp[offset]);
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedX(int x, int y, int z, float[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((float)tmp[offset]);
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedX(int x, int y, int z, double[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = tmp[offset];
        offset++;
      }
    }
    catch (Exception e) {
      throw_get("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedY(int x, int y, int z, byte[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((byte)(int)tmp[offset]);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedY(int x, int y, int z, short[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((short)(int)tmp[offset]);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedY(int x, int y, int z, float[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((float)tmp[offset]);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedY(int x, int y, int z, double[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        buffer[i] = tmp[offset];
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_get("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedZ(int x, int y, int z, byte[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((byte)(int)((double[])(double[])this.data[k])[offset]);
        k++;
      }
    }
    catch (Exception e) {
      throw_get("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedZ(int x, int y, int z, short[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((short)(int)((double[])(double[])this.data[k])[offset]);
        k++;
      }
    }
    catch (Exception e) {
      throw_get("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedZ(int x, int y, int z, float[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((float)((double[])(double[])this.data[k])[offset]);
        k++;
      }
    }
    catch (Exception e) {
      throw_get("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedZ(int x, int y, int z, double[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        buffer[i] = ((double[])(double[])this.data[k])[offset];
        k++;
      }
    }
    catch (Exception e) {
      throw_get("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXY(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((byte)(int)tmp[offset]);
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXY(int x, int y, int z, short[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((short)(int)tmp[offset]);
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXY(int x, int y, int z, float[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((float)tmp[offset]);
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXY(int x, int y, int z, double[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = tmp[offset];
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[z])[offset]);
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[z])[offset]);
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((float)((double[])(double[])this.data[z])[offset]);
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((double[])(double[])this.data[z])[offset];
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedYZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[z])[offset]);
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedYZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[z])[offset]);
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedYZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((float)((double[])(double[])this.data[z])[offset]);
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedYZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          buffer[i][j] = ((double[])(double[])this.data[z])[offset];
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXYZ(int x, int y, int z, byte[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            buffer[i][j][k] = ((byte)(int)tmp[offset]);
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXYZ(int x, int y, int z, short[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            buffer[i][j][k] = ((short)(int)tmp[offset]);
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXYZ(int x, int y, int z, float[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            buffer[i][j][k] = ((float)tmp[offset]);
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBoundedXYZ(int x, int y, int z, double[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            buffer[i][j][k] = tmp[offset];
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void getBlockX(int x, int y, int z, byte[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = x + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = ((byte)(int)tmp[(xp + yp)]);
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockX(int x, int y, int z, short[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = x + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = ((short)(int)tmp[(xp + yp)]);
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockX(int x, int y, int z, float[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = x + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = ((float)tmp[(xp + yp)]);
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockX(int x, int y, int z, double[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = x + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = tmp[(xp + yp)];
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockY(int x, int y, int z, byte[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = y + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = ((byte)(int)tmp[(xp + yp * this.nx)]);
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockY(int x, int y, int z, short[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = y + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = ((short)(int)tmp[(xp + yp * this.nx)]);
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockY(int x, int y, int z, float[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = y + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = ((float)tmp[(xp + yp * this.nx)]);
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockY(int x, int y, int z, double[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = y + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = tmp[(xp + yp * this.nx)];
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockZ(int x, int y, int z, byte[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      for (int i = 0; i < leni; i++) {
        int zp = z + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((byte)(int)((double[])(double[])this.data[zp])[xyp]);
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockZ(int x, int y, int z, short[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      for (int i = 0; i < leni; i++) {
        int zp = z + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((short)(int)((double[])(double[])this.data[zp])[xyp]);
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockZ(int x, int y, int z, float[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      for (int i = 0; i < leni; i++) {
        int zp = z + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((float)((double[])(double[])this.data[zp])[xyp]);
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockZ(int x, int y, int z, double[] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      for (int i = 0; i < leni; i++) {
        int zp = z + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((double[])(double[])this.data[zp])[xyp];
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXY(int x, int y, int z, byte[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = y + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((byte)(int)tmp[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXY(int x, int y, int z, short[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = y + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((short)(int)tmp[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXY(int x, int y, int z, float[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = y + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((float)tmp[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXY(int x, int y, int z, double[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = y + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = tmp[(xp + yp)];
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXZ(int x, int y, int z, byte[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[zp])[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXZ(int x, int y, int z, short[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[zp])[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXZ(int x, int y, int z, float[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((float)((double[])(double[])this.data[zp])[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXZ(int x, int y, int z, double[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = x + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((double[])(double[])this.data[zp])[(xp + yp)];
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockYZ(int x, int y, int z, byte[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = y + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[zp])[(xp + yp * this.nx)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockYZ(int x, int y, int z, short[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = y + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[zp])[(xp + yp * this.nx)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockYZ(int x, int y, int z, float[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = y + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((float)((double[])(double[])this.data[zp])[(xp + yp * this.nx)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockYZ(int x, int y, int z, double[][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      for (int j = 0; j < lenj; j++) {
        int zp = z + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = y + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((double[])(double[])this.data[zp])[(xp + yp * this.nx)];
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXYZ(int x, int y, int z, byte[][][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; k++) {
        int zp = z + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = y + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = x + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = ((byte)(int)tmp[(xp + yp)]);
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXYZ(int x, int y, int z, short[][][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; k++) {
        int zp = z + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = y + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = x + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = ((short)(int)tmp[(xp + yp)]);
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXYZ(int x, int y, int z, float[][][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; k++) {
        int zp = z + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = y + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = x + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = ((float)tmp[(xp + yp)]);
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getBlockXYZ(int x, int y, int z, double[][][] buffer, byte boundaryConditions)
  {
    int xperiod = 0;
    int yperiod = 0;
    int zperiod = 0;
    switch (boundaryConditions) {
    case 2:
      xperiod = this.nx <= 1 ? 1 : 2 * this.nx - 2;
      yperiod = this.ny <= 1 ? 1 : 2 * this.ny - 2;
      zperiod = this.nz <= 1 ? 1 : 2 * this.nz - 2;
      break;
    case 3:
      xperiod = this.nx;
      yperiod = this.ny;
      zperiod = this.nz;
      break;
    default:
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }

    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      for (int k = 0; k < lenk; k++) {
        int zp = z + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = y + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = x + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = tmp[(xp + yp)];
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodX(int x, int y, int z, byte[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = xs + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = ((byte)(int)tmp[(xp + yp)]);
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodX(int x, int y, int z, short[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = xs + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = ((short)(int)tmp[(xp + yp)]);
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodX(int x, int y, int z, float[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = xs + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = ((float)tmp[(xp + yp)]);
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodX(int x, int y, int z, double[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int xp = xs + i;
        while (xp < 0)
          xp += xperiod;
        while (xp >= this.nx) {
          xp = xperiod - xp;
          xp = xp < 0 ? -xp : xp;
        }
        buffer[i] = tmp[(xp + yp)];
      }
    }
    catch (Exception e) {
      throw_get("X", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodY(int x, int y, int z, byte[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = ys + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = ((byte)(int)tmp[(xp + yp * this.nx)]);
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodY(int x, int y, int z, short[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = ys + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = ((short)(int)tmp[(xp + yp * this.nx)]);
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodY(int x, int y, int z, float[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = ys + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = ((float)tmp[(xp + yp * this.nx)]);
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodY(int x, int y, int z, double[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      double[] tmp = (double[])this.data[zp];
      for (int i = 0; i < leni; i++) {
        int yp = ys + i;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        buffer[i] = tmp[(xp + yp * this.nx)];
      }
    }
    catch (Exception e) {
      throw_get("Y", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodZ(int x, int y, int z, byte[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      int zs = z - leni / 2;
      for (int i = 0; i < leni; i++) {
        int zp = zs + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((byte)(int)((double[])(double[])this.data[zp])[xyp]);
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodZ(int x, int y, int z, short[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      int zs = z - leni / 2;
      for (int i = 0; i < leni; i++) {
        int zp = zs + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((short)(int)((double[])(double[])this.data[zp])[xyp]);
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodZ(int x, int y, int z, float[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      int zs = z - leni / 2;
      for (int i = 0; i < leni; i++) {
        int zp = zs + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((float)((double[])(double[])this.data[zp])[xyp]);
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodZ(int x, int y, int z, double[] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      int xyp = xp + yp * this.nx;
      int zs = z - leni / 2;
      for (int i = 0; i < leni; i++) {
        int zp = zs + i;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        buffer[i] = ((double[])(double[])this.data[zp])[xyp];
      }
    }
    catch (Exception e) {
      throw_get("Z", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXY(int x, int y, int z, byte[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = ys + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((byte)(int)tmp[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXY(int x, int y, int z, short[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = ys + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((short)(int)tmp[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXY(int x, int y, int z, float[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = ys + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((float)tmp[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXY(int x, int y, int z, double[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int zp = z;
      while (zp < 0)
        zp += zperiod;
      while (zp >= this.nz) {
        zp = zperiod - zp;
        zp = zp < 0 ? -zp : zp;
      }
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      double[] tmp = (double[])this.data[zp];
      for (int j = 0; j < lenj; j++) {
        int yp = ys + j;
        while (yp < 0)
          yp += yperiod;
        while (yp >= this.ny) {
          yp = yperiod - yp;
          yp = yp < 0 ? -yp : yp;
        }
        yp *= this.nx;
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = tmp[(xp + yp)];
        }
      }
    }
    catch (Exception e) {
      throw_get("XY", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXZ(int x, int y, int z, byte[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[zp])[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXZ(int x, int y, int z, short[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[zp])[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXZ(int x, int y, int z, float[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((float)((double[])(double[])this.data[zp])[(xp + yp)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXZ(int x, int y, int z, double[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int yp = y;
      while (yp < 0)
        yp += yperiod;
      while (yp >= this.ny) {
        yp = yperiod - yp;
        yp = yp < 0 ? -yp : yp;
      }
      yp *= this.nx;
      int xs = x - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - yp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int xp = xs + i;
          while (xp < 0)
            xp += xperiod;
          while (xp >= this.nx) {
            xp = xperiod - xp;
            xp = xp < 0 ? -xp : xp;
          }
          buffer[i][j] = ((double[])(double[])this.data[zp])[(xp + yp)];
        }
      }
    }
    catch (Exception e) {
      throw_get("XZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodYZ(int x, int y, int z, byte[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = ys + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((byte)(int)((double[])(double[])this.data[zp])[(xp + yp * this.nx)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodYZ(int x, int y, int z, short[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = ys + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((short)(int)((double[])(double[])this.data[zp])[(xp + yp * this.nx)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodYZ(int x, int y, int z, float[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = ys + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((float)((double[])(double[])this.data[zp])[(xp + yp * this.nx)]);
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodYZ(int x, int y, int z, double[][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int xp = x;
      while (xp < 0)
        xp += xperiod;
      while (xp >= this.nx) {
        xp = xperiod - xp;
        xp = xp < 0 ? -xp : xp;
      }
      int ys = y - leni / 2;
      int zs = z - lenj / 2;
      for (int j = 0; j < lenj; j++) {
        int zp = zs + j;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        for (int i = 0; i < leni; i++) {
          int yp = ys + i;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          buffer[i][j] = ((double[])(double[])this.data[zp])[(xp + yp * this.nx)];
        }
      }
    }
    catch (Exception e) {
      throw_get("YZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXYZ(int x, int y, int z, byte[][][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      int zs = z - lenk / 2;
      for (int k = 0; k < lenk; k++) {
        int zp = zs + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = ys + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = xs + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = ((byte)(int)tmp[(xp + yp)]);
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXYZ(int x, int y, int z, short[][][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      int zs = z - lenk / 2;
      for (int k = 0; k < lenk; k++) {
        int zp = zs + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = ys + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = xs + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = ((short)(int)tmp[(xp + yp)]);
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXYZ(int x, int y, int z, float[][][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      int zs = z - lenk / 2;
      for (int k = 0; k < lenk; k++) {
        int zp = zs + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = ys + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = xs + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = ((float)tmp[(xp + yp)]);
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void getNeighborhoodXYZ(int x, int y, int z, double[][][] buffer, byte boundaryConditions)
  {
    int xperiod = boundaryConditions == 2 ? 2 * this.nx - 2 : this.nx <= 1 ? 1 : this.nx;
    int yperiod = boundaryConditions == 2 ? 2 * this.ny - 2 : this.ny <= 1 ? 1 : this.ny;
    int zperiod = boundaryConditions == 2 ? 2 * this.nz - 2 : this.nz <= 1 ? 1 : this.nz;
    try
    {
      int offset = x + y * this.nx;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      int xs = x - leni / 2;
      int ys = y - lenj / 2;
      int zs = z - lenk / 2;
      for (int k = 0; k < lenk; k++) {
        int zp = zs + k;
        while (zp < 0)
          zp += zperiod;
        while (zp >= this.nz) {
          zp = zperiod - zp;
          zp = zp < 0 ? -zp : zp;
        }
        double[] tmp = (double[])this.data[zp];
        for (int j = 0; j < lenj; j++) {
          int yp = ys + j;
          while (yp < 0)
            yp += yperiod;
          while (yp >= this.ny) {
            yp = yperiod - yp;
            yp = yp < 0 ? -yp : yp;
          }
          yp *= this.nx;
          for (int i = 0; i < leni; i++) {
            int xp = xs + i;
            while (xp < 0)
              xp += xperiod;
            while (xp >= this.nx) {
              xp = xperiod - xp;
              xp = xp < 0 ? -xp : xp;
            }
            buffer[i][j][k] = tmp[(xp + yp)];
          }
        }
      }
    }
    catch (Exception e) {
      throw_get("XYZ", "Mirror or periodic boundaray conditions", buffer, x, y, z);
    }
  }

  public void putBoundedX(int x, int y, int z, byte[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = (buffer[i] & 0xFF);
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedX(int x, int y, int z, short[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = (buffer[i] & 0xFFFF);
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedX(int x, int y, int z, float[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = buffer[i];
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedX(int x, int y, int z, double[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int offset = x + iinf + y * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = buffer[i];
        offset++;
      }
    }
    catch (Exception e) {
      throw_put("X", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedY(int x, int y, int z, byte[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = (buffer[i] & 0xFF);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedY(int x, int y, int z, short[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = (buffer[i] & 0xFFFF);
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedY(int x, int y, int z, float[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = buffer[i];
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedY(int x, int y, int z, double[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int offset = x + (y + iinf) * this.nx;
      int k = z;
      int leni = buffer.length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      double[] tmp = (double[])this.data[z];
      for (int i = iinf; i < isup; i++) {
        tmp[offset] = buffer[i];
        offset += this.nx;
      }
    }
    catch (Exception e) {
      throw_put("Y", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedZ(int x, int y, int z, byte[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        ((double[])this.data[k])[offset] = (buffer[i] & 0xFF);
        k++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedZ(int x, int y, int z, short[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        ((double[])this.data[k])[offset] = (buffer[i] & 0xFFFF);
        k++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedZ(int x, int y, int z, float[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        ((double[])this.data[k])[offset] = buffer[i];
        k++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedZ(int x, int y, int z, double[] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = z < 0 ? -z : 0;
      int k = z + iinf;
      int offset = x + y * this.nx;
      int leni = buffer.length;
      if (x < 0) return;
      if (y < 0) return;
      if (z + leni < 0) return;
      int isup = z + leni >= this.nz ? this.nz - z : leni;
      for (int i = iinf; i < isup; i++) {
        ((double[])this.data[k])[offset] = buffer[i];
        k++;
      }
    }
    catch (Exception e) {
      throw_put("Z", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXY(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          tmp[offset] = (buffer[i][j] & 0xFF);
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_put("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXY(int x, int y, int z, short[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          tmp[offset] = (buffer[i][j] & 0xFFFF);
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_put("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXY(int x, int y, int z, float[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          tmp[offset] = buffer[i][j];
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_put("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXY(int x, int y, int z, double[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int offset = 0;
      int k = z;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      double[] tmp = (double[])this.data[z];
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + (y + j) * this.nx;
        for (int i = iinf; i < isup; i++) {
          tmp[offset] = buffer[i][j];
          offset++;
        }
      }
    }
    catch (Exception e) {
      throw_put("XY", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = (buffer[i][j] & 0xFF);
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = (buffer[i][j] & 0xFFFF);
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = buffer[i][j];
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x + leni < 0) return;
      if (y < 0) return;
      if (z + lenj < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + iinf + y * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = buffer[i][j];
          offset++;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("YZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedYZ(int x, int y, int z, byte[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = (buffer[i][j] & 0xFF);
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedYZ(int x, int y, int z, short[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = (buffer[i][j] & 0xFFFF);
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedYZ(int x, int y, int z, float[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = buffer[i][j];
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedYZ(int x, int y, int z, double[][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = y < 0 ? -y : 0;
      int jinf = z < 0 ? -z : 0;
      int k = z + jinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      if (x < 0) return;
      if (y + leni < 0) return;
      if (z + lenj < 0) return;
      int isup = y + leni >= this.ny ? this.ny - y : leni;
      int jsup = z + lenj >= this.nz ? this.nz - z : lenj;
      for (int j = jinf; j < jsup; j++) {
        offset = x + (y + iinf) * this.nx;
        for (int i = iinf; i < isup; i++) {
          ((double[])this.data[k])[offset] = buffer[i][j];
          offset += this.nx;
        }
        k++;
      }
    }
    catch (Exception e) {
      throw_put("XZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXYZ(int x, int y, int z, byte[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            tmp[offset] = (buffer[i][j][k] & 0xFF);
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_put("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXYZ(int x, int y, int z, short[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            tmp[offset] = (buffer[i][j][k] & 0xFFFF);
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_put("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXYZ(int x, int y, int z, float[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            tmp[offset] = buffer[i][j][k];
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_put("XYZ", "Bounded check", buffer, x, y, z);
    }
  }

  public void putBoundedXYZ(int x, int y, int z, double[][][] buffer)
  {
    try
    {
      if (x >= this.nx) return;
      if (y >= this.ny) return;
      if (z >= this.nz) return;
      int iinf = x < 0 ? -x : 0;
      int jinf = y < 0 ? -y : 0;
      int kinf = z < 0 ? -z : 0;
      int ko = z + kinf;
      int offset = 0;
      int leni = buffer.length;
      int lenj = buffer[0].length;
      int lenk = buffer[0][0].length;
      if (x + leni < 0) return;
      if (y + lenj < 0) return;
      if (z + lenk < 0) return;
      int isup = x + leni >= this.nx ? this.nx - x : leni;
      int jsup = y + lenj >= this.ny ? this.ny - y : lenj;
      int ksup = z + lenk >= this.nz ? this.nz - z : lenk;
      for (int k = kinf; k < ksup; k++) {
        double[] tmp = (double[])this.data[ko];
        for (int j = jinf; j < jsup; j++) {
          offset = x + iinf + (y + j) * this.nx;
          for (int i = iinf; i < isup; i++) {
            tmp[offset] = buffer[i][j][k];
            offset++;
          }
        }
        ko++;
      }
    }
    catch (Exception e) {
      throw_put("XYZ", "Bounded check", buffer, x, y, z);
    }
  }
}