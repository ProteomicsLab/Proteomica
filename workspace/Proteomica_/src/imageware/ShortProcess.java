package imageware;

import ij.ImageStack;
import java.awt.Image;

public class ShortProcess extends ShortPointwise
  implements Process
{
  protected ShortProcess(int nx, int ny, int nz)
  {
    super(nx, ny, nz); } 
  protected ShortProcess(Image image, int mode) { super(image, mode); } 
  protected ShortProcess(ImageStack stack, int mode) {
    super(stack, mode); } 
  protected ShortProcess(ImageStack stack, byte chan) { super(stack, chan); } 
  protected ShortProcess(byte[] array, int mode) {
    super(array, mode); } 
  protected ShortProcess(byte[][] array, int mode) { super(array, mode); } 
  protected ShortProcess(byte[][][] array, int mode) { super(array, mode); } 
  protected ShortProcess(short[] array, int mode) { super(array, mode); } 
  protected ShortProcess(short[][] array, int mode) { super(array, mode); } 
  protected ShortProcess(short[][][] array, int mode) { super(array, mode); } 
  protected ShortProcess(float[] array, int mode) { super(array, mode); } 
  protected ShortProcess(float[][] array, int mode) { super(array, mode); } 
  protected ShortProcess(float[][][] array, int mode) { super(array, mode); } 
  protected ShortProcess(double[] array, int mode) { super(array, mode); } 
  protected ShortProcess(double[][] array, int mode) { super(array, mode); } 
  protected ShortProcess(double[][][] array, int mode) { super(array, mode); }


  public void smoothGaussian(double sigma)
  {
    smoothGaussian(sigma, sigma, sigma);
  }

  public void smoothGaussian(double sigmaX, double sigmaY, double sigmaZ)
  {
    int n = 3;
    double N = n;
    double step = 1.0D / (this.nx + this.ny + this.nz);
    double[] poles = new double[n];

    if ((this.nx > 1) && (sigmaX > 0.0D)) {
      double s2 = sigmaX * sigmaX;
      double alpha = 1.0D + N / s2 - Math.sqrt(N * N + 2.0D * N * s2) / s2;
      double tmp96_95 = (poles[2] = alpha); poles[1] = tmp96_95; poles[0] = tmp96_95;
      double[] line = new double[this.nx];
      for (int z = 0; z < this.nz; z++) {
        for (int y = 0; y < this.ny; y++) {
          getX(0, y, z, line);
          putX(0, y, z, Convolver.convolveIIR(line, poles));
        }
      }
    }

    if ((this.ny > 1) && (sigmaY > 0.0D)) {
      double s2 = sigmaY * sigmaY;
      double alpha = 1.0D + N / s2 - Math.sqrt(N * N + 2.0D * N * s2) / s2;
      double tmp233_232 = (poles[2] = alpha); poles[1] = tmp233_232; poles[0] = tmp233_232;
      double[] line = new double[this.ny];
      for (int x = 0; x < this.nx; x++) {
        for (int z = 0; z < this.nz; z++) {
          getY(x, 0, z, line);
          putY(x, 0, z, Convolver.convolveIIR(line, poles));
        }
      }
    }

    if ((this.nz > 1) && (sigmaZ > 0.0D)) {
      double s2 = sigmaZ * sigmaZ;
      double alpha = 1.0D + N / s2 - Math.sqrt(N * N + 2.0D * N * s2) / s2;
      double tmp373_372 = (poles[2] = alpha); poles[1] = tmp373_372; poles[0] = tmp373_372;
      double[] line = new double[this.nz];
      for (int y = 0; y < this.ny; y++)
        for (int x = 0; x < this.nx; x++) {
          getZ(x, y, 0, line);
          putZ(x, y, 0, Convolver.convolveIIR(line, poles));
        }
    }
  }

  public void max(ImageWare imageware)
  {
    if (!isSameSize(imageware)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to get the maximum because the two operands are not the same size.\n[" + this.nx + "," + this.ny + "," + "," + this.nz + "] != " + "[" + imageware.getSizeX() + "," + imageware.getSizeY() + "," + imageware.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    switch (imageware.getType()) {
    case 1:
      for (int z = 0; z < this.nz; z++) {
        byte[] tmp = ((ByteSet)imageware).getSliceByte(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] < (short)tmp[k])
            ((short[])this.data[z])[k] = ((short)tmp[k]);
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] < tmp[k])
            ((short[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] < (short)(int)tmp[k])
            ((short[])this.data[z])[k] = ((short)(int)tmp[k]);
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] < (short)(int)tmp[k])
            ((short[])this.data[z])[k] = ((short)(int)tmp[k]);
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }

  public void min(ImageWare imageware)
  {
    if (!isSameSize(imageware)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to get the minimum because the two operands are not the same size.\n[" + this.nx + "," + this.ny + "," + "," + this.nz + "] != " + "[" + imageware.getSizeX() + "," + imageware.getSizeY() + "," + imageware.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    switch (imageware.getType()) {
    case 1:
      for (int z = 0; z < this.nz; z++) {
        byte[] tmp = ((ByteSet)imageware).getSliceByte(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] > (short)tmp[k])
            ((short[])this.data[z])[k] = ((short)tmp[k]);
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] > tmp[k])
            ((short[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] > (short)(int)tmp[k])
            ((short[])this.data[z])[k] = ((short)(int)tmp[k]);
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((short[])(short[])this.data[z])[k] > (short)(int)tmp[k])
            ((short[])this.data[z])[k] = ((short)(int)tmp[k]);
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }

  public void add(ImageWare imageware)
  {
    if (!isSameSize(imageware)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to add because the two operands are not the same size.\n[" + this.nx + "," + this.ny + "," + "," + this.nz + "] != " + "[" + imageware.getSizeX() + "," + imageware.getSizeY() + "," + imageware.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    switch (imageware.getType()) {
    case 1:
      for (int z = 0; z < this.nz; z++) {
        byte[] tmp = ((ByteSet)imageware).getSliceByte(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp205_203 = k;
          short[] tmp205_200 = ((short[])this.data[z]); tmp205_200[tmp205_203] = ((short)(tmp205_200[tmp205_203] + (short)tmp[k]));
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp275_273 = k;
          short[] tmp275_270 = ((short[])this.data[z]); tmp275_270[tmp275_273] = ((short)(tmp275_270[tmp275_273] + tmp[k]));
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp344_342 = k;
          short[] tmp344_339 = ((short[])this.data[z]); tmp344_339[tmp344_342] = ((short)(tmp344_339[tmp344_342] + (short)(int)tmp[k]));
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp415_413 = k;
          short[] tmp415_410 = ((short[])this.data[z]); tmp415_410[tmp415_413] = ((short)(tmp415_410[tmp415_413] + (short)(int)tmp[k]));
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }

  public void multiply(ImageWare imageware)
  {
    if (!isSameSize(imageware)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to multiply because the two operands are not the same size.\n[" + this.nx + "," + this.ny + "," + "," + this.nz + "] != " + "[" + imageware.getSizeX() + "," + imageware.getSizeY() + "," + imageware.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    switch (imageware.getType()) {
    case 1:
      for (int z = 0; z < this.nz; z++) {
        byte[] tmp = ((ByteSet)imageware).getSliceByte(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp205_203 = k;
          short[] tmp205_200 = ((short[])this.data[z]); tmp205_200[tmp205_203] = ((short)(tmp205_200[tmp205_203] * (short)tmp[k]));
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp275_273 = k;
          short[] tmp275_270 = ((short[])this.data[z]); tmp275_270[tmp275_273] = ((short)(tmp275_270[tmp275_273] * tmp[k]));
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp344_342 = k;
          short[] tmp344_339 = ((short[])this.data[z]); tmp344_339[tmp344_342] = ((short)(tmp344_339[tmp344_342] * (short)(int)tmp[k]));
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp415_413 = k;
          short[] tmp415_410 = ((short[])this.data[z]); tmp415_410[tmp415_413] = ((short)(tmp415_410[tmp415_413] * (short)(int)tmp[k]));
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }

  public void subtract(ImageWare imageware)
  {
    if (!isSameSize(imageware)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to subtract because the two operands are not the same size.\n[" + this.nx + "," + this.ny + "," + "," + this.nz + "] != " + "[" + imageware.getSizeX() + "," + imageware.getSizeY() + "," + imageware.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    switch (imageware.getType()) {
    case 1:
      for (int z = 0; z < this.nz; z++) {
        byte[] tmp = ((ByteSet)imageware).getSliceByte(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp205_203 = k;
          short[] tmp205_200 = ((short[])this.data[z]); tmp205_200[tmp205_203] = ((short)(tmp205_200[tmp205_203] - (short)tmp[k]));
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp275_273 = k;
          short[] tmp275_270 = ((short[])this.data[z]); tmp275_270[tmp275_273] = ((short)(tmp275_270[tmp275_273] - tmp[k]));
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp344_342 = k;
          short[] tmp344_339 = ((short[])this.data[z]); tmp344_339[tmp344_342] = ((short)(tmp344_339[tmp344_342] - (short)(int)tmp[k]));
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp415_413 = k;
          short[] tmp415_410 = ((short[])this.data[z]); tmp415_410[tmp415_413] = ((short)(tmp415_410[tmp415_413] - (short)(int)tmp[k]));
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }

  public void divide(ImageWare imageware)
  {
    if (!isSameSize(imageware)) {
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to divide because the two operands are not the same size.\n[" + this.nx + "," + this.ny + "," + "," + this.nz + "] != " + "[" + imageware.getSizeX() + "," + imageware.getSizeY() + "," + imageware.getSizeZ() + "].\n" + "-------------------------------------------------------\n");
    }

    switch (imageware.getType()) {
    case 1:
      for (int z = 0; z < this.nz; z++) {
        byte[] tmp = ((ByteSet)imageware).getSliceByte(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp205_203 = k;
          short[] tmp205_200 = ((short[])this.data[z]); tmp205_200[tmp205_203] = ((short)(tmp205_200[tmp205_203] / (short)tmp[k]));
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp275_273 = k;
          short[] tmp275_270 = ((short[])this.data[z]); tmp275_270[tmp275_273] = ((short)(tmp275_270[tmp275_273] / tmp[k]));
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp344_342 = k;
          short[] tmp344_339 = ((short[])this.data[z]); tmp344_339[tmp344_342] = ((short)(tmp344_339[tmp344_342] / (short)(int)tmp[k]));
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++)
        {
          int tmp415_413 = k;
          short[] tmp415_410 = ((short[])this.data[z]); tmp415_410[tmp415_413] = ((short)(tmp415_410[tmp415_413] / (short)(int)tmp[k]));
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }
}