package imageware;

import ij.ImageStack;
import java.awt.Image;

public class FloatProcess extends FloatPointwise
  implements Process
{
  protected FloatProcess(int nx, int ny, int nz)
  {
    super(nx, ny, nz); } 
  protected FloatProcess(Image image, int mode) { super(image, mode); } 
  protected FloatProcess(ImageStack stack, int mode) {
    super(stack, mode); } 
  protected FloatProcess(ImageStack stack, byte chan) { super(stack, chan); } 
  protected FloatProcess(byte[] array, int mode) {
    super(array, mode); } 
  protected FloatProcess(byte[][] array, int mode) { super(array, mode); } 
  protected FloatProcess(byte[][][] array, int mode) { super(array, mode); } 
  protected FloatProcess(short[] array, int mode) { super(array, mode); } 
  protected FloatProcess(short[][] array, int mode) { super(array, mode); } 
  protected FloatProcess(short[][][] array, int mode) { super(array, mode); } 
  protected FloatProcess(float[] array, int mode) { super(array, mode); } 
  protected FloatProcess(float[][] array, int mode) { super(array, mode); } 
  protected FloatProcess(float[][][] array, int mode) { super(array, mode); } 
  protected FloatProcess(double[] array, int mode) { super(array, mode); } 
  protected FloatProcess(double[][] array, int mode) { super(array, mode); } 
  protected FloatProcess(double[][][] array, int mode) { super(array, mode); }


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
          if (((float[])(float[])this.data[z])[k] < tmp[k])
            ((float[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((float[])(float[])this.data[z])[k] < tmp[k])
            ((float[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((float[])(float[])this.data[z])[k] < tmp[k])
            ((float[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((float[])(float[])this.data[z])[k] < (float)tmp[k])
            ((float[])this.data[z])[k] = ((float)tmp[k]);
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
          if (((float[])(float[])this.data[z])[k] > tmp[k])
            ((float[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((float[])(float[])this.data[z])[k] > tmp[k])
            ((float[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((float[])(float[])this.data[z])[k] > tmp[k])
            ((float[])this.data[z])[k] = tmp[k];
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          if (((float[])(float[])this.data[z])[k] > (float)tmp[k])
            ((float[])this.data[z])[k] = ((float)tmp[k]);
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
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] += tmp[k];
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] += tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] += tmp[k];
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] += (float)tmp[k];
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
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] *= tmp[k];
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] *= tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] *= tmp[k];
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] *= (float)tmp[k];
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
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] -= tmp[k];
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] -= tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] -= tmp[k];
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] -= (float)tmp[k];
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
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] /= tmp[k];
        }
      }
      break;
    case 2:
      for (int z = 0; z < this.nz; z++) {
        short[] tmp = ((ShortSet)imageware).getSliceShort(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] /= tmp[k];
        }
      }
      break;
    case 3:
      for (int z = 0; z < this.nz; z++) {
        float[] tmp = ((FloatSet)imageware).getSliceFloat(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] /= tmp[k];
        }
      }
      break;
    case 4:
      for (int z = 0; z < this.nz; z++) {
        double[] tmp = ((DoubleSet)imageware).getSliceDouble(z);
        for (int k = 0; k < this.nxy; k++) {
          ((float[])this.data[z])[k] /= (float)tmp[k];
        }
      }
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + imageware.getType() + "].\n" + "-------------------------------------------------------\n");
    }
  }
}