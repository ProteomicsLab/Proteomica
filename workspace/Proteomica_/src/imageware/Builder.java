package imageware;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import java.awt.Image;

public class Builder
{
  public static ImageWare wrapOnFocus()
  {
    ImagePlus imp = WindowManager.getCurrentImage();
    return wrap(imp.getStack());
  }

  public static ImageWare wrap(ImagePlus imp)
  {
    return wrap(imp.getStack());
  }

  public static ImageWare wrap(ImageStack stack)
  {
    if (stack == null)
      throw_null();
    ImageProcessor ip = stack.getProcessor(1);
    if ((ip instanceof ByteProcessor)) {
      return new ByteSet(stack, 2);
    }
    if ((ip instanceof ShortProcessor)) {
      return new ShortSet(stack, 2);
    }
    if ((ip instanceof FloatProcessor)) {
      return new FloatSet(stack, 2);
    }

    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to wrap this ImageStack object.\nSupport only the 8-bit, 16-bit and 32-bits type.\n-------------------------------------------------------\n");
  }

  public static ImageWare create(int nx, int ny, int nz, int type)
  {
    switch (type) {
    case 1:
      return new ByteSet(nx, ny, nz);
    case 2:
      return new ShortSet(nx, ny, nz);
    case 3:
      return new FloatSet(nx, ny, nz);
    case 4:
      return new DoubleSet(nx, ny, nz);
    }
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + ".\n" + "-------------------------------------------------------\n");
  }

  public static ImageWare create(Image image)
  {
    if (image == null)
      throw_null();
    return new ByteSet(image, 1);
  }

  public static ImageWare createOnFocus()
  {
    return create(WindowManager.getCurrentImage());
  }

  public static ImageWare create(ImageStack stack)
  {
    if (stack == null)
      throw_null();
    ImageWare wrapped = wrap(stack);
    return wrapped.duplicate();
  }

  public static ImageWare create(ImagePlus imp)
  {
    if (imp == null)
      throw_null();
    ImageWare wrapped = wrap(imp);
    return wrapped.duplicate();
  }

  public static ImageWare[] createColors(ImagePlus imp)
  {
    if (imp == null)
      throw_null();
    return createColors(imp.getStack());
  }

  public static ImageWare[] createColors(ImageStack stack)
  {
    if (stack == null)
      throw_null();
    ImageWare[] color = new ImageWare[3];
    color[0] = new ByteSet(stack, 0);
    color[1] = new ByteSet(stack, 1);
    color[2] = new ByteSet(stack, 2);
    return color;
  }

  public static ImageWare createColorChannel(ImagePlus imp, byte channel)
  {
    if (imp == null)
      throw_null();
    return createColorChannel(imp.getStack(), channel);
  }

  public static ImageWare createColorChannel(ImageStack stack, byte channel)
  {
    if (stack == null)
      throw_null();
    return new ByteSet(stack, channel);
  }

  public static ImageWare create(byte[] object)
  {
    if (object == null)
      throw_null();
    return new ByteSet(object, 1);
  }

  public static ImageWare create(short[] object)
  {
    if (object == null)
      throw_null();
    return new ShortSet(object, 1);
  }

  public static ImageWare create(float[] object)
  {
    if (object == null)
      throw_null();
    return new FloatSet(object, 1);
  }

  public static ImageWare create(double[] object)
  {
    if (object == null)
      throw_null();
    return new DoubleSet(object, 1);
  }

  public static ImageWare create(byte[][] object)
  {
    if (object == null)
      throw_null();
    return new ByteSet(object, 1);
  }

  public static ImageWare create(short[][] object)
  {
    if (object == null)
      throw_null();
    return new ByteSet(object, 1);
  }

  public static ImageWare create(float[][] object)
  {
    if (object == null)
      throw_null();
    return new FloatSet(object, 1);
  }

  public static ImageWare create(double[][] object)
  {
    if (object == null)
      throw_null();
    return new DoubleSet(object, 1);
  }

  public static ImageWare create(byte[][][] object)
  {
    if (object == null)
      throw_null();
    return new ByteSet(object, 1);
  }

  public static ImageWare create(short[][][] object)
  {
    if (object == null)
      throw_null();
    return new ShortSet(object, 1);
  }

  public static ImageWare create(float[][][] object)
  {
    if (object == null)
      throw_null();
    return new FloatSet(object, 1);
  }

  public static ImageWare create(double[][][] object)
  {
    if (object == null)
      throw_null();
    return new DoubleSet(object, 1);
  }

  public static ImageWare create(Image image, int type)
  {
    if (image == null)
      throw_null();
    switch (type) {
    case 1:
      return new ByteSet(image, 1);
    case 2:
      return new ShortSet(image, 1);
    case 3:
      return new FloatSet(image, 1);
    case 4:
      return new DoubleSet(image, 1);
    }
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + ".\n" + "-------------------------------------------------------\n");
  }

  public static ImageWare createOnFocus(int type)
  {
    return create(WindowManager.getCurrentImage(), type);
  }

  public static ImageWare create(ImagePlus imp, int type)
  {
    if (imp == null)
      throw_null();
    switch (type) {
    case 1:
      return new ByteSet(imp.getStack(), 1);
    case 2:
      return new ShortSet(imp.getStack(), 1);
    case 3:
      return new FloatSet(imp.getStack(), 1);
    case 4:
      return new DoubleSet(imp.getStack(), 1);
    }
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + ".\n" + "-------------------------------------------------------\n");
  }

  public static ImageWare create(ImageStack object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare wrapped = wrap(object);
    return wrapped.convert(type);
  }

  public static ImageWare[] createColor(ImagePlus imp, int type)
  {
    if (imp == null)
      throw_null();
    return createColor(imp.getStack(), type);
  }

  public static ImageWare[] createColor(ImageStack stack, int type)
  {
    if (stack == null)
      throw_null();
    ImageWare[] color = new ImageWare[3];
    switch (type) {
    case 1:
      color[0] = new ByteSet(stack, 0);
      color[1] = new ByteSet(stack, 1);
      color[2] = new ByteSet(stack, 2);
      break;
    case 2:
      color[0] = new ShortSet(stack, 0);
      color[1] = new ShortSet(stack, 1);
      color[2] = new ShortSet(stack, 2);
      break;
    case 3:
      color[0] = new FloatSet(stack, 0);
      color[1] = new FloatSet(stack, 1);
      color[2] = new FloatSet(stack, 2);
      break;
    case 4:
      color[0] = new DoubleSet(stack, 0);
      color[1] = new DoubleSet(stack, 1);
      color[2] = new DoubleSet(stack, 2);
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + ".\n" + "-------------------------------------------------------\n");
    }

    return color;
  }

  public static ImageWare createColorChannel(ImagePlus imp, byte channel, int type)
  {
    if (imp == null)
      throw_null();
    return createColorChannel(imp.getStack(), channel, type);
  }

  public static ImageWare createColorChannel(ImageStack stack, byte channel, int type)
  {
    if (stack == null)
      throw_null();
    ImageWare out = null;
    switch (type) {
    case 1:
      out = new ByteSet(stack, channel);
      break;
    case 2:
      out = new ShortSet(stack, channel);
      break;
    case 3:
      out = new FloatSet(stack, channel);
      break;
    case 4:
      out = new DoubleSet(stack, channel);
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnknown type " + type + ".\n" + "-------------------------------------------------------\n");
    }

    return out;
  }

  public static ImageWare create(byte[] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, 1, 1, type);
    out.putX(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(short[] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, 1, 1, type);
    out.putX(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(float[] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, 1, 1, type);
    out.putX(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(double[] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, 1, 1, type);
    out.putX(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(byte[][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, 1, type);
    out.putXY(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(short[][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, 1, type);
    out.putXY(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(float[][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, 1, type);
    out.putXY(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(double[][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, 1, type);
    out.putXY(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(byte[][][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, object[0][0].length, type);
    out.putXYZ(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(short[][][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, object[0][0].length, type);
    out.putXYZ(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(float[][][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, object[0][0].length, type);
    out.putXYZ(0, 0, 0, object);
    return out;
  }

  public static ImageWare create(double[][][] object, int type)
  {
    if (object == null)
      throw_null();
    ImageWare out = createType(object.length, object[0].length, object[0][0].length, type);
    out.putXYZ(0, 0, 0, object);
    return out;
  }

  private static ImageWare createType(int nx, int ny, int nz, int type)
  {
    ImageWare out = null;
    switch (type) {
    case 1:
      out = new ByteSet(nx, ny, nz);
      break;
    case 2:
      out = new ShortSet(nx, ny, nz);
      break;
    case 3:
      out = new FloatSet(nx, ny, nz);
      break;
    case 4:
      out = new DoubleSet(nx, ny, nz);
      break;
    default:
      throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to create this object.\n-------------------------------------------------------\n");
    }

    return out;
  }

  private static void throw_null() {
    throw new ArrayStoreException("\n-------------------------------------------------------\nError in imageware package\nUnable to wrap the ImagePlus.\nThe object parameter is null.\n-------------------------------------------------------\n");
  }
}