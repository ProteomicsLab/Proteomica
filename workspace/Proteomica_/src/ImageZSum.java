 
import ij.*; 
import ij.gui.GenericDialog; 
import ij.process.*;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.filter.*; 
import ij.measure.Measurements;
import java.lang.*; 
import java.awt.*; 
import java.awt.event.*; 
import java.util.Arrays;

/** This plugin performs a z-projection of the input stack. Type of
    output image is same as type of input image.

    @author Patrick Kelly <phkelly@ucsd.edu> 
    @modified by Bhushan N Kharbikar for implementation of Proteomica 
    */

public class ImageZSum implements PlugIn {
    public static final int AVG_METHOD = 0; 
    public static final int MAX_METHOD = 1;
    public static final int MIN_METHOD = 2;
    public static final int SUM_METHOD = 3;
	public static final int SD_METHOD = 4;
	public static final int MEDIAN_METHOD = 5;
	public static final String[] METHODS = 
		{"Average Intensity", "Max Intensity", "Min Intensity", "Sum Slices", "Standard Deviation", "Median"}; 
    private static final String METHOD_KEY = "zproject.method";
    private int method = (int)Prefs.get(METHOD_KEY, AVG_METHOD);

    private static final int BYTE_TYPE  = 0; 
    private static final int SHORT_TYPE = 1; 
    private static final int FLOAT_TYPE = 2;
    
    public static final String lutMessage =
    	"Stacks with inverter LUTs may not project correctly.\n"
    	+"To create a standard LUT, invert the stack (Edit/Invert)\n"
    	+"and invert the LUT (Image/Lookup Tables/Invert LUT)."; 

    /** Image to hold z-projection. */
    private ImagePlus projImage = null; 

    /** Image stack to project. */
    private ImagePlus imp = null; 

    /** Projection starts from this slice. */
    private int startSlice = 1;
    /** Projection ends at this slice. */
    private int stopSlice = 1;
    /** Project all time points? */
    private boolean allTimeFrames = true;
    
    private String color = "";
    private boolean isHyperstack;
    private int increment = 1;
    private int sliceCount;

    public ImageZSum() {
    }

    /** Construction of ZProjector with image to be projected. */
    public ImageZSum(ImagePlus imp) {
		setImage(imp); 
    }

    /** Explicitly set image to be projected. This is useful if
	ZProjection_ object is to be used not as a plugin but as a
	stand alone processing object.  */
    public void setImage(ImagePlus imp) {
    	this.imp = imp; 
		startSlice = 1; 
		stopSlice = imp.getStackSize(); 
    }

    public void setStartSlice(int slice) {
		if(imp==null || slice < 1 || slice > imp.getStackSize())
	    	return; 
		startSlice = slice; 
    }

    public void setStopSlice(int slice) {
		if(imp==null || slice < 1 || slice > imp.getStackSize())
	    	return; 
		stopSlice = slice; 
    }

	public void setMethod(int projMethod){
		method = projMethod;
	}
    
    /** Retrieve results of most recent projection operation.*/
    public ImagePlus getProjection() {
		return projImage; 
    }

    public void run(String arg) {
		imp = IJ.getImage();
		int stackSize = imp.getStackSize();
		if (imp==null) {
	    	IJ.noImage(); 
	    	return; 
		}

		//  Make sure input image is a stack.
		if(stackSize==1) {
	    	IJ.error("Z Project", "Stack required"); 
	    	return; 
		}
	
		//  Check for inverting LUT.
		if (imp.getProcessor().isInvertedLut()) {
	    	if (!IJ.showMessageWithCancel("ZProjection", lutMessage))
	    		return; 
		}

		// Set default bounds.
		int channels = imp.getNChannels();
		int frames = imp.getNFrames();
		int slices = imp.getNSlices();
		isHyperstack = imp.isHyperStack()||( ij.macro.Interpreter.isBatchMode()&&((frames>1&&frames<stackSize)||(slices>1&&slices<stackSize)));
		boolean simpleComposite = channels==stackSize;
		if (simpleComposite)
			isHyperstack = false;
		startSlice = 1; 
		if (isHyperstack) {
			int nSlices = imp.getNSlices();
			if (nSlices>1)
				stopSlice = nSlices;
			else
				stopSlice = imp.getNFrames();
		} else
			stopSlice  = stackSize;
			
		// Build control dialog
		GenericDialog gd = buildControlDialog(startSlice,stopSlice); 
		//gd.showDialog(); 
		//if(gd.wasCanceled()) return; 

		if (!imp.lock()) return;   // exit if in use
		long tstart = System.currentTimeMillis();
		setStartSlice((int)startSlice); //modified by Bhushan N Kharbikar
		setStopSlice((int)stopSlice);  // modified by Bhushan N Kharbikar
		//method = gd.getNextChoiceIndex();
		//Prefs.set(METHOD_KEY, method);
		Prefs.set(METHOD_KEY, 3);
		if (isHyperstack) {
			//allTimeFrames = imp.getNFrames()>1&&imp.getNSlices()>1?gd.getNextBoolean():false;
			allTimeFrames = false;
			doHyperStackProjection(allTimeFrames);
		} else if (imp.getType()==ImagePlus.COLOR_RGB)
			doRGBProjection();
		else 
			doProjection(); 

		if (arg.equals("") && projImage!=null) {
			long tstop = System.currentTimeMillis();
			projImage.setCalibration(imp.getCalibration());
			if (simpleComposite) IJ.run(projImage, "Grays", "");
			projImage.show("ZProjector: " +IJ.d2s((tstop-tstart)/1000.0,2)+" seconds");
		}

		imp.unlock();
		IJ.register(ImageZSum.class);
		return; 
    }
    
    public void doRGBProjection() {
		doRGBProjection(imp.getStack());
    }

    private void doRGBProjection(ImageStack stack) {
        ImageStack[] channels = ChannelSplitter.splitRGB(stack, true);
        ImagePlus red = new ImagePlus("Red", channels[0]);
        ImagePlus green = new ImagePlus("Green", channels[1]);
        ImagePlus blue = new ImagePlus("Blue", channels[2]);
        imp.unlock();
        ImagePlus saveImp = imp;
        imp = red;
		color = "(red)"; doProjection();
		ImagePlus red2 = projImage;
        imp = green;
		color = "(green)"; doProjection();
		ImagePlus green2 = projImage;
        imp = blue;
		color = "(blue)"; doProjection();
		ImagePlus blue2 = projImage;
        int w = red2.getWidth(), h = red2.getHeight(), d = red2.getStackSize();
        if (method==SD_METHOD) {
        	ImageProcessor r = red2.getProcessor();
        	ImageProcessor g = green2.getProcessor();
        	ImageProcessor b = blue2.getProcessor();
        	double max = 0;
        	double rmax = r.getStatistics().max; if (rmax>max) max=rmax;
        	double gmax = g.getStatistics().max; if (gmax>max) max=gmax;
        	double bmax = b.getStatistics().max; if (bmax>max) max=bmax;
        	double scale = 255/max;
        	r.multiply(scale); g.multiply(scale); b.multiply(scale);
        	red2.setProcessor(r.convertToByte(false));
        	green2.setProcessor(g.convertToByte(false));
        	blue2.setProcessor(b.convertToByte(false));
        }
        RGBStackMerge merge = new RGBStackMerge();
        ImageStack stack2 = merge.mergeStacks(w, h, d, red2.getStack(), green2.getStack(), blue2.getStack(), true);
        imp = saveImp;
        projImage = new ImagePlus(makeTitle(), stack2);
    }

    /** Builds dialog to query users for projection parameters.
	@param start starting slice to display
	@param stop last slice */
    protected GenericDialog buildControlDialog(int start, int stop) {
		GenericDialog gd = new GenericDialog("ZProjection",IJ.getInstance()); 
		gd.addNumericField("Start slice:",startSlice,0/*digits*/); 
		gd.addNumericField("Stop slice:",stopSlice,0/*digits*/);
		gd.addChoice("Projection type", METHODS, METHODS[method]); 
		if (isHyperstack && imp.getNFrames()>1&& imp.getNSlices()>1)
			gd.addCheckbox("All time frames", allTimeFrames); 
		return gd; 
    }

    /** Performs actual projection using specified method. */
    public void doProjection() {
		if (imp==null)
			return;
		sliceCount = 0;
		if (method<AVG_METHOD || method>MEDIAN_METHOD)
			method = AVG_METHOD;
    	for (int slice=startSlice; slice<=stopSlice; slice+=increment)
    		sliceCount++;
		if (method==MEDIAN_METHOD) {
			projImage = doMedianProjection();
			return;
		} 
		
		// Create new float processor for projected pixels.
		FloatProcessor fp = new FloatProcessor(imp.getWidth(),imp.getHeight()); 
		ImageStack stack = imp.getStack();
		RayFunction rayFunc = getRayFunction(method, fp);
		if (IJ.debugMode==true) {
	    	IJ.log("\nProjecting stack from: "+startSlice
		     	+" to: "+stopSlice); 
		}

		// Determine type of input image. Explicit determination of
		// processor type is required for subsequent pixel
		// manipulation.  This approach is more efficient than the
		// more general use of ImageProcessor's getPixelValue and
		// putPixel methods.
		int ptype; 
		if(stack.getProcessor(1) instanceof ByteProcessor) ptype = BYTE_TYPE; 
		else if(stack.getProcessor(1) instanceof ShortProcessor) ptype = SHORT_TYPE; 
		else if(stack.getProcessor(1) instanceof FloatProcessor) ptype = FLOAT_TYPE; 
		else {
	    	IJ.error("Z Project", "Non-RGB stack required"); 
	    	return; 
		}

		// Do the projection.
		for(int n=startSlice; n<=stopSlice; n+=increment) {
	    	IJ.showStatus("ZProjection " + color +": " + n + "/" + stopSlice);
	    	IJ.showProgress(n-startSlice, stopSlice-startSlice);
	    	projectSlice(stack.getPixels(n), rayFunc, ptype);
		}

		// Finish up projection.
		if (method==SUM_METHOD) {
			fp.resetMinAndMax();
			projImage = new ImagePlus(makeTitle(),fp); 
		} else if (method==SD_METHOD) {
			rayFunc.postProcess();
			fp.resetMinAndMax();
			projImage = new ImagePlus(makeTitle(), fp); 
		} else {
			rayFunc.postProcess(); 
			projImage = makeOutputImage(imp, fp, ptype);
		}

		if(projImage==null)
	    	IJ.error("Z Project", "Error computing projection.");
    }

	public void doHyperStackProjection(boolean allTimeFrames) {
		int start = startSlice;
		int stop = stopSlice;
		int firstFrame = 1;
		int lastFrame = imp.getNFrames();
		if (!allTimeFrames)
			firstFrame = lastFrame = imp.getFrame();
		ImageStack stack = new ImageStack(imp.getWidth(), imp.getHeight());
		int channels = imp.getNChannels();
		int slices = imp.getNSlices();
		if (slices==1) {
			slices = imp.getNFrames();
			firstFrame = lastFrame = 1;
		}
		int frames = lastFrame-firstFrame+1;
		increment = channels;
		boolean rgb = imp.getBitDepth()==24;
		for (int frame=firstFrame; frame<=lastFrame; frame++) {
			for (int channel=1; channel<=channels; channel++) {
				startSlice = (frame-1)*channels*slices + (start-1)*channels + channel;
				stopSlice = (frame-1)*channels*slices + (stop-1)*channels + channel;
				if (rgb)
					doHSRGBProjection(imp);
				else
					doProjection();
				stack.addSlice(null, projImage.getProcessor());
			}
		}
        projImage = new ImagePlus(makeTitle(), stack);
        projImage.setDimensions(channels, 1, frames);
        if (channels>1) {
           	projImage = new CompositeImage(projImage, 0);
        	((CompositeImage)projImage).copyLuts(imp);
      		if (method==SUM_METHOD || method==SD_METHOD)
        			((CompositeImage)projImage).resetDisplayRanges();
        }
        if (frames>1)
        	projImage.setOpenAsHyperStack(true);
        IJ.showProgress(1, 1);
	}
	
	private void doHSRGBProjection(ImagePlus rgbImp) {
		ImageStack stack = rgbImp.getStack();
		ImageStack stack2 = new ImageStack(stack.getWidth(), stack.getHeight());
		for (int i=startSlice; i<=stopSlice; i++)
			stack2.addSlice(null, stack.getProcessor(i));
		startSlice = 1;
		stopSlice = stack2.getSize();
		doRGBProjection(stack2);
	}

 	private RayFunction getRayFunction(int method, FloatProcessor fp) {
 		switch (method) {
 			case AVG_METHOD: case SUM_METHOD:
	    		return new AverageIntensity(fp, sliceCount); 
			case MAX_METHOD:
	    		return new MaxIntensity(fp);
	    	case MIN_METHOD:
	    		return new MinIntensity(fp); 
			case SD_METHOD:
	    		return new StandardDeviation(fp, sliceCount); 
			default:
	    		IJ.error("Z Project", "Unknown method.");
	    		return null;
	    }
	}

    /** Generate output image whose type is same as input image. */
    private ImagePlus makeOutputImage(ImagePlus imp, FloatProcessor fp, int ptype) {
		int width = imp.getWidth(); 
		int height = imp.getHeight(); 
		float[] pixels = (float[])fp.getPixels(); 
		ImageProcessor oip=null; 

		// Create output image consistent w/ type of input image.
		int size = pixels.length;
		switch (ptype) {
			case BYTE_TYPE:
				oip = imp.getProcessor().createProcessor(width,height);
				byte[] pixels8 = (byte[])oip.getPixels(); 
				for(int i=0; i<size; i++)
					pixels8[i] = (byte)pixels[i];
				break;
			case SHORT_TYPE:
				oip = imp.getProcessor().createProcessor(width,height);
				short[] pixels16 = (short[])oip.getPixels(); 
				for(int i=0; i<size; i++)
					pixels16[i] = (short)pixels[i];
				break;
			case FLOAT_TYPE:
				oip = new FloatProcessor(width, height, pixels, null);
				break;
		}
	
		// Adjust for display.
	    // Calling this on non-ByteProcessors ensures image
	    // processor is set up to correctly display image.
	    oip.resetMinAndMax(); 

		// Create new image plus object. Don't use
		// ImagePlus.createImagePlus here because there may be
		// attributes of input image that are not appropriate for
		// projection.
		return new ImagePlus(makeTitle(), oip); 
    }

    /** Handles mechanics of projection by selecting appropriate pixel
	array type. We do this rather than using more general
	ImageProcessor getPixelValue() and putPixel() methods because
	direct manipulation of pixel arrays is much more efficient.  */
	private void projectSlice(Object pixelArray, RayFunction rayFunc, int ptype) {
		switch(ptype) {
			case BYTE_TYPE:
	    		rayFunc.projectSlice((byte[])pixelArray); 
	    		break; 
			case SHORT_TYPE:
	    		rayFunc.projectSlice((short[])pixelArray); 
	    		break; 
			case FLOAT_TYPE:
	    		rayFunc.projectSlice((float[])pixelArray); 
	    		break; 
		}
    }
    
    String makeTitle() {
    	String prefix = "AVG_";
 		switch (method) {
 			case SUM_METHOD: prefix = "SUM_"; break;
			case MAX_METHOD: prefix = "MAX_"; break;
	    	case MIN_METHOD: prefix = "MIN_"; break;
			case SD_METHOD:  prefix = "STD_"; break;
			case MEDIAN_METHOD:  prefix = "MED_"; break;
	    }
    	return WindowManager.makeUniqueName(prefix+imp.getTitle());
    }

	ImagePlus doMedianProjection() {
		IJ.showStatus("Calculating median...");
		ImageStack stack = imp.getStack();
		ImageProcessor[] slices = new ImageProcessor[sliceCount];
		int index = 0;
		for (int slice=startSlice; slice<=stopSlice; slice+=increment)
			slices[index++] = stack.getProcessor(slice);
		ImageProcessor ip2 = slices[0].duplicate();
		ip2 = ip2.convertToFloat();
		float[] values = new float[sliceCount];
		int width = ip2.getWidth();
		int height = ip2.getHeight();
		int inc = Math.max(height/30, 1);
		for (int y=0; y<height; y++) {
			if (y%inc==0) IJ.showProgress(y, height-1);
			for (int x=0; x<width; x++) {
				for (int i=0; i<sliceCount; i++)
				values[i] = slices[i].getPixelValue(x, y);
				ip2.putPixelValue(x, y, median(values));
			}
		}
		if (imp.getBitDepth()==8)
			ip2 = ip2.convertToByte(false);
		IJ.showProgress(1, 1);
		return new ImagePlus(makeTitle(), ip2);
	}

	float median(float[] a) {
		Arrays.sort(a);
		int middle = a.length/2;
		if ((a.length&1)==0) //even
			return (a[middle-1] + a[middle])/2f;
		else
			return a[middle];
	}

     /** Abstract class that specifies structure of ray
	function. Preprocessing should be done in derived class
	constructors.
	*/
    abstract class RayFunction {
		/** Do actual slice projection for specific data types. */
		public abstract void projectSlice(byte[] pixels);
		public abstract void projectSlice(short[] pixels);
		public abstract void projectSlice(float[] pixels);
		
		/** Perform any necessary post processing operations, e.g.
	    	averging values. */
		public void postProcess() {}

    } // end RayFunction


    /** Compute average intensity projection. */
    class AverageIntensity extends RayFunction {
     	private float[] fpixels;
 		private int num, len; 

		/** Constructor requires number of slices to be
	    	projected. This is used to determine average at each
	    	pixel. */
		public AverageIntensity(FloatProcessor fp, int num) {
			fpixels = (float[])fp.getPixels();
			len = fpixels.length;
	    	this.num = num;
		}

		public void projectSlice(byte[] pixels) {
	    	for(int i=0; i<len; i++)
				fpixels[i] += (pixels[i]&0xff); 
		}

		public void projectSlice(short[] pixels) {
	    	for(int i=0; i<len; i++)
				fpixels[i] += pixels[i]&0xffff;
		}

		public void projectSlice(float[] pixels) {
	    	for(int i=0; i<len; i++)
				fpixels[i] += pixels[i]; 
		}

		public void postProcess() {
			float fnum = num;
	    	for(int i=0; i<len; i++)
				fpixels[i] /= fnum;
		}

    } // end AverageIntensity


     /** Compute max intensity projection. */
    class MaxIntensity extends RayFunction {
    	private float[] fpixels;
 		private int len; 

		/** Simple constructor since no preprocessing is necessary. */
		public MaxIntensity(FloatProcessor fp) {
			fpixels = (float[])fp.getPixels();
			len = fpixels.length;
			for (int i=0; i<len; i++)
				fpixels[i] = -Float.MAX_VALUE;
		}

		public void projectSlice(byte[] pixels) {
	    	for(int i=0; i<len; i++) {
				if((pixels[i]&0xff)>fpixels[i])
		    		fpixels[i] = (pixels[i]&0xff); 
	    	}
		}

		public void projectSlice(short[] pixels) {
	    	for(int i=0; i<len; i++) {
				if((pixels[i]&0xffff)>fpixels[i])
		    		fpixels[i] = pixels[i]&0xffff;
	    	}
		}

		public void projectSlice(float[] pixels) {
	    	for(int i=0; i<len; i++) {
				if(pixels[i]>fpixels[i])
		    		fpixels[i] = pixels[i]; 
	    	}
		}
		
    } // end MaxIntensity

     /** Compute min intensity projection. */
    class MinIntensity extends RayFunction {
    	private float[] fpixels;
 		private int len; 

		/** Simple constructor since no preprocessing is necessary. */
		public MinIntensity(FloatProcessor fp) {
			fpixels = (float[])fp.getPixels();
			len = fpixels.length;
			for (int i=0; i<len; i++)
				fpixels[i] = Float.MAX_VALUE;
		}

		public void projectSlice(byte[] pixels) {
	    	for(int i=0; i<len; i++) {
				if((pixels[i]&0xff)<fpixels[i])
		    		fpixels[i] = (pixels[i]&0xff); 
	    	}
		}

		public void projectSlice(short[] pixels) {
	    	for(int i=0; i<len; i++) {
				if((pixels[i]&0xffff)<fpixels[i])
		    		fpixels[i] = pixels[i]&0xffff;
	    	}
		}

		public void projectSlice(float[] pixels) {
	    	for(int i=0; i<len; i++) {
				if(pixels[i]<fpixels[i])
		    		fpixels[i] = pixels[i]; 
	    	}
		}
		
    } // end MaxIntensity


    /** Compute standard deviation projection. */
    class StandardDeviation extends RayFunction {
    	private float[] result;
    	private double[] sum, sum2;
		private int num,len; 

		public StandardDeviation(FloatProcessor fp, int num) {
			result = (float[])fp.getPixels();
			len = result.length;
		    this.num = num;
			sum = new double[len];
			sum2 = new double[len];
		}
	
		public void projectSlice(byte[] pixels) {
			int v;
		    for(int i=0; i<len; i++) {
		    	v = pixels[i]&0xff;
				sum[i] += v;
				sum2[i] += v*v;
			} 
		}
	
		public void projectSlice(short[] pixels) {
			double v;
		    for(int i=0; i<len; i++) {
		    	v = pixels[i]&0xffff;
				sum[i] += v;
				sum2[i] += v*v;
			} 
		}
	
		public void projectSlice(float[] pixels) {
			double v;
		    for(int i=0; i<len; i++) {
		    	v = pixels[i];
				sum[i] += v;
				sum2[i] += v*v;
			} 
		}
	
		public void postProcess() {
			double stdDev;
			double n = num;
		    for(int i=0; i<len; i++) {
				if (num>1) {
					stdDev = (n*sum2[i]-sum[i]*sum[i])/n;
					if (stdDev>0.0)
						result[i] = (float)Math.sqrt(stdDev/(n-1.0));
					else
						result[i] = 0f;
				} else
					result[i] = 0f;
			}
		}

    } // end StandardDeviation

}  // end ZProjection


