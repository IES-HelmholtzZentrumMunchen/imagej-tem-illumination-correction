import ij.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.*;

import java.lang.Math;

import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.FFT;


/**
 * Compute spectral coefficients for texture analysis.
 * @author Julien Pontabry
 */
public class Spectral_Texture_Analysis implements PlugInFilter {
	/**
	 * Structure for FFT coordinates (polar coordinates).
	 * @author Julien Pontabry
	 */
	private class FFTCoordinates
	{
		/**
		 * Radius of the polar system.
		 */
		public double r;
		
		/**
		 * Angle of the polar system.
		 */
		public double theta;
		
		/**
		 * Default constructor (coordinates 0,0).
		 */
		public FFTCoordinates()
		{
			this(0.0, 0.0);
		}
		
		/**
		 * Constructor
		 * @param _r Radius of the polar system.
		 * @param _theta Angle of the polar system.
		 */
		public FFTCoordinates(double _r, double _theta)
		{
			r     = _r;
			theta = _theta;
		}
	}
	
	/**
	 * Structure for cartesian coordinates.
	 * @author Julien Pontabry
	 *
	 */
	private class CartesianCoordinates
	{
		/**
		 * X coordinate.
		 */
		public double x = 0.0;
		
		/**
		 * Y coordinate.
		 */
		public double y = 0.0;
		
		/**
		 * Default constructor (coordinates 0,0).
		 */
		public CartesianCoordinates()
		{
			this(0.0, 0.0);
		}
		
		/**
		 * Constructor
		 * @param _x X coordinate.
		 * @param _y Y coordinate.
		 */
		public CartesianCoordinates(double _x, double _y)
		{
			x = _x;
			y = _y;
		}
	}
	
	/**
	 * FHT object of the input image.
	 */
	private FHT m_frequencyImage;
	
	/**
	 * Number of samples for each coefficient to be computed.
	 */
	private int m_numberOfSamples;
	
	/**
	 * Calibration of the image to analyze.
	 */
	private Calibration m_calibration;
	
	/**
	 * Center of the image in the image coordinate system.
	 */
	private double m_center;
	
	/**
	 * Height of the image.
	 */
	private int m_height;
	
	/**
	 * Radius scale.
	 */
	private double[] m_rScale;
	
	/**
	 * Angle scale.
	 */
	private double[] m_thetaScale;
	
	
	/**
	 * Setup method
	 * @param arg Arguments of the plugin.
	 * @param imp Image on which the plugin will be processed.
	 */
	public int setup(String arg, ImagePlus imp)
	{
		// Check if there is an image
		if(imp == null)
		{
			IJ.error("Spectral texture analysis", "There is no image");
			return DONE;
		}
		
		// Set the number of samples to compute for coefficients
		m_numberOfSamples = 1000;
		
		// Get calibration
		m_calibration = imp.getCalibration();
		

		return DOES_8G+DOES_16+DOES_32+NO_CHANGES;
	}
	
	/**
	 * Actually run the plugin processing.
	 * @param ip Image processor on which processing is done.
	 */
	public void run(ImageProcessor ip)
	{
		/////////////////////////////////////////////////////////////
		// Initialization
		//
		
		// Get the mask
		byte[] mask = ip.getMaskArray();
		
		// Crop using the ROI
		ip = ip.crop();
		
		// Remove pixels outside mask
		if(mask != null)
		{
			for(int index = 0; index < ip.getWidth()*ip.getHeight(); index++)
			{
				if(mask[index] == 0)
				{
					ip.set(index, 0);
				}
			}
		}
		
		// Compute FHT transform
		m_frequencyImage = new FHT(this.pad(ip));
		m_frequencyImage.transform();
		ImageProcessor ps = this.computeRawPowerSpectrum(m_frequencyImage);
		
//		Fast_FourierTransform fft = new Fast_FourierTransform();
		
		// Compute useful values
		m_height = m_frequencyImage.getHeight();
		m_center = (double)m_frequencyImage.getWidth() / 2.0;
		
		double rMin = 1.0, rMax = m_center/2, rStep = (rMax - rMin) / (double)m_numberOfSamples;
		double thetaMin = 0.0, thetaMax = Math.PI, thetaStep = (thetaMax - thetaMin) / (double)m_numberOfSamples;
		
		// Compute radius scale
		m_rScale = new double[m_numberOfSamples+1];
		double currentR = rMin;
		
		for(int i = 0; i < m_rScale.length; i++)
		{
			m_rScale[i] = currentR;
			currentR   += rStep;
		}
		
		// Compute angle scale
		m_thetaScale = new double[m_numberOfSamples+1];
		double currentTheta = thetaMin;
		
		for(int i = 0; i < m_thetaScale.length; i++)
		{
			m_thetaScale[i] = currentTheta;
			currentTheta   += thetaStep;
		}
		
		
		/////////////////////////////////////////////////////////////
		// Compute coefficients
		//
		
		// Compute r coefficient
		double[] rCoefficient = new double[m_rScale.length];
		
		for(int i = 0; i < rCoefficient.length; i++)
		{
			rCoefficient[i] = this.computeRCoefficient(m_rScale[i], ps);
		}
		
		// Compute theta coefficient
		double[] thetaCoefficient = new double[m_thetaScale.length];
		
		for(int i = 0; i < thetaCoefficient.length; i++)
		{
			thetaCoefficient[i] = this.computeThetaCoefficient(m_thetaScale[i], ps);
		}
		
		// Convert radians into degrees
		for(int i = 0; i < m_thetaScale.length; i++)
		{
			m_thetaScale[i] = m_thetaScale[i] * 180.0 / Math.PI;
		}
		
		// Convert pixels into image's unit
		if(m_calibration.scaled())
		{
			for(int i = 0; i < m_rScale.length; i++)
			{
				m_rScale[i] = ((double)ip.getWidth()/m_rScale[i]) * m_calibration.pixelWidth;
			}
		}
		
		// Build graphics
		Plot rPlot = new Plot("Radius coefficient", "Radius ("+ m_calibration.getUnit() +"/c)", "", m_rScale, rCoefficient);
		Plot thetaPlot = new Plot("Angle coefficient", "Angle (degrees)", "", m_thetaScale, thetaCoefficient);
		
		// Display graphics
		rPlot.show();
		thetaPlot.show();
	}
	
	/**
	 * Try to padd image if needed (in order to compute the FHT).
	 * @param ip Image processor.
	 * @return The padded image processor if needed.
	 */
	private ImageProcessor pad(ImageProcessor ip)
	{
        int  originalWidth = ip.getWidth();
        int originalHeight = ip.getHeight();
        int           maxN = Math.max(originalWidth, originalHeight);
        
        int i = 2;
        while(i<maxN) i *= 2;
        
        // We do not need to pad the image
        if (i == maxN && originalWidth == originalHeight) {
            return ip;
        }
        
        maxN = i;
        
        // Image need to be padded
        ImageStatistics stats = ImageStatistics.getStatistics(ip, ImageStatistics.MEAN, null);
        ImageProcessor    ip2 = ip.createProcessor(maxN, maxN);
        ip2.setValue(stats.mean);
        ip2.fill();
        ip2.insert(ip, 0, 0);
        Undo.reset();

        return ip2;
	}
	
	/**
	 * Compute the raw power spectrum of the FHT.
	 * @param frequencyImage The FHT.
	 * @return The raw power spectrum.
	 */
	private ImageProcessor computeRawPowerSpectrum(FHT frequencyImage)
	{
		int maxN = frequencyImage.getWidth();
		int base;
		float  r, scale;
		float min = Float.MAX_VALUE;
  		float max = Float.MIN_VALUE;
   		float[] fps = new float[maxN*maxN];
		float[] fht = (float[])frequencyImage.getPixels();

  		for (int row=0; row<maxN; row++) {
			FHTps(row, maxN, fht, fps);
			base = row * maxN;
			for (int col=0; col<maxN; col++) {
				r = fps[base+col];
				if (r<min) min = r;
				if (r>max) max = r;
			}
		}

		max = (float)Math.log(max);
		min = (float)Math.log(min);
		if (Float.isNaN(min) || max-min>50)
			min = max - 50; //display range not more than approx e^50
		scale = (float)(253.999/(max-min));

		for (int row=0; row<maxN; row++) {
			base = row*maxN;
			for (int col=0; col<maxN; col++) {
				r = fps[base+col];
				r = ((float)Math.log(r)-min)*scale;
				if (Float.isNaN(r) || r<0)
					r = 0f;
			}
		}

		ImageProcessor ip2 = new FloatProcessor(maxN, maxN, fps, null);
		swapQuadrants(ip2);
		
		return ip2;
	}
	
	/** Power Spectrum of one row from 2D Hartley Transform. */
 	private void FHTps(int row, int maxN, float[] fht, float[] ps) {
 		int base = row*maxN;
		int l;
		for (int c=0; c<maxN; c++) {
			l = ((maxN-row)%maxN) * maxN + (maxN-c)%maxN;
			ps[base+c] = (sqr(fht[base+c]) + sqr(fht[l]))/2f;
 		}
	}
 	
	/**	Swap quadrants 1 and 3 and 2 and 4 of the specified ImageProcessor 
	so the power spectrum origin is at the center of the image.
	<pre>
	    2 1
	    3 4
	</pre>
	 */
 	private void swapQuadrants(ImageProcessor ip) {
 		//IJ.log("swap");
 		ImageProcessor t1, t2;
 		int size = ip.getWidth()/2;
 		ip.setRoi(size,0,size,size);
 		t1 = ip.crop();
 		ip.setRoi(0,size,size,size);
 		t2 = ip.crop();
 		ip.insert(t1,0,size);
 		ip.insert(t2,size,0);
 		ip.setRoi(0,0,size,size);
 		t1 = ip.crop();
 		ip.setRoi(size,size,size,size);
 		t2 = ip.crop();
 		ip.insert(t1,size,size);
 		ip.insert(t2,0,0);
 		ip.resetRoi();
 	}
 	
	private float sqr(float x) {
		return x*x;
	}
	
	/**
	 * Compute the coefficient for the given radius (r) value.
	 * @param r Radius value.
	 * @param ip Image processor on which the processing is done.
	 * @return The coefficient of the given radius value.
	 */
	private double computeRCoefficient(double r, ImageProcessor ip)
	{
		double coefficient = 0.0;
		
		// Integrate over all theta values
		for(int i = 0; i < m_thetaScale.length; i++)
		{
			// Compute cartesian coordinates from FFT coordinates
			CartesianCoordinates cc = getCartesianCoordinates(r, m_thetaScale[i]);
			
			// Add the current value of FFT
			coefficient += ip.getInterpolatedValue(cc.x, cc.y);
		}
		
		return coefficient;
	}
	
	/**
	 * Compute the coefficient for the given angle (theta) value.
	 * @param theta Angle value.
	 * @param ip Image processor on which the processing is done.
	 * @return The coefficient of the given angle value.
	 */
	private double computeThetaCoefficient(double theta, ImageProcessor ip)
	{
		double coefficient = 0.0;
		
		// Integrate over all r values
		for(int i = 0; i < m_rScale.length; i++)
		{	
			// Compute cartesian coordinates from FFT coordinates
			CartesianCoordinates cc = getCartesianCoordinates(m_rScale[i], theta);
			
			// Add the current value of FFT
			coefficient += ip.getInterpolatedValue(cc.x, cc.y);
		}
		
		return coefficient;
	}
	
	/**
	 * Convert cartesian coordinates into FFT coordinates (polar system).
	 * @param x X coordinates.
	 * @param y Y coordinates.
	 * @return Structure containing polar coordinates.
	 */
	private FFTCoordinates getFFTCoordinates(int x, int y)
	{	
		FFTCoordinates fc = new FFTCoordinates();
		             fc.r = Math.sqrt((x-m_center)*(x-m_center) + (y-m_center)*(y-m_center));
		         fc.theta = Math.atan2(y-m_center, x-m_center);
		
		if(fc.r < 1.0)
		{
			fc.r = 1.0;
		}
		
		fc.theta = fc.theta * 180.0 / Math.PI;
		
		if(fc.theta < 0)
		{
			fc.theta = 360.0 + fc.theta;
		}
		
		return fc;
	}
	
	/**
	 * Convert FFT coordinates (polar system) into cartesian coordinates.
	 * @param r Radius
	 * @param theta Angle
	 * @return Structure containing cartesian coordinates
	 */
	private CartesianCoordinates getCartesianCoordinates(double r, double theta)
	{
		CartesianCoordinates cc = new CartesianCoordinates(
				(int)Math.round( r * Math.cos(theta) + m_center ),
				m_height - (int)Math.round( r * Math.sin(theta) + m_center )
			);
		
		return cc;
	}
	
	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Spectral_Texture_Analysis.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

//		// open the Clown sample
//		ImagePlus image = IJ.openImage("http://imagej.net/images/clown.jpg");
//		image.show();

//		// run the plugin
//		IJ.runPlugIn(clazz.getName(), "");
	}
}
