import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.plugin.filter.GaussianBlur;

import MyJama.*;

/**
 * Correct illumination bias of TEM images.
 *
 * This plug-in is an implementation of Tasdizen's method (see Tasdizen et al. (2008).
 * Non-uniform illumination correction in transmission electron microscopy. In: MICCAI
 * Workshop on Microscopic Image Analysis with Applications in Biology, 5-6).
 *
 * @author Julien Pontabry
 */
public class TEM_Illumination_Bias_Correction implements PlugInFilter {
    private boolean m_displayBias;
    private float m_mu2 = 0.001f;

	/**
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
        GenericDialog gui = new GenericDialog("TEM illumination bias correction");
        gui.addCheckbox("Display illumination bias", false);
        gui.addNumericField("Weighting parameter", m_mu2, 3);
        gui.showDialog();

        if (gui.wasCanceled()) {
            return DONE;
        }

        m_displayBias = gui.getNextBoolean();
        m_mu2 = (float)gui.getNextNumber();

		return DOES_8G+DOES_16+DOES_32;
	}

	/**
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
		/////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//
		
		IJ.showStatus("Initialization...");
		
		// Some modeling constants
		int     polynomeDegree = 2;
		int numberOfParameters = (polynomeDegree+1) * (polynomeDegree+2) / 2 - 1;
		
		// Image information
		int          width = ip.getWidth();
		int         height = ip.getHeight();
		int numberOfPixels = width * height;
		
		// Blur image (in order to suppress noise - we assume this is in very high frequencies)
		ImagePlus  originalImage = new ImagePlus("Original",ip);
		ImagePlus   currentImage = originalImage.duplicate();
		ImageProcessor currentIp = currentImage.getProcessor();
		
		GaussianBlur blur = new GaussianBlur();
		blur.blurGaussian(currentIp, 10, 10, 0.02);

        // Initialize progress bar
        int currentProgress = 0;
        int     maxProgress = numberOfPixels * 7;
        IJ.showProgress(currentProgress++, maxProgress);
		
		
		/////////////////////////////////////////////////////////////////////////////////
		// Build matrices
		//
		
		// Monomial matrix data
		float[][] monomialMatrix_data = new float[2*numberOfPixels][numberOfParameters];
		
		// Image gradient matrix data
		float[][] gradientMatrix_data = new float[2*numberOfPixels][1];

        // Weights of pixels
        float[] weights_data = new float[2*numberOfPixels];

        IJ.showStatus("Computing data matrices...");
		
		// For each pixel, compute the monomials coefficients and the image gradient
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				int currentParameter = 0;
				
				// Compute monomial matrix data
				for(int i = 1; i <= polynomeDegree; i++)
				{
					for(int j = 0; j <= i; j++)
					{
						// Compute derivative of polynomial along direction x 
						if(i-j > 0)
						{
							monomialMatrix_data[x*height+y][currentParameter] = (i-j) * (float)Math.pow(x,i-j-1) * (float)Math.pow(y,j);
						}
						else
						{
							monomialMatrix_data[x*height+y][currentParameter] = 0.0f;
						}
						
						// Compute derivative of polynomial along direction y
						if(j > 0)
						{
							monomialMatrix_data[x*height+y + numberOfPixels][currentParameter] = j * (float)Math.pow(x,i-j) * (float)Math.pow(y,j-1);
						}
						else
						{
							monomialMatrix_data[x*height+y + numberOfPixels][currentParameter] = 0.0f;
						}
						
						currentParameter++;
					} // for j
				} // for i
				
				// Compute gradient of the image log of direction x
				if(x == 0)
				{
					gradientMatrix_data[x*height+y][0] = (float)Math.log(currentIp.getf(2,y)) - (float)Math.log(currentIp.getf(1,y));
				}
				else if(x == width-1)
				{
					gradientMatrix_data[x*height+y][0] = (float)Math.log(currentIp.getf(width-1,y)) - (float)Math.log(currentIp.getf(width-2,y));
				}
				else
				{
					gradientMatrix_data[x*height+y][0] = 0.5f * ( (float)Math.log(currentIp.getf(x+1,y)) - (float)Math.log(currentIp.getf(x-1,y)) );
				}
				
				// Compute gradient of the image log of direction y
				if(y == 0)
				{
					gradientMatrix_data[x*height+y + numberOfPixels][0] = (float)Math.log(currentIp.getf(x,2)) - (float)Math.log(currentIp.getf(x,1));
				}
				else if(y == height-1)
				{
					gradientMatrix_data[x*height+y + numberOfPixels][0] = (float)Math.log(currentIp.getf(x,height-1)) - (float)Math.log(currentIp.getf(x,height-2));
				}
				else
				{
					gradientMatrix_data[x*height+y + numberOfPixels][0] = 0.5f * ( (float)Math.log(currentIp.getf(x,y+1)) - (float)Math.log(currentIp.getf(x,y-1)) );
				}

                // Compute weights
                weights_data[x*height+y]                  = (float) Math.exp(-Math.sqrt(gradientMatrix_data[x + height + y][0] * gradientMatrix_data[x + height + y][0] + gradientMatrix_data[x * height + y + numberOfPixels][0] * gradientMatrix_data[x * height + y + numberOfPixels][0]) / m_mu2);
                weights_data[x*height+y + numberOfPixels] = weights_data[x*height+y];
				
				// Update progress bar
				IJ.showProgress(currentProgress++, maxProgress);
			} // for y
		} // for x


        /////////////////////////////////////////////////////////////////////////////////
        // Estimate modeling parameters by linear fitting
        //

        IJ.showStatus("Estimating illumination modeling...");

		// Setup monomial and image gradient matrices
		Matrix monomialMatrix = new Matrix(monomialMatrix_data);
		Matrix gradientMatrix = new Matrix(gradientMatrix_data);

        monomialMatrix_data = null;
        gradientMatrix_data = null;

        // Pre-compute weighting
        Matrix weightedMonomialMatrix = new Matrix(monomialMatrix.getRowDimension(), monomialMatrix.getColumnDimension());

        for (int i = 0; i < weightedMonomialMatrix.getRowDimension(); i++) {
            for (int j = 0; j < weightedMonomialMatrix.getColumnDimension(); j++) {
                weightedMonomialMatrix.set(i,j, monomialMatrix.get(i,j) * weights_data[i]);
            }

            IJ.showProgress(currentProgress++, maxProgress);
        }

        Matrix weightedGradientMatrix = new Matrix(gradientMatrix.getRowDimension(), gradientMatrix.getColumnDimension());

        for (int i = 0; i < weightedGradientMatrix.getRowDimension(); i++) {
            weightedGradientMatrix.set(i,0, gradientMatrix.get(i,0) * weights_data[i]);

            IJ.showProgress(currentProgress++, maxProgress);
        }

        gradientMatrix = null;
        weights_data   = null;
	
		// Least-square fitting (parameters = (M^T . W . M)^-1 . M^T . W . g
//        Matrix parameters = monomialMatrix.transpose().times(weightedMonomialMatrix).inverse().times(monomialMatrix.transpose()).times(weightedGradientMatrix);
        monomialMatrix = monomialMatrix.transpose();
        weightedMonomialMatrix = monomialMatrix.times(weightedMonomialMatrix);
        weightedGradientMatrix = monomialMatrix.times(weightedGradientMatrix);
        Matrix parameters = weightedMonomialMatrix.inverse().times(weightedGradientMatrix);

        currentProgress = currentProgress + numberOfPixels;
        IJ.showProgress(currentProgress, maxProgress);

        monomialMatrix = null;
		
		
		/////////////////////////////////////////////////////////////////////////////////
		// Compute the illumination image given the modeling
		//
		
		IJ.showStatus("Fixing illumination bias...");

        ImagePlus               bias = NewImage.createFloatImage("Illumination bias of "+ originalImage.getTitle(), width, height, 1, NewImage.FILL_BLACK);
        ImageProcessor biasProcessor = bias.getProcessor();
		
		// For each pixel, compute the illumination modeling
		for(int x = 0; x < width; x++)
		{
			for(int y = 0; y < height; y++)
			{
				double sum = 0;
				int currentParameter = 0;
				
				// Compute the polynomial value
				for(int i = 1; i <= polynomeDegree; i++)
				{
					for(int j = 0; j <= i; j++)
					{
						sum += parameters.get(currentParameter,0) * Math.pow(x,i-j) * Math.pow(y,j);
						currentParameter++;
					}
				}
				
				// Set the illumination model
                biasProcessor.setf(x,y, (float)Math.exp(sum));
				ip.setf(x,y, ip.getf(x,y) / biasProcessor.getf(x,y) );
				
				// Update progress bar
				IJ.showProgress(currentProgress++, maxProgress);
			} // for y
		} // for x

        IJ.showProgress(maxProgress, maxProgress);

        if (m_displayBias) {
            bias.show();
        }
        else { // !m_displayBias
            bias.close();
        }
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
		Class<?> clazz = TEM_Illumination_Bias_Correction.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();
	}
}
