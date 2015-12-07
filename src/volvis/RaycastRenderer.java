/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     
    short interpolation (double x1, double x2, double alpha) {
        // Perform linear interpolation corresponding to the equation used
        // in the lecture slides
        return (short)Math.round((1 - alpha) * x1 + alpha * x2); 
    }
    
    double biInterpolation (double x00, double x01, double x10, double x11, 
            double alpha, double beta) {
        // Performing two linear interpolations for a 2D figure (interpolation
        // for the top line and bottom line of a 2D square)
        double y0 = interpolation (x00, x01, alpha); 
        double y1 = interpolation (x10, x11, alpha); 
        // Perform linear interpolation between the two earlier calculated
        // values
        return interpolation (y0, y1, beta);
    }
    
    short getVoxel(double[] coord) {
        // computing the bottom left corners to the corresponding 
        // input coordinates
        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        // Determine if any of the input coordinates are outside the image
        // boundaries. 
        // Also determine if computed bottom left corner coordinates plus 2 
        // are outside the image boundaries. Plus 2 is taken to overcome the
        // start at 0 and the plus 1 used in a later stage of this method. 
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 
                || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ() 
                || (z + 2) > volume.getDimZ() || (y + 2) > volume.getDimY() 
                || (x + 2) > volume.getDimX()) {
            return 0;
        }
        
        // Determine relative place of the input coordinates compared to the 
        // bottom left corner coordinates, used to compute alpha, beta and 
        // gamma. 
        double alpha = coord[0] - x; 
        double beta = coord[1] - y; 
        double gamma = coord[2] - z; 
        
        // Setting an initial 3D array of voxel values of each corner of a 
        // 3D cube around the input coordinates 
        int[][][] s = new int[2][2][2]; 
        // Setting an initial array where the values of the bi-linear
        // interpolation are stored. 
        double[] temp = new double[2];
        
        // A 3 part loop for determining the different voxel values of the 
        // eight corners, and also directly determine the value of bi-linear
        // interpolation in the outer loop
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k <2; k++) {
                    // Determine voxel value
                    s[i][j][k] = volume.getVoxel(x + i, y + j, z + k);
                }
            }
            // Bi-linear interpolation for the plane corresponding to the i 
            // value for the x coordinate
            temp[i] = biInterpolation (s[i][0][0], s[i][0][1], s[i][1][0], 
                    s[i][1][1], alpha, beta);
        }

        // Linear interpolation between the two bi-linear interpolation values
        return interpolation(temp[0], temp[1], gamma);
        //return volume.getVoxel(x, y, z);
    }

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < (image.getHeight() - 1); j++) {
            for (int i = 0; i < (image.getWidth() - 1); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];
                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }
       
    void mip(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
 
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < (image.getHeight() - 1); j++) {
            for (int i = 0; i < (image.getWidth() - 1); i++) {
                // Using the original pixel coordinates set in the skeleton code,
                // But now the starting point is not the middle of the image, but 
                // is set just behind the image using: 
                // imageCenter * viewVec[x];
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0]  - imageCenter * viewVec[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] - imageCenter * viewVec[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] - imageCenter * viewVec[2];

                // Creating a boolean which is used for the while-loop to see
                // when the viewvector leaves the image boundaries
                boolean cont = true;
                int maxVal = getVoxel(pixelCoord);

                // A while-loop used to "walk" across the view vector
                while(cont) {
                    // Each iteration, a step of X is taken across the view
                    // vector and the corresponding pixel coordinates of the 
                    // new voxel is determined. 
                    pixelCoord[0] = pixelCoord[0] + 5 * viewVec[0];
                    pixelCoord[1] = pixelCoord[1] + 5 * viewVec[1];
                    pixelCoord[2] = pixelCoord[2] + 5 * viewVec[2];
                    
                    // If the view vector leaves the image boundaries, the while
                    // loop is terminated. The boundaries are set to correspond
                    // ......
                    if (pixelCoord[0] > (volumeCenter[0] + Math.sqrt(2) * imageCenter) 
                            || pixelCoord[0] < (volumeCenter[0] - Math.sqrt(2) * imageCenter)
                        || pixelCoord[1] > (volumeCenter[1] + Math.sqrt(2) * imageCenter) 
                            || pixelCoord[1] < (volumeCenter[1] - Math.sqrt(2) * imageCenter) 
                        || pixelCoord[2] > (volumeCenter[2] + Math.sqrt(2) * imageCenter) 
                            || pixelCoord[2] < (volumeCenter[2] - Math.sqrt(2) * imageCenter)
                            ) {
                        cont = false; 
                    }
                    
                    // Get the Voxel value for the new pixel coordinates
                    int val = getVoxel(pixelCoord);
                    // Determine if the voxel value of the new coordinates
                    // is larger than the already largest value. If so, the 
                    // maxValue is changed to this value, otherwise it remains
                    // untouched
                    maxVal = val > maxVal ? val : maxVal; 
                }
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxVal/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void compositing(double[] viewMatrix) {
        
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;
        
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
 
        // sample on a plane through the origin of the volume data
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < (image.getHeight() - 1); j++) {
            for (int i = 0; i < (image.getWidth() - 1); i++) {
                // Using the original pixel coordinates set in the skeleton code,
                // But now the starting point is not the middle of the image, but 
                // is set just behind the image using: 
                // imageCenter * viewVec[x];
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0]  - imageCenter * viewVec[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] - imageCenter * viewVec[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] - imageCenter * viewVec[2];
                
                int val = getVoxel(pixelCoord);
                voxelColor = tFunc.getColor(val);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                double c_a = voxelColor.a <= 1.0 ? voxelColor.a : 1.0;
                double c_r = voxelColor.r <= 1.0 ? voxelColor.r : 1.0;
                double c_g = voxelColor.g <= 1.0 ? voxelColor.g : 1.0;
                double c_b = voxelColor.b <= 1.0 ? voxelColor.b : 1.0;
                
                // Creating a boolean which is used for the while-loop to see
                // when the viewvector leaves the image boundaries
//                boolean cont = true;
                
//                // A while-loop used to "walk" across the view vector
//                while(cont) {
//                    // Each iteration, a step of X is taken across the view
//                    // vector and the corresponding pixel coordinates of the 
//                    // new voxel is determined. 
//                    pixelCoord[0] = pixelCoord[0] + 10 * viewVec[0];
//                    pixelCoord[1] = pixelCoord[1] + 10 * viewVec[1];
//                    pixelCoord[2] = pixelCoord[2] + 10 * viewVec[2];
//                    
//                    // If the view vector leaves the image boundaries, the while
//                    // loop is terminated. The boundaries are set to correspond
//                    // ......
//                    if (pixelCoord[0] > (volumeCenter[0] + Math.sqrt(2) * imageCenter) 
//                            || pixelCoord[0] < (volumeCenter[0] - Math.sqrt(2) * imageCenter)
//                        || pixelCoord[1] > (volumeCenter[1] + Math.sqrt(2) * imageCenter) 
//                            || pixelCoord[1] < (volumeCenter[1] - Math.sqrt(2) * imageCenter) 
//                        || pixelCoord[2] > (volumeCenter[2] + Math.sqrt(2) * imageCenter) 
//                            || pixelCoord[2] < (volumeCenter[2] - Math.sqrt(2) * imageCenter)
//                            ) {
//                        cont = false; 
//                    } else {
//                    
//                        // Get the voxelColor of the new pixel coordinates       
//                        val = getVoxel(pixelCoord); 
//                        voxelColor = tFunc.getColor(val);
//
//                        // Determine the intensity corresponding to the new pixel
//                        // coordinates
//                        double tau = voxelColor.a <= 1.0 ? voxelColor.a : 1.0;
//
//                        // Calculating the individual ARGB values of the new pixel
//                        // coordinates
//                        c_a = c_a * (1 - tau);
//
//                        c_r = voxelColor.r * tau + (1 - tau) * c_r; 
//                        c_g = voxelColor.g * tau + (1 - tau) * c_g; 
//                        c_b = voxelColor.b * tau + (1 - tau) * c_b; 
//                    }
//                }
                int scaling = 1; 
                for (double k = 0.0; k < 4 * imageCenter; k = k + scaling) {
                    pixelCoord[0] = pixelCoord[0] + scaling * viewVec[0];
                    pixelCoord[1] = pixelCoord[1] + scaling * viewVec[1];
                    pixelCoord[2] = pixelCoord[2] + scaling * viewVec[2];
                    
                    if (pixelCoord[0] < 0 || pixelCoord[0] > (volume.getDimX() - 1) 
                            || pixelCoord[1] < 0 || pixelCoord[1] > (volume.getDimY() - 1)
                            || pixelCoord[2] < 0 || pixelCoord[2] > (volume.getDimZ() - 1)) {
                        
                    } else {
                        val = getVoxel(pixelCoord); 
                        voxelColor = tFunc.getColor(val);

                        // Determine the intensity corresponding to the new pixel
                        // coordinates
                        double tau = voxelColor.a <= 1.0 ? voxelColor.a : 1.0;

                        // Calculating the individual ARGB values of the new pixel
                        // coordinates
                        c_a = c_a * (1 - tau);

                        c_r = voxelColor.r * tau + (1 - tau) * c_r; 
                        c_g = voxelColor.g * tau + (1 - tau) * c_g; 
                        c_b = voxelColor.b * tau + (1 - tau) * c_b; 
                    }
                }
                               
                int c_alpha = (1 - c_a) <= 1.0 ? (int) Math.floor((1 - c_a) * 255) : 255;
                int c_red = c_r <= 1.0 ? (int) Math.floor(c_r * 255) : 255;
                int c_green = c_g <= 1.0 ? (int) Math.floor(c_g * 255) : 255;
                int c_blue = c_b <= 1.0 ? (int) Math.floor(c_b * 255) : 255;
                // Determine the pixel colorand set the image to that color
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }
    
    void TF2D(double[] viewMatrix) {
        
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;
        
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
 
        // sample on a plane through the origin of the volume data
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        double opacity = tfEditor2D.triangleWidget.color.a;
        double radius = tfEditor2D.triangleWidget.radius;
        double bIntensity = tfEditor2D.triangleWidget.baseIntensity;
           
        double colorR = tfEditor2D.triangleWidget.color.r;
        double colorG = tfEditor2D.triangleWidget.color.g;
        double colorB = tfEditor2D.triangleWidget.color.b;        

        for (int j = 0; j < (image.getHeight() - 1); j++) {
            for (int i = 0; i < (image.getWidth() - 1); i++) {
                // Using the original pixel coordinates set in the skeleton code,
                // But now the starting point is not the middle of the image, but 
                // is set just in front of the image using: 
                // imageCenter * viewVec[x];
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] - imageCenter * viewVec[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] - imageCenter * viewVec[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] - imageCenter * viewVec[2];
                
 
//                 Creating a boolean which is used for the while-loop to see
//                 when the viewvector leaves the image boundaries
//                boolean cont = true;
                double c_a = 1;
                double scaling = 1;
                
                for (double k = 0.0; k < 2 * imageCenter; k = k + scaling) {
                    
                    if (pixelCoord[0] < 0 || pixelCoord[0] > (volume.getDimX() - 1) 
                            || pixelCoord[1] < 0 || pixelCoord[1] > (volume.getDimY() - 1)
                            || pixelCoord[2] < 0 || pixelCoord[2] > (volume.getDimZ() - 1)) {
                        
                    } 
                    else {
                        int val = getVoxel(pixelCoord);
                        
                        double alpha; 
                    
                        int pixelCoordX = (int) Math.floor(pixelCoord[0]);
                        int pixelCoordY = (int) Math.floor(pixelCoord[1]);
                        int pixelCoordZ = (int) Math.floor(pixelCoord[2]);
                        
//                        float mag = -1;
                        float mag = gradients.getGradient(pixelCoordX, pixelCoordY, pixelCoordZ).mag;
                        
                        if (mag == 0 && val == bIntensity) {
                            alpha = opacity * 1; 
                        } 
                        else if (mag > 0 && (val - radius * mag) <= bIntensity && bIntensity <= (val + radius * mag)) {
                            alpha = opacity * (1 - ((1 / radius) * Math.abs((bIntensity - val) / mag)));
                        } 
                        else {
                            alpha = opacity * 0;
                        }
                        
                        c_a = c_a * (1 - alpha);
                    }
                    
                    pixelCoord[0] = pixelCoord[0] + scaling * viewVec[0];
                    pixelCoord[1] = pixelCoord[1] + scaling * viewVec[1];
                    pixelCoord[2] = pixelCoord[2] + scaling * viewVec[2];
                }
                               
//                // A while-loop used to "walk" across the view vector
//                while(cont) {
//
//                    // Get the voxelColor of the new pixel coordinates       
//                    int val = getVoxel(pixelCoord);
//                    
//                    double alpha; 
//                    
//                    int pixelCoordX = (int) Math.floor(pixelCoord[0]);
//                    int pixelCoordY = (int) Math.floor(pixelCoord[1]);
//                    int pixelCoordZ = (int) Math.floor(pixelCoord[2]);
//                    float mag = -1;
//                    if (pixelCoordX < 0 || pixelCoordX + 1 > volume.getDimX() 
//                            || pixelCoordY < 0 || pixelCoordY + 1 > volume.getDimY()
//                            || pixelCoordZ < 0 || pixelCoordZ + 1 > volume.getDimZ()) {
//                    }
//                    else{
//                        mag = gradients.getGradient(pixelCoordX, pixelCoordY, pixelCoordZ).mag;
//                    }
//
//                    if (mag == 0 && val == bIntensity) {
//                        alpha = opacity * 1.0; 
//                        
//                    } else if (mag > 0 && (val - radius * mag) <= bIntensity && bIntensity <= (val + radius * mag)) {
//                        alpha = opacity * (1 - (1 / radius) * Math.abs((bIntensity - val) / mag));
//                    } else {
//                        alpha = opacity * 0.0;
//                    }
//
//                    c_a = c_a * (1 - alpha);
//                    
//                    // Each iteration, a step of X is taken across the view
//                    // vector and the corresponding pixel coordinates of the 
//                    // new voxel is determined. 
//                    pixelCoord[0] = pixelCoord[0] + 8 * viewVec[0];
//                    pixelCoord[1] = pixelCoord[1] + 8 * viewVec[1];
//                    pixelCoord[2] = pixelCoord[2] + 8 * viewVec[2];
//                    
//                    // If the view vector leaves the image boundaries, the while
//                    // loop is terminated. The boundaries are set to correspond
//                    // ......
//                    if (pixelCoord[0] > (volumeCenter[0] + Math.sqrt(2) * imageCenter) 
//                            || pixelCoord[0] < (volumeCenter[0] - Math.sqrt(2) * imageCenter)
//                        || pixelCoord[1] > (volumeCenter[1] + Math.sqrt(2) * imageCenter) 
//                            || pixelCoord[1] < (volumeCenter[1] - Math.sqrt(2) * imageCenter) 
//                        || pixelCoord[2] > (volumeCenter[2] + Math.sqrt(2) * imageCenter) 
//                            || pixelCoord[2] < (volumeCenter[2] - Math.sqrt(2) * imageCenter)
//                            ) {
//                        cont = false; 
//                    }
//                }
                               
                int c_alpha = (1 - c_a) <= 1.0 ? (int) Math.floor((1 - c_a) * 255) : 255;
                
                int c_red = colorR <= 1.0 ? (int) Math.floor(colorR * 255) : 255;
                int c_green = colorG <= 1.0 ? (int) Math.floor(colorG * 255) : 255;
                int c_blue = colorB <= 1.0 ? (int) Math.floor(colorB * 255) : 255;
                
                // Determine the pixel colorand set the image to that color
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }
    
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        switch(renderOption){
            case 1: slicer(viewMatrix);
                break;
            case 2: mip(viewMatrix);
                break; 
            case 3: compositing(viewMatrix); 
                break;
            case 4: TF2D(viewMatrix);
                break;
            default: slicer(viewMatrix);
;               break; 
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
    public int renderOption;

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
