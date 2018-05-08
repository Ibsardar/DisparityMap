////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       llfd.cpp
//  Date:           03/31/2018
//  Description:    Main implementation for LLFD (low level feature detection) class.
//
////////////////////////////////////////////////////////////////////////////////
//
//  Honor Pledge:
//
//  I pledge that I have neither given nor received any help on this project.
//
//  ibsardar
//
////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2018 Copyright Holder All Rights Reserved.
//
////////////////////////////////////////////////////////////////////////////////

#include "pnmio.h"
#include "llfd.h"
#include "llip.h"
#include "matrix.h"

#include "math.h"
#include <cmath>  //need this for correct abs function
#include <iostream>
#include <stdexcept>
#include <algorithm>

// for quicker PI calculations
#define PI 3.14159265

namespace llfd {

    //
    //  get_mat_window
    //
    Matrix get_mat_window(Image& img, int i, int w, int h, double replacer) {

        // helper variables
        const int irow = i / img.width; // the row or y coordinate in the image of the i'th pixel
        const int icol = i % img.width; // the column or x coordinate in the image of the i'th pixel
        const int halfw = w / 2;
        const int halfh = h / 2;
        const int imgw = img.width;
        const int imgh = img.height;

        // prepare
        w += (w % 2 ? 0 : 1); //if even, add 1
        h += (h % 2 ? 0 : 1); //if even, add 1
        Matrix window(h,w); //window

        // gather values for window
        for (int r=0, y=irow-halfh; y<=irow+halfh; r++, y++) {
            for (int c=0, x=icol-halfw; x<=icol+halfw; c++, x++) {
                // does this pixel exist?
                if (x>=0 && x<imgw && y>=0 && y<imgh) {
                    window(r,c) = img.data[x+y*imgw]; //YES!
                } else {
                    window(r,c) = replacer; //NO!
                }
            }
        }

        // return a matrix of intensities centered at the i'th pixel
        return window;
    }

    //
    //  get_mat_window (for matrix input)
    //
    Matrix get_mat_window(Matrix& m, int w, int h, int cc, int cr,
            int ww, int wh, double replacer) {

        // helper variables
        const int halfw = ww / 2;
        const int halfh = wh / 2;

        // prepare
        ww += (ww % 2 ? 0 : 1); //if even, add 1
        wh += (wh % 2 ? 0 : 1); //if even, add 1
        Matrix window(wh,ww); //window

        // gather values for window
        for (int r=0, j=cr-halfh; r<wh; r++, j++) {
            for (int c=0, i=cc-halfw; c<ww; c++, i++) {
                // does this pixel exist?
                if (i>=0 && i<w && j>=0 && j<h) {
                    window(r,c) = m(j,i); //YES!
                } else {
                    window(r,c) = replacer; //NO!
                }
            }
        }

        // return a matrix of intensities
        return window;
    }

    //
    //  trace
    //
    double trace(Matrix m, unsigned int n) {
        double t = 0;
        for (int i=0; i<n; i++)
            t += m(i,i);
        return t;
    }

    //
    //  det2d
    //
    double det2d(Matrix m) {
        return m(0,0) * m(1,1) - m(1,0) * m(0,1);
    }

    //
    //  normalize_gradients
    //
    void normalize_gradients(Matrix& a, Matrix& b, unsigned int w, unsigned int h) {
        // find min, max
        double min_val=999999999, max_val=-999999999;
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                if (min_val > a(r,c))
                    min_val = a(r,c);
                if (max_val < a(r,c))
                    max_val = a(r,c);
                if (min_val > b(r,c))
                    min_val = b(r,c);
                if (max_val < b(r,c))
                    max_val = b(r,c);
            }
        }

        // normalize
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                a(r,c) = 255 * (a(r,c) - min_val) / (max_val - min_val);
                b(r,c) = 255 * (b(r,c) - min_val) / (max_val - min_val);
            }
        }
    }

    //
    //  nearest_max
    //
    std::vector<int> nearest_max(int c, int r, Matrix& m, unsigned int w, unsigned int h) {
        std::vector<int> next(2);
        next[0] = c;
        next[1] = r;
        double maxval = m(r,c);
        for (int j=r-1; j<=r+1; j++) {
            for (int i=c-1; i<=c+1; i++) {
                if (i>=0 && i<w && j>=0 && j<h) {
                    if (maxval < m(j,i)) {
                        maxval = m(j,i);
                        next[0] = i;
                        next[1] = j;
                    }
                }
            }
        }
        if (next[0]==c && next[1]==r)
            return next;
        else
            return nearest_max(next[0],next[1],m,w,h);
    }

    bool is_angle_between(double radians, double start, double finish) {
        double fuel = 2*PI;
        while(radians < start) radians += 2*PI;
        fuel -= fmod(radians-start,2*PI);
        while(finish < radians) finish += 2*PI;
        fuel -= fmod(finish-radians,2*PI);
        return fuel >= 0;
    }

    //
    //  round2
    //
    double round2(double a) {
        return round(a*100)/100;
    }

    //
    //  min_eigen_vec_2d
    //
    std::vector<double> min_eigen_vec_2d(Matrix& m) {
        // use Oliver knill's method (from Harvard)
        // source: http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html

        // Eigenvalues:
        // e1 = T/2 + sqrt(T*T/4-D)
        // e2 = T/2 - sqrt(T*T/4-D)
        double T = llfd::trace(m,2);
        double D = llfd::det2d(m);
        double e1 = T/2 + sqrt(T*T/4-D);
        double e2 = T/2 - sqrt(T*T/4-D);

        // Eigenvectors:
        std::vector<double> min_vec(2);
        if (e1 >= e2) {
            if (m(1,0)!=0) {
                min_vec[0]=e1-m(1,1);
                min_vec[1]=m(1,0);
            } else if (m(0,1)!=0) {
                min_vec[0]=m(0,1);
                min_vec[1]=e1-m(0,0);
            } else if (m(1,0)==0 && m(0,1)==0) {
                min_vec[0]=1;
                min_vec[1]=0;
            }
        } else {
            if (m(1,0)!=0) {
                min_vec[0]=e2-m(1,1);
                min_vec[1]=m(1,0);
            } else if (m(0,1)!=0) {
                min_vec[0]=m(0,1);
                min_vec[1]=e2-m(0,0);
            } else if (m(1,0)==0 && m(0,1)==0) {
                min_vec[0]=0;
                min_vec[1]=1;
            }
        }

        return min_vec;
    }

    //
    //  detect_corners
    //
    std::vector<std::vector<int> > detect_corners(Image& img, unsigned int N, double sigma,
            double alpha, unsigned int threshold, bool colored) {
        //  -----------------------------------------------------------------------
        //  Colors a red dot on all detected corners on the input image:
        //  1) Initialize 5 matrices, MIx, MIy, MIxx, MIyy, and MIxy to 0s (same size as input image)
        //  2) Calculate spatial derivative masks given an NxN window
        //     a) create 1st-Order-Partial Gaussian Masks for x and y directions
        //        - mask Mx will be a NxN matrix, mask My will be Mx transposed
        //        - Mx and My will be normalized between 0 and 1
        //        - x partial formula = -1*(x/(2*PI*sigma^4))*exp(-1*(x^2+y^2)/(2*sigma^2))
        //        - y partial formula = -1*(y/(2*PI*sigma^4))*exp(-1*(x^2+y^2)/(2*sigma^2))
        //        - Ix = Mx convolved with a pixel
        //        - Iy = My convolved with a pixel
        //        - MIx = matrix, same size of the input image, containing x-dir gradient values
        //        - MIy = matrix, same size of the input image, containing y-dir gradient values
        //        - MIxz = matrix, same size of the input image, containing Ix * Ix calculations
        //        - MIyy = matrix, same size of the input image, containing Iy * Iy calculations
        //        - MIxy = matrix, same size of the input image, containing Ix * Iy calculations
        //        - When convolving, if the mask reaches out of the image, pretend a zero exists
        //        - Gradient direction, g_theta = arctan(Iy / Ix)
        //        - Gradient magnitude, g_mag = sqrt(Ix*Ix + Iy*Iy)
        //     b) Calculate MIx and MIy by convolving masks Mx and My with each pixel
        //  3) Initialize a Score matrix with 0s, same size as the input image
        //  4) Create a NxN Gaussian mask matrix and normalize it from 0 to 1
        //  5) For each pixel in the input image:
        //     a) Initialize a 2x2 M matrix with 0s along with Ixx, Iyy, and Ixy doubles (also set to 0)
        //     b) For each pixel in an NxN window (centered at current pixel) of both gradient matrices:
        //        - Add a gaussed pixel intensity to Ixx from MIxx
        //        - Add a gaussed pixel intensity to Iyy from MIyy
        //        - Add a gaussed pixel intensity to Ixy from MIxy
        //     c) Construct the M (structure tensor) matrix:
        //                                  [ Ixx  Ixy ]
        //                                  [ Ixy  Iyy ]
        //     d) Calculate the corner-ness score, R:
        //        - R = determinant(M) - alpha*trace(M)*trace(M)
        //        - where alpha is a constant between 0.04 and 0.06,
        //        - determinant(M) = Ixx*Iyy - Ixy*Ixy,
        //        - and trace(M) = Ixx + Iyy
        //     e) set the corresponding point in the Score matrix equal to R
        //  6) Create an empty corners list
        //  7) For each score in the Score matrix greater than or equal to threshold:
        //     a) push the following score's local maximum to a list (don't repeat values)
        //  8) Convert the input image into RGB format
        //  9) Loop through the corners list and paint a red dot on each corner
        //  10) Return list of corners (vector of int[2] arrays)
        //  -----------------------------------------------------------------------

        // helper variables
        int n = (N % 2 ? N : N+1);
        int rad = n / 2; //radius of window
        int pixels = img.height * img.width;

        // 1.
        Matrix MIx(img.height, img.width);
        Matrix MIy(img.height, img.width);
        Matrix MIxx(img.height, img.width);
        Matrix MIyy(img.height, img.width);
        Matrix MIxy(img.height, img.width);

        // 2.
        Matrix Mx(n, n);

        // a.
        for (int r=-rad; r<n-rad; r++) {
            for (int c=-rad; c<n-rad; c++) {
                // x = c, y = 0
                double s2 = sigma*sigma;
                Mx(r+rad,c+rad) = -1*(c/(2*PI*s2*s2))*exp(-1*(c*c+r*r)/(2*s2)); // 1st x-partial of Gaussian
            }
        }
        Mx = Mx * (-1 * Matrix::createIdentity(n)); // horizontal flip
        Matrix My = Mx.transpose();

        // b.
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                int i = c + (r * img.width);
                Matrix window = llip::get_window(img, rad, i, true);
                // convolve
                double Ix = 0.0, Iy = 0.0;
                for (int wr=0; wr<n; wr++) {
                    for (int wc=0; wc<n; wc++) {
                        Ix += Mx(wr,wc) * window(wr,wc);
                        Iy += My(wr,wc) * window(wr,wc);
                    }
                }
                MIx(r,c) = Ix;
                MIy(r,c) = Iy;
                MIxx(r,c) = Ix * Ix;
                MIyy(r,c) = Iy * Iy;
                MIxy(r,c) = Ix * Iy;
            }
        }

        // this commented-out snippet alters the image data to be its x-gradient
        /*llfd::normalize_gradients(MIy,MIx,img.width,img.height);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = MIx(r,c);
            }
        }
        return std::vector< std::vector<int> >();*/

        // 3.
        Matrix Scores(img.height, img.width);

        // 4.
        Matrix G(n,n);
        for (int gr=-rad; gr<=rad; gr++) {
            for (int gc=-rad; gc<=rad; gc++) {
                G(gr+rad,gc+rad) = (1/(2*PI*sigma))*exp((gc*gc+gr*gr)/(-2*sigma));
            }
        }
        llip::normalize_matrix2(G,n,n);

        // 5.
        for (int i=0, ir=0, ic=0; i<pixels; i++) {
            ir = i / img.width;
            ic = i % img.width;

            // a.
            Matrix M(2,2);
            double Ixx = 0; // these are actually smoothed versions of Ixx, Iyy, and Ixy
            double Iyy = 0;
            double Ixy = 0;

            // b.
            for (int r=ir-rad, gr=0; r<=ir+rad; r++, gr++) {
                for (int c=ic-rad, gc=0; c<=ic+rad; c++, gc++) {
                    if (c>=0 && c<img.width && r>=0 && r<img.height) {
                        Ixx += MIxx(r,c) * G(gr,gc);
                        Iyy += MIyy(r,c) * G(gr,gc);
                        Ixy += MIxy(r,c) * G(gr,gc);
                    }
                }
            }

            // c.
            M(0,0) = Ixx;
            M(1,0) = Ixy;
            M(0,1) = Ixy;
            M(1,1) = Iyy;

            // d.
            double R = llfd::det2d(M) - alpha * llfd::trace(M,2) * llfd::trace(M,2);

            // e.
            Scores(ir,ic) = R;
        }

        // 6.
        std::vector<std::vector<int> > corners; // x,y coordinates of corners

        // 7.
        int rad2 = rad*2+5;
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                if (Scores(r,c) > threshold) {
                    // get local max score in NxN square
                    std::vector<int> max_corner = llfd::nearest_max(c,r,Scores,img.width,img.height);
                    // add corner if it doesn't already exist
                    bool exists = false;
                    for (int i=0;i<corners.size();i++) {
                        if (corners[i][0]==max_corner[0] && corners[i][1]==max_corner[1]) {
                            exists = true;
                            break;
                        }
                    }
                    if (!exists)
                        corners.push_back(max_corner);
                }
            }
        }

        // check colored?
        if (colored) {
            // colored red on image
            // 8.
            pnmio converter;
            converter.convert_gray_to_rgb(img);

            // 9.
            for (int i=0; i<corners.size(); i++) {
                int index = 3 * (corners[i][0] + corners[i][1] * img.width);
                img.data[index] = 255; // R
                img.data[index+1] = 0; // G
                img.data[index+2] = 0; // B
            }
        } else {
            // binary black on white
            fill(img.data.begin(), img.data.end(), 0); // Black
            for (int i=0; i<corners.size(); i++) {
                img.data[corners[i][0] + corners[i][1] * img.width] = 255; // White
            }
        }

        // 10.
        return corners;
    }

    //
    //  detect_edges
    //
    void detect_edges(Image& img, int N, double sigma, double alpha, double low_th, double high_th) {
        //  -----------------------------------------------------------------------
        //  1. Blur image using Gaussian blur
        //  2. Make Gaussian masks for obtaining the gradients
        //  3. Gather x and y gradients (in matrices same size as image)
        //  4. Gather gradient magnitudes and directions (in matrices same size as image)
        //  5. Apply non-max suppression on each pixel:
        //     a) Obtain front and back pixel values from 8-neighborhood of current pixel:
        //          i   loop through neighborhood incrementing by 45 degrees
        //          ii  once past the gradient angle, interpolate intensities between current & last neighbor
        //          iii then flip through 180 degrees and do ii again
        //     b) If g_mag of current pixel <= g_mag of pixel behind or pixel in front:
        //          -> set the current pixel value to 0
        //  6. Gather strengths of edge pixels based on low and high thresholds
        //  7. Preform hysteresis; filter Gsup based on if weak edges have surrounding strong edges
        //  8. Alter the image data using the final Canny matrix (normalized of course)
        //  -----------------------------------------------------------------------

        // helper variables
        int n = (N % 2 ? N : N+1);
        int rad = n / 2; //radius of window
        int pixels = img.height * img.width;

        // 1.
        llip::gaussian_averaging_filter(img, sigma);

        // 2.
        Matrix Mx(n,n);
        for (int r=-rad; r<n-rad; r++) {
            for (int c=-rad; c<n-rad; c++) {
                double s2 = sigma*sigma;
                Mx(r+rad,c+rad) = -1*(c/(2*PI*s2*s2))*exp(-1*(c*c+r*r)/(2*s2));
            }
        }
        Mx = Mx * (-1 * Matrix::createIdentity(n)); // horizontal flip
        Matrix My = Mx.transpose();

        // 3. & 4.
        Matrix MIx(img.height,img.width);
        Matrix MIy(img.height,img.width);
        Matrix Gmag(img.height,img.width);
        Matrix Gdir(img.height,img.width);
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                int i = c + (r * img.width);
                Matrix window = llip::get_window(img,rad,i,true);
                double Ix = 0.0, Iy = 0.0;
                for (int wr=0; wr<n; wr++) {
                    for (int wc=0; wc<n; wc++) {
                        Ix += Mx(wr,wc) * window(wr,wc);
                        Iy += My(wr,wc) * window(wr,wc);
                    }
                }
                MIx(r,c) = Ix;
                MIy(r,c) = Iy;
                Gmag(r,c) = sqrt(Ix*Ix+Iy*Iy);
                Gdir(r,c) = atan2(Iy,Ix);
            }
        }

        /* this commented-out snippet alters the image data to be its x,y averaged gradient
        llfd::normalize_gradients(MIy,MIx,img.width,img.height);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = (MIx(r,c)+MIy(r,c))/2;
            }
        }
        return;//*/

        /* this commented-out snippet alters the image data to be its gradient magnitude
        llip::normalize_matrix(Gmag,img.width,img.height,0,255);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = Gmag(r,c);
            }
        }
        return;//*/

        /* this commented-out snippet alters the image data to be its gradient direction
        llip::normalize_matrix(Gdir,img.width,img.height,0,255);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = Gdir(r,c);
            }
        }
        return;//*/

        // 5.
        Matrix Gsup = Gmag;
        for (int i=0; i<pixels; i++) {
            int ir = i / img.width;
            int ic = i % img.width;
            Matrix w = llfd::get_mat_window(Gmag,img.width,img.height,ic,ir);
            int ang = std::abs(Gdir(ir,ic) * 180 / PI); //degrees
            ang %= 180;
            double ifront=-1, iback=-1;
            double small=0, large=0;

            // if 0->44
            if (ang >= 0 && ang < 45) {
                if (ang%45 >= 23) {
                    large = double(ang%45)/45;
                    small = double(45-ang%45)/45;
                    ifront = w(2,1)*small + w(2,0)*large;
                    iback = w(0,1)*small + w(0,2)*large;
                } else {
                    small = double(ang%45)/45;
                    large = double(45-ang%45)/45;
                    ifront = w(2,1)*large + w(2,0)*small;
                    iback = w(0,1)*large + w(0,2)*small;
                }
            }
            // if 45->89
            else if (ang >= 45 && ang < 90) {
                if (ang%45 >= 23) {
                    large = double(ang%45)/45;
                    small = double(45-ang%45)/45;
                    ifront = w(2,0)*small + w(1,0)*large;
                    iback = w(0,2)*small + w(1,2)*large;
                } else {
                    small = double(ang%45)/45;
                    large = double(45-ang%45)/45;
                    ifront = w(2,0)*large + w(1,0)*small;
                    iback = w(0,2)*large + w(1,2)*small;
                }
            }
            // if 90->134
            else if (ang >= 90 && ang < 135) {
                if (ang%45 >= 23) {
                    large = double(ang%45)/45;
                    small = double(45-ang%45)/45;
                    ifront = w(1,0)*small + w(0,0)*large;
                    iback = w(1,2)*small + w(2,2)*large;
                } else {
                    small = double(ang%45)/45;
                    large = double(45-ang%45)/45;
                    ifront = w(1,0)*large + w(0,0)*small;
                    iback = w(1,2)*large + w(2,2)*small;
                }
            }
            // if 135->179
            else if (ang >= 135 && ang < 180) {
                if (ang%45 >= 23) {
                    large = double(ang%45)/45;
                    small = double(45-ang%45)/45;
                    ifront = w(0,0)*small + w(0,1)*large;
                    iback = w(2,2)*small + w(2,1)*large;
                } else {
                    small = double(ang%45)/45;
                    large = double(45-ang%45)/45;
                    ifront = w(0,0)*large + w(0,1)*small;
                    iback = w(2,2)*large + w(2,1)*small;
                }
            }

            if (iback >= w(1,1) || ifront >= w(1,1))
                Gsup(ir,ic) = 0; //suppress non-edge pixel
            else
                Gsup(ir,ic) = Gmag(ir,ic); //preserve edge pixel
        }

        /* this commented-out snippet alters the image data to show line-thinned gradients
        llip::normalize_matrix(Gsup,img.width,img.height,0,255);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = Gsup(r,c);
            }
        }
        return;//*/

        // 6.
        Matrix Strengths(img.height,img.width); // 0=no edge, 1=weak edge, 2=strong edge
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                if (Gsup(r,c) >= high_th)
                    Strengths(r,c) = 2;
                else if (Gsup(r,c))
                    Strengths(r,c) = 1;
                else
                    Strengths(r,c) = 0;
            }
        }

        // 7.
        Matrix Canny(img.height,img.width); // a filtered Gsup
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                if (Strengths(r,c)==2)
                    Canny(r,c) = Gsup(r,c);
                else if (Strengths(r,c)==1) {
                    Matrix w = llfd::get_mat_window(Strengths,img.width,img.height,c,r);
                    bool is_edge = false;
                    for (int wr=0;wr<3;wr++){
                        for (int wc=0;wc<3;wc++){
                            if (wc==1 && wr==1)
                                continue;
                            if (w(wr,wc)==2)
                                is_edge = true;
                            if (is_edge)
                                break;
                        }
                        if (is_edge)
                            break;
                    }
                    if (is_edge)
                        Canny(r,c) = Gsup(r,c);
                }
            }
        }

        // 8.
        llip::normalize_matrix(Canny,img.width,img.height,0,255);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = Canny(r,c);
            }
        }
    }

    //
    //  detect_and_fit_lines
    //
    std::vector<std::vector<double> > detect_and_fit_lines(Image& img, int N, double sigma, double theta_th,
            int grid_size, int line_th) {
        //  -----------------------------------------------------------------------
        //  0. Apply Canny edge detection to the input image before running this function
        //  1. Get the gradient direction of the input image (and a matrix version of the image)
        //  2. Initialize the Accumulator as the same size as the image, filled with zeros
        //  3. For each pixel in the image gradient image at x,y:
        //     a) Use the formula: rho = x*cos(theta) + y*sin(theta) to get a rho value for every theta value
        //        - The hough space matrix should have rho on its y-axis and theta on its x-axis
        //        - The number of theta values can be equivalent to img.width
        //        - Theta begins at zero and is incremented by 2*PI / img.width for each pixel
        //        - Skip when theta is not perpendicular enough to the gradient direction of a pixel
        //     b) Round the rho and theta values to that of the Accumulator's axis points
        //     c) Increment the bin at position (theta,rho) for each theta
        //  4. Obtain the local maxima from the Accumulator
        //  5. Back-project the image from the hough space to x,y Cartesian space using the following:
        //     - if theta=0:    x=rho
        //     - if theta=/=0:  y=(-cos(theta)/sin(theta))*x + (r/sin(theta))
        //     a) scan the x and y variables to get lines
        //  6. Line fit using eigenvectors from the total least squares fitting method
        //  7. Alter the image data to show red lines over the original image
        //  8. Return list of parameters of best fit lines
        //  -----------------------------------------------------------------------

        // helper variables
        int n = (N % 2 ? N : N+1);
        int rad = n / 2; //radius of window
        int pixels = img.height * img.width;
        theta_th = theta_th * PI / 180; //convert from degs to rads

        // 0. (better to do outside this function)
        //llfd::detect_edges(img, N, sigma, c_alpha, c_low_th, c_high_th);

        // 1.
        Matrix Mx(n,n);
        for (int r=-rad; r<n-rad; r++) {
            for (int c=-rad; c<n-rad; c++) {
                double s2 = sigma*sigma;
                Mx(r+rad,c+rad) = -1*(c/(2*PI*s2*s2))*exp(-1*(c*c+r*r)/(2*s2));
            }
        }
        Mx = Mx * (-1 * Matrix::createIdentity(n)); // horizontal flip
        Matrix My = Mx.transpose();
        Matrix MIx(img.height,img.width);
        Matrix MIy(img.height,img.width);
        Matrix Gmag(img.height,img.width);
        Matrix Gdir(img.height,img.width);
        Matrix Canny(img.height,img.width);
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                int i = c + (r * img.width);
                Matrix window = llip::get_window(img,rad,i,true);
                double Ix = 0.0, Iy = 0.0;
                for (int wr=0; wr<n; wr++) {
                    for (int wc=0; wc<n; wc++) {
                        Ix += Mx(wr,wc) * window(wr,wc);
                        Iy += My(wr,wc) * window(wr,wc);
                    }
                }
                MIx(r,c) = Ix;
                MIy(r,c) = Iy;
                Gmag(r,c) = sqrt(Ix*Ix+Iy*Iy);
                Gdir(r,c) = atan2(Iy,Ix);
                Canny(r,c) = img.data[i];
            }
        }

        // 2.
        // rows : rho
        // cols : theta
        Matrix Accumulator(grid_size,grid_size);

        // 3.
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                // ignore black pixels
                if (!Canny(r,c)) continue;
                // map to accumulator space
                for (int t=0; t<grid_size; t++) {
                    if (llfd::is_angle_between(Gdir(r,c)+PI/2, 2*PI/grid_size*t-theta_th, 2*PI/grid_size*t+theta_th))
                        continue;
                    int p = round(c * cos(2*PI/grid_size*t) + r * sin(2*PI/grid_size*t) + grid_size/2);
                    if (p >= 0 && p < grid_size)
                        Accumulator(p,t)++;
                }
            }
        }

        // 4.
        std::vector< std::vector<int> > hough_points;
        for (int p=0; p<grid_size; p++) {
            for (int t=0; t<grid_size; t++) {
                if (Accumulator(p,t) > line_th) {
                    // get local max score in NxN square
                    std::vector<int> point = llfd::nearest_max(t,p,Accumulator,grid_size,grid_size);
                    // add point if it doesn't already exist
                    bool exists = false;
                    for (int i=0;i<hough_points.size();i++) {
                        if (hough_points[i][0]==point[0] && hough_points[i][1]==point[1]) {
                            exists = true;
                            break;
                        }
                    }
                    if (!exists)
                        hough_points.push_back(point);
                }
            }
        }

        /*// For converting the image into hough space output + maxima
        img.height = grid_size;
        img.width = grid_size;
        img.data = std::vector<unsigned int>(grid_size*grid_size);
        llip::normalize_matrix(Accumulator,grid_size,grid_size,0,255);
        for (int r=0;r<grid_size;r++){
            for (int c=0;c<grid_size;c++){
                img.data[c+r*grid_size] = Accumulator(r,c);
            }
        }
        pnmio converter0;
        converter0.convert_gray_to_rgb(img);
        for (int i=0;i<hough_points.size();i++){
            int c = hough_points[i][0];
            int r = hough_points[i][1];
            img.data[(c+r*grid_size)*3] = 255;
            img.data[(c+r*grid_size)*3+1] = 0;
            img.data[(c+r*grid_size)*3+2] = 0;
        }
        return std::vector<std::vector<double> >();//*/

        // 5.
        // map hough points to Cartesian points via slope-intercept form (t,p -> m,b -> x,y)
        Matrix Hough(img.height,img.width);
        std::vector<std::vector<double> > params;
        for (int i=0; i<hough_points.size(); i++) {
            // map accumulator values to actual hough space values:
            //   theta = theta_inc * theta value from hough space point
            //   rho   = rho_inc   * hough_points[i][1]
            double theta = hough_points[i][0] * 2*PI / grid_size;
            double p = hough_points[i][1];
            double avgx = 0;
            double avgy = 0;
            std::vector<int> edge_x_pxs;
            std::vector<int> edge_y_pxs;

            for (int x=0; x<img.width; x++) {
                int y = round((p - grid_size/2 - x*cos(theta)) / sin(theta));
                if (y >= 0 && y < img.height) {
                    Hough(y,x) = 255; //color pixel white
                    edge_x_pxs.push_back(x);
                    edge_y_pxs.push_back(y);
                    avgx += x;
                    avgy += y;
                }
            }

            for (int y=0; y<img.height; y++) {
                int x = round((p - grid_size/2 - y*sin(theta)) / cos(theta));
                if (x >= 0 && x < img.height)
                    Hough(y,x) = 255; //color pixel white
            }

            // WARNING: MAY SLOW DOWN THIS LOOP A LOT
            // 6. find best fit line from edge pixels collected
            int edge_px_amt = edge_x_pxs.size();
            avgx /= edge_px_amt;
            avgy /= edge_px_amt;
            Matrix mat_A(edge_px_amt,2);
            for (int r=0; r<edge_px_amt; r++) {
                mat_A(r,0) = edge_x_pxs[r]-avgx; //x
                mat_A(r,1) = edge_y_pxs[r]-avgy; //y
            }
            Matrix ATA = mat_A.transpose() * mat_A;
            std::vector<double> min_eig_vec = llfd::min_eigen_vec_2d(ATA);
            params.push_back(min_eig_vec);//*/
        }

        /*// 7. (alternate - white lines on black bg)
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                img.data[c+r*img.width] = Hough(r,c);
            }
        }
        return params;//*/


        // 7.
        pnmio converter;
        converter.convert_gray_to_rgb(img);
        for (int r=0;r<img.height;r++){
            for (int c=0;c<img.width;c++){
                if (Hough(r,c) == 255) {
                    img.data[3*(c+r*img.width)] = 255; // R
                    img.data[3*(c+r*img.width)+1] = 0; // G
                    img.data[3*(c+r*img.width)+2] = 0; // B
                }
            }
        }

        // 8.
        return params;
    }
}

