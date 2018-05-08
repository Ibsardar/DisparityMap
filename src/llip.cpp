////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       llip.cpp
//  Date:           03/14/2018
//  Description:    Main implementation for LLIP (low level image processing) class.
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
#include "llip.h"
#include "matrix.h"

#include "math.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>

// for quicker PI calculations
#define PI 3.14159265

namespace llip {

    //
    //  get_window
    //
    Matrix get_window(Image& img, int r, int i, bool nearest) {

        // helper variables
        const int n = r+r+1; // window size
        const int irow = i / img.width; // the row or y coordinate in the image of the i'th pixel
        const int icol = i % img.width; // the column or x coordinate in the image of the i'th pixel

        // create and fill window matrix
        Matrix m(n,n);
        for (int y=irow-r, mrow=0; y<=irow+r; y++, mrow++) {
            for (int x=icol-r, mcol=0; x<=icol+r; x++, mcol++) {
                if (x>=0 && x<(int)img.width && y>=0 && y<(int)img.height) { // if within image bounds
                    int intensity = img.data.at(x+(y*(int)img.width));
                    m(mrow,mcol) = intensity;
                } else {
                    if (!nearest) {
                        m(mrow,mcol) = -1;
                    } else {

                        // *** TOO LAZY TO MAKE A SEPERATE FUNCTION FOR THIS SORRY ***
                        // -Find distances and intensities of all existing pixels in this window from the current pixel
                        // vector of Manhattan distances and a vector of their corresponding intensity:
                        std::vector<int> dists, intensities;
                        for (int y2=irow-r, mrow2=0; y2<=irow+r; y2++, mrow2++) {
                            for (int x2=icol-r; x2<=icol+r; x2++) {
                                if (x2>=0 && x2<(int)img.width && y2>=0 && y2<(int)img.height) { // if within image bounds
                                    int intensity = img.data.at(x2+(y2*(int)img.width));
                                    dists.push_back(abs(x2-x) + abs(y2-y));
                                    intensities.push_back(intensity);
                                }
                            }
                        }

                        // -Find nearest colored pixel intensity
                        int nearest_index = 0;
                        int nearest_distance = 999999;
                        for (int ni=0; ni<dists.size(); ni++) {
                            if (nearest_distance > dists[ni]) {
                                nearest_distance = dists[ni];
                                nearest_index = ni;
                            }
                        }

                        // -Set nearest colored pixel intensity
                        m(mrow,mcol) = intensities[nearest_index];
                    }
                }
            }
        }

        // return window
        return m;
    }

    //
    //  normalize_matrix
    //
    void normalize_matrix(Matrix& m, unsigned int w, unsigned int h, double smin, double smax) {

        // get min and max values
        double min_val=999999999, max_val=-999999999;
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                if (min_val > m(r,c))
                    min_val = m(r,c);
                if (max_val < m(r,c))
                    max_val = m(r,c);
            }
        }

        // scale to range: smin to smax (ex: 0 to 1)
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                m(r,c) = smin + (smax - smin) * (m(r,c) - min_val) / (max_val - min_val);
            }
        }

        /*** (faster if only between 0 an 1):

        // get sum of matrix values
        double sum=0;
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                sum += m(r,c);
            }
        }

        // normalize each element
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                m(r,c) = m(r,c) / sum;
            }
        }
        ***/
    }

    //
    //  normalize_matrix2
    //
    void normalize_matrix2(Matrix& m, unsigned int w, unsigned int h) {

        // get sum of matrix values
        double sum=0;
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                sum += m(r,c);
            }
        }

        // normalize each element
        for (unsigned int r=0; r<h; r++) {
            for (unsigned int c=0; c<w; c++) {
                m(r,c) = m(r,c) / sum;
            }
        }
    }

    //
    //  exists_in
    //
    bool exists_in(double value, Matrix& m, unsigned int w, unsigned int h) {
        for (unsigned int r=0; r<h; r++)
            for (unsigned int c=0; c<w; c++)
                if (m(r,c) == value)
                    return true;
        return false;
    }

    //
    //  mat_to_vec
    //
    std::vector<double> mat_to_vec(Matrix& m, unsigned int w, unsigned int h) {
        std::vector<double> vec;
        for (unsigned int r=0; r<h; r++)
            for (unsigned int c=0; c<w; c++)
                vec.push_back(m(r,c));
        return vec;
    }

    //
    //  histogram_equalization
    //
    void histogram_equalization(Image& img) {
        //  -----------------------------------------------------------------------
        //  Manipulates the input gray-scale image by applying histogram normalization:
        //  1) Record frequency of each unique pixel intensity from 0 to 255
        //  2) Record probability of each unique pixel intensity (frequency / total # of pixels)
        //  3) Record the cumulative probability for each unique pixel intensity from 0 to 255 (previous P + current P)
        //  4) Normalize the CP's to 255 since that is our maximum range ( floor(CP * 255) )
        //  5) Match the normalized CP's to the input gray-scale image
        //  -----------------------------------------------------------------------

        // helper variables
        const unsigned int pixels = img.width * img.height;
        const unsigned int max_intensity = 256; // technically it's 255, but this makes the below easier

        // 1.
        std::vector<unsigned int> freqs(max_intensity, 0);
        for (unsigned int i=0; i<pixels; i++) {
            unsigned int intensity = img.data[i];
            freqs[intensity]++;
        }

        // 2.
        std::vector<double> probs(max_intensity, 0.0);
        for (unsigned int i=0; i<max_intensity; i++)
            probs[i] = (double) freqs[i] / pixels;

        // 3.
        std::vector<double> cprobs(max_intensity, 0.0);
        cprobs[0] = probs[0];
        for (unsigned int i=1; i<max_intensity; i++)
            cprobs[i] = cprobs[i-1] + probs[i];

        // 4.
        std::vector<unsigned int> norms(max_intensity, 0);
        for (unsigned int i=0; i<max_intensity; i++)
            norms[i] = floor(cprobs[i] * (max_intensity-1));

        // 5.
        for (unsigned int i=0; i<pixels; i++) {
            unsigned int intensity = img.data[i];
            img.data[i] = norms[intensity];
        }
    }

    //
    //  logarithmic_mapping
    //
    void logarithmic_mapping(Image& img, double factor) {
        //  -----------------------------------------------------------------------
        //  Manipulates the input gray-scale image by mapping the intensities using a log function:
        //  1) For each pixel's intensity value 'i', change it to: factor * log (base 10) of i+1
        //  2) Re-normalize the data to fit in the 0 to 255 range
        //  -----------------------------------------------------------------------

        // helper variables
        const unsigned int pixels = img.width * img.height;

        // 1.
        std::vector<double> tmp_data;
        for (unsigned int i=0; i<pixels; i++)
            tmp_data.push_back(factor * log10(double(img.data[i]+1)));

        // 2.
        for (unsigned int i=0; i<pixels; i++)
            img.data[i] = floor(tmp_data[i] * 255 / log10(255+1));
            /*
        double imax=0, imin=999999;
        for (unsigned int i=0; i<pixels; i++) {
            if (imax < tmp_data[i])
                imax = tmp_data[i];
            if (imin > tmp_data[i])
                imin = tmp_data[i];
        }
        for (unsigned int i=0; i<pixels; i++)
            img.data[i] = round((tmp_data[i]-imin) / (imax-imin) * 255);
            */
    }

    //
    //  centered_rotation
    //
    void centered_rotation(Image& img, double angle) {
        //  -----------------------------------------------------------------------
        //  Rotates the input gray-scale image around its center:
        //  1) Create a 2xN matrix holding all rotated pixel coordinates, where N=# of pixels
        //  2) Calculate the 2D rotation matrix using the negative of the input angle (convert angle to radians first)
        //  3) Calculate the original 2xN set of pixel coordinates from the rotated pixel coordinates:
        //     - translate leftwards width/2 and upwards height/2
        //     - rotate
        //     - translate rightwards width/2 and downwards height/2
        //  4) Map all the rotated intensities to a new Image data vector:
        //     - rotated Image should initialize to black (intensity=0)
        //  5) Set the the new vector to the original Image data vector
        //  -----------------------------------------------------------------------

        // helper variables
        const unsigned int pixels = img.width * img.height;

        // 1.
        Matrix new_coords_mat(2, pixels);
        for (unsigned int c=0; c<pixels; c++) {
            new_coords_mat(0,c) = c % img.width; // x-coordinate of pixel
            new_coords_mat(1,c) = c / img.width; // y-coordinate of pixel
        }

        // 2.
        double rads = angle * PI / 180 * -1;
        Matrix rot_mat(2, 2);
        rot_mat(0,0) = cos(rads); rot_mat(1,0) = -sin(rads);
        rot_mat(0,1) = sin(rads); rot_mat(1,1) = cos(rads);

        // 3.
        unsigned int shiftx = img.width/2;
        unsigned int shifty = img.height/2;
        Matrix orig_coords_mat = new_coords_mat;
        for (unsigned int c=0; c<pixels; c++) {
            orig_coords_mat(0,c) -= shiftx;
            orig_coords_mat(1,c) -= shifty;
        }
        orig_coords_mat = rot_mat * orig_coords_mat;
        for (unsigned int c=0; c<pixels; c++) {
            orig_coords_mat(0,c) += shiftx;
            orig_coords_mat(1,c) += shifty;
        }

        // 4.
        std::vector<unsigned int> new_data(pixels, 0);
        for (unsigned int c=0; c<pixels; c++) {
            int x = round(orig_coords_mat(0,c));
            int y = round(orig_coords_mat(1,c));
            if (x >= 0 && x < img.width &&
                y >= 0 && y < img.height) {
                new_data[c] = img.data[x+(y*img.width)];
            }
        }

        // 5.)
        img.data = new_data;
    }

    //
    //  gaussian_averaging_filter
    //
    void gaussian_averaging_filter(Image& img, double sigma) {
        //  -----------------------------------------------------------------------
        //  Blurs the input gray-scale image:
        //  1) Create a Gaussian mask using the formula: (1/(2*PI*sigma))*e^((x*x+y*y)/(-2*sigma))
        //     a) determine the size of the mask square solely based on sigma (std-deviation from the center)
        //     b) fill the mask using the above formula while making sure middle is the center
        //     c) normalize the mask
        //  2) Create a vector of the same size as the # of pixels
        //  3) For each pixel:
        //     - create a window centered on the pixel
        //     - re-normalize mask if there are out-of-bounds pixels (set those areas to 0 on mask)
        //     - multiply then add up the mask on the window
        //     - set the corresponding pixel in the new vector to the result of the mask on this pixel
        //  4) Set the the new vector to the original Image data vector
        //  -----------------------------------------------------------------------

        // helper variables
        const unsigned int pixels = img.width * img.height;

        // 1.a)
        const unsigned int mask_size = ((int)round(5*sigma)%2==1 ? round(5*sigma) : round(5*sigma)+1); // always odd
        const unsigned int mask_radius = (mask_size-1)/2;

        // 1.b)
        Matrix mask(mask_size, mask_size);
        for (unsigned int i=0; i<mask_size; i++) {
            for (unsigned int j=0; j<mask_size; j++) {
                unsigned int x=i-mask_radius, y=j-mask_radius;
                mask(i,j) = (1/(2*PI*sigma))*exp((x*x+y*y)/(-2*sigma)); // 2D Gaussian filter function
            }
        }

        // 1.c)
        normalize_matrix2(mask, mask_size, mask_size);

        // 2.)
        std::vector<unsigned int> new_data(pixels, 0);

        // 3.)
        for (unsigned int p=0; p<pixels; p++) {
            Matrix window = get_window(img, mask_radius, p, false);
            Matrix mask_to_apply = mask;
            double gaussed_pixel = 0;

            // * comment this part out if the boolean above is set to true instead of false * <-SLOW METHOD
            // re-normalize mask to match only in-bounds pixels <-FAST METHOD
            if (exists_in(-1, window, mask_size, mask_size)) {
                for (unsigned int r=0; r<mask_size; r++)
                    for (unsigned int c=0; c<mask_size; c++)
                        if (window(r,c) < 0)
                            mask_to_apply(r,c) = 0;
            }
            normalize_matrix2(mask_to_apply, mask_size, mask_size);

            // apply weights to pixel
            for (unsigned int r=0; r<mask_size; r++)
                for (unsigned int c=0; c<mask_size; c++)
                    if (window(r,c) >= 0)
                        gaussed_pixel += mask_to_apply(r,c) * window(r,c);
            new_data[p] = gaussed_pixel;
        }

        // 4.)
        img.data = new_data;
    }

    //
    //  median_neighborhood_filter
    //
    void median_neighborhood_filter(Image& img, unsigned int radius) {
        //  -----------------------------------------------------------------------
        //  Smooths the input gray-scale image:
        //  1) Create vector of the same size as the # of pixels
        //  2) For each pixel:
        //     - get a r x r window list at the pixel and convert it into a vector (use nearest in-bounds pixel)
        //     - filter its values to exclude -1s (if needed) and sort them in ascending order
        //     - set the corresponding pixel in the new vector to the median of the sorted list
        //  3) Set the the new vector to the original Image data vector
        //  -----------------------------------------------------------------------

        // helper variables
        const unsigned int pixels = img.width * img.height;
        const unsigned int filter_size = 1 + 2 * radius;

        // 1.)
        std::vector<unsigned int> new_data(pixels, 0);

        // 2.)
        for (unsigned int p=0; p<pixels; p++) {

            // get sorted elements in 3x3 window
            Matrix window = get_window(img, radius, p, true);
            std::vector<double> data = mat_to_vec(window, filter_size, filter_size);
            /*
            // un-comment this if the boolean above is set to false instead of true
            for (unsigned int i=0; i<data.size(); i++)
                if (data[i] < 0)
                    data.erase(data.begin()+i);
            */
            std::sort(data.begin(), data.end());

            // get median
            unsigned int median = 0;
            // if size is odd:
            if (data.size() % 2 == 1)
                median = data[data.size()/2+1];
            // if size is even:
            else
                median = round((data[data.size()/2] + data[data.size()/2+1]) / 2);

            // set new pixel
            new_data[p] = median;
        }

        // 3.)
        img.data = new_data;
    }

}
