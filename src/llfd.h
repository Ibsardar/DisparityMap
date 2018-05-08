////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       llfd.h
//  Date:           03/31/2018
//  Description:    Header for LLFD (low level feature detection) class.
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

#ifndef _LLFD_H_
#define _LLFD_H_

#include "pnmio.h"
#include "matrix.h"

/**
*   @namespace llfd
*   @see pnmio::Image @struct
*   @see Matrix @class
*   @see llip @namespace
*
*   Contains multiple low level feature detection functions.
*   Note: utilizes pnmio's Image structure
*   Note: utilizes Matrix class
*   Note: utilizes llip namespace
*/
namespace llfd {

    /**
    *   * Helper function *
    *   Returns a window of values from a gray-scale image around a certain point
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          int                     index of row-major, image vector
    *   @param          int                     desired window width (always odd)
    *   @param          int                     desired window height (always odd)
    *   @param          double                  desired filler for non-existent cells
    *   @return         Matrix                  window of intensity values centered at selected pixel
    */
    Matrix get_mat_window(Image&, int, int, int, double rep=0);

    /**
    *   * Helper function *
    *   Returns a window of values from a matrix around a certain point
    *
    *   @param          Matrix &                matrix of values
    *   @param          int                     matrix width
    *   @param          int                     matrix height
    *   @param          int                     window center column
    *   @param          int                     window center row
    *   @param          int                     desired window width (always odd)
    *   @param          int                     desired window height (always odd)
    *   @param          double                  desired filler for non-existent cells
    *   @return         Matrix                  window of intensity values centered at selected pixel
    */
    Matrix get_mat_window(Matrix&, int, int, int, int, int ww=3, int wh=3, double rep=0);

    /**
    *   * Helper function *
    *   Returns the trace of an NxN matrix
    *
    *   @param          Matrix                  input square matrix
    *   @param          u-int                   size of matrix (NxN)
    *   @return         double                  trace of matrix
    */
    double trace(Matrix, unsigned int);

    /**
    *   * Helper function *
    *   Returns the determinant of a 2x2 matrix
    *
    *   @param          Matrix                  input square matrix
    *   @return         double                  determinant of matrix
    */
    double det2d(Matrix);

    /**
    *   * Helper function *
    *   Normalizes 2 gradient matrices (between 0 and 255)
    *
    *   @param          Matrix &                a Matrix object to normalize
    *   @param          Matrix &                a Matrix object to normalize
    *   @param          u-int                   width of matrix
    *   @param          u-int                   height of matrix
    */
    void normalize_gradients(Matrix&, Matrix&, unsigned int, unsigned int);

    /**
    *   * Helper function *
    *   Returns closest local maximum point (uses 3x3 neighbor)
    *
    *   @param          int                     column
    *   @param          int                     row
    *   @param          Matrix &                Matrix of values
    *   @param          u-int                   width of matrix
    *   @param          u-int                   height of matrix
    *   @return         vector                  col,row coordinates of local max
    */
    std::vector<int> nearest_max(int, int, Matrix&, unsigned int, unsigned int);

    /**
    *   * Helper function *
    *   Detects whether some angle is between two others
    *
    *   @param          double                  angle to check
    *   @param          double                  start angle
    *   @param          double                  end angle
    *   @return         true                    angle is between start and end
    *   @return         false                   angle is not between start and end
    */
    bool is_angle_between(double, double, double);

    /**
    *   * Helper function *
    *   Returns double roughly rounded to 2 decimal places
    *
    *   @param          double                  decimal to round
    *   @return         double                  rounded decimal output
    */
    double round2(double);

    /**
    *   * Helper function *
    *   Returns the eigenvector of the minimum eigenvalue of the input matrix
    *
    *   @param          Matrix &                input 2x2 matrix
    *   @return         vector                  eigenvector of min eigenvalue
    */
    std::vector<double> min_eigen_vec_2d(Matrix&);

    /**
    *   Returns a list of corners detected by the Harris Corner Detector.
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          u-int                   N (window size)
    *   @param          double                  sigma; constant affecting the Gaussian blur effect
    *   @param          double                  alpha; constant affecting the image gradient effect
    *   @param          u-int                   threshold (limits how much "cornerness" allowed)
    *   @param          boolean                 whether corners should be red on image or white on black
    *   @return         vector                  list of corner coordinates in input image
    */
    std::vector<std::vector<int> > detect_corners(Image&, unsigned int, double, double, unsigned int, bool=true);

    /**
    *   Using the Canny method, converts the image so we can see the edges
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          int                     N (window size)
    *   @param          double                  sigma; constant affecting the Gaussian blur effect
    *   @param          double                  alpha; constant affecting the image gradient effect
    *   @param          double                  low threshold for defining weak edges
    *   @param          double                  high threshold for defining strong edges
    */
    void detect_edges(Image&, int, double, double, double, double);

    /**
    *   Using the Hough space and total squares minimization method, detects lines and fits them by
    *   drawing them red on the input image
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          int                     N (window size)
    *   @param          double                  sigma; constant affecting the Gaussian blur effect
    *   @param          double                  Delta threshold for theta (higher value means allow more lines)
    *   @param          int                     size of accumulator (shouldn't be too large or too small)
    *   @param          int                     line threshold; higher value means fewer line captures
    *   @return         vector                  list of optimal parameters for best fitting line
    */
    std::vector<std::vector<double> > detect_and_fit_lines(Image&, int, double, double, int, int);

}

#endif   // !defined _LLFD_H_

