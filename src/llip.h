////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       llip.h
//  Date:           03/14/2018
//  Description:    Header for LLIP (low level image processing) class.
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

#ifndef _LLIP_H_
#define _LLIP_H_

#include "pnmio.h"
#include "matrix.h"

/**
*   @namespace llip
*   @see pnmio::Image @struct
*   @see Matrix @class
*
*   Contains multiple low level image processing functions.
*   Note: utilizes pnmio's Image structure
*   Note: utilizes Matrix class
*/
namespace llip {

    /**
    *   * Helper function *
    *   Returns a (1+2R)x(1+2R) window of elements from a gray-scale image at index I
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          int (R)                 "radius" of the window from the index i
    *   @param          int (I)                 index of window (center of window matrix)
    *   @param          bool=false              replace out-of-bound pixels with -1
    *   @param          bool=true               replace out-of-bound pixels with closest pixel value
    *   @return         Matrix                  square window matrix
    */
    Matrix get_window(Image&, int, int, bool=false);

    /**
    *   * Helper function *
    *   Normalizes a matrix (more like scales it)
    *
    *   @param          Matrix &                Matrix object to normalize
    *   @param          u-int                   width of matrix
    *   @param          u-int                   height of matrix
    *   @param          double                  lowest possible normalized value
    *   @param          double                  highest possible normalized value
    */
    void normalize_matrix(Matrix&, unsigned int, unsigned int, double smin=0, double smax=1);

    /**
    *   * Helper function *
    *   Normalizes a matrix (sum of all values = to 1)
    *
    *   @param          Matrix &                Matrix object to normalize
    *   @param          u-int                   width of matrix
    *   @param          u-int                   height of matrix
    */
    void normalize_matrix2(Matrix&, unsigned int, unsigned int);

    /**
    *   * Helper function *
    *   Returns true if a value is found in a Matrix object, else false
    *
    *   @param          double                  value to check for
    *   @param          Matrix &                Matrix object to search
    *   @param          u-int                   width of matrix
    *   @param          u-int                   height of matrix
    *   @return         true                    value exists in Matrix
    *   @return         false                   value does not exist in Matrix
    */
    bool exists_in(double, Matrix&, unsigned int, unsigned int);

    /**
    *   * Helper function *
    *   Returns a vector of doubles, sorted in row-major format from the input Matrix
    *
    *   @param          Matrix &                Matrix object to search
    *   @param          u-int                   width of matrix
    *   @param          u-int                   height of matrix
    *   @return         vector                  row-major vector of doubles
    */
    std::vector<double> mat_to_vec(Matrix&, unsigned int, unsigned int);

    /**
    *   Manipulates the input gray-scale image by applying histogram equalization
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    */
    void histogram_equalization(Image&);

    /**
    *   Manipulates the input gray-scale image by mapping the intensities using a log function
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          double                  factor that scales the logarithmic mapped output
    */
    void logarithmic_mapping(Image&, double);

    /**
    *   Rotates the input gray-scale image around its center (crops parts that are out of bounds)
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          double                  amount of rotation in degrees
    */
    void centered_rotation(Image&, double);

    /**
    *   Blurs the input gray-scale image using a 2D Gaussian averaging filter
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          double                  sigma value - also determines size of window/mask
    */
    void gaussian_averaging_filter(Image&, double);

    /**
    *   Smooths the input gray-scale image using a 2D median filter
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          u-int                   radius of square window/mask
    */
    void median_neighborhood_filter(Image&, unsigned int=1);

}

#endif   // !defined _LLIP_H_
