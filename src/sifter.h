////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       sifter.h
//  Date:           04/14/2018
//  Description:    Header for SIFT'er (Scale Invariant Feature Transform) class.
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

#ifndef _SIFTER_H_
#define _SIFTER_H_

#include "pnmio.h"
#include "matrix.h"

#define SQRT2 1.41421356237 //quick square root of 2

/**
*   ...
*/
struct Keypoint {
    int x = 0;
    int y = 0;
    double mag = 0.0;
    double dir = 0.0;
    int bin = 0; // dir = 0-9 means 1st bin, 10-19 2nd bin, 20-29 etc...
};

/**
*   @namespace sifter
*   @see pnmio::Image @struct
*   @see Matrix @class
*
*   Contains SIFT related functions.
*   Note: utilizes pnmio's Image structure
*   Note: utilizes Matrix class
*/
namespace sifter {

    /**
    *   *Helper Function*
    *
    */
    double dot(Matrix&, Matrix&);

    /**
    *   *Helper Function*
    *
    */
    Matrix img_to_mat(Image&);

    /**
    *   *Helper Function*
    *
    */
    Image mat_to_img(Matrix&, int, int, int type=0, int mxv=255);

    /**
    *   *Helper Function*
    *
    */
    Matrix shrink_octave(Matrix&,int,int);

    /**
    *   *Helper Function*
    *
    */
    bool exists_in(std::vector<int>, std::vector<std::vector<int> >);

    /*/**
    *   *Helper Function*
    *
    * /
    std::vector<std::vector<int> > find_min_max(Matrix&, int, int);*/

    /**
    *   *Helper Function*
    *
    */
    Matrix quick_DOG_blur(Matrix&,int,int,int,double,double k=SQRT2);

    /**
    *   *Helper Function*
    *
    */
    std::vector<std::vector<Matrix> > get_DOGs(Image&, std::vector<std::vector<Matrix> >&);

    /**
    *   Returns octaves of matrices from the DOG of the scale space of an image
    *   (k value always begins at sqrt of 2)
    *   (DoG = Difference of Gaussian)
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          int                     blur window size (always odd)
    *   @param          double                  sigma value for Gaussian
    *   @return         vector                  scale space octave list
    */
    std::vector<std::vector<Matrix> > get_scale_space(Image&,int,double);

    /**
    *   Returns a list of sift interest points from a list of DoG octaves
    *
    *   @param          Image &                 gray-scale pnmio Image structure
    *   @param          vector &                scale space octave list
    *   @param          vector &                DoG octave list
    *   @param          int                     Corner filtering threshold (window size=3,sigma=1,alpha=0.04)
    *   @return         vector                  sift interest point-painted octaves of s filtered DoGs
    */
    std::vector<std::vector<Matrix> > get_sift_pts(Image&, std::vector<std::vector<Matrix> >&,
                                                   std::vector<std::vector<Matrix> >&, int);

    /**
    *
    */
    std::vector<std::vector<Keypoint> > get_keypoints(Image&, std::vector<std::vector<Matrix> >&,
                                                      std::vector<std::vector<Matrix> >&);

}

#endif   // !defined _SIFTER_H_


