////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       depthmapper.h
//  Date:           04/29/2018
//  Description:    Header for depth mapper class.
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

#ifndef _DEPTHMAPPER_H_
#define _DEPTHMAPPER_H_

#include "pnmio.h"
#include "matrix.h"

/**
*   @namespace depthmapper
*   @see pnmio::Image @struct
*   @see Matrix @class
*
*   Contains multiple depth map generating related functions.
*   Note: utilizes pnmio's Image structure
*   Note: utilizes Matrix class
*/
namespace depthmapper {

    /**
    *   Converts the given left & right images and the 8 matching points on each image into a fundamental matrix
    */
    Matrix create_fundamental(Image&, Image&, std::vector<std::vector<int> >, std::vector<std::vector<int> >);

    /**
    *   Computes the disparity between the left and right rectified images
    */
    Matrix get_disparity(Image&, Image&, int ds=3, int ws=5, double maxs=999);

}

#endif   // !defined _DEPTHMAPPER_H_
