
////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       llfd.cpp
//  Date:           04/29/2018
//  Description:    Main implementation for depth mapper class.
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
#include "depthmapper.h"
#include "sifter.h"
#include "llfd.h"
#include "llip.h"
#include "matrix.h"

#include "math.h"
#include <iostream>
#include <stdexcept>
#include <algorithm>

// for quicker PI calculations
#define PI 3.14159265

namespace depthmapper {

    Matrix create_fundamental(Image& left, Image& right, std::vector<std::vector<int> > left_matches,
                              std::vector<std::vector<int> > right_matches) {

        // convert to mats
        Matrix lm = sifter::img_to_mat(left);
        Matrix rm = sifter::img_to_mat(right);

        // translation matrix
        Matrix T(3,3);
        //...ran out of time (spent too much time on SIFT)

    }

    Matrix get_disparity(Image& left, Image& right, int ds, int ws, double maxs) {
        // ds = descriptor matrix size
        // ws = epipolar line window size
        // maxs = max SSD allowed when matching points
        std::cout<<"Calculating depth (actually it's disparity) map..."<<std::endl;
        Matrix lm = sifter::img_to_mat(left);
        Matrix rm = sifter::img_to_mat(right);
        int w = left.width;
        int h = left.height;
        Matrix depth(h,w);

        // looping through the left, look for corresponding points in the right
        // (along the epipolar line - horiz since they are both rectified)
        // for each pixel in left
        for (int r=0; r<h; r++) {
            for (int c=0; c<w; c++) {
                Matrix left_descriptor = llfd::get_mat_window(lm,w,h,c,r,ds,ds,lm(r,c));
                std::vector<Matrix> right_descriptors;
                int w_rad = ws/2;
                // for each pixel in right (same row only), gather a ds x ds window
                int start_coord = -1;
                for (int c2=c-w_rad; c2<=c+w_rad; c2++) {
                    if (c2 >= 0 && c2 < w) {
                        right_descriptors.push_back(llfd::get_mat_window(rm,w,h,c2,r,ds,ds,rm(r,c)));
                        // record starting coord:
                        if (start_coord == -1)
                            start_coord = c2;
                    }
                }
                // compute sum of squared differences (for each right descriptor gathered)
                std::vector<double> SSDs;
                for (int i=0; i<right_descriptors.size(); i++) {
                    Matrix& m = right_descriptors[i];
                    double ssd = 0;
                    for (int j=0; j<ds; j++) {
                        for (int k=0; k<ds; k++) {
                            if (m(j,k) != -1 && left_descriptor(j,k) != -1) {
                                double diff = left_descriptor(j,k) - m(j,k);
                                ssd += diff*diff;
                            }
                        }
                    }
                    SSDs.push_back(ssd);
                }
                // pick best match (if below maxs threshold)
                double smallest = 99999999;
                int smallest_c = -1;
                for (int i=0; i<SSDs.size(); i++) {
                    if (smallest > SSDs[i]) {
                        smallest = SSDs[i];
                        smallest_c = start_coord+i;
                    }
                }
                if (smallest < maxs) {
                    depth(r,smallest_c) = (c - smallest_c)*(c - smallest_c); //disparity
                } else {
                    depth(r,smallest_c) = 0;
                }
            }
            std::cout<<"Row "<<r<<" completed."<<std::endl;
        }

        // normalize depth map to be able to see depth gradient
        llip::normalize_matrix(depth, w, h, 0, 255);

        std::cout<<"Depth (actually it's disparity) map completed."<<std::endl;
        return depth;
    }

}
