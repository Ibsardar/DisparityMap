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

// best recommended SIFT scale space parameters
#define SCALES  5 // blur levels
#define OCTAVES 4 // resolutions
#define DOGCNT  4 // 4 Difference of Gaussian images (minimum recommended for SIFT)

namespace sifter {

    //
    //  dot
    //
    double dot(Matrix& a, Matrix& b) {
        return (a*b)(0,0);
    }

    //
    //  img_to_mat
    //
    Matrix img_to_mat(Image& img) {
        int n = img.height * img.width;
        Matrix m(img.height, img.width);
        for (int i=0; i<n; i++){
            int r = i / img.width;
            int c = i % img.width;
            m(r,c) = img.data[i];
        }
        return m;
    }

    //
    //  mat_to_img
    //
    Image mat_to_img(Matrix& m, int w, int h, int type, int mxv) {
        Image img;
        img.width = w;
        img.height = h;
        img.type = type;
        img.max_val = mxv;
        for (int r=0; r<h; r++)
            for (int c=0; c<w; c++)
                img.data.push_back(m(r,c));
        return img;
    }

    //
    //  shrink_octave
    //
    Matrix shrink_octave(Matrix& m, int w, int h) {
        // shrinks matrix by half via nearest neighbor
        // TODO (MAYBE):
        // - get 3x3 window around each pixel gathered from m
        // - take the avg of that window and set it at m2 instead of N.N. method
        // - ignore out-of-bounds pixels
        int w2 = w/2;
        int h2 = h/2;
        Matrix m2(h2,w2);
        for (int r=0;r<h2;r++){
            for (int c=0;c<w2;c++){
                m2(r,c) = m(r*2,c*2);
            }
        }
        return m2;
    }

    //
    //  exists_in
    //
    bool exists_in(std::vector<int> pt, std::vector<std::vector<int> > pts) {
        for (int i=0; i<pts.size(); i++)
            if (pts[i][0] == pt[0] && pts[i][1] == pt[1])
                return true;
        return false;
    }

    /*//
    //  find_min_max
    //
    std::vector<std::vector<int> > find_min_max(Matrix& m, int w, int h) {
        std::vector<std::vector<int> > pts;
        for (int r=0;r<h;r++){
            for (int c=0;c<w;c++){
                std::vector<int> ptmax = llfd::nearest_max(c,r,m,w,h);
                std::vector<int> ptmin = llfd::nearest_min(c,r,m,w,h);
                if (!exists_in(ptmax, pts))
                    pts.push_back(ptmax);
                if (!exists_in(ptmin, pts))
                    pts.push_back(ptmin);
            }
        }
        return pts;
    }*/

    ///JUNK***
    //
    //  quick_DOG_blur
    //
    Matrix quick_DOG_blur(Matrix& m, int w, int h, int n, double sigma, double k) {
        // calculate Difference of Gaussian
        Matrix a(h,w);
        Matrix b(h,w);
        int radius=n/2;


        /****** IF n is large, utilizing separability will be FASTER ******

        Matrix ga(1,n); // row mask
        Matrix gb(1,n); // row mask
        // create 2 row & 2 col masks (utilizing separability of Gaussian)
        for (int i=-radius; i<=radius; i++){
            ga(0,i+radius) = (1/(2*PI*sigma*sigma))*exp((i*i)/(-2*sigma*sigma));
            gb(0,i+radius) = (1/(2*PI*sigma*k*sigma*k))*exp((i*i)/(-2*sigma*k*sigma*k));
        }
        llip::normalize_matrix2(ga,n,1);
        llip::normalize_matrix2(gb,n,1);
        Matrix gaT = ga.transpose(); // col mask
        Matrix gbT = gb.transpose(); // col mask

        // Gaussian blur a and b in x direction
        Matrix atmp = a; // intermediate image of a (convolved only in x-dir)
        Matrix btmp = b; // intermediate image of b (convolved only in x-dir)
        for (int y=0;y<h;y++){
            for (int x=0;x<w;x++){
                Matrix win_row = llfd::get_mat_window(m,w,h,x,y,n,1,-1);

                // handle edge cases
                Matrix mask1 = gaT;
                Matrix mask2 = gbT;
                bool edge = false;
                for (int i=0;i<n;i++) {
                    if (win_row(0,i) == -1) {
                        edge = true;
                        mask1(i,0) = 0;
                        mask2(i,0) = 0;
                        win_row(0,i) = 0;
                    }
                }
                if (edge) {
                    llip::normalize_matrix2(mask1,1,n);
                    llip::normalize_matrix2(mask2,1,n);
                    std::cout<<mask1<<std::endl;
                }

                atmp(y,x) = dot(win_row,mask1);
                btmp(y,x) = dot(win_row,mask2);
            }
        }

        // Gaussian blur a and b in y direction
        for (int y=0;y<h;y++){
            for (int x=0;x<w;x++){
                Matrix wina_col = llfd::get_mat_window(atmp,w,h,x,y,1,n,-1);
                Matrix winb_col = llfd::get_mat_window(btmp,w,h,x,y,1,n,0);

                // handle edge cases
                Matrix mask1 = ga;
                Matrix mask2 = gb;
                bool edge = false;
                for (int i=0;i<n;i++) {
                    if (wina_col(i,0) == -1) {
                        edge = true;
                        mask1(0,i) = 0;
                        mask2(0,i) = 0;
                        wina_col(i,0) = 0;
                        winb_col(i,0) = 0;
                    }
                }
                if (edge) {
                    llip::normalize_matrix2(mask1,n,1);
                    llip::normalize_matrix2(mask2,n,1);
                }

                a(y,x) = dot(wina_col,mask1);
                b(y,x) = dot(winb_col,mask2);
            }
        }
        ****************************************************/

        Matrix ga(n,n);
        Matrix gb(n,n);
        // create 2 masks
        for (int r=-radius;r<=radius;r++){
            for (int c=-radius;c<=radius;c++){
                ga(r+radius,c+radius) = (1/(2*PI*sigma*sigma))*exp((c*c+r*r)/(-2*sigma*sigma));
                gb(r+radius,c+radius) = (1/(2*PI*sigma*k*sigma*k))*exp((c*c+r*r)/(-2*sigma*k*sigma*k));
            }
        }
        llip::normalize_matrix2(ga,n,n);
        llip::normalize_matrix2(gb,n,n);
        // Gaussian blur a and b
        for (int y=0;y<h;y++){
            for (int x=0;x<w;x++){
                Matrix win = llfd::get_mat_window(m,w,h,x,y,n,n,-1);

                // handle edge cases
                Matrix mask1 = ga;
                Matrix mask2 = gb;
                bool edge = false;
                for (int i=0;i<n;i++) {
                    for (int j=0;j<n;j++) {
                        if (win(j,i) == -1) {
                            edge = true;
                            mask1(j,i) = 0;
                            mask2(j,i) = 0;
                            win(j,i) = 0;
                        }
                    }
                }
                if (edge) {
                    llip::normalize_matrix2(mask1,n,n);
                    llip::normalize_matrix2(mask2,n,n);
                }

                for (int wr=0, r=y-radius; wr<n; wr++, r++){
                    for (int wc=0, c=x-radius; wc<n; wc++, c++){
                        if (c>=0 && c<w && r>=0 && r<h) {
                            a(y,x) += win(wr,wc) * mask1(wr,wc);
                            b(y,x) += win(wr,wc) * mask2(wr,wc);
                        }
                    }
                }
            }
        }//*/

        // return Difference of Gaussian convolution
        return b-a;
    }

    //
    //  gaussian_blur
    //
    Matrix gaussian_blur(Matrix& m, int w, int h, int n, double sigma) {
        Matrix out(h,w);
        int radius=n/2;
        Matrix mask(n,n);
        // create mask
        for (int r=-radius;r<=radius;r++)
            for (int c=-radius;c<=radius;c++)
                mask(r+radius,c+radius) = (1/(2*PI*sigma*sigma))*exp((c*c+r*r)/(-2*sigma*sigma));
        llip::normalize_matrix2(mask,n,n);
        // Gaussian blur out matrix
        for (int y=0;y<h;y++){
            for (int x=0;x<w;x++){
                Matrix win = llfd::get_mat_window(m,w,h,x,y,n,n,-1);

                // handle edge cases
                Matrix mask_tmp = mask;
                bool edge = false;
                for (int i=0;i<n;i++) {
                    for (int j=0;j<n;j++) {
                        if (win(j,i) == -1) {
                            edge = true;
                            mask_tmp(j,i) = 0;
                            win(j,i) = 0;
                        }
                    }
                }
                if (edge)
                    llip::normalize_matrix2(mask_tmp,n,n);

                for (int wr=0, r=y-radius; wr<n; wr++, r++)
                    for (int wc=0, c=x-radius; wc<n; wc++, c++)
                        if (c>=0 && c<w && r>=0 && r<h)
                            out(y,x) += win(wr,wc) * mask_tmp(wr,wc);
            }
        }

        // return blurred image matrix
        return out;
    }

    //
    //  get_scale_space
    //
    std::vector<std::vector<Matrix> > get_scale_space(Image& img, int N, double sigma) {

        std::cout << "Creating scale space: "<<OCTAVES<<" octaves and "<<SCALES<<" blur levels...\n";

        // helper variables
        int n = (N % 2 ? N : N+1); //always odd
        int rad = n / 2; //radius of window
        int pixels = img.height * img.width;

        // convert img to matrix
        Matrix I = img_to_mat(img);

        // create scale space
        Matrix tmp = I;
        double k = SQRT2;
        double power = 0;
        int tmpw = img.width;
        int tmph = img.height;
        std::vector<std::vector<Matrix> > scale_space;
        for (int r=0; r<OCTAVES; r++) {
            // add octave
            scale_space.push_back(std::vector<Matrix>());
            int subpow = power;
            for (int c=0; c<SCALES; c++) {
                // calculate sub k's
                double ksig = pow(k,subpow) * sigma;
                double kpow = pow(k,subpow+1);
                // blur image
                //scale_space[r].push_back(quick_DOG_blur(tmp,tmpw,tmph,n,ksig,kpow)); //old method
                scale_space[r].push_back(gaussian_blur(tmp,tmpw,tmph,n,ksig));// new method
                // increase sub power
                subpow++;
            }
            // shrink octave
            tmp = shrink_octave(tmp,tmpw,tmph);
            tmpw /= 2;
            tmph /= 2;
            // increment power
            power += 2;
        }
        std::cout << "Scale space created."<<std::endl;
        return scale_space;
    }

    //
    //  get_DOGs
    //
    std::vector<std::vector<Matrix> > get_DOGs(Image& img, std::vector<std::vector<Matrix> >& scalespace) {
        std::cout << "Calculating difference of Gaussians...\n";
        std::vector<std::vector<Matrix> > dogs(OCTAVES);
        for (int i=0; i<OCTAVES; i++) {
            int w = img.width  / int(pow(2,i));
            int h = img.height / int(pow(2,i));
            for (int j=0; j<SCALES-1; j++) {
                Matrix dog = scalespace[i][j+1] - scalespace[i][j];
                llip::normalize_matrix(dog,w,h,0,255);
                dogs[i].push_back(dog);
            }
        }
        std::cout << "Difference of Gaussians calculated."<<std::endl;
        return dogs;
    }

    //
    //  get_sift_pts
    //
    std::vector<std::vector<Matrix> > get_sift_pts(Image& img, std::vector<std::vector<Matrix> >& scalespace,
                                                   std::vector<std::vector<Matrix> >& dogs, int corner_thresh) {
        // find maxima & minima in each octave (using all blur levels except the first and last)
        //std::vector<std::vector<int> > interests; //old method
        std::vector<std::vector<Matrix> > interests; // new method
        for (int i=0; i<OCTAVES; i++) {
            std::cout<<"Calculating min/max for octave "<<i<<"..."<<std::endl;
            interests.push_back(std::vector<Matrix>(DOGCNT-2));
            int w = img.width  / int(pow(2,i));
            int h = img.height / int(pow(2,i));
            // loop through DoG's; skip first and last DoG
            for (int j=1; j<DOGCNT-1; j++) {
                // init interest matrix
                interests[i][j-1] = Matrix(h,w);
                // loop through image
                for (int r=0;r<h;r++){
                    for (int c=0;c<w;c++){
                        // mark current point if it's the largest/smallest
                        // among its neighbors in all 3 dimensions
                        double curr = dogs[i][j](r,c);
                        bool is_max = true, is_min = true, stop = false;
                        for (int x=-1; x<=1; x++) {
                            for (int y=-1; y<=1; y++) {
                                for (int z=-1; z<=1; z++) {
                                    if (x!=0 || y!=0 || z!=0) {
                                        int ir = r+y, ic = c+x;
                                        if (ir >= 0 && ir < h && ic >= 0 && ic < w) {// check out of bounds
                                            double nbor = dogs[i][j+z](ir,ic);
                                            if (nbor > curr) is_max = false;
                                            if (nbor < curr) is_min = false;
                                            if (!is_max && !is_min) stop = true;
                                        }
                                    }
                                    if (stop) break;
                                }
                                if (stop) break;
                            }
                            if (stop) break;
                        }
                        if ((is_max || is_min) && !(is_max && is_min)) {
                            /* OLD method
                            std::vector<int> pt(2);
                            pt[0] = c; pt[1] = r;
                            if (!sifter::exists_in(pt, interests))
                                interests.push_back(pt);
                            */

                            // NEW method
                            interests[i][j-1](r,c) = 255;
                        }
                    }
                }
            }
        }
        std::cout<<"Maxima/minima gathered."<<std::endl;

        //...use harris corner detection to ignore edges
        int sublvls = interests[0].size();
        for (int i=0; i<OCTAVES; i++) {
            std::cout<<"Filtering edges from octave "<<i<<"..."<<std::endl;
            int w = img.width  / int(pow(2,i));
            int h = img.height / int(pow(2,i));

            // use least blurred image in current octave to match corners
            Image tmp_img = mat_to_img(scalespace[i][0], w, h);
            llfd::detect_corners(tmp_img, 3, 1, 0.04, corner_thresh, false);
            Matrix corners = img_to_mat(tmp_img);

            for (int j=0; j<sublvls; j++)
                for (int r=0;r<h;r++)
                    for (int c=0;c<w;c++)
                        if (interests[i][j](r,c) > 0 && corners(r,c) == 0)
                            interests[i][j](r,c) = 0;

            std::cout<<"Octave "<<i<<"'s DoG images' edges filtered."<<std::endl;
        }

        return interests;
    }

    //
    //  get_keypoints
    //
    std::vector<std::vector<Keypoint> > get_keypoints(Image& img, std::vector<std::vector<Matrix> >& pts,
                                                      std::vector<std::vector<Matrix> >& scalespace) {

        std::cout<<"Calculating Key Points..."<<std::endl;
        std::vector<std::vector<Keypoint> > kpts(OCTAVES);
        // thru octaves
        for (int i=0; i<OCTAVES; i++) {
            int w = img.width  / int(pow(2,i));
            int h = img.height / int(pow(2,i));
            Matrix L = scalespace[i][0]; // use first blur level image
            // thru interest point matrices (ignore interest points on edges of image)
            for (int r=1;r<h-1;r++) {
                for (int c=1;c<w-1;c++) {
                    if (pts[i][0](r,c)>0 || pts[i][1](r,c)>0) { // use both DoG interest point matrices
                        // gather all neighbor gradient magnitude, direction into a histogram
                        std::vector<Keypoint> hist;
                        Matrix win = llfd::get_mat_window(L,w,h,c,r,5,5,L(r,c));
                        for (int wr=1;wr<4;wr++) {
                            for (int wc=1;wc<4;wc++) {
                                double xdif = win(wr,wc+1) - win(wr,wc-1);
                                double ydif = win(wr+1,wc) - win(wr-1,wc);
                                Keypoint kpt;
                                kpt.x = c+wc-2;
                                kpt.y = r+wr-2;
                                kpt.mag = sqrt(xdif*xdif + ydif*ydif);
                                kpt.dir = atan2(ydif, xdif);
                                kpt.bin = round(kpt.dir / 10); //36 bins (indexes 0 to 35)
                                hist.push_back(kpt);
                            }
                        }
                        // find max mag kpt in hist
                        double max_mag = -99999999;
                        for (int j=0; j<hist.size(); j++)
                            if (max_mag < hist[j].mag)
                                max_mag = hist[j].mag;
                        // add a Keypoint for every point with mag > 80% of max_mag
                        int min_dir = max_mag * 0.8;
                        for (int j=0; j<hist.size(); j++)
                            if (hist[j].dir >= min_dir)
                                kpts[i].push_back(hist[j]);
                    }
                }
            }
        }
        std::cout<<"Key Points Calculated."<<std::endl;
        return kpts;
    }

}//end sifter name space
