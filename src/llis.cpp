////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       llis.cpp
//  Date:           04/21/2018
//  Description:    Main implementation for LLIS (low level image segmentation) class.
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
#include "llis.h"
#include "llip.h"
#include "matrix.h"

#include "math.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

// for quicker PI calculations
#define PI 3.14159265

namespace llis {

    //
    //  img_to_mat
    //
    Matrix img_to_mat(Image& img) {
        Matrix m(img.height, img.width);
        int len = img.data.size();
        for (int i=0; i<len; i++) {
            m(i/img.width, i%img.width) = img.data[i];
        }
        return m;
    }

    //
    //  mat_to_img
    //
    Image mat_to_img(Matrix& m, int w, int h) {
        Image img;
        img.type = 0;
        img.width = w;
        img.height = h;
        img.max_val = 255;
        for (int r=0; r<h; r++) {
            for (int c=0; c<w; c++) {
                img.data.push_back(m(r,c));
            }
        }
        return img;
    }

    //
    //  write_descriptors
    //
    void write_descriptors(std::vector<Descriptor>& props, std::string fname, int thresh) {

        // open
        std::ofstream f;
        f.open (fname);

        // write
        f << "Given:\n";
        f << "Smoothed Histogram Threshold: " << thresh << "\n\n";
        f << "--- BEGIN COMPONENT DESCRIPTION LIST ----------------\n\n";
        for (int i=0; i<props.size(); i++) {
            f << "~~~ Component #" << i << ":\n";
            f << "Intensity: " << props[i].intensity << "\n";
            f << "Area: " << props[i].area << "\n";
            f << "Centroid: " << props[i].centroid[0] << "," << props[i].centroid[1] << "\n";
            f << "Central Moment: " << props[i].central_moments[0] << "," << props[i].central_moments[1] << "," << props[i].central_moments[2] << "\n";
            f << "Scale Invariant Central Moment: " << props[i].si_central_moments[0] << "," << props[i].si_central_moments[1] << "," << props[i].si_central_moments[2] << "\n";
            f << "Perimeter: " << props[i].perimeter << "\n";
            f << "Compactness: " << props[i].compactness << "\n";
            f << "Elongation: " << props[i].elongation << "\n\n";
        }
        f << "--- END COMPONENT DESCRIPTION LIST ----------------\n";

        // close
        f.close();

        // notify
        if (f.bad())
            std::cout<<"ERROR: Segment descriptors could not be written properly."<<std::endl;
        else
            std::cout<<"\""<<fname<<"\" saved successfully."<<std::endl;
    }

    //
    //  guess_valley
    //
    int guess_valley(std::vector<double>& hist, int thresh, int prev_thresh) {
        // get means of both sides
        int new_thresh=0;
        double left_sum=0, left_weighted_sum=0, left_mean=0;
        double right_sum=0, right_weighted_sum=0, right_mean=0;
        for (int i=0; i<thresh; i++){
            left_sum          += hist[i];
            left_weighted_sum += hist[i] * i;
        }
        for (int i=thresh; i<hist.size(); i++){
            right_sum          += hist[i];
            right_weighted_sum += hist[i] * i;
        }
        left_mean = left_weighted_sum / left_sum;
        right_mean = right_weighted_sum / right_sum;

        // recurse with the new threshold if it's different, else, valley found
        new_thresh = round((left_mean+right_mean)/2);
        if (new_thresh != thresh && new_thresh != prev_thresh)
            return guess_valley(hist, new_thresh, thresh);
        else
            return thresh;
    }

    //
    //  scan_equivs
    //
    void scan_equivs(std::vector< std::vector<int> >& eq, std::vector<int>& eq2, int curr, int lbl) {
        eq2[curr] = lbl;
        for (int i=0;i<eq[curr].size();i++)
            if (!eq2[eq[curr][i]-1])
                scan_equivs(eq,eq2,eq[curr][i]-1,lbl);

        // look for curr in the rest of the class items
        for (int j=0;j<eq.size();j++)
            if (j != curr)
                if (!eq2[j])
                    if (std::find(eq[j].begin(), eq[j].end(), curr+1) != eq[j].end())
                        scan_equivs(eq, eq2, j, lbl);
    }

    //
    //  simplify_equivs
    //
    std::vector<int> simplify_equivs(std::vector< std::vector<int> >& eq) {
        int label = 1;
        std::vector<int> new_eq(eq.size(), 0);
        for (int i=0;i<eq.size();i++) {
            // skip if this index has already been calculated (not 0)
            if (new_eq[i])
                continue;

            // recurse through current class items
            scan_equivs(eq, new_eq, i, label);

            // next label
            label++;
        }
        return new_eq;
    }

    //
    //  label_connected_components
    //
    Matrix label_connected_components(Matrix& m, int w, int h) {
        // 1st pass (4-neighbor raster scan method)
        Matrix lbls(h,w);
        std::vector< std::vector<int> > equivs;
        int lbl = 0;
        for (int r=0; r<h; r++)
            for (int c=0; c<w; c++)
                if (m(r,c) > 0) {
                    // get top and left values, respectively
                    int duo[2] = {-1,-1}; //-1 means DNE
                    if (r>0) duo[0] = m(r-1,c); //top
                    if (c>0) duo[1] = m(r,c-1); //left

                    // mimic top
                    if (duo[0] > 0)
                        lbls(r,c) = lbls(r-1,c);

                    // mimic left
                    if (duo[1] > 0)
                        lbls(r,c) = lbls(r,c-1);

                    // mimic left, add top to equiv of left
                    if (duo[0] > 0 && duo[1] > 0) {
                        int top = lbls(r-1,c);
                        int left = lbls(r,c-1);
                        if (top != left) {
                            std::vector<int>& v = equivs[left-1];
                            if (std::find(v.begin(), v.end(), top) == v.end())
                                v.push_back(top);
                        }
                    }

                    // new label, add empty equiv
                    if (duo[0] <= 0 && duo[1] <= 0) {
                        lbl++;
                        lbls(r,c) = lbl;
                        equivs.push_back(std::vector<int>());
                    }
                }
        // simplify equivalence class list into a direct 1:1 mapping
        std::vector<int> codes = simplify_equivs(equivs);

        // 2nd pass (replace equivalent classes)
        for (int r=0; r<h; r++)
            for (int c=0; c<w; c++)
                if (m(r,c))
                    lbls(r,c) = codes[lbls(r,c)-1];
                else
                    lbls(r,c) = 0;

        return lbls;
    }

    //
    //  neighbor_diffs
    //
    int neighbor_diffs(Matrix& m, int w, int h, int c, int r) {
        int edges=0;
        if (c+1 < w && m(r,c+1) != m(r,c))
            edges++;
        if (c-1 >= 0 && m(r,c-1) != m(r,c))
            edges++;
        if (r+1 < h && m(r+1,c) != m(r,c))
            edges++;
        if (r-1 >= 0 && m(r-1,c) != m(r,c))
            edges++;
        return edges;
    }

    //
    //  segment_image
    //
    std::vector<Image> segment_image(Image& img, std::string desc_txt_fname) {

        // helper vars
        int pixels = img.width * img.height;

        // return var
        std::vector<Image> outputs;

        //
        //  Part 1.
        //

        // 1. Get a matrix of the image in the space: x=color, y=(unsorted), value=coordinate(x,y)
        std::vector< std::vector< std::pair<int,int> > > clr_coords(256);
        Matrix img_mat = img_to_mat(img);
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                std::pair<int, int> coord(c,r);
                clr_coords[(int)img_mat(r,c)].push_back(coord);
            }
        }

        // 2. Get image histogram vector from the matrix
        std::vector<double> hist(256);
        for (int i=0; i<clr_coords.size(); i++) {
            hist[i] = clr_coords[i].size();
        }

        // output an image representing the histogram
        Matrix out_hist_mat(100,256);
        std::vector<double> hist_scaled = hist;
        for (int c=0; c<256; c++) {
            // scale to range of 0-99
            double min_val=999999999, max_val=-999999999;
            for (int i=0; i<hist.size(); i++) {
                    if (min_val > hist[i])
                        min_val = hist[i];
                    if (max_val < hist[i])
                        max_val = hist[i];
            }
            for (int i=0; i<hist.size(); i++)
                hist_scaled[i] = round(99*(hist[i]-min_val)/(max_val-min_val));
            // make bar chart below
            for (int r=100-1; r>=0; r--) {
                if (100-r-1 <= (int)hist_scaled[c])
                    out_hist_mat(r,c) = 255;
                else
                    out_hist_mat(r,c) = 0;
            }
        }
        outputs.push_back(mat_to_img(out_hist_mat, 256, 100));

        // 3. Smooth histogram using this mask:
        double mask[5] = {1.0/9.0, 2.0/9.0, 3.0/9.0, 2.0/9.0, 1.0/9.0};
        std::vector<double> hist_smooth(256);
        for (int i=0; i<hist.size(); i++) {
            int oob_val; // out of bounds guess value for hist
            if (i < 2)
                oob_val = hist[0];
            else if (i > hist.size()-3)
                oob_val = hist[hist.size()-1];

            for (int j=i-2, k=0; k<5; j++, k++) {
                if (j>=0 && j<hist.size())
                    hist_smooth[i] += mask[k]*hist[j];
                else
                    hist_smooth[i] += mask[k]*oob_val;
            }
        }

        // output an image representing the smoothed histogram
        for (int c=0; c<256; c++) {
            // scale to range of 0-99
            double min_val=999999999, max_val=-999999999;
            for (int i=0; i<hist.size(); i++) {
                    if (min_val > hist[i])
                        min_val = hist[i];
                    if (max_val < hist[i])
                        max_val = hist[i];
            }
            for (int i=0; i<hist.size(); i++)
                hist_scaled[i] = round(99*(hist_smooth[i]-min_val)/(max_val-min_val));
            // make bar chart below
            for (int r=100-1; r>=0; r--) {
                if (100-r-1 <= (int)hist_scaled[c])
                    out_hist_mat(r,c) = 255;
                else
                    out_hist_mat(r,c) = 0;
            }
        }
        Image out_hist_img = mat_to_img(out_hist_mat, 256, 100);
        pnmio mngr;
        mngr.convert_gray_to_rgb(out_hist_img);

        // 4. Estimate the threshold index (valley of histogram)
        //    Initial guess should be the avg intensity of the image
        int sum = 0;
        for (int i=0; i<pixels; i++) sum += img.data[i];
        int threshold = guess_valley(hist_smooth, round((double)sum/pixels));
        std::cout<<"\nThreshold guess (at histogram valley): "<<threshold<<"\n"<<std::endl;

        // draw red line at threshold
        for (int i=0;i<out_hist_img.data.size();i++) {
            int x = i%(out_hist_img.width*3)/3;
            if (x == threshold) {
                if (i%3 == 0)
                    out_hist_img.data[i] = 255;
                else
                    out_hist_img.data[i] = 0;
            }
        }
        outputs.push_back(out_hist_img);

        // 5. Alter the colors in the matrix using the threshold:
        //    - everything before thresh (0 to N-1) will be black;
        //    - everything after thresh (N to 255) will be white
        for (int r=0; r<img.height; r++) {
            for (int c=0; c<img.width; c++) {
                if (img_mat(r,c) < threshold)
                    img_mat(r,c) = 0;
                else
                    img_mat(r,c) = 255;
            }
        }
        outputs.push_back(mat_to_img(img_mat, img.width, img.height));

        //
        //  Part 2. (img_mat is now binary)
        //

        // label connected components from 1 to N
        Matrix lbl_mat = label_connected_components(img_mat, img.width, img.height);
        Matrix lbl_mat_norm = lbl_mat;
        llip::normalize_matrix(lbl_mat_norm, img.width, img.height, 50, 255);
        for (int r=0; r<img.height; r++)
            for (int c=0; c<img.width; c++)
                if (lbl_mat_norm(r,c) <= 50)
                    lbl_mat_norm(r,c) = 0;
        // round down normalized labels
        for (int r=0; r<img.height; r++)
            for (int c=0; c<img.width; c++)
                lbl_mat_norm(r,c) = ceil(lbl_mat_norm(r,c));
        outputs.push_back(mat_to_img(lbl_mat_norm, img.width, img.height));

        //
        //  Part 3. (calculate properties using lbl_mat_norm)
        //

        // NOTE: Properties for each component will be saved in a text file

        // get sorted list of labels
        std::vector<int> comps;
        for (int r=0; r<img.height; r++)
            for (int c=0; c<img.width; c++)
                if (std::find(comps.begin(), comps.end(), lbl_mat_norm(r,c)) == comps.end())
                    comps.push_back(lbl_mat_norm(r,c));
        std::sort(comps.begin(), comps.end());

        // list of properties for each component
        std::vector<Descriptor> props(comps.size());

        // ID Property
        for (int i=0;i<props.size();i++)
            props[i].intensity = comps[i];

        std::cout<<"Identification Property: Component Intensity:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].intensity<<"\n";
        std::cout<<std::endl;

        // Property #1: Area
        for (int i=0, area=0; i<props.size(); i++, area=0) {
            for (int r=0; r<img.height; r++)
                for (int c=0; c<img.width; c++)
                    if (lbl_mat_norm(r,c) == comps[i])
                        area++;
            props[i].area = area;
        }

        std::cout<<"Property #1: Areas:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].area<<"\n";
        std::cout<<std::endl;

        // Property #2: Centroid
        for (int i=0, cx=0, cy=0; i<props.size(); i++, cx=0, cy=0) {
            for (int r=0; r<img.height; r++)
                for (int c=0; c<img.width; c++)
                    if (lbl_mat_norm(r,c) == comps[i]) {
                        cx += c;
                        cy += r;
                    }
            props[i].centroid[0] = (double)cx / props[i].area;
            props[i].centroid[1] = (double)cy / props[i].area;
        }

        std::cout<<"Property #2: Moments:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].centroid[0]<<","<<props[i].centroid[1]<<"\n";
        std::cout<<std::endl;

        // Property #3: Central Moment
        for (int i=0, cxx=0, cyy=0, cxy=0; i<props.size(); i++, cxx=0, cyy=0, cxy=0) {
            for (int r=0; r<img.height; r++)
                for (int c=0; c<img.width; c++)
                    if (lbl_mat_norm(r,c) == comps[i]) {
                        cxx += c*c;
                        cyy += r*r;
                        cxy += c*r;
                    }
            props[i].central_moments[0] = (double)cxx / props[i].area;
            props[i].central_moments[1] = (double)cxy / props[i].area;
            props[i].central_moments[2] = (double)cyy / props[i].area;
        }

        std::cout<<"Property #3: Central Moments:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].central_moments[0]<<","<<props[i].central_moments[1]<<","<<props[i].central_moments[2]<<"\n";
        std::cout<<std::endl;

        // Mediating Property (to make elongation calculation simpler)
        for (int i=0; i<props.size(); i++) {
            props[i].si_central_moments[0] = props[i].central_moments[0] / (props[i].area * props[i].area);
            props[i].si_central_moments[1] = props[i].central_moments[1] / (props[i].area * props[i].area);
            props[i].si_central_moments[2] = props[i].central_moments[2] / (props[i].area * props[i].area);
        }

        std::cout<<"Mediating Property: Scale Invariant Central Moments:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].si_central_moments[0]<<","<<props[i].si_central_moments[1]<<","<<props[i].si_central_moments[2]<<"\n";
        std::cout<<std::endl;

        // Property #4: Perimeter
        for (int i=0, edges=0; i<props.size(); i++, edges=0) {
            for (int r=0; r<img.height; r++)
                for (int c=0; c<img.width; c++)
                    if (lbl_mat_norm(r,c) == comps[i])
                        edges += neighbor_diffs(lbl_mat_norm,img.width,img.height,c,r);
            props[i].perimeter = edges;
        }

        std::cout<<"Property #4: Perimeter:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].perimeter<<"\n";
        std::cout<<std::endl;

        // Property #5: Compactness
        for (int i=0; i<props.size(); i++)
            props[i].compactness = (props[i].perimeter * props[i].perimeter) / 4*PI*props[i].area;

        std::cout<<"Property #5: Compactness:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].compactness<<"\n";
        std::cout<<std::endl;

        // Property #6: Elongation
        for (int i=0;i<props.size();i++) {
            double si_xx = props[i].si_central_moments[0];
            double si_xy = props[i].si_central_moments[1];
            double si_yy = props[i].si_central_moments[2];
            double thing = sqrt(((si_xx-si_yy)*(si_xx-si_yy))+4*(si_xy*si_xy));
            props[i].elongation = sqrt(thing/(si_xx+si_yy+thing));
        }

        std::cout<<"Property #6: Elongation:\n";
        for (int i=0;i<props.size();i++)
            std::cout<<i<<": "<<props[i].elongation<<"\n";
        std::cout<<std::endl;

        write_descriptors(props, desc_txt_fname, threshold);

        return outputs;
    }

}
