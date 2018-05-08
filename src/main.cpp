/*
*   Author: Ibrahim Sardar
*   03/31/2018
*/

#include "pnmio.h"  /* REQUIRES C++11 */
#include "llfd.h"   /* LOW LEVEL FEATURE DETECTION FUNCTIONS HERE */
#include "llip.h"   /* LOW LEVEL PROCESSING FUNCTIONS HERE */
#include "sifter.h" /* SIFT (Scale Invariant Feature Transform) FUCTIONS HERE */
#include "depthmapper.h" /* DEPTH MAP RELATED FUNCTIONS HERE */


#include <iostream>
#include <sstream>
#include <exception>
#include <string>
#include "math.h"

// helpers
void say(std::string s) {std::cout<<s<<std::endl;}
void save(Image img, std::string original, pnmio image_manager) {
    // ask for new filename for processed image
    say("Enter a name for the processed version of " + original);
    std::string newfilename;
    std::getline(std::cin, newfilename);
    try {
        image_manager.store_image(newfilename, img);
        say("\""+newfilename+"\" saved successfully.");
    } catch (std::exception &e) {
        say("\""+newfilename+"\" failed to save.");
    }
}
void quicksave(Image img, std::string name, pnmio mngr) {
    try {
        mngr.store_image(name, img);
        say("\""+name+"\" saved successfully.");
    } catch (std::exception &e) {
        say("\""+name+"\" failed to save.");
    }
}

int main() {

    // print instructions
    say("Welcome to the Stereo Image Depth Map Generator."
        "This program expects at least 2 stereo images."
        "Type 'q' or 'quit' to quit.");

    // menu
    while (true) {

        // ask for left stereo image
        say("\nEnter a filename (in this exe's directory) of the LEFT stereo image:");
        std::string filename_left;
        std::getline(std::cin, filename_left);

        // check quit flag
        if (filename_left == "q" || filename_left == "quit")
            break;

        // ask for right stereo image
        say("\nEnter a filename (in this exe's directory) of the RIGHT stereo image:");
        std::string filename_right;
        std::getline(std::cin, filename_right);

        // check quit flag
        if (filename_right == "q" || filename_right == "quit")
            break;

        //junk
        //filename_left = "umbl.pnm";
        //filename_right = "umbr.pnm";

        // read image
        pnmio manager;
        Image input_left;
        Image input_right;
        try {
            input_left = manager.load_image(filename_left);
            input_right = manager.load_image(filename_right);
        } catch (std::exception &e) {
            say(e.what()); break;
        };

        // notify the user the image is being converted to gray-scale if it is RGB; warn/notify otherwise
        if (input_left.type == 1) {
            manager.convert_rgb_to_gray(input_left);
            say("(LEFT) RGB image detected; converted to gray-scale for compatibility.");
        } else if (input_left.type == 0) {
            say("(LEFT) Gray-scale image detected.");
        } else {
            say("(LEFT) Unknown image-type detected; some features may behave unexpectedly.");
        }
        if (input_right.type == 1) {
            manager.convert_rgb_to_gray(input_right);
            say("(RIGHT) RGB image detected; converted to gray-scale for compatibility.");
        } else if (input_right.type == 0) {
            say("(RIGHT) Gray-scale image detected.");
        } else {
            say("(RIGHT) Unknown image-type detected; some features may behave unexpectedly.");
        }

        // ask for processing type
        say("\nEnter a processing type: (type the index)\n"
            "0: Depth Map (actually it's disparity at the moment)\n"
            "1: SIFT point detection (multiple outputs) - VERY slow on large images!");
        std::string ptype;
        std::getline(std::cin, ptype);
        say("");

        // check quit flag
        if (ptype == "quit" || ptype == "QUIT")
            break;

        // process image and save it as a new image file
        switch (std::stoi(ptype)) {
        case 0: {

            // calculate & save disparity depth map (assumes left and right are rectified)
            Matrix depth_mat = depthmapper::get_disparity(input_left, input_right, 5, 50);
            Image depthmap = sifter::mat_to_img(depth_mat, input_left.width, input_left.height);
            quicksave(depthmap, "depthmap.pnm", manager);

            break;
        } case 1: {
            // save original data
            std::vector<unsigned int> original_left = input_left.data;
            std::vector<unsigned int> original_right = input_right.data;

            // create a scale space (4 levels of 5 sub levels)
            std::vector<std::vector<Matrix> > ss_left = sifter::get_scale_space(input_left,5,0.707107);
            std::vector<std::vector<Matrix> > ss_right = sifter::get_scale_space(input_right,5,0.707107);

            // get difference of Gaussian images (reduces number of sub levels to 4)
            std::vector<std::vector<Matrix> > DOGs_left = sifter::get_DOGs(input_left, ss_left);
            std::vector<std::vector<Matrix> > DOGs_right = sifter::get_DOGs(input_right, ss_right);

            // gather interest points (reduces number of sub levels to 2)
            std::vector<std::vector<Matrix> > interests_left = sifter::get_sift_pts(input_left, ss_left, DOGs_left, 0);
            std::vector<std::vector<Matrix> > interests_right = sifter::get_sift_pts(input_right, ss_right, DOGs_right, 0);

            // gather key points for each octave (reduces number of sub levels to 1...so 5 octaves, 1 set of key points)
            std::vector<std::vector<Keypoint> > KPs_left = sifter::get_keypoints(input_left, interests_left, ss_left);
            std::vector<std::vector<Keypoint> > KPs_right = sifter::get_keypoints(input_right, interests_right, ss_right);

            // calculate feature vectors for each key point
            ///...ran out of time :(

            //*// save octave images (left)
            int OCTAVES = 4, SCALES = 5;
            for (int i=0; i<OCTAVES; i++) {
                for (int j=0; j<SCALES; j++) {
                    Image im;
                    im.max_val=input_left.max_val;
                    im.type=input_left.type;
                    im.width = input_left.width / int((pow(2,i)));
                    im.height = input_left.height / int((pow(2,i)));
                    llip::normalize_matrix(ss_left[i][j],im.width,im.height,0,255);
                    for (int r=0;r<im.height;r++){
                        for (int c=0;c<im.width;c++){
                            im.data.push_back(int(round(ss_left[i][j](r,c))));
                        }
                    }
                    std::stringstream s;
                    s << "_ss_oct" << i << "_blur" << j << ".pnm";
                    quicksave(im, s.str(), manager);
                }
            }//*/

            //*// save DoG images (left)
            OCTAVES = 4, SCALES = 4;
            for (int i=0; i<OCTAVES; i++) {
                for (int j=0; j<SCALES; j++) {
                    Image im;
                    im.max_val=input_left.max_val;
                    im.type=input_left.type;
                    im.width = input_left.width / int((pow(2,i)));
                    im.height = input_left.height / int((pow(2,i)));
                    //llip::normalize_matrix(DOGs_left[i][j],im.width,im.height,0,255);
                    for (int r=0;r<im.height;r++){
                        for (int c=0;c<im.width;c++){
                            im.data.push_back(int(round(DOGs_left[i][j](r,c))));
                        }
                    }
                    std::stringstream s;
                    s << "_dg_oct" << i << "_dog" << j << ".pnm";
                    quicksave(im, s.str(), manager);
                }
            }//*/

            //*// save interest point images (left)
            OCTAVES = 4, SCALES = 2;
            for (int i=0; i<OCTAVES; i++) {
                for (int j=0; j<SCALES; j++) {
                    Image im;
                    im.max_val=input_left.max_val;
                    im.type=input_left.type;
                    im.width = input_left.width / int((pow(2,i)));
                    im.height = input_left.height / int((pow(2,i)));
                    for (int r=0;r<im.height;r++){
                        for (int c=0;c<im.width;c++){
                            im.data.push_back(int(round(interests_left[i][j](r,c))));
                        }
                    }
                    std::stringstream s;
                    s << "_ip_oct" << i << "_dog" << j << ".pnm";
                    quicksave(im, s.str(), manager);
                }
            }//*/

            /*
            // reset the data
            pnmio::convert_rgb_to_gray(input_left);
            pnmio::convert_rgb_to_gray(input_right);
            input_left.data = original_left;
            input_right.data = original_right;

            // gather sift features from both images
            std::vector<sifter::Feature> sifts_left = sifter::get_sift_features(input_left, interests_left);
            std::vector<sifter::Feature> sifts_right = sifter::get_sift_features(input_right, interests_right);

            // gather matching sift points from both images
            std::vector<sifter::Feature> sifts_matched = sifter::get_matching_sifts(sifts_left, sifts_right);
*/
            // save depth map image
            //save(input, filename, manager);
            break;
        } default:
            say("Input could not be understood.");
            break;
        }//end switch
    }

    // end program
    return 0;
}
