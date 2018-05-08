////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       pnmio.cpp
//  Date:           02/04/2018
//  Description:    Main implementation for PNM IO (read/write) class.
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

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <sstream>

//
//  bin_to_int
//
unsigned int pnmio::bin_to_uint(unsigned char c) {

    // place binary character byte into a bitset
    std::bitset<8> binary(c);
    return (unsigned int)(binary.to_ulong());
}

unsigned char pnmio::uint_to_bin(unsigned int i)
{
	// place binary character byte into a bitset
	std::bitset<8> binary(i);
	return (unsigned char)(binary.to_ulong());
}

//
//  char_to_int
//
unsigned int pnmio::str_to_uint(std::string s) {

    // use a string-stream to get appropriate int value
    std::stringstream ss(s);
    unsigned int i;
    ss >> i;
    return i;
}

std::string pnmio::uint_to_str(unsigned int i) {

	// use a string-stream to get appropriate string value
	std::stringstream ss;
	ss << i;
	return ss.str();
}

std::string pnmio::str_until_space(std::ifstream &ifs) {

    // only stop getting chars if eof, newline, or space is encountered
    std::string out = "";
    while (ifs.peek() && ifs.peek() != '\n' && ifs.peek() != ' ') {
        out.push_back(ifs.get());
    }
    return out;
}

//
//  load_image
//
Image pnmio::load_image(std::string fp, bool report, bool has_comment) {

    // open binary mode
    std::ifstream ifs(fp.c_str(), std::ios::binary|std::ios::in);

    // abort if failed to open
    if (!ifs) throw read_failure_exception(fp);

    // create empty Image structure
    Image img;

    // create empty (temporary) magic number array
    unsigned char mn[2];

    // create empty (temporary) dimensions array
    unsigned int dim[2];

    // create nil (temporary) max rgb magnitude
    unsigned int max_value = 0;

    // create a string that is used to ignore lines
    std::string ignore_text = "";

    // get magick number
    mn[0] = ifs.get();
    mn[1] = ifs.get();

    // assert format: 'P5' or 'P6'
    if (mn[0] != 'P' && (mn[1] == '5' || mn[1] == '6')) // Gray-scale or RGB only
        throw read_failure_exception(fp);

    // next line
    ifs.get();

    // skip comment lines
    while (ifs.peek() == '#')
        std::getline(ifs, ignore_text);

    // get image dimensions
    dim[0] = str_to_uint(str_until_space(ifs)); //width
    ifs.get(); // skip space character
    dim[1] = str_to_uint(str_until_space(ifs)); // height

    // next line
    ifs.get();

    // get max rgb/gray-scale magnitude
    max_value = str_to_uint(str_until_space(ifs));

    // next line
    ifs.get();

    // calculate amount of expected data values
    unsigned int amt = dim[0] * dim[1];
    if (mn[1] == '6') amt *= 3; // 3x more data for RGB than Gray-scale

    // gather data
    int cntr = 0;
    while (cntr < amt) {
        cntr++;
        img.data.push_back(bin_to_uint(ifs.get())); // gather data
    }

    // close the filestream
    ifs.close();

	// exception if something went wrong
	if (ifs.bad()) throw read_failure_exception(fp);

    // print read report
    if (report) {
        std::cout << "\n--- PNM Read Report:" << '\n';
        std::cout << "\tType: " << mn[0] << mn[1] << '\n';
        std::cout << "\tWidth: " << dim[0] << '\n';
        std::cout << "\tHeight: " << dim[1] << '\n';
        std::cout << "\tMax Data Point Value: " << max_value << '\n';
        std::cout << "\tData Point Value Matrix:";
        int col_elems = (mn[1] == '6') ? dim[0]*3 : dim[0];
        for (unsigned int i=0; i<amt; i++) {
            int row = i / col_elems; // first row is = to 0
            int col = i % col_elems; // first column is = to 0

            if (col == 0) std::cout << '\n' << '\t' << '\t';
            std::cout << img.data[i];
            if (col != col_elems - 1) std::cout << ',';
        }
        std::cout << "\n--------------------" << std::endl;
    }

    // construct and return the Image structure
    img.type = (mn[1] == '6') ? 1 : 0; //type
    img.width = dim[0]; //width
    img.height = dim[1]; //height
    img.max_val = max_value; //max data point value
    return img;
}

//
// store_image
//
void pnmio::store_image(std::string fp, Image & img, bool report) {

	// open binary mode
	std::ofstream ofs(fp.c_str(), std::ios::binary | std::ios::out);

	// abort if failed to open
	if (!ofs) throw write_failure_exception(fp);

	// write magic number
	ofs.put('P');
	if (img.type) ofs.put('6');
	else ofs.put('5');
	ofs.put('\n');

	// write comment
	std::string comment = "# Written by pnmio.cpp PNM Writer";
	ofs << comment;
	ofs.put('\n');

	// write width + space + height
	ofs << uint_to_str(img.width);
	ofs.put(' ');
	ofs << uint_to_str(img.height);
	ofs.put('\n');

	// write maximum data point value
	ofs << uint_to_str(img.max_val);
	ofs.put('\n');

	// write data point matrix
	for (unsigned int i = 0; i < img.data.size(); i++)
		ofs.put(uint_to_bin(img.data[i]));

	// print write report
	if (report) {
		std::cout << "\n--- PNM Write Report:" << '\n';
		std::cout << "\tType: P" << (img.type ? '6' : '5') << '\n';
		std::cout << "\tWidth: " << img.width << '\n';
		std::cout << "\tHeight: " << img.height << '\n';
		std::cout << "\tMax Data Point Value: " << img.max_val << '\n';
		std::cout << "\tData Point Value Matrix:";
		int col_elems = (img.type ? img.width * 3 : img.width);
		for (unsigned int i = 0; i<img.data.size(); i++) {
			int row = i / col_elems; // first row is = to 0
			int col = i % col_elems; // first column is = to 0

			if (col == 0) std::cout << '\n' << '\t' << '\t';
			std::cout << img.data[i];
			if (col != col_elems - 1) std::cout << ',';
		}
		std::cout << "\n--------------------" << std::endl;
	}

	// close the filestream
	ofs.close();

	// exception if something went wrong
	if (ofs.bad()) throw write_failure_exception(fp);
}

//
//  convert_gray_to_rgb
//
void pnmio::convert_gray_to_rgb(Image & img) {

    // ignore images which are already RGB
    if (img.type == 1)
        return;

    // helper variables
    const unsigned int pixels = img.width * img.height;

    // convert gray to rgb
    std::vector<unsigned int> rgb_data;
    for (unsigned int i=0; i<pixels; i++) {
        rgb_data.push_back(img.data[i]);
        rgb_data.push_back(img.data[i]);
        rgb_data.push_back(img.data[i]);
    }

    // set new data
    img.data = rgb_data;
    img.type = 1;
}

//
//  convert_rgb_to_gray
//
void pnmio::convert_rgb_to_gray(Image & img) {

    // ignore images which are already Gray-scale
    if (img.type == 0)
        return;

    // helper variables
    const unsigned int values = img.width * img.height * 3;

    // convert rgb to gray
    std::vector<unsigned int> gray_data;
    for (unsigned int i=0; i<values; i+=3) {
        gray_data.push_back((img.data[i]+img.data[i]+img.data[i]) / 3); // avg of 3 values
    }

    // set new data
    img.data = gray_data;
    img.type = 0;
}
