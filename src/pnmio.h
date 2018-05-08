////////////////////////////////////////////////////////////////////////////////
//
//  Author:         Ibrahim Sardar
//  Class:          CSCI 557
//  Filename:       pnmio.h
//  Date:           02/04/2018
//  Description:    Header for PNM IO (read/write) class.
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

#ifndef _PNMIO_H_
#define _PNMIO_H_

#include <string>
#include <exception>
#include <fstream>
#include <vector>

/**
*   @struct Image
*
*   Image structure that represents data held in a pnm image file.
*/
struct Image {
    unsigned int type; // an integer value specifying the image type (Gray = 0 and RGB = 1)
    unsigned int width; // an integer value specifying the width of the image
    unsigned int height; // an integer value specifying the height of the image
    unsigned int max_val; // an integer value specifying the maximum data point value used
    std::vector<unsigned int> data; // an array of integers representing data point values (row-major format)
};

/**
*   @class pnmio
*
*   Assists with reading and writing of pnm image files.
*   Note: use of 'noexcept' requires C++11
*/
class pnmio {

    /**
     *  @class read_failure_exception
     *
     *  Exception thrown to indicate a failure
     *  while reading the PNM file.
     */
    class read_failure_exception : public std::exception {
    public:
        /// Initialization constructor
        read_failure_exception(std::string fp) : std::exception(), filepath_(fp) {}

        /// message of exception
        const char* what() const noexcept {
            std::string str = "pnmio failed to read from \"" + filepath_ + "\"";
            return str.c_str();
        }
    private:
        /// file directory + file name trying to read from
        std::string filepath_;
    };

    /**
     *  @class write_failure_exception
     *
     *  Exception thrown to indicate a failure
     *  while writing the PNM file.
     */
    class write_failure_exception : public std::exception {
    public:
        /// Initialization constructor
        write_failure_exception(std::string fp) : std::exception(), filepath_(fp) {}

        /// message of exception
        const char* what() const noexcept {
            std::string str = "pnmio failed to write to \"" + filepath_ + "\"";
            return str.c_str();
        }
    private:
        /// file directory + file name trying to write to
        std::string filepath_;
    };

    public:
        /**
        *   Reads in a PNM image file and outputs an Image structure.
        *
        *   @param          string                  PNM image path + \filename (image to be loaded)
        *   @param          boolean (optional)      indicates that a report should be written after the read
        *   @param          boolean (optional)      indicates a comment on the second line
        *   @return         Image                   PNM image structure
        */
        Image load_image(std::string, bool=false, bool=true);

        /**
        *   Writes a PNM image file using data from an Image structure.
        *
        *   @param          string                  PNM image path + \filename
        *   @param          Image &                 PNM image structure to be saved
        *   @param          boolean (optional)      indicates that a report should be written after the write
        */
        void store_image(std::string, Image &, bool=false);

        /**
        *   Converts a Gray-scale Image structure into an RGB Image structure
        *
        *   @param          Image &                 PNM image structure to be manipulated
        */
        void convert_gray_to_rgb(Image &);

        /**
        *   Converts an RGB Image structure into a Gray-scale Image structure
        *
        *   @param          Image &                 PNM image structure to be manipulated
        */
        void convert_rgb_to_gray(Image &);

    private:
        /**
        *   Converts a 1 byte binary value stored in a character into a 4 byte uint.
        *
        *   @param          char                    binary value (from PNM file)
        *   @return         unsigned int            converted integer value
        */
        unsigned int bin_to_uint(unsigned char);

		/**
		*   Converts a 4 byte uint into a 1 byte binary value stored in a character.
		*
		*   @param          unsigned int            integer value
		*   @return         char                    converted binary value (for PNM file)
		*/
		unsigned char uint_to_bin(unsigned int);

        /**
        *   Converts a string into its apparent 4 byte uint.
        *
        *   @param          string                  string value
        *   @return         unsigned int            converted integer value
        */
        unsigned int str_to_uint(std::string);

		/**
		*   Converts a 4 byte uint to its apparent string.
		*
		*   @param          unsigned int            integer value
		*   @return         string                  converted string value
		*/
		std::string uint_to_str(unsigned int);

        /**
        *   Returns a string of all characters from an ifstream until a space,
        *   newline, or end of file is encountered.
        *
        *   @param          ifstream &              input file stream
        *   @return         string                  all characters until a space is encountered
        */
        std::string str_until_space(std::ifstream &);
};

#endif   // !defined _PNMIO_H_
