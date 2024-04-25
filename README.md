# Radon Transform in C++

This project implements the Radon Transform in C++. The Radon Transform is a mathematical transformation used in medical imaging, particularly in computed tomography (CT) scans, to reconstruct an image from its projections.

## Table of Contents

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)

## Introduction

The Radon Transform is a technique used in medical imaging and other fields to extract information about an object by measuring how it absorbs or reflects radiation from different angles. This project provides a C++ implementation of the Radon Transform algorithm.

## Dependencies

This project has the following dependencies:

- C++ compiler (supporting C++11 or later)
- [Optional] OpenCV library for image input/output

## Installation

To build the project, follow these steps:

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/radon-transform-cpp.git

2. Navigate to the project directory:

   ```bash
   cd radon-transform-cpp

3. Compile the source code using your preferred C++ compiler:

   ```bash
   g++ -std=c++11 -o radon_transform radon_transform.cpp

## Usage
To use the Radon Transform implementation, follow these steps:

1. Ensure you have compiled the source code (see Installation section).
2. Run the compiled executable with input parameters:

   ```bash
   ./radon_transform input_image.jpg output_image.jpg

Replace input_image.jpg with the path to your input image and output_image.jpg with the desired output file name.
3. The program will compute the Radon Transform of the input image and save the result to the specified output file.

## Contact

For questions or feedback, please contact Bastien Youssfi.

