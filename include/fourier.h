/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FOURIER_H
#define DY4_FOURIER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// declaration of a function prototypes
void DFT(const std::vector<float> &,std::vector<std::complex<float>> &);

// you should add your own IDFT
// time-permitting you can build your own function for FFT
//std::vector<float> downsample(std::vector<float> &, int)


void computeVectorMagnitude(const std::vector<std::complex<float>> &,
	std::vector<float> &);

// provide the prototype to estimate PSD
// ...
void IDFT(const std::vector<std::complex<float>> &, std::vector<std::complex<float>> &);

void FFT_recursive(const std::vector<std::complex<float>> &, std::vector<std::complex<float>> &);

void FFT_improved(const std::vector<std::complex<float>> &, std::vector<std::complex<float>> &, const std::vector<std::complex<float>> &, const unsigned char);

void FFT_optimized(const std::vector<std::complex<float>> &, std::vector<std::complex<float>> &, const std::vector<std::complex<float>> &);

void compute_twiddles(std::vector<std::complex<float>> &);

void estimatePSD(std::vector<std::vector<std::complex<float>>> &, const std::vector<float> &,int , float);

#endif // DY4_FOURIER_H
