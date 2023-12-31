/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (unsigned int m = 0; m < Xf.size(); m++) {
		for (unsigned int k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
 	for (unsigned int i = 0; i < Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i])/Xf.size();
	}
}

// add your own code to estimate the PSD

//////////////////////////////////////////////////////

// added IDFT

void IDFT(const std::vector<std::complex<float>> &Xf, std::vector<std::complex<float>> &x) {
	x.resize(Xf.size(), static_cast<std::complex<float>>(0));
	for (unsigned int k = 0; k < x.size(); k++) {
		for (unsigned int m = 0; m < x.size(); m++) {
			std::complex<float> expval(0, 2*PI*(k*m) / Xf.size());
			x[k] += Xf[m] * std::exp(expval);
		}
		x[k] /= Xf.size();
	}
}

// added FFT

unsigned int swap_bits(unsigned int x, unsigned char i, unsigned char j) {

  unsigned char bit_i = (x >> i) & 0x1L;
  unsigned char bit_j = (x >> j) & 0x1L;

  unsigned int val = x;
  val = bit_i ? (val | (0x1L << j)) : (val & ~(0x1L << j));
  val = bit_j ? (val | (0x1L << i)) : (val & ~(0x1L << i));

  return val;
}

unsigned int bit_reversal(unsigned int x, unsigned char bit_size) {

  unsigned int val = x;

  for (int i=0; i < int(bit_size/2); i++)
    val = swap_bits(val, i, bit_size-1-i);

  return val;
}

void compute_twiddles(std::vector<std::complex<float>> &twiddles) {
  for (int k=0; k<(int)twiddles.size(); k++) {
      std::complex<float> expval(0.0, -2*PI*float(k)/ NFFT);
      twiddles[k] = std::exp(expval);
  }
}

void FFT_recursive(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf) {

  if (x.size() > 1) {
    // declare vectors and allocate space for the even and odd halves
    std::vector<std::complex<float>> xe(int(x.size()/2)), xo(int(x.size()/2));
    std::vector<std::complex<float>> Xfe(int(x.size()/2)), Xfo(int(x.size()/2));

    // split into even and odd halves
    for (int k=0; k<(int)x.size(); k++)
      if ((k%2) == 0) xe[k/2] = x[k];
      else xo[k/2] = x[k];

    // call recursively FFT of half size for even and odd halves respectively
    FFT_recursive(xe, Xfe);
    FFT_recursive(xo, Xfo);

    // merge the results from the odd/even FFTs (each of half the size)
    for (int k=0; k<(int)xe.size(); k++) {
        std::complex<float> expval(0.0, -2*PI*float(k)/ x.size());
        std::complex<float> twiddle = std::exp(expval);
        Xf[k]           = Xfe[k] + twiddle * Xfo[k];
        Xf[k+xe.size()] = Xfe[k] - twiddle * Xfo[k];
    }
  } else {
    // end of recursion - copy time domain samples to frequency bins (default values)
    Xf[0] = x[0];
  }
}

void FFT_improved(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf, \
  const std::vector<std::complex<float>> &twiddles, \
  const unsigned char recursion_level) {

  if (x.size() > 1) {
    int half_size = int(x.size()/2);
    std::vector<std::complex<float>> xe(half_size), xo(half_size);
    std::vector<std::complex<float>> Xfe(half_size), Xfo(half_size);

    for (int k=0; k<half_size; k++) {
      xe[k] = x[k*2];
      xo[k] = x[k*2+1];
    }

    FFT_improved(xe, Xfe, twiddles, recursion_level+1);
    FFT_improved(xo, Xfo, twiddles, recursion_level+1);

    for (int k=0; k<half_size; k++) {
        Xf[k]           = Xfe[k] + twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
        Xf[k+half_size] = Xfe[k] - twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
    }
  } else {
    Xf[0] = x[0];
  }
}

void FFT_optimized(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf, \
  const std::vector<std::complex<float>> &twiddles) {

  unsigned char no_levels = (unsigned char)std::log2((float)x.size());
  for (unsigned int i=0; i<x.size(); i++) {
    Xf[i] = x[bit_reversal(i, no_levels)];
  }

  unsigned int step_size = 1;

  std::complex<float> tmp;
  for (unsigned char l=0; l<no_levels; l++) {
    for (unsigned int p=0; p<x.size(); p+=2*step_size) {
      for (unsigned int k=p; k<p+step_size; k++) {
        tmp             = Xf[k] + twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
        Xf[k+step_size] = Xf[k] - twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
        Xf[k]           = tmp;
      }
    }
    step_size *= 2;
  }
}


void estimatePSD(std::vector<std::vector<std::complex<float>>> &graph, const std::vector<float> &samples, const int freq_bins, const float Fs){
	float df =0.0;
	std::cout << "Debug" << "\n";

	df = Fs/float(freq_bins);
	std::cout << "Df " << df << "\n";

	std::vector<std::complex<float>> freq; //
	std::cout << "Debug 2"<< "\n";

	freq.clear();
	int size = int(Fs/2);
	freq.resize(int(size/df),0.0);
	float start =0;

	for (int i =0;i<int(size/df);i++){
		std::cout << "Debug 7"<< start << "\n";

		freq[i] = start;
		start += df;
	}
	//std::cout << "Debug 3"<< "\n";

	std::vector<float> hann;
	hann.clear();
	hann.resize(freq_bins,0.0);

	for(int i =0;i<int(hann.size());i++){

		hann[i] = pow(sin((i*PI)/freq_bins),2);
	}


	int no_segments =0;
	no_segments = int(floor(samples.size()/float(freq_bins)));


	std::vector<std::complex<float>> Xf;
	Xf.clear();
	std::vector<std::complex<float>> psd_list;
	psd_list.clear();

	for(int k =0;k<no_segments;k++){

		// Starting and Ending iterators
		auto start = samples.begin()+ k*freq_bins;
		auto end = samples.begin() + (k+1)*freq_bins ;

		// To store the sliced vector
		int s = end - start;
		std::vector<float> window(s);

		// Copy vector using copy function()
		copy(start, end+1, window.begin());

		std::vector<float> windowed_samples;
		windowed_samples.resize(int(window.size()),0.0);

		for(int j=0;j < int(window.size());j++){
			windowed_samples[j] = window[j]*hann[j];
		}
		DFT(windowed_samples, Xf);

		auto start2 = Xf.begin();
		auto end2 = Xf.begin() + int(freq_bins/2);

		// To store the sliced vector

		std::vector<std::complex<float>> Xfpos(int(freq_bins/2));
		//Xfpos.clear();

		// Copy vector using copy function()
		copy(start2, end2+1, Xfpos.begin());

		std::vector<std::complex<float>> psd_seg;
		psd_seg.resize(Xfpos.size(),0.0);

		for(int m =0;m<int(Xfpos.size());m++){
			psd_seg[m] = 10*log10(2*((1/(Fs*freq_bins/2)) * pow(abs(Xfpos[m]),2)));
			std::cout << "PSD SEG VALUE: " << psd_seg[m].real() << "\n";

		}
		psd_list.insert( psd_list.end(), psd_seg.begin(), psd_seg.end() );

	}

	std::vector<std::complex<float>> psd_est;
	psd_est.clear();
	psd_est.resize(int(freq_bins/2),0.0);
	for(int k =0;k<int(freq_bins/2);k++){
		for(int l =0;l<no_segments;l++){
			psd_est[k] += psd_list[k + l*int(freq_bins/2)];
		//	std::cout << "PSD VALUE: " << psd_est[k] << "\n";

		}
		psd_est[k] = std::complex<float>(psd_est[k] / float(no_segments));
	}
	graph.push_back(freq);
	graph.push_back(psd_est);

}

// add your own code to estimate the PSD
