/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>


// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void block_processing0(const std::vector<float> &, const std::vector<float> &, std::vector<float> &, std::vector<float> &, int);//void blocks_processing_with_state(std::vector<float> &,const std::vector<float> &,const std::vector<float> &,std::vector<float> &);
void my_fmdemod(std::vector<float> &, std::vector<float> &, float &, float &, std::vector<float> &);
std::vector<float> slicing(const std::vector<float>&, int, int);
void bandpass(float, float, float, unsigned short int, std::vector<float> &);
void mixer(const std::vector<float>, const std::vector<float>,std::vector<float> &);
void fmPll(std::vector<float> &, std::vector<float> &, float &, float &, float , \
  float &, float &, float &, float , float , float , float);
void recombine_mono_stereo(const std::vector<float> ,const std::vector<float>, std::vector<float> &);
void allpass(const std::vector<float>, std::vector<float>&, std::vector<float>&);
void block_processing1(const std::vector<float> , const std::vector<float> , std::vector<float> &, std::vector<float> &, int , int);
void Squaring_nonlinearity(std::vector<float> &);
void impulseResponseRRC(float, float, std::vector<float> &);
void findSamplingPoint(const std::vector<float>, int &);
void CDR(const std::vector<float>, std::vector<int> &, int &, const int, int &lastSymbol);
void differentialDecoder(const std::vector<int> &, std::vector<int> &);
void Pcheck(const std::vector<int>, std::vector<int> &);
void ApplicationLayer(const std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<int> &);

#endif // DY4_FILTER_H
