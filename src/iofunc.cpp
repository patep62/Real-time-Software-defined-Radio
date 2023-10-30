/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "iofunc.h"

// some basic functions for printing information from vectors
// or to read from/write to binary files in 32-bit float format
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

void printComplexVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (unsigned int i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

// assumes data in the raw binary file is in 32-bit float format
void readBinData(const std::string in_fname, std::vector<float> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(float);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(float));
	fdin.close();
}
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &i_data, std::vector<float> &q_data){

	std::vector<char> raw_data(num_samples);
std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
int counter = 0;

for(unsigned int i = 0; i < num_samples; i+=2){

	i_data[counter] = float(((unsigned char)raw_data[i] - 128)/128.0);
	q_data[counter] = float(((unsigned char)raw_data[i+1] - 128)/128.0);
	counter++;
	//std::cout << "ran";

}
}




// assumes data in the raw binary file is 32-bit float format
void writeBinData(const std::string out_fname, const std::vector<float> &bin_data)
{
	std::cout << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (unsigned int i=0; i<bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char*>(&bin_data[i]),\
								sizeof(bin_data[i]));
	}
	fdout.close();
}


// function to split an audio data where the left channel is in even samples
// and the right channel is in odd samples


// function to write audio data to a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the python script that can read this type of files
// and then reformat them to .wav files to be run on third-party players
