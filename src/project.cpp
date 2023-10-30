/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <chrono>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

#define QUEUE_ELEMS 8

//Producer Thread
void RF_thread(std::queue<std::vector<float>> &my_queue_MonoStereo, std::queue<std::vector<float>> &my_queue_RDS, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode){

	//Default mode 0 paramaters
	int num_taps = 151;
	int rf_Fs = 2.4e6;
	int rf_Fc = 100e3;
	int rf_decim = 10;

	int mono_Fs = 48e3;
	int mono_Fc = 16e3;
	int mono_decim = 5;
	int mono_upscale = 1;

	//Parameter selection for modes 1,2,3.
	if(mode == 1){

		rf_Fs = 1.44e6;
		rf_decim = 6;

	}

	else if(mode == 2){

		mono_Fs = 44.1e3;
		mono_decim = 800;
		mono_upscale = 147;

	}

	else if(mode == 3){

		rf_Fs = 1.152e6;
		rf_decim = 4;

		mono_Fs = 44.1e3;
		mono_decim = 320;
		mono_upscale = 49;

	}

	int block_size = 1024 * rf_decim * (mono_decim/mono_upscale) * 2;

	//IQ Extraction Initializations
	std::vector<float> rf_coeff(num_taps, 0.0);
	impulseResponseLPF(rf_Fs, rf_Fc, num_taps, rf_coeff);
	std::vector<float> I_statelpf(num_taps-1, 0.0);
	std::vector<float> Q_statelpf(num_taps-1, 0.0);
	std::vector<float> i_data(block_size/2, 0.0);
	std::vector<float> q_data(block_size/2, 0.0);

	//FM Demodulation Initializations
	std::vector<float> i_ds(block_size/(rf_decim*2), 0.0);
	std::vector<float> q_ds(block_size/(rf_decim*2), 0.0);
	std::vector<float> fm_demod(block_size/(rf_decim*2), 0.0);
	float I_state = 0;
	float Q_state = 0;

	for(unsigned int block_id = 0; ; block_id++){

		readStdinBlockData(block_size, block_id, i_data, q_data);

		if((std::cin.rdstate()) != 0){

			std::cerr << "End of input stream reached" << std::endl;
			exit(1);

		}

		//Extract IQ data
		block_processing0(rf_coeff, i_data, i_ds, I_statelpf, rf_decim);
		block_processing0(rf_coeff, q_data, q_ds, Q_statelpf, rf_decim);

		//Demodulate IQ Data
		my_fmdemod(i_ds, q_ds, I_state, Q_state, fm_demod);

		//Thread Synchronization
		std::unique_lock<std::mutex> my_lock(my_mutex);
		while(my_queue_MonoStereo.size() >= QUEUE_ELEMS || my_queue_RDS.size() >= QUEUE_ELEMS){

			my_cvar.wait(my_lock);

		}

		my_queue_RDS.push(fm_demod);
		my_queue_MonoStereo.push(fm_demod);
		my_cvar.notify_all();
		my_lock.unlock();

	}
}

//Consumer Thread 1
void MonoStereo_thread(std::queue<std::vector<float>> &my_queue, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode){

	//Default mode 0 paramaters
	int num_taps = 151;
	int rf_Fs = 2.4e6;
	int rf_Fc = 100e3;
	int rf_decim = 10;

	int mono_Fs = 48e3;
	int mono_Fc = 16e3;
	int mono_decim = 5;
	int mono_upscale = 1;

	//Parameter selection for modes 1,2,3.
	if(mode == 1){

		rf_Fs = 1.44e6;
		rf_decim = 6;

	}

	else if(mode == 2){
		mono_Fs = 44.1e3;
		mono_decim = 800;
		mono_upscale = 147;

	}

	else if(mode == 3){

		rf_Fs = 1.152e6;
		rf_decim = 4;

		mono_Fs = 44.1e3;
		mono_decim = 320;
		mono_upscale = 49;

	}

	int block_size = 1024 * rf_decim * (mono_decim/mono_upscale) * 2;
	int IF_sample_rate = rf_Fs/rf_decim;

	std::vector<float> audio_coeff(num_taps*mono_upscale, 0.0);
	if((mono_upscale*mono_Fs)/(mono_decim*2) < (mono_Fs/2)){
		impulseResponseLPF(mono_Fs*mono_upscale, (mono_upscale*mono_Fs)/(mono_decim*2), num_taps*mono_upscale, audio_coeff);


	}else{
		impulseResponseLPF(mono_Fs*mono_upscale, (mono_Fs/2), num_taps*mono_upscale, audio_coeff);
	}
	std::vector<float> state(num_taps - 1, 0.0);
	std::vector<float> audio_block(block_size/(rf_decim*mono_decim/mono_upscale*2), 0.0);

	int fe_st = 54e3;
  int fb_st = 22e3;

  int fe_carrier = 18.5e3;
  int fb_carrier = 19.5e3;
  std::vector<float> stereo_coeff(num_taps, 0.0);
  std::vector<float> carrier_coeff(num_taps, 0.0);

  std::vector<float> stereo_state(num_taps - 1, 0.0);
  std::vector<float> stereo_block(block_size/(rf_decim*2), 0.0);

  std::vector<float> carrier_state(num_taps - 1, 0.0);
  std::vector<float> carrier_block(block_size/(rf_decim*2), 0.0);

  std::vector<float> ncoOut(block_size/(rf_decim*2), 0.0);
  std::vector<float> mixer_output(block_size/(rf_decim*2), 0.0);

  std::vector<float> audio_stereo_state(num_taps - 1, 0.0);
  std::vector<float> audio_stereo_block(block_size/(rf_decim*mono_decim/mono_upscale*2), 0.0);

  std::vector<float> LR_channel(block_size/(rf_decim*mono_decim/mono_upscale), 0.0);

	ncoOut[0] = 1.0;
	float integrator = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	float phaseEst = 0.0;
	float trigOffset = 0.0;
	int phase =0;

	// Stereo Coeffecents
	bandpass(rf_Fs/rf_decim,fe_st,fb_st,num_taps,stereo_coeff); // Stereo coeff at 240 K
	bandpass(rf_Fs/rf_decim,fe_carrier,fb_carrier,num_taps,carrier_coeff); // Stereo Carrier at 240 K

	//All Pass delay blocks
	std::vector<float> delay_state_block((num_taps-1)/2, 0.0);
	std::vector<float> output_fm_block(block_size/(rf_decim*2), 0.0);
	int index =0;
	for(unsigned int block_id = 0; ; block_id++){

		std::unique_lock<std::mutex> my_lock(my_mutex);
		while(my_queue.empty()){
			my_cvar.wait(my_lock);
		}

		std::vector<float> fm_demod = my_queue.front();
		my_queue.pop();
		my_cvar.notify_one();
		my_lock.unlock();

		allpass(fm_demod, delay_state_block, output_fm_block);

		block_processing1(audio_coeff, output_fm_block, audio_block, state, mono_decim,mono_upscale);

		block_processing1(stereo_coeff, fm_demod, stereo_block, stereo_state, 1,1); // stereo channel extraction

		block_processing1(carrier_coeff, fm_demod, carrier_block, carrier_state, 1,1); // carrier extraction

		if(block_id > 0){
			ncoOut[0] = ncoOut[ncoOut.size()-1];
		}

		fmPll(carrier_block, ncoOut, integrator, phaseEst, 0.0, feedbackI, feedbackQ, trigOffset, 19e3, rf_Fs/rf_decim, 2.0, 0.01);

		// Mixing
		mixer(stereo_block, ncoOut, mixer_output);

		// Digital filtering with decimator of 5
		block_processing1(audio_coeff, mixer_output, audio_stereo_block, audio_stereo_state, mono_decim, mono_upscale);

		recombine_mono_stereo(audio_block,audio_stereo_block, LR_channel);

		std::vector<short int> audio_data(LR_channel.size());

		for(unsigned int k = 0; k < LR_channel.size(); k++){

			if(std::isnan(LR_channel[k])) audio_data[k] = 0;

			else audio_data[k] = static_cast<short int>(LR_channel[k]*16384);

		}

		fwrite(&audio_data[0], sizeof(short int), LR_channel.size(), stdout);

		if((std::cin.rdstate()) != 0 && my_queue.empty()){

			std::cerr << "End of input stream reached" << std::endl;
			exit(1);

		}
	}
}

//Consumer Thread 2
void RDS_thread(std::queue<std::vector<float>> &my_queue_RDS, std::mutex &my_mutex, std::condition_variable &my_cvar, int mode){

	//Default mode 0 paramaters
	int num_taps = 151;
	int rf_Fs = 2.4e6;
	int rf_Fc = 100e3;
	int rf_decim = 10;

	int mono_decim = 5;
	int mono_upscale = 1;

	int SPS = 20;
	int RDS_decim = 96;
	int RDS_upscale = 19;

	//Parameter selection for modes 0,1.

	if(mode == 5){

		RDS_decim = 1920;
		RDS_upscale = 703;
		SPS = 37;

	}

	//Filter Parameters
	int IF_sample_rate = rf_Fs/rf_decim;
	int block_size = 1024 * rf_decim * (mono_decim/mono_upscale) * 2;

	//RDS Channel Initializations
	std::vector<float> RDS_channel_coeff(num_taps, 0.0);
	std::vector<float> RDS_channel(block_size/(rf_decim*2), 0.0);
	std::vector<float> RDS_channel_state(num_taps-1, 0.0);
	bandpass(IF_sample_rate, 60e3, 54e3, num_taps, RDS_channel_coeff);

	//RDS Carrier Initializations
	std::vector<float> RDS_carrier_coeff(num_taps, 0.0);
	std::vector<float> RDS_carrier(block_size/(rf_decim*2), 0.0);
	std::vector<float> RDS_carrier_state(num_taps-1, 0.0);
	bandpass(IF_sample_rate, 114.5e3, 113.5e3, num_taps, RDS_carrier_coeff);

	//Mixer Initialization
	std::vector<float> mixer_output(block_size/(rf_decim*2), 0.0);

	//RDS Demodulation Initializations
	std::vector<float> RDS_demod_coeff(num_taps, 0.0);
	impulseResponseLPF(IF_sample_rate, 3e3, num_taps, RDS_demod_coeff);

	//RDS Baseband Initialization
	std::vector<float> RDS_baseband(block_size/(rf_decim*2), 0.0);
	std::vector<float> RDS_baseband_state(num_taps-1, 0.0);

	//Root Raised Cosine Filter Initializations
	std::vector<float> RRC_coeff(num_taps, 0.0);
	std::vector<float> RDS_RRC_Output(block_size/(RDS_decim/RDS_upscale*rf_decim*2), 0.0);
	std::vector<float> RDS_RRC_State(num_taps-1, 0.0);
	impulseResponseRRC(2375*SPS, num_taps, RRC_coeff);

	//Rational Resamler Initializations
	std::vector<float> RDS_resampled(block_size/(RDS_decim/RDS_upscale*rf_decim*2), 0.0);
	std::vector<float> RDS_resample_coeff(num_taps*RDS_upscale, 0.0);
	std::vector<float> RDS_resample_state(num_taps*RDS_upscale-1, 0.0);
	impulseResponseLPF(IF_sample_rate*RDS_upscale, 16e3, num_taps*RDS_upscale, RDS_resample_coeff);

	//Delay block Initializations
	std::vector<float> state_block((num_taps-1)/2, 0.0);
	std::vector<float> output_block(block_size/(rf_decim*2), 0.0);

	//Pll Initializations
	std::vector<float> ncoOut(block_size/(rf_decim*2), 0.0);
	ncoOut[0] = 1.0;
	float integrator = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	float phaseEst = 0.0;
	float trigOffset = 0.0;

	//CDR Initializations
	int sampling_point = 0;
	std::vector<int> bitstream;
	std::vector<int> output_bitstream;
	int lastSymbol = 0;

	//Application Layer Initializations
	std::vector<int> PI_Code(16,0);
	std::vector<int> program_type(5,0);
	std::vector<int> program_service_name(16,0);

	for(unsigned int block_id = 0; ; block_id++){

		//Thread Synchronization
		std::unique_lock<std::mutex> my_lock(my_mutex);
		while(my_queue_RDS.empty()){
			my_cvar.wait(my_lock);
		}

		std::vector<float> fm_demod = my_queue_RDS.front();
		my_queue_RDS.pop();
		my_cvar.notify_one();
		my_lock.unlock();

		//When running in mode 3, all of the RDS thread needs to be commented out, or else the program will not run.
		//We think this had something to do with the additional load created by the RDS thread.

		//RDS CHANNEL EXTRACTION AND RECOVERY
		block_processing0(RDS_channel_coeff, fm_demod, RDS_channel, RDS_channel_state, 1);
		allpass(RDS_channel, state_block, output_block); //Send to mixer.
		Squaring_nonlinearity(RDS_channel);
		block_processing0(RDS_carrier_coeff, RDS_channel, RDS_carrier, RDS_carrier_state, 1);

		//Pll State Saving
		if(block_id > 0){
			ncoOut[0] = ncoOut[ncoOut.size()-1];
		}

		//Pll Processing
		fmPll(RDS_carrier, ncoOut, integrator, phaseEst, 0.0, feedbackI, feedbackQ, trigOffset, 114e3, IF_sample_rate, 0.5, 0.003);

		//Mixing
		mixer(ncoOut, output_block, mixer_output);

		//RDS DEMODULATION
		block_processing0(RDS_demod_coeff, mixer_output, RDS_baseband, RDS_baseband_state, 1);

		//Resampling
		block_processing1(RDS_resample_coeff, RDS_baseband, RDS_resampled, RDS_resample_state, RDS_decim, RDS_upscale); //RDS_baseband = input, RDS_resampled = output.

		//Root Raised cosine filtering.
		block_processing0(RRC_coeff, RDS_resampled, RDS_RRC_Output, RDS_RRC_State, 1);


		//Clock and data recovery.
		if(block_id == 0){

			//Analyze first block of data to derive the sampling point.
			findSamplingPoint(RDS_RRC_Output, sampling_point);
			CDR(RDS_RRC_Output, bitstream, sampling_point, SPS, lastSymbol);

			//Modulo The sampling point with the SPS and pass as state to next block.
			sampling_point = sampling_point % SPS;

		}

		else{

			CDR(RDS_RRC_Output, bitstream, sampling_point, SPS, lastSymbol);

		}

		differentialDecoder(bitstream, output_bitstream);

		//Frame Synchronization is not working properly so we left it commented out.
		//std::vector<int> synched_bitstream(output_bitstream.size(), 0);
		//Pcheck(output_bitstream, synched_bitstream);

		//Final Application Layer to extract RDS data
		//ApplicationLayer(synched_bitstream, PI_Code, program_type, program_service_name);

		if((std::cin.rdstate()) != 0 && my_queue_RDS.empty()){

			std::cerr << "End of input stream reached" << std::endl;
			exit(1);

		}
	}
}

int main(int argc, char *argv[])
{
	int mode = 0;

	if(argc < 2){
		std::cerr << "Operating in default mode 0" << std::endl;
	} else if(argc == 2){
		mode = atoi(argv[1]);
		if (mode > 3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);
		}
	} else{
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << "<mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
		exit(1);
	}

	std::cerr << "Operating in mode " << mode << std::endl;

	std::queue<std::vector<float>> my_queue_MonoStereo;
	std::queue<std::vector<float>> my_queue_RDS;
	std::mutex my_mutex;
	std::condition_variable my_cvar;

	std::thread ta = std::thread(RF_thread, std::ref(my_queue_MonoStereo), std::ref(my_queue_RDS), std::ref(my_mutex), std::ref(my_cvar), mode);
	std::thread tb = std::thread(MonoStereo_thread, std::ref(my_queue_MonoStereo), std::ref(my_mutex), std::ref(my_cvar), mode);
  std::thread tc = std::thread(RDS_thread, std::ref(my_queue_RDS), std::ref(my_mutex), std::ref(my_cvar), mode);

	ta.join();
	tb.join();
	tc.join();

	return 0;
}
