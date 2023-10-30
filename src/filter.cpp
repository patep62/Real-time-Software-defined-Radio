/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"

void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	float Norm_cutoff =0.0f;
	Norm_cutoff = Fc/(Fs/2);

	for(int i =0;i<num_taps;i++)
	{
		if(i == ((num_taps-1)/2)){
			h[i] = Norm_cutoff;
		}else{

			h[i] = Norm_cutoff* (sin(PI*Norm_cutoff*(i-(num_taps-1)/2)))/(PI*Norm_cutoff*(i-(num_taps-1)/2));
		}
		h[i] = h[i]*pow(sin((i*PI)/num_taps),2);
	}

}

void impulseResponseRRC(float Fs, float num_taps, std::vector<float> &h){

	float T_symbol = 1/2375.0;
	float beta = 0.90;

	for(int k = 0; k < num_taps; k++){

		float t = (float)((k-(num_taps/2))/Fs);

		if(t == 0.0){
			h[k] = 1.0 + beta*((4/PI)-1);
		}
		else if(t == (-T_symbol/(4*beta)) || (t == T_symbol/(4*beta))){

			h[k] = (beta/sqrt(2))*(((1+2/PI)*(std::sin(PI/(4*beta)))) + ((1-2/PI)*(std::cos(PI/(4*beta)))));

		}
		else{

			h[k] = (std::sin(PI*t*(1-beta)/T_symbol) + 4*beta*(t/T_symbol)*std::cos(PI*t*(1+beta)/T_symbol))/(PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);

		}
	}
}

void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.clear();
	y.resize(x.size()+h.size()-1, 0.0);

	for( int n =0;n<y.size();n++){
		for( int k =0;k<h.size();k++){
			if((n-k)>=0 && (n-k)<x.size()){
				y[n] += h[k] * x[n-k];

			}

		}



	}

}
void block_processing0(const std::vector<float> &h, const std::vector<float> &xb, std::vector<float> &yb, std::vector<float> &state, int decim){

	unsigned int counter = 0;

	for(int n = 0; n < xb.size(); n+=decim){
		yb[counter] = 0;
		for(int k = 0; k < h.size(); k++){

			if ((n-k >= 0) && ((n-k) < (xb.size()))){

				yb[counter] += h[k] * xb[n-k];

			}

			else{

				yb[counter] += h[k] * state[state.size() + (n-k)];
			}
		}

		counter++;

	}

	state = slicing(xb, xb.size() - h.size(), xb.size()-1);

}
void block_processing1(const std::vector<float> h, const std::vector<float> xb, std::vector<float> &yb, std::vector<float> &state, int decim, int u){

    int n;
    unsigned int counter = 0;
    int phase;
    for(n = 0; n < xb.size()*u; n+=decim){
        yb[counter] = 0;
        phase = n%u;
        for(int k = phase; k < h.size(); k+=u){
            if ((n/u)-(k/u) >= 0){
                yb[counter] += u * h[k] * xb[(n/u)-(k/u)];
            }else{
                yb[counter] += u * state[state.size() + (n/u)-(k/u)] * h[k];
            }
        }
        counter++;
    }

    state = slicing(xb, xb.size() - h.size(), xb.size()-1);

}

void Squaring_nonlinearity(std::vector<float> &RDS_carrier){

	for(int i = 0; i < RDS_carrier.size(); i++){

		RDS_carrier[i] = RDS_carrier[i] * RDS_carrier[i];

	}
}

void my_fmdemod(std::vector<float> &I, std::vector<float> &Q, float &I_state, float &Q_state, std::vector<float> &fm_demod){

	float I_derivative;
	float Q_derivative;

	float I_previous = I_state;
	float Q_previous = Q_state;

	for(int i = 0; i < I.size(); i++){

		I_derivative = I[i] - I_previous;
		Q_derivative = Q[i] - Q_previous;

		if(I[i] == 0 && Q[i] == 0){
			fm_demod[i] = 0;
		}

		else{

			fm_demod[i] = ((I[i] * Q_derivative) - (Q[i] * I_derivative)) / (pow(I[i],2) + pow(Q[i],2));

		}

		I_previous = I[i];
		Q_previous = Q[i];

	}

	I_state = I_previous;
	Q_state = Q_previous;
}

std::vector<float>slicing(const std::vector<float>& arr, int X, int Y){

    // Starting and Ending iterators
    auto start = arr.begin() + X;
    auto end = arr.begin() + Y + 1;

    // To store the sliced vector
    std::vector<float> result(Y - X + 1);

    // Copy vector using copy function()
    std::copy(start, end, result.begin());

    // Return the final sliced vector
    return result;
}

std::vector<float>downsample(std::vector<float> &arr, int decim){

	int counter = 0;
	std::vector<float> arr_ds(arr.size()/decim, 0.0);

	for(unsigned int i = 0; i < arr.size(); i += decim){

		arr_ds[counter] = arr[i];
		counter++;

	}

	return arr_ds;
}

void bandpass(float Fs, float fe, float fb, unsigned short int num_taps, std::vector<float> &h){

	float Norm_c = ((fe+fb)/2)/(Fs/2);
	float Norm_p = (fe-fb)/(Fs/2);


	for(int i =0;i<num_taps-1;i++){
		if(i == ((num_taps-1)/2)){
			h[i] = Norm_p; // avoid division by 0
		}
		else{
			h[i]=Norm_p*((sin(PI*(Norm_p/2)*(i-(num_taps-1)/2)))/(PI*(Norm_p/2)*(i-(num_taps-1)/2)));
		}
		h[i] = h[i]*cos(i*PI*Norm_c);
		h[i] = h[i]*pow(sin((i*PI)/num_taps),2);

	}

}

void mixer(const std::vector<float> arr, const std::vector<float> arr2, std::vector<float> &result){
	for(unsigned int i =0;i<arr2.size();i++){
		result[i] = arr[i]*arr2[i]*2;
	}
}



void fmPll(std::vector<float> &pllIn, std::vector<float> &ncoOut, float &integrator, float &phaseEst, float phaseAdjust, float &feedbackI, float &feedbackQ, float &trigOffset, float freq, float Fs, float ncoScale, float normBandwidth){

	float Cp = 2.666;
	float Ci = 3.555;

	float Kp = (normBandwidth*Cp);
	float Ki = (normBandwidth*normBandwidth*Ci);

	float errorI = 0;
	float errorQ = 0;
	float errorD = 0;
	float trigArg = 0;

	for(int k = 0; k < pllIn.size(); k++){

		errorI = pllIn[k] * (+feedbackI);
		errorQ = pllIn[k] * (-feedbackQ);

		errorD = std::atan2(errorQ, errorI);

		integrator += (Ki*(errorD));
		phaseEst += (Kp*errorD) + integrator;

		trigOffset++;
		trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;

		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);

		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);

	}
}

void recombine_mono_stereo(const std::vector<float>mono,const std::vector<float>stereo, std::vector<float> &LR_channel){

	for(unsigned int i =0;i<mono.size();i++){

		LR_channel[i*2] = (mono[i]+stereo[i])/2; // left
		LR_channel[(i*2)+1] = (mono[i]-stereo[i])/2; // right channel

	// L0R0L1R1

	}
}

void allpass(const std::vector<float>input_block, std::vector<float>&state_block, std::vector<float>&output_block){
	output_block.clear();
	unsigned int size = state_block.size();
	output_block.insert(output_block.begin(), state_block.begin(), state_block.end());
	output_block.insert(output_block.end(), input_block.begin(), input_block.end()-size);

	state_block.clear();
	state_block.insert(state_block.begin(),input_block.end()-size,input_block.end());

}

void findSamplingPoint(const std::vector<float>RDS_RRC_Output, int &sampling_point){

	int start = 0;

	for(int i = 0; i < RDS_RRC_Output.size(); i++){

		if(RDS_RRC_Output[i] >= 0.25){
			start = i;
			break;
		}
	}

	float current = RDS_RRC_Output[start+1];
	float previous = RDS_RRC_Output[start];

	if(previous >= 0){

		if(current >= previous){

			for(int i = start+1; i < RDS_RRC_Output.size(); i++){

				current = RDS_RRC_Output[i];
				if(current < previous){
					sampling_point = i-1;
					break;
				}

				else{
					previous = current;
				}
			}
		}

		else{

			for(int i = start+1; i < RDS_RRC_Output.size(); i++){

				current = RDS_RRC_Output[i];
				if(current > previous){
					sampling_point = i-1;
					break;
				}
				else{
					previous = current;
				}
			}
		}
	}

	else{

		if(current < previous){

			for(int i = start+1; i < RDS_RRC_Output.size(); i++){

				current = RDS_RRC_Output[i];
				if(current > previous){
					sampling_point = i-1;
					break;
				}

				else{
					previous = current;
				}
			}
		}

		else{

			for(int i = start+1; i < RDS_RRC_Output.size(); i++){

				current = RDS_RRC_Output[i];

				if(current < previous){
					sampling_point = i-1;
					break;
				}

				else{
					previous = current;
				}
			}
		}
	}
}

void CDR(const std::vector<float>RDS_RRC_Output, std::vector<int> &bitstream, int &sampling_point, const int SPS, int &lastSymbol){

	//Very Initial Version of CDR. Error checking for Invalid bitstreams has not been implemented.

	float next;
	float current;
	int i;
	bitstream.clear();

	for(i = sampling_point; i < RDS_RRC_Output.size(); i+=SPS*2){

		current = RDS_RRC_Output[i];
		next = RDS_RRC_Output[i+SPS];
		if(current >= 0){
			if(next >= 0){
				//INVALID BIT
			}
			else{
				bitstream.push_back(1);
			}
		}
		else{
			if(next >= 0){
				bitstream.push_back(0);
			}
			else{
				//INVALID BIT
			}
		}
	}
	if(current >= 0){
		lastSymbol = 1;
	}
	else{
		lastSymbol = 0;
	}
}

void differentialDecoder(const std::vector<int> &bitstream, std::vector<int> &output){

	output.push_back(bitstream[0]);
	for(int i = 1; i < bitstream.size(); i++){

		output.push_back(bitstream[i]^bitstream[i-1]);

	}
}

void Pcheck(const std::vector<int> input, std::vector<int> &output){
    std::vector<int>newinput(26,0);

    std::vector<int> blockA{1,1,1,1,0,1,1,0,0,0};


    std::vector<int> blockB{1,1,1,1,0,1,0,1,0,0};

    std::vector<int> blockC{1,0,0,1,0,1,1,1,0,0};


    std::vector<int> blockCPrime{1,1,1,1,0,0,1,1,0,0};

    std::vector<int> blockD{1,0,0,1,0,1,1,0,0,0};

    std::vector<std::vector<int>> vect{
                    {1,0,0,0,0,0,0,0,0,0},
                    {0,1,0,0,0,0,0,0,0,0},
                    {0,0,1,0,0,0,0,0,0,0},
                    {0,0,0,1,0,0,0,0,0,0},
                    {0,0,0,0,1,0,0,0,0,0},
                    {0,0,0,0,0,1,0,0,0,0},
                    {0,0,0,0,0,0,1,0,0,0},
                    {0,0,0,0,0,0,0,1,0,0},
                    {0,0,0,0,0,0,0,0,1,0},
                    {0,0,0,0,0,0,0,0,0,1},
                    {1,0,1,1,0,1,1,1,0,0},
                    {0,1,0,1,1,0,1,1,1,0},
                    {0,0,1,0,1,1,0,1,1,1},
                    {1,0,1,0,0,0,0,1,1,1},
                    {1,1,1,0,0,1,1,1,1,1},
                    {1,1,0,0,0,1,0,0,1,1},
                    {1,1,0,1,0,1,0,1,0,1},
                    {1,1,0,1,1,1,0,1,1,0},
                    {0,1,1,0,1,1,1,0,1,1},
                    {1,0,0,0,0,0,0,0,0,1},
                    {1,1,1,1,0,1,1,1,0,0},
                    {0,1,1,1,1,0,1,1,1,0},
                    {0,0,1,1,1,1,0,1,1,1},
                    {1,0,1,0,1,0,0,1,1,1},
                    {1,1,1,0,0,0,1,1,1,1},
                    {1,1,0,0,0,1,1,0,1,1}
                };

for(int s=0; s<10; s++){
		// Starting and Ending iterators
		auto start = input.begin() + s;
		auto end = input.begin() + s + 26;

		// To store the sliced vector
		std::vector<int> newinput(s -(s + 26));

		// Copy vector using copy function()
		std::copy(start, end, newinput.begin());
    for(int i=0; i<1; i++){
            for(int j=0; j<10; j++){
                    for(int k=0; k<26; k++){
                            output[j] += (newinput[i] * vect[k][j]);
                    }
            }
    }
    if(std::equal(output.begin(), output.end(),blockA.end())){
         std::cerr << "Match with Block A";
    }
    if(std::equal(output.begin(), output.end(),blockB.end())){
         std::cerr << "Match with Block B";
    }

    if(std::equal(output.begin(), output.end(),blockC.end())){
         std::cerr << "Match with Block C";
    }

    if(std::equal(output.begin(), output.end(),blockCPrime.end())){
         std::cerr << "Match with Block C Prime";
    }

    if(std::equal(output.begin(), output.end(),blockD.end())){
         std::cerr << "Match with Block D";
    }
}
}

void ApplicationLayer(const std::vector<int> &bitstream, std::vector<int> &PI_Code, std::vector<int> &program_type, std::vector<int> &program_service_name){

	int counter = 0;
	bool flag = true;
	for(int i = 0; i < bitstream.size(); i+= 104){

		//PI Code Extraction
		for(int a = i; a < 15+i; a++){
			PI_Code[counter] = bitstream[a];
			counter++;
		}
		counter = 0;

		//Program Type Extraction
		for(int b = 5+26+i; b < 10+26+i; b++){
			program_type[counter] = bitstream[b];
			counter++;
		}
		counter = 0;

		//Program Service Name Extraction
		//Check if the group is of type A
		for(int c = 26+i; c < 26+4+i; c++){
			if(bitstream[c] != 0){
				flag = false;
				break;
			}
		}

		//If we are on the correct group, extract the data
		if(flag){
			for(int d = 78+i; d < 78+15+i; d++){
				program_service_name[counter] = bitstream[d];
				counter++;
			}
		}
	}
}
