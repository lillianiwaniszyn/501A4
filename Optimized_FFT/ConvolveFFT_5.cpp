
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <malloc.h>
#include <cstdint>
#include <iostream>
#include <cstdio>
#include <ctime>
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr


using namespace std;


char *format;
int chunkSize;
char *chunkID;
char *subChunk1ID;
int size;
int16_t blockAlign;
int16_t bitsPerSample;
int16_t audioFormat;
int16_t numChannels;
int sampleRate;
char *subChunk2ID;
int dataSize;
short* fileData;
int subChunk1Size;
int byteRate;

float minRead = 0;
float maxRead = 0;

//function definitions
void wavWrite(char *fileName, int numSamples, float *signal);
float* wavRead(char *fileName, float *signal, int *Thesize);
void convolve(float x[], int N, float h[], int M, float y[], int P);
void fft(float data[], int nn, int isign);
void FFTScale (float signal[], int N);



int main(int argc, char* args[])
{
	std::clock_t start;
    double duration;
    start = std::clock();
	if (argc!= 4){ //check if we hve 3 command line arguments
		
		printf("Please enter an input file, and IR file and an output file name." );
		return 0;
	}

	char *inputFileName = args[1];
	char *IRFileName = args[2];
	char *outputFileName = args[3];
	float *inFileSignal;
	int inFileSignalSize;
	float *IRFileSignal;
	int IRFileSignalSize;
	inFileSignal = wavRead(inputFileName, inFileSignal, &size );
	inFileSignalSize = size;
	IRFileSignal = wavRead(IRFileName, IRFileSignal, &size);
	IRFileSignalSize = size;
	int outFileSignalSize = inFileSignalSize + IRFileSignalSize - 1;
	float *outFileSignal = new float[outFileSignalSize];
	printf("Convolving...");
	convolve(inFileSignal, inFileSignalSize, IRFileSignal, IRFileSignalSize, outFileSignal, outFileSignalSize);
	//scale output below
	float min = 0, max = 0;
	int i = 0;

	for(i = 0; i < outFileSignalSize; i++)
	{
		if(outFileSignal[i] > max)
			max = outFileSignal[i];
		if(outFileSignal[i] < min)
			min = outFileSignal[i];
	}

	min = min * -1;
	if(min > max)
		max = min;
	for(i = 0; i < outFileSignalSize; i++)
	{
		outFileSignal[i] = outFileSignal[i] / max;
	}
	wavWrite(outputFileName, outFileSignalSize, outFileSignal);
	printf("Written to file: %s\n", outputFileName);
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("The duration of the program is: ");
	std::cout<< duration <<" seconds";
	return 0;
}

void wavWrite(char *fileName, int numSamples, float *signal)
{
	ofstream outFile( fileName, ios::out | ios::binary);
	//  Calculate the total number of bytes for the data chunk  
	int chunkSize = numChannels * numSamples * (bitsPerSample / 8);
	chunkID = "RIFF";
	outFile.write( chunkID, 4);
	outFile.write( (char*) &chunkSize, 4);
	format = "WAVE";
	outFile.write( format, 4);
	outFile.write( subChunk1ID, 4);
	subChunk1Size = 16;
	outFile.write( (char*) &subChunk1Size, 4);
	outFile.write( (char*) &audioFormat, 2);
	outFile.write( (char*) &numChannels, 2);
	outFile.write( (char*) &sampleRate, 4);
	outFile.write( (char*) &byteRate, 4);
	outFile.write( (char*) &blockAlign, 2);
	outFile.write( (char*) &bitsPerSample, 2);
	outFile.write( subChunk2ID, 4);
	dataSize = numSamples * 2;
	outFile.write( (char*)&dataSize, 4);
	int16_t val;
	for(int i = 0; i < numSamples; i++)
	{
		val = (int16_t)(signal[i] * (32767));
		outFile.write((char*)&val, 2);
	}
	outFile.close();
}


float* wavRead(char *fileName, float *signal, int *Thesize)
{
	ifstream inputFile( fileName, ios::in | ios::binary);
	inputFile.seekg(ios::beg);
	chunkID = new char[4];
	inputFile.read( chunkID, 4);
	inputFile.read( (char*) &chunkSize, 4);
	format = new char[4];
	inputFile.read( format, 4);
	subChunk1ID = new char[4];
	inputFile.read( subChunk1ID, 4);
	inputFile.read( (char*) &subChunk1Size, 4);
	inputFile.read( (char*) &audioFormat, 2);
	inputFile.read( (char*) &numChannels, 2);
	inputFile.read( (char*) &sampleRate, 4);
	inputFile.read( (char*) &byteRate, 4);
	inputFile.read( (char*) &blockAlign, 2);
	inputFile.read( (char*) &bitsPerSample, 2);

	if(subChunk1Size == 18)
	{
		char *garbage;
		garbage = new char[2];
		inputFile.read( garbage, 2);
	}

	subChunk2ID = new char[4];
	inputFile.read( subChunk2ID, 4);
	//DataSize
	inputFile.read( (char*)&dataSize, 4);
	//GetData
	*Thesize = dataSize / 2;
	int size = dataSize / 2;
	fileData = new short[size];
	for(int j = 0 ; j < size; j++)
	{
		inputFile.read((char*) &fileData[j], 2);
	}

	short val;
	signal = new float[size];
	for(int i = 0; i < size; i++)
	{
		val = fileData[i];
		signal[i] = (val * 1.0) / (32767); //changed this because strength reduction
		if(signal[i] < -1.0)
			signal[i] = -1.0;

	}
	inputFile.close();
	return signal;
}


//uses FFT algorithm to convolve
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
	int myArraySize = 1;
	int i = 0;
	// For FFT we need array size of a power of 2
	while (myArraySize < P) {
		myArraySize += myArraySize;
	}
	float *paddedInput = new float[2 * myArraySize];
	for (i = 0; i < (N * 2); i+=2) {
		paddedInput[i] = x[i/2];
		paddedInput[i+1] = 0;
	}
	for (; i <myArraySize-1; i+=2) {
		paddedInput[i] = 0;
		paddedInput[i+1] = 0; //changed code here. unrolled loop once.
	}
	if (i == myArraySize-1){
		paddedInput[myArraySize-1] = myArraySize-1;
	}
	
	float *paddedImpulseResponse = new float[2 * myArraySize];
	for (i = 0; i < (M * 2); i+=2) {
		paddedImpulseResponse[i] = h[i/2];
		paddedImpulseResponse[i+1] = 0;
	}
	for (; i < myArraySize-1; i+=2) {
		paddedImpulseResponse[i] = 0;
		paddedImpulseResponse[i+1] = 0; //unrolled second loop
	}
	if (i == myArraySize-1){
		paddedImpulseResponse[myArraySize-1] = myArraySize-1;
	}
		

	
	float *paddedOutput = new float[2 * myArraySize];
	for (i = 0; i < myArraySize; i++) { 
		paddedOutput[i] = 0.0; // used proper format for constraints
	}
	fft((paddedInput - 1), myArraySize, 1);
	fft((paddedImpulseResponse - 1), myArraySize, 1);
	for (i = 0; i < (myArraySize * 2); i+=2) {
		paddedOutput[i] = (paddedInput[i] * paddedImpulseResponse[i]) - (paddedInput[i+1] * paddedImpulseResponse[i+1]);
		paddedOutput[i+1] = (paddedInput[i+1] * paddedImpulseResponse[i]) + (paddedInput[i] * paddedImpulseResponse[i+1]);
	}
	fft((paddedOutput - 1), myArraySize, -1);
	
	// FFT scaling.. we need to scale as per class notes
	FFTScale(paddedOutput, myArraySize);
	
	// removing padding
	for (i = 0; i < P; i++) {
		y[i] = paddedOutput[i*2];
	}
}
//fft algorithm... sourced from https://ftp.samba.org/pub/unpacked/junkcode/speech/myfft.c
void fft(float data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
    }

    mmax = 2;
    while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax); //changed just this
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
    }
}

//the output from either the FFTs or the IFFT (but not both) will have to be scaled by dividing by each data point by N
//this algorithm is taken from class Notes
void FFTScale (float x[], int N)
{
	int k;
	int i;
	for (k = 0, i = 0; k < N; k++, i+=2) {
		x[i] /= (float)N;
		x[i+1] /= (float)N;
	}
}




