
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
		val = (int16_t)(signal[i] * (pow(2,15) - 1));
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
		signal[i] = (val * 1.0) / (pow(2,15) - 1);
		if(signal[i] < -1.0)
			signal[i] = -1.0;

	}
	inputFile.close();
	return signal;
}


/*****************************************************************************
*
*    Function:     convolve
*
*    Description:  Convolves two signals, producing an output signal.
*                  The convolution is done in the time domain using the
*                  "Input Side Algorithm" (see Smith, p. 112-115).
*
*    Parameters:   x[] is the signal to be convolved
*                  N is the number of samples in the vector x[]
*                  h[] is the impulse response, which is convolved with x[]
*                  M is the number of samples in the vector h[]
*                  y[] is the output signal, the result of the convolution
*                  P is the number of samples in the vector y[].  P must
*                       equal N + M - 1
*
*	Source:			Leonard Manzara's convolution demo program
*****************************************************************************/
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
	int n, m;
	 
	 /*  Make sure the output buffer is the right size: P = N + M - 1  */
	if (P != (N + M - 1)) {
		printf("Output signal vector is the wrong size\n");
		printf("It is %-d, but should be %-d\n", P, (N + M - 1));
		printf("Aborting convolution\n");
		return;
	}
	
	/*  Clear the output buffer y[] to all zero values  */  
	for (n = 0; n < P; n++)
		y[n] = 0.0;
		
	  /*  Do the convolution  */
  /*  Outer loop:  process each input value x[n] in turn  */	
	for (n = 0; n < N; n++) {
		/*  Inner loop:  process x[n] with each sample of h[]  */
		for (m = 0; m < M; m++)
			y[n+m] += x[n] * h[m];
	}
}




