#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <ctime>

using namespace std;


char ChunkID[4];
int myChunkSize;
char myFormat[2];
char mySubChunk1ID[4];
int mySubChunk1Size;
int16_t audioFormat;
int16_t myChannels;
int sampleRate;	
int byteRate;
int16_t blockAlign;
int16_t bitsPerSample;
char subChunk2ID[4];
int subChunk2Size;
short* fileData;

int size;

float *wavRead(char *filename);
void convolve(float x[], int N, float h[], int M, float y[], int P);
void wavWrite(char *filename, float *signal, int signalSize);


int main(int argc, char **argv)
{
	std::clock_t start;
    double duration;

    start = std::clock();
	if (argc!= 4){ //check if we hve 3 command line arguments
		
		printf("Please enter an input file, and IR file and an output file name." );
	}
	else{
	char *inputName = argv[1]; //first arg is input file
	char *impulseResponseName = argv[2]; //IR file name
	char *outputName = argv[3]; //output file name
	
	float *inputSignal = NULL;
	float *IRSignal = NULL;
	float *outputSignal = NULL;
	int inputSignalSize;
	int IRSignalSize;
	int outputSignalSize;
	
	inputSignal = wavRead(inputName);
	inputSignalSize = size;
	printf("\n");
	IRSignal = wavRead(impulseResponseName);
	IRSignalSize = size;
	outputSignalSize = inputSignalSize + IRSignalSize - 1;
	outputSignal = new float[outputSignalSize];	
	convolve(inputSignal, inputSignalSize, IRSignal, IRSignalSize, outputSignal, outputSignalSize);
	wavWrite(outputName, outputSignal, outputSignalSize);
	}
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("The duration of the program is: ");
	std::cout<< duration <<" seconds";
	return 0;
	
}



	
float *wavRead (char *filename) 
{
	float *signal =NULL;
	ifstream inFile( filename, ios::in | ios::binary);
	inFile.seekg(ios::beg);
	inFile.read(ChunkID, 4);
	inFile.read((char*) &myChunkSize, 4); // read the ChunkSize
	inFile.read(myFormat, 4);
	inFile.read(mySubChunk1ID, 4);
	inFile.read((char*) &mySubChunk1Size, 4);
	inFile.read((char*) &audioFormat, 2);
	inFile.read((char*) &myChannels, 2);
	inFile.read((char*) &sampleRate, 4);
	inFile.read((char*) &byteRate, 4);
	inFile.read((char*) &blockAlign, 2);
	inFile.read((char*) &bitsPerSample, 2);
	if (mySubChunk1Size == 18) {
		inFile.seekg(2, ios::cur);
	}
	inFile.read(subChunk2ID, 4);
	inFile.read((char*)&subChunk2Size, 4);
	size = subChunk2Size / 2;
	
	short *data = new short[size];
	for (int i = 0; i < size; i++) {
		inFile.read((char *) &data[i], 2);
	}
	inFile.close();
	short sample;
	signal = new float[size];
	printf("Size: %d\n", size);
	for (int i = 0; i < size; i++) {
		sample = data[i];
		signal[i] = (sample * 1.0) / (pow(2.0, 15.0) - 1);
		if (signal[i] < -1.0)
			signal[i] = -1.0;
	}
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



void wavWrite(char *filename, float *signal, int sigSize)
{

	ofstream outFile(filename, ios::out | ios::binary);
	
	// File corrupted without hardcoded values
	char *ChunkID = "RIFF";
	char *format = "WAVE";
	// PCM = 18 was unnecessary
	mySubChunk1Size = 16;
	subChunk2Size = myChannels * sigSize * (bitsPerSample / 8);
	myChunkSize = subChunk2Size + 36;
	outFile.write(ChunkID, 4);
	outFile.write((char*) &myChunkSize, 4);
	outFile.write(format, 4);
	outFile.write(mySubChunk1ID, 4);
	outFile.write((char*) &mySubChunk1Size, 4);
	outFile.write((char*) &audioFormat, 2);
	outFile.write((char*) &myChannels, 2);
	outFile.write((char*) &sampleRate, 4);
	outFile.write((char*) &byteRate, 4);
	outFile.write((char*) &blockAlign, 2);
	outFile.write((char*) &bitsPerSample, 2);
	outFile.write(subChunk2ID, 4);
	outFile.write((char*)&subChunk2Size, 4);
	int16_t sample;
	
	// converting float to int between -2^15 to 2^15 - 1
	for(int i = 0; i < sigSize; i++)
	{
		sample = (int16_t)(signal[i] * (pow(2.0, 15.0) - 1));
		outFile.write((char*)&sample, 2);
	}
	outFile.close();
}
