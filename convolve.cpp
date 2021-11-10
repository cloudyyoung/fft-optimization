#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <cstring>
#include <string>
#include <complex>
#include <valarray>

// Constants
#define PI              	3.14159265358979
#define NEG_DOUBLE_PI       -2 * PI
#define TONE_FREQUENCY		440 // Frequency of tone to be created (A = 440 Hz)
#define SAMPLE_RATE     	44100.0 //  Standard sample rate in Hz
#define BITS_PER_SAMPLE		16 //  Standard sample size in bits
#define BYTES_PER_SAMPLE	(BITS_PER_SAMPLE / 8) // Standard sample size in bytes
#define MAX_SHORT_VALUE		32768 // Rescaling factor to convert between 16-bit shorts and doubles between -1 and 1
#define MONOPHONIC			1
#define STEREOPHONIC		2
#define FMT_OFFSET			12 // Offset of the fmt chunk in the WAV header

using namespace std;

typedef complex<double> Complex;
typedef valarray<Complex> ComplexArray;

class WaveFile {
public:
    string filename;
    bool succeed;

    int numberOfSample;

    // FMT
    int subchunk1ID;
    int subchunck1Size; // fmtSize

    int audioFormat;
    int channels;
    int sampleRate;
    int byteRate; // bytesPerSecond

    int blockAlign;
    int bitsPerSample;

    // DATA
    int subchunk2ID;
    int subchunck2Size;
    double* array;
    int arraySize;

    void read(string filename) {
        this->filename = filename;
        this->succeed = false;

        this->inputFileStream = fopen(filename.c_str(), "rb");
        if (this->inputFileStream == NULL) {
            return;
        }

        this->readHeader();

        if (this->numberOfSample <= 0 || this->bitsPerSample != BITS_PER_SAMPLE || this->sampleRate != SAMPLE_RATE) {
            return;
        }

        this->readData();

        fclose(this->inputFileStream);
        this->succeed = true;
    }

    void write(string filename) {
        this->filename = filename;
        this->succeed = false;

        this->inputFileStream = fopen(filename.c_str(), "wb");
        if (this->inputFileStream == NULL) {
            return;
        }

        this->writeHeader();
        this->writeData();

        fclose(this->inputFileStream);
        this->succeed = true;
    }

private:
    FILE* inputFileStream;

    void readHeader() {
        unsigned char buffer[64];
        fread(buffer, sizeof(unsigned char), 12, this->inputFileStream);

        this->subchunk1ID = nextIntLSB(); // subchunck1ID
        this->subchunck1Size = nextIntLSB(); // subchunck1Size / fmtSize

        this->audioFormat = nextShortLSB(); // audioFormat = 1
        this->channels = nextShortLSB();

        this->sampleRate = nextIntLSB();
        this->byteRate = nextIntLSB(); // bytesPerSecond

        this->blockAlign = nextShortLSB(); // blockAlign / frameSize
        this->bitsPerSample = nextShortLSB();

        fread(buffer, sizeof(unsigned char), this->subchunck1Size - 16, this->inputFileStream);

        this->subchunk2ID = nextIntLSB(); // subchunk2ID
        this->subchunck2Size = nextIntLSB();

        this->numberOfSample = this->subchunck2Size / (BYTES_PER_SAMPLE * (this->channels));

        cout << "sampleRate: " << this->sampleRate << ", byteRate: " << this->byteRate << endl;
        cout << "subchunk1 size: " << this->subchunck1Size << ", subchunk2 size: " << this->subchunck2Size << endl;
        cout << "channels: " << this->channels << endl;
        cout << "bits per sample: " << this->bitsPerSample << ", bytes per sample: " << BYTES_PER_SAMPLE << endl;
        cout << "Number of samples: " << this->numberOfSample << endl;
        cout << endl;
    }

    void readData() {
        int arraySize = this->numberOfSample * this->channels;

        cout << "array size: " << arraySize << endl;

        short* intArray = new short[arraySize];
        int count = fread(intArray, BYTES_PER_SAMPLE, this->numberOfSample, this->inputFileStream);

        int largest = 0;
        for (int t = 0; t < arraySize; t++) {
            if (intArray[t] > largest) {
                largest = intArray[t];
            }
        }
        // cout << "largest int: " << largest << endl;

        double* array = new double[arraySize];
        for (int t = 0; t < arraySize; t++) {
            array[t] = ((double)intArray[t]) / largest;
        }

        delete[] intArray;
        this->array = array;
        this->arraySize = arraySize;

        cout << endl;
    }

    void writeHeader() {
        int dataChunkSize = this->numberOfSample * BYTES_PER_SAMPLE;
        cout << "subchunk2 size: " << dataChunkSize << endl;

        // RIFF
        fputs("RIFF", this->inputFileStream);
        nextIntLSB(36 + this->subchunck2Size);
        fputs("WAVE", this->inputFileStream);

        // FMT
        fputs("fmt ", this->inputFileStream); // Subchunk1ID
        nextIntLSB(this->subchunck1Size); // Subchunk1Size
        nextShortLSB(this->audioFormat); // AudioFormat
        nextShortLSB(this->channels);
        nextIntLSB(this->sampleRate);
        nextIntLSB(this->byteRate); // bytesPerSecond
        nextShortLSB(this->channels * BYTES_PER_SAMPLE); // BlockAlign / frameSize
        nextShortLSB(BITS_PER_SAMPLE);

        // DATA
        fputs("data", this->inputFileStream); // Subchunk2ID
        nextIntLSB(this->subchunck2Size);
    }

    void writeData() {
        double largest = 1;
        int arraySize = this->arraySize;

        for (int t = 0; t < arraySize; t++) {
            if (this->array[t] > largest) {
                largest = this->array[t];
            }
        }

        short* intArray = new short[arraySize];
        for (int t = 0; t < arraySize; t++) {
            intArray[t] = (short)((this->array[t] / largest) * MAX_SHORT_VALUE);
        }

        fwrite(intArray, sizeof(short), arraySize, this->inputFileStream);
        delete[] intArray;
    }

    //reads an integer from the provided stream in little-endian form
    int nextIntLSB() {
        unsigned char array[4];
        fread(array, sizeof(unsigned char), 4, this->inputFileStream);

        int data = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);
        return data;
    }

    void nextIntLSB(int data) {
        unsigned char array[4];
        array[3] = (unsigned char)((data >> 24) & 0xFF);
        array[2] = (unsigned char)((data >> 16) & 0xFF);
        array[1] = (unsigned char)((data >> 8) & 0xFF);
        array[0] = (unsigned char)(data & 0xFF);

        fwrite(array, sizeof(unsigned char), 4, this->inputFileStream);
    }

    //reads a short integer from the provided stream in little-endian form
    short int nextShortLSB() {
        unsigned char array[2];
        fread(array, sizeof(unsigned char), 2, this->inputFileStream);

        int data = array[0] | (array[1] << 8);
        return data;
    }

    void nextShortLSB(short int data) {
        unsigned char array[2];
        array[1] = (unsigned char)((data >> 8) & 0xFF);
        array[0] = (unsigned char)(data & 0xFF);

        fwrite(array, sizeof(unsigned char), 2, this->inputFileStream);
    }
};

// https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2
unsigned long upper_power_of_two(unsigned long v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

// Cooley–Tukey FFT (in-place, divide-and-conquer)
void fft(ComplexArray& x) {
    const size_t n = x.size();
    if (n <= 1) return;

    int n_over_2 = int(n / 2);

    // divide
    ComplexArray even = x[slice(0, n_over_2, 2)];
    ComplexArray odd = x[slice(1, n_over_2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < n_over_2; ++k) {
        Complex t = polar(1.0, NEG_DOUBLE_PI * k / n) * odd[k];
        Complex even_k = even[k];
        
        x[k] = even_k + t;
        x[k + n_over_2] = even_k - t;
    }
}

// inverse fft (in-place)
void ifft(ComplexArray& x) {
    // conjugate the complex numbers
    x = x.apply(conj);

    // forward fft
    fft(x);

    // conjugate the complex numbers again
    x = x.apply(conj);

    // scale the numbers
    x /= x.size();
}

WaveFile convolution(WaveFile input, WaveFile IR) {
    WaveFile output;

    // Copy configurations    
    output.channels = IR.channels;
    output.sampleRate = input.sampleRate;
    output.byteRate = input.byteRate;
    output.bitsPerSample = input.bitsPerSample;

    // Fixed subchunk1 size and audio format
    output.subchunck1Size = 16;
    output.audioFormat = 1;

    // New subchunk2 size
    int outputSample = input.numberOfSample + IR.numberOfSample - 1; // m + n - 1
    output.numberOfSample = outputSample;
    output.subchunck2Size = outputSample * BYTES_PER_SAMPLE * output.channels;

    // New subchunk2 array data
    int outputArraySize = outputSample * output.channels;
    output.arraySize = outputArraySize;
    output.array = new double[outputArraySize];

    int complexArraySize = upper_power_of_two(outputSample);
    cout << "complex array per channel size: " << complexArraySize << endl;
    
    ComplexArray inputComplexArray;
    ComplexArray IRComplexArray;
    ComplexArray outputComplexArray;
    
    cout << "complex arrays build" << endl;
    // For each channel, FFT the input and output and multiply, and copy to output
    for (int r = 0; r < output.channels; r++) {
        cout << "IR complex array channel #" << r << endl;

        inputComplexArray.resize(complexArraySize, 0);
        IRComplexArray.resize(complexArraySize, 0);
        outputComplexArray.resize(complexArraySize, 0);

        // FFT input
        for (int t = 0; t < input.arraySize; t++) {
            inputComplexArray[t] = input.array[t];
        }
        fft(inputComplexArray);
        
        // FFT IR
        int IRArrayIndex = r;
        for (int t = 0; t < IR.numberOfSample; t++) {
            IRComplexArray[t] = IR.array[IRArrayIndex];
            IRArrayIndex += output.channels;
        }
        fft(IRComplexArray);

        // Multiplication
        outputComplexArray = inputComplexArray * IRComplexArray;
        ifft(outputComplexArray);
        cout << "output complex array ifft" << endl;

        // Copy real to output intertwined
        int outputIndex = r;
        for (int t = 0; t < outputSample; t++) {
            output.array[outputIndex] = outputComplexArray[t].real();
            outputIndex += output.channels;
        }
        cout << "output complex array to real" << endl << endl;
    }
    cout << "complex arrays fft" << endl;

    return output;
}


int main(int argc, char** argv) {
    if (argc != 4) {
        cout << "Format: " << "./convolve.o <inputFileName> <IRFileName> <outputFileName>" << endl;
        exit(-1);
    }


    string inputFileName(argv[1]);
    string IRFilename(argv[2]);
    string outputFilename(argv[3]);

    // cout << inputFileName << " " << IRFilename << " " << outputFilename << endl;

    WaveFile input;
    input.read(inputFileName);

    WaveFile IR;
    IR.read(IRFilename);

    WaveFile output = convolution(input, IR);
    output.write(outputFilename);

    cout << "output subchunk2 size: " << output.subchunck2Size << endl;
}