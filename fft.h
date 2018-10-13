//   fft.h - declaration of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#ifndef _FFT_H_
#define _FFT_H_

//   Include complex numbers header
#include "comp.h"

class CFFT
{
public:
	//   FORWARD FOURIER TRANSFORM
	//     Input  - input data
	//     Output - transform result
	//     N      - length of both input data and result
	static bool Forward(const comp *const Input, comp *const Output, const unsigned int N);

	//   FORWARD FOURIER TRANSFORM, INPLACE VERSION
	//     Data - both input data and output
	//     N    - length of input data
	static bool Forward(comp *const Data, const unsigned int N);

	//   INVERSE FOURIER TRANSFORM
	//     Input  - input data
	//     Output - transform result
	//     N      - length of both input data and result
	//     Scale  - if to scale result
	static bool Inverse(const comp *const Input, comp *const Output, const unsigned int N, const bool Scale = true);

	//   INVERSE FOURIER TRANSFORM, INPLACE VERSION
	//     Data  - both input data and output
	//     N     - length of both input data and result
	//     Scale - if to scale result
	static bool Inverse(comp *const Data, const unsigned int N, const bool Scale = true);

protected:
	//   Rearrange function and its inplace version
	static void Rearrange(const comp *const Input, comp *const Output, const unsigned int N);
	static void Rearrange(comp *const Data, const unsigned int N);

	//   FFT implementation
	static void Perform(comp *const Data, const unsigned int N, const bool Inverse = false);

	//   Scaling of inverse FFT result
	static void Scale(comp *const Data, const unsigned int N);

//   Conjugate comp array
    static void Conjugate(comp *const Data, const unsigned int N);
};

void fourier_transform(comp *array, unsigned int size);
void inverse_fourier_transform(comp *array, unsigned int size);

#endif
