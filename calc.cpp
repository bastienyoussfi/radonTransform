#include "stdafx.h"
#include "Calc.h"
#include "LocalMatrixStudiable.h"
#include "num_recipes/fourier.h"
#include <vector>

#define PI 3.14159265359

//------------------------------------------------------------------------------------------------------------------------
// PreProcessing
//
// This function pre-processes the images in order to make the colocalization easier/better, it calls different methods according to what the user chooses
//
// inputStudiable : an interface to the source studiable (the surface from which a profile will be extracted)
// resultStudiable : return value, the processed image
//------------------------------------------------------------------------------------------------------------------------
void Calc::PreProcessing(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, int sigma)
{

    long z = inputStudiable.GetZMax();

    // Apply a test filter to the image
    if (testCalc) {
        for (int i = 0; i < inputStudiable.Height(); i++)
        {
            for (int j = 0; j < inputStudiable.Width(); j++)
            {
                long val = z - inputStudiable(i, j);
                resultStudiable(i, j) = val;
            }
        }
    }

    // Apply a Gaussian filter to the image
    if (gaussianFilterCalc) {
        // Gaussian kernel creation
        double** kernel = new double* [5];
        for (int k = 0; k < 5; k++) {
            kernel[k] = new double[5];
        }
        gaussianKernel(kernel, 5, sigma);
        //gaussian_kernel(kernel, 5, 1, 0, 1);

        convoluteIMatrix(resultStudiable, inputStudiable, kernel, 5);

        // Frees the kernel memory ---
        for (int l = 0; l < 5; ++l) {
            delete[] kernel[l];
        }

        delete[] kernel;
        // ---------------------------

        return;
    }

    // Apply a mean filter to the image
    if (meanFilterCalc) {
        double** kernel = new double* [7];
        for (int k = 0; k < 7; k++) {
            kernel[k] = new double[7];
        }
        meanKernel(kernel, 7);

        convoluteIMatrix(resultStudiable, inputStudiable, kernel, 7);

        // Frees the kernel memory ---
        for (int l = 0; l < 7; ++l) {
            delete[] kernel[l];
        }

        delete[] kernel;
        // ---------------------------

        return;
    }

    // Apply a ZNCC filter to the image.
    if (znccCalc) {
        return;
    }

    // Apply an outline detection filter to the image.
    if (outlineDetectionCalc) {
        return;
    }

    if (!testCalc && !gaussianFilterCalc && !meanFilterCalc && !znccCalc && !outlineDetectionCalc) {
        for (int i = 0; i < inputStudiable.Height(); i++)
        {
            for (int j = 0; j < inputStudiable.Width(); j++)
            {
                resultStudiable(i, j) = inputStudiable(i, j);
            }
        }
    }
    resultStudiable.SaveValues();
}


//------------------------------------------------------------------------------------------------------------------------
// Creates a kernel for a gaussian filter
//------------------------------------------------------------------------------------------------------------------------
void Calc::gaussianKernel(double** kernel, long n, double sig)
{
    double sigma = sig;
    double sum = 0.0;
    long m = (long)n / 2;

    for (int x = -m; x <= m; x++) {
        for (int y = -m; y <= m; y++) {
            kernel[x+m][y+m] = exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * PI * sigma * sigma);
            sum += kernel[x+m][y+m];
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            kernel[i][j] /= sum;
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Creates a kernel for a mean filter
//------------------------------------------------------------------------------------------------------------------------
void Calc::meanKernel(double** kernel, long n)
{

    long m = (long)(n / 2);

    for (int x = -m; x <= m; x++) {
        for (int y = -m; y <= m; y++) {
            kernel[x + m][y + m] = 1;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            kernel[i][j] /= n*n;
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Regular convolution function.
//------------------------------------------------------------------------------------------------------------------------
void Calc::convoluteIMatrix(CLocalMatrixStudiable& resultStudiable, const CLocalMatrixStudiable& inputStudiable, double** kernel, long n)
{
    double convolute = 0; // This holds the convolution results for an index.
    long edges = (long)(n/2);
    int x, y; // Used for input matrix index
    long M = inputStudiable.Height();
    long N = inputStudiable.Width();

    // Fill output matrix: rows and columns are i and j respectively
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            x = i;
            y = j;

            // Kernel rows and columns are k and l respectively
            for (int k = -edges; k <= edges; k++)
            {
                for (int l = -edges; l <= edges; l++)
                {
                    // Convolute here.
                    if (x + k >= 0 && x + k < N && y + l >= 0 && y + l < M)
                    {
                        convolute += kernel[k + edges][l + edges] * inputStudiable(y + l, x + k);
                    }
                    y++; // Move right.
                }
                x++; // Move down.
                y = j; // Restart column position
            }
            resultStudiable(j, i) = (long) convolute; // Add result to output matrix.
            convolute = 0; // Needed before we move on to the next index.
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Set every boolean to false so there is no pre-processing of the image.
//------------------------------------------------------------------------------------------------------------------------
void Calc::reset() 
{
    testCalc = false;
    gaussianFilterCalc = false;
    meanFilterCalc = false;
    znccCalc = false;
    outlineDetectionCalc = false;
}


//------------------------------------------------------------------------------------------------------------------------
// Rotate
// 
// Rotate an image of theta dregrees
//
// inputStudiable : an interface to the source studiable
// resultStudiable : return value, the rotated image
//------------------------------------------------------------------------------------------------------------------------
void Calc::Rotate(double theta, const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable)
{
    int a, b, c, d;
    double prevX, prevY, intPrevX, intPrevY, dIntPrevX, dIntPrevY;
    double res = 0;
    double nm = LONG_MIN;
    if (inputStudiable.Width() % 2 == 0) 
    {
        a = (int) floor(inputStudiable.Width() / 2);
        b = (int) floor(inputStudiable.Height() / 2);
    }
    else
    {
        a = (int) floor(inputStudiable.Width() / 2);
        b = (int) floor(inputStudiable.Height() / 2);
    }
    if (inputStudiable.Width() % 2 == 0)
    {
        c = (int) floor(resultStudiable.Width() / 2);
        d = (int) floor(resultStudiable.Height() / 2);
    }
    else
    {
        c = (int) floor(resultStudiable.Width() / 2) + 1;
        d = (int) floor(resultStudiable.Height() / 2) + 1;
    }
    
    for (int i = -c; i < c-1; i++)
    {
        for (int j = -d; j < d-1; j++)
        {
            res = 0;
            prevX = cos(theta) * i - sin(theta) * j;
            prevY = sin(theta) * i + cos(theta) * j;
            intPrevX = floor(prevX);
            intPrevY = floor(prevY);
            dIntPrevX = prevX - intPrevX;
            dIntPrevY = prevY - intPrevY;
            if (abs(intPrevX) < (double) a-1 && abs(intPrevY) < (double) b-1)
            {
                if (dIntPrevX <= 0.5)
                {
                    if (dIntPrevY <= 0.5) 
                    {
                        res = res + (dIntPrevX + 0.5) * (dIntPrevY + 0.5) * inputStudiable((long) intPrevX + a, (long) intPrevY+b);
                        res = res + (dIntPrevX + 0.5) * (0.5 - dIntPrevY) * inputStudiable((long) intPrevX + a, (long) intPrevY - 1 + b);
                        res = res + (0.5 - dIntPrevX) * (0.5 - dIntPrevY) * inputStudiable((long) intPrevX + a - 1, (long) intPrevY - 1 + b);
                        res = res + (0.5 - dIntPrevX) * (dIntPrevY + 0.5) * inputStudiable((long) intPrevX + a - 1, (long) intPrevY + b);
                    }
                    else
                    {
                        res = res + (dIntPrevX + 0.5) * (1.5 - dIntPrevY) * inputStudiable((long) intPrevX + a, (long) intPrevY + b);
                        res = res + (0.5 - dIntPrevX) * (1.5 - dIntPrevY) * inputStudiable((long) intPrevX - 1 + a, (long) intPrevY + b);
                        res = res + (dIntPrevX + 0.5) * (dIntPrevY - 0.5) * inputStudiable((long) intPrevX + a, (long) intPrevY + 1 + b);
                        res = res + (0.5 - dIntPrevX) * (dIntPrevY - 0.5) * inputStudiable((long) intPrevX - 1 + a, (long) intPrevY + 1 + b);
                    }
                }
                else 
                {
                    if (dIntPrevY <= 0.5) 
                    {
                        res = res + (1.5 - dIntPrevX) * (dIntPrevY + 0.5) * inputStudiable((long) intPrevX + a, (long) intPrevY + b);
                        res = res + (dIntPrevX - 0.5) * (dIntPrevY + 0.5) * inputStudiable((long) intPrevX + 1 + a, (long) intPrevY + b);
                        res = res + (dIntPrevX + 0.5) * (0.5 - dIntPrevY) * inputStudiable((long)intPrevX + 1 + a, (long) intPrevY - 1 + b);
                        res = res + (1.5 - dIntPrevX) * (0.5 - dIntPrevY) * inputStudiable((long) intPrevX + a, (long) intPrevY - 1 + b);
                    }
                    else 
                    {
                        res = res + (1.5 - dIntPrevX) * (1.5 - dIntPrevY) * inputStudiable((long) intPrevX + a, (long) intPrevY + b);
                        res = res + (dIntPrevX - 0.5) * (1.5 - dIntPrevY) * inputStudiable((long) intPrevX + 1 + a, (long)intPrevY + b);
                        res = res + (dIntPrevX - 0.5) * (dIntPrevY - 0.5) * inputStudiable((long) intPrevX + 1 + a, (long) intPrevY + 1 + b);
                        res = res + (1.5 - dIntPrevX) * (dIntPrevY - 0.5) * inputStudiable((long) intPrevX + a, (long) intPrevY + 1 + b);
                    }
                }
                resultStudiable(i+c, j+d) = (long) res;
            }
            else
            {
                resultStudiable(i + c, j + d) = (long) nm;
            }
        }
    }
    resultStudiable.SaveValues();
}


//------------------------------------------------------------------------------------------------------------------------
// Translate
// 
// Translate an image of dx and dy
//
// inputStudiable : an interface to the source studiable
// resultStudiable : return value, the translated image
//------------------------------------------------------------------------------------------------------------------------
void Calc::Translate(double dx, double dy, IDispatch* pDispStudiable)
{
    IMatrixStudiablePtr pStudiable(pDispStudiable);
    pStudiable->AxisOffsetInLengthUnit[X_AXIS] = pStudiable->AxisOffsetInLengthUnit[X_AXIS] + dx;
    pStudiable->AxisOffsetInLengthUnit[Y_AXIS] = pStudiable->AxisOffsetInLengthUnit[Y_AXIS] - dy;
}


//------------------------------------------------------------------------------------------------------------------------
// Scale
// 
// Scale an image with an alpha factor
//
// inputStudiable : an interface to the source studiable
// resultStudiable : return value, the scaled image
//------------------------------------------------------------------------------------------------------------------------
void Calc::Scale(double alpha, IDispatch* pDispStudiable)
{
    IMatrixStudiablePtr pStudiable(pDispStudiable);
    //double x = pStudiable->AxisSpacing[X_AXIS];
    pStudiable->AxisSpacing[X_AXIS] = pStudiable->AxisSpacing[X_AXIS] * alpha;
    //double y = pStudiable->AxisSpacing[X_AXIS];
    pStudiable->AxisSpacing[Y_AXIS] = pStudiable->AxisSpacing[Y_AXIS] * alpha;
}


//------------------------------------------------------------------------------------------------------------------------
// Chirp-z transform horizontal 1D
// 
// Chirp-z transform 1D written in C++
//------------------------------------------------------------------------------------------------------------------------
/**/void Calc::Chirpz1DHorizontal(double A, double W, long M, const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, long index) 
{
    Complex A_comp = Complex(A);
    Complex W_comp = Complex(W);

    Int signFFT = 1;
    Int signIFFT = -1;

    long N = inputStudiable.Width();

    // We choose L, the smallest integer greater than or equal to N + M - 1
    long L = (long) floor(pow(2, ceil(log2(M + N - 1))));

    VecDoub res(N);

    for (int i = 0; i < N; i++)
    {
        res[i] = inputStudiable(index, i);
    }

    // We compute the L yn values from xn
    VecComplex Y(L);
    for (int i = 0; i < L; i++) 
    {
        if (i < N) 
        {
            Y[i] = pow(A_comp, -i) * pow(W_comp, (i*i)/ 2) * res[i];
        }
        else
        {
            Y[i] = 0;
        }
    }

    // We compute the fft of Y
    four1(Y, signFFT);

    // We compute the L vn values
    VecComplex V(L);
    for (int i = 0; i < L; i++) 
    {
        if (i <= M - 1) 
        {
            V[i] = pow(W_comp, (-i*i)/2);
        }
        else if (L - N + 1 <= i && i < L) 
        {
            V[i] = pow(W_comp, -((L-i)*(L-i))/2);
        }
        else
        {
            V[i] = 0;
        }
    }

    // We compute the fft of V
    four1(V, signFFT);

    // We multiply Y and V
    VecDoub GReal(L);
    VecDoub GImag(L);
    VecComplex G(L);
    for (int i = 0; i < L; i++)
    {
        GReal[i] = (Y[i].real() * V[i].real()) - (Y[i].imag() * V[i].imag());
        GImag[i] = (Y[i].real() * V[i].imag()) + (V[i].imag() * Y[i].imag());
    }

    for (int i = 0; i < L; i++)
    {
        G[i] = Complex(GReal[i], GImag[i]);
    }

    // And we compute its ifft
    VecComplex g(L);
    four1(G, signIFFT); // TODO: add the normalization factor
    for (int i = 0; i < L; i++)
    {
        g[i] = Complex(G[i].real()*2/N, G[i].imag()*2/N);
    }

    // Finally we multiply by the factor to get X
    VecComplex X(L);
    for (int i = 0; i < M; i++)
    {
        X[i] = g[i]*pow(W_comp, (i * i) / 2);
    }

    // We put the result in resultStudiable
    double realRes;
    for (int i = 0; i < N; i++)
    {
        realRes = X[i].real();
        resultStudiable(index, i) = (long)realRes; // TODO: Faire la conversion entier/réel pour que Mountains ait la bonne valeur
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Chirp-z transform vertical 1D
// 
// Chirp-z transform 1D written in C++
//------------------------------------------------------------------------------------------------------------------------
/**/void Calc::Chirpz1DVertical(double A, double W, long M, const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, long index)
{
    Complex A_comp = Complex(A);
    Complex W_comp = Complex(W);

    Int signFFT = 1;
    Int signIFFT = -1;

    long N = inputStudiable.Width();

    // We choose L, the smallest integer greater than or equal to N + M - 1
    long L = (long)floor(pow(2, ceil(log2(M + N - 1))));

    VecDoub res(N);

    for (int i = 0; i < N; i++)
    {
        res[i] = inputStudiable(i, index);
    }

    // We compute the L yn values from xn
    VecComplex Y(L);
    for (int i = 0; i < L; i++)
    {
        if (i < N)
        {
            Y[i] = pow(A_comp, -i) * pow(W_comp, (i * i) / 2) * res[i];
        }
        else
        {
            Y[i] = 0;
        }
    }

    // We compute the fft of Y
    four1(Y, signFFT);

    // We compute the L vn values
    VecComplex V(L);
    for (int i = 0; i < L; i++)
    {
        if (i <= M - 1)
        {
            V[i] = pow(W_comp, (-i * i) / 2);
        }
        else if (L - N + 1 <= i && i < L)
        {
            V[i] = pow(W_comp, -((L - i) * (L - i)) / 2);
        }
        else
        {
            V[i] = 0;
        }
    }

    // We compute the fft of V
    four1(V, signFFT);

    // We multiply Y and V
    VecDoub GReal(L);
    VecDoub GImag(L);
    VecComplex G(L);
    for (int i = 0; i < L; i++)
    {
        GReal[i] = (Y[i].real() * V[i].real()) - (Y[i].imag() * V[i].imag());
        GImag[i] = (Y[i].real() * V[i].imag()) + (V[i].imag() * Y[i].imag());
    }

    for (int i = 0; i < L; i++)
    {
        G[i] = Complex(GReal[i], GImag[i]);
    }

    // And we compute its ifft
    VecComplex g(L);
    four1(G, signIFFT); // TODO: add the normalization factor
    for (int i = 0; i < L; i++)
    {
        g[i] = Complex(G[i].real() * 2 / N, G[i].imag() * 2 / N);
    }

    // Finally we multiply by the factor to get X
    VecComplex X(L);
    for (int i = 0; i < M; i++)
    {
        X[i] = g[i] * pow(W_comp, (i * i) / 2);
    }

    // We put the result in resultStudiable
    double realRes;
    for (int i = 0; i < N; i++)
    {
        realRes = X[i].real();
        resultStudiable(i, index) = (long)realRes; // TODO: Faire la conversion entier/réel pour que Mountains ait la bonne valeur
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Chirp-z transform 2D
// 
// Chirp-z transform 2D written in C++
//------------------------------------------------------------------------------------------------------------------------
void Calc::Chirpz2D(double A, double W, long M, const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, long index)
{
    long n = inputStudiable.Width();
    long m = inputStudiable.Height();

    IMatrixStudiablePtr pTemp;

    IStudiableFactoryPtr factorytemp;

    HRESULT hrtemp = factorytemp.CreateInstance(__uuidof(CoStudiableFactory));
    if (FAILED(hrtemp))
        return;

    pTemp = factorytemp->CreateMatrixStudiable(kDSStudiableSurface, n, m, 1, 1);
    if (pTemp == nullptr)
        return;

    CLocalMatrixStudiable temp(pTemp, TRUE, TRUE);

    // We iterate the 1D horizontal transform over all the lines of the image
    for (int i = 0; i < n; i++) {
        Chirpz1DHorizontal(A, W, M, inputStudiable, temp, i);
    }
    // We iterate the 1D vertical transform over all the columns of the image
    for (int j = 0; j < m; j++) {
        Chirpz1DVertical(A, W, M, temp, resultStudiable, j);
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Radon transform
// 
// Radon transform transform 2D written in C++, the PI version (180°)
//------------------------------------------------------------------------------------------------------------------------
/*void Calc::Radon(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, double thetaoffset)
{
    double deltax, deltay, costheta, sintheta, deltarho, rhooffset, alpha, beta, sum, thetarad, thetazero, x;
    long M, N, xmin, xmax, ymin, ymax, mmin, mmax, nmin, nmax, rhomax, rhomin, theta, offset, diag, diagmin, diagmax;

    M = inputStudiable.Width();
    N = inputStudiable.Height();

    deltax = 1;
    deltay = 1;
    deltarho = 1;
    offset = 0;
    thetazero = thetaoffset;
    diag = (long) floor(sqrt(2) * M);
    diagmin = -diag / 2;
    diagmax = diag / 2;
    x = 0;

    if (inputStudiable.Width() % 2 == 0)
    {
        xmax = inputStudiable.Width() / 2;
        ymax = inputStudiable.Height() / 2;
        xmin = -xmax + 1;
        ymin = -ymax + 1;
        offset = 1;
    }
    else
    {
        xmax = inputStudiable.Width() / 2;
        ymax = inputStudiable.Height() / 2;
        xmin = -xmax;
        ymin = -ymax;
    }

    for (theta = (long) thetazero; theta < 180 + thetazero; theta++)
    {
        thetarad = theta * PI / 180;
        costheta = cos(thetarad);
        sintheta = sin(thetarad);
        rhooffset = xmin * (costheta + sintheta);
        if (sintheta > 1 / sqrt(2))
        {
            rhomax = sqrt(2) * xmax / sintheta;
            rhomin = -rhomax + offset;
            alpha = -costheta / sintheta;
            for (long rho = -diag / 2; rho < diag / 2; rho += (long)deltarho)
            {
                if (rho >= rhomin && rho <= rhomax)
                {
                    beta = (rho - rhooffset) / (deltax * sintheta);
                    if (alpha > 0)
                    {
                        mmin = (long)max(0, ceil(-(beta + 0.5) / alpha));
                        mmax = (long)min(M - 1, floor((N - 0.5 - beta) / alpha));
                    }
                    else if (alpha < 0)
                    {
                        mmin = (long)max(0, ceil((N - 0.5 - beta) / alpha));
                        mmax = (long)min(M - 1, floor(-(beta + 0.5) / alpha));
                    }
                    else
                    {
                        mmin = 0;
                        mmax = xmax - xmin;
                    }
                    sum = 0;
                    for (int m = mmin; m <= mmax; m++)
                    {
                        if (inputStudiable((long)round(alpha * m + beta), m) != LONG_MIN)
                        {
                            sum += inputStudiable((long)round(alpha * m + beta), m);
                        }
                    }
                    resultStudiable(rho - diagmin, theta - thetazero) = (long)(deltax * sum / sintheta);
                }
                else
                {
                    resultStudiable(rho - diagmin, theta - thetazero) = 0;
                }
            }
        }
        else
        {
            rhomax = sqrt(2) * abs(xmax / costheta);
            rhomin = -rhomax + offset;
            alpha = -sintheta / costheta;
            for (long rho = -diag / 2; rho < diag / 2; rho += (long)deltarho)
            {
                if (rho >= rhomin && rho <= rhomax)
                {
                    beta = (rho - rhooffset) / (deltax * costheta);
                    if (alpha > 0)
                    {
                        nmin = (long)max(0, ceil(-(beta + 0.5) / alpha));
                        nmax = (long)min(M - 1, floor((N - 0.5 - beta) / alpha));
                    }
                    else if (alpha < 0)
                    {
                        nmin = (long)max(0, ceil((N - 0.5 - beta) / alpha));
                        nmax = (long)min(M - 1, floor(-(beta + 0.5) / alpha));
                    }
                    else
                    {
                        nmin = 0;
                        nmax = xmax - xmin;
                    }
                    sum = 0;
                    for (int n = nmin; n < nmax; n++)
                    {
                        if ((long)round(alpha * n + beta) >= 0)
                        {
                            if (inputStudiable(n, (long)round(alpha * n + beta)) != LONG_MIN)
                            {
                                sum += inputStudiable(n, (long)round(alpha * n + beta));
                            }
                        }
                    }
                    resultStudiable(rho - diagmin, theta - thetazero) = (long)(deltax * sum / abs(costheta));
                }
                else
                {
                    resultStudiable(rho - diagmin, theta - thetazero) = 0;
                }
            }
        }
    }
}*/


//------------------------------------------------------------------------------------------------------------------------
// Radon transform
// 
// Radon transform transform 2D written in C++, the 2PI version (360°)
//------------------------------------------------------------------------------------------------------------------------
void Calc::RadonV2(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, double thetaoffset, long coeffNorm)
{
    double deltax, deltay, costheta, sintheta, deltarho, rhooffset, alpha, beta, sum, thetarad, thetazero, x;
    long M, N, xmin, xmax, ymin, ymax, mmin, mmax, nmin, nmax, rhomax, rhomin, theta, offset, diag, diagmin, diagmax;

    M = inputStudiable.Width();
    N = inputStudiable.Height();

    if (false)
    {
        coeffNorm = 1;
    }

    deltax = 1;
    deltay = 1;
    deltarho = 1;
    offset = 0;
    thetazero = thetaoffset;
    diag = (long)floor(sqrt(2) * M);
    diagmin = -diag / 2;
    diagmax = diag / 2;
    x = 0;

    if (inputStudiable.Width() % 2 == 0)
    {
        xmax = inputStudiable.Width() / 2;
        ymax = inputStudiable.Height() / 2;
        xmin = -xmax + 1;
        ymin = -ymax + 1;
        offset = 1;
    }
    else
    {
        xmax = inputStudiable.Width() / 2;
        ymax = inputStudiable.Height() / 2;
        xmin = -xmax;
        ymin = -ymax;
    }

    for (theta = (long)thetazero; theta < 180 + thetazero; theta++)
    {
        thetarad = theta * PI / 180;
        if (theta%180 == 0)
        {
            costheta = 1;
            sintheta = 0;
        }
        else if ((theta+90)%180 == 0)
        {
            costheta = 0;
            sintheta = 1;
        }
        else
        { 
            costheta = cos(thetarad);
            sintheta = sin(thetarad);
        }
        rhooffset = xmin * (costheta + sintheta);
        if (sintheta > 1 / sqrt(2))
        {
            rhomax = sintheta!=1?(long)(sqrt(2)*abs(xmax/sintheta)):xmax;
            rhomin = -rhomax+offset;
            alpha = sintheta!=0?-costheta/sintheta:0;
            for (long rho = -diag/2; rho < diag/2; rho += (long)deltarho)
            {
                if (rho >= rhomin && rho <= rhomax)
                {
                    beta = (rho - rhooffset) / (deltax * sintheta);
                    if (alpha > 0)
                    {
                        mmin = (long)max(0, ceil(-(beta + 0.5) / alpha));
                        mmax = (long)min(M - 1, floor((N - 0.5 - beta) / alpha));
                    }
                    else if (alpha < 0)
                    {
                        mmin = (long)max(0, ceil((N - 0.5 - beta) / alpha));
                        mmax = (long)min(M - 1, floor(-(beta + 0.5) / alpha));
                    }
                    else
                    {
                        mmin = 0;
                        mmax = xmax - xmin;
                    }
                    sum = 0;
                    for (int m = mmin; m <= mmax; m++)
                    {
                        if ((long)round(alpha * m + beta) >= 0)
                        {
                            if (inputStudiable((long)round(alpha * m + beta), m) != LONG_MIN)
                            {
                                sum += inputStudiable((long)round(alpha * m + beta), m);//coeffNorm;
                            }
                        }
                    }
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = (long)(deltax * sum / sintheta);
                }
                else
                {
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = 0;
                }
            }
        }
        else
        {
            rhomax = costheta==0||costheta==1?xmax:(long)(sqrt(2) * abs(xmax / costheta));
            rhomin = -rhomax + offset;
            alpha = costheta!=0?-sintheta/costheta:0;
            for (long rho = -diag/2; rho < diag/2; rho += (long)deltarho)
            {
                if (rho >= rhomin && rho <= rhomax)
                {
                    beta = (rho - rhooffset) / (deltax * costheta);
                    if (alpha > 0)
                    {
                        nmin = (long)max(0, ceil(-(beta + 0.5) / alpha));
                        nmax = (long)min(M - 1, floor((N - 0.5 - beta) / alpha));
                    }
                    else if (alpha < 0)
                    {
                        nmin = (long)max(0, ceil((N - 0.5 - beta) / alpha));
                        nmax = (long)min(M - 1, floor(-(beta + 0.5) / alpha));
                    }
                    else
                    {
                        nmin = 0;
                        nmax = xmax - xmin;
                    }
                    sum = 0;
                    for (int n = nmin; n < nmax; n++)
                    {
                        if ((long)round(alpha * n + beta) >= 0)
                        {
                            if (inputStudiable(n, (long)round(alpha * n + beta)) != LONG_MIN)
                            {
                                sum += inputStudiable(n, (long)round(alpha * n + beta));//coeffNorm;
                            }
                        }
                    }
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = (long)(deltax * sum / abs(costheta));
                }
                else
                {
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = 0;
                }
            }
        }
    }
    for (theta = (long)thetazero + 180; theta < thetazero + 360; theta++)
    {
        for (long rho = -diag / 2; rho < diag / 2; rho += (long)deltarho)
        {
            resultStudiable(- rho - diagmin - 1, theta - (long)thetazero) = resultStudiable(rho - diagmin, theta - (long)thetazero - 180);
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Radon transform
// 
// Radon transform transform 2D written in C++, the 2PI version (360°), dealing with rectangle shaped studiables
//------------------------------------------------------------------------------------------------------------------------
void Calc::RadonV3(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, double thetaoffset)
{
    double deltax, deltay, costheta, sintheta, deltarho, rhooffset, alpha, beta, sum, thetarad, thetazero, x;
    long M, N, xmin, xmax, ymin, ymax, mmin, mmax, nmin, nmax, rhomax, rhomin, theta, offset, diag, diagmin, diagmax;

    M = inputStudiable.Width();
    N = inputStudiable.Height();

    deltax = 1;
    deltay = 1;
    deltarho = 1;
    offset = 0;
    thetazero = thetaoffset;
    diag = (long)floor(sqrt(2) * M);
    diagmin = -diag / 2;
    diagmax = diag / 2;
    x = 0;

    if (inputStudiable.Width() % 2 == 0)
    {
        xmax = inputStudiable.Width() / 2;
        ymax = inputStudiable.Height() / 2;
        xmin = -xmax + 1;
        ymin = -ymax + 1;
        offset = 1;
    }
    else
    {
        xmax = inputStudiable.Width() / 2;
        ymax = inputStudiable.Height() / 2;
        xmin = -xmax;
        ymin = -ymax;
    }

    for (theta = (long)thetazero; theta < 180 + thetazero; theta++)
    {
        thetarad = theta * PI / 180;
        if (theta % 180 == 0)
        {
            costheta = 1;
            sintheta = 0;
        }
        else if ((theta + 90) % 180 == 0)
        {
            costheta = 0;
            sintheta = 1;
        }
        else
        {
            costheta = cos(thetarad);
            sintheta = sin(thetarad);
        }
        rhooffset = xmin * (costheta + sintheta);
        if (sintheta > 1 / sqrt(2))
        {
            rhomax = sintheta != 1 ? (long)(sqrt(2) * abs(xmax / sintheta)) : xmax;
            rhomin = -rhomax + offset;
            alpha = sintheta != 0 ? -costheta / sintheta : 0;
            for (long rho = -diag / 2; rho < diag / 2; rho += (long)deltarho)
            {
                if (rho >= rhomin && rho <= rhomax)
                {
                    beta = (rho - rhooffset) / (deltax * sintheta);
                    if (alpha > 0)
                    {
                        mmin = (long)max(0, ceil(-(beta + 0.5) / alpha));
                        mmax = (long)min(M - 1, floor((N - 0.5 - beta) / alpha));
                    }
                    else if (alpha < 0)
                    {
                        mmin = (long)max(0, ceil((N - 0.5 - beta) / alpha));
                        mmax = (long)min(M - 1, floor(-(beta + 0.5) / alpha));
                    }
                    else
                    {
                        mmin = 0;
                        mmax = xmax - xmin;
                    }
                    sum = 0;
                    for (int m = mmin; m <= mmax; m++)
                    {
                        if ((long)round(alpha * m + beta) >= 0)
                        {
                            if (inputStudiable((long)round(alpha * m + beta), m) != LONG_MIN)
                            {
                                sum += inputStudiable((long)round(alpha * m + beta), m);
                            }
                        }
                    }
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = (long)(deltax * sum / sintheta);
                }
                else
                {
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = 0;
                }
            }
        }
        else
        {
            rhomax = costheta == 0 || costheta == 1 ? xmax : (long)(sqrt(2) * abs(xmax / costheta));
            rhomin = -rhomax + offset;
            alpha = costheta != 0 ? -sintheta / costheta : 0;
            for (long rho = -diag / 2; rho < diag / 2; rho += (long)deltarho)
            {
                if (rho >= rhomin && rho <= rhomax)
                {
                    beta = (rho - rhooffset) / (deltax * costheta);
                    if (alpha > 0)
                    {
                        nmin = (long)max(0, ceil(-(beta + 0.5) / alpha));
                        nmax = (long)min(M - 1, floor((N - 0.5 - beta) / alpha));
                    }
                    else if (alpha < 0)
                    {
                        nmin = (long)max(0, ceil((N - 0.5 - beta) / alpha));
                        nmax = (long)min(M - 1, floor(-(beta + 0.5) / alpha));
                    }
                    else
                    {
                        nmin = 0;
                        nmax = xmax - xmin;
                    }
                    sum = 0;
                    for (int n = nmin; n < nmax; n++)
                    {
                        if ((long)round(alpha * n + beta) >= 0)
                        {
                            if (inputStudiable(n, (long)round(alpha * n + beta)) != LONG_MIN)
                            {
                                sum += inputStudiable(n, (long)round(alpha * n + beta));
                            }
                        }
                    }
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = (long)(deltax * sum / abs(costheta));
                }
                else
                {
                    resultStudiable(rho - diagmin, theta - (long)thetazero) = 0;
                }
            }
        }
    }
    for (theta = (long)thetazero + 180; theta < thetazero + 360; theta++)
    {
        for (long rho = -diag / 2; rho < diag / 2; rho += (long)deltarho)
        {
            resultStudiable(-rho - diagmin - 1, theta - (long)thetazero) = resultStudiable(rho - diagmin, theta - (long)thetazero - 180);
        }
    }
}


void Calc::matchingOperator(const CLocalMatrixStudiable& inputStudiable1, const CLocalMatrixStudiable& inputStudiable2)
{

}


//------------------------------------------------------------------------------------------------------------------------
// HRT Histogram
// 
// Computes the HRT descriptor of the Radon transform
// inputStudiable : A radon transform
// resultStudiable : The matrix of frequencies computed on the Radon transform for the angle parameter (size has to be rhomax - 180)
//------------------------------------------------------------------------------------------------------------------------
void Calc::histogram(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable, long numBins)
{
    long N = numBins;
    long M = inputStudiable.Height();
    long zMax = inputStudiable.GetZMaxLocal();
    long zMin = inputStudiable.GetZMinLocal();
    double range = zMax - zMin;
    double rangeBin = range / N;
long offset = 0;
long delta = (zMax - zMin) / N;
long deltasurdeux = delta / 2;

if (zMin > 0)
{
    offset = -zMin;
}
else
{
    offset = zMin;
}

for (int j = 0; j < 180; j++)
{
    for (int i = 0; i < M; i += 1)
    {
        long a = inputStudiable(i, j);
        if (/*inputStudiable(i, j) != 0 && */inputStudiable(i, j) != LONG_MIN)
        {
            double intVal1 = inputStudiable(i, j) - offset;
            long index = intVal1 / rangeBin;
            if (index > 0 && index < N)
            {
                resultStudiable(index, j) += 1;
            }
        }
    }
}

/*for (int j = 0; j < 180; j++)
{
    for (int i = 0; i < M; i += 1)
    {
        if (resultStudiable(i, j) == 0)
        {
            resultStudiable(i, j) = LONG_MIN;
        }
    }
}*/
}


//------------------------------------------------------------------------------------------------------------------------
// Angle correlation computation
// 
// Compute the angle correlation between two studiables using their descriptors
// inputStudiable1 : The first studiable
// inputStudiable2 : The second studiable
// resultStudiable : The angle correlation matrix
//------------------------------------------------------------------------------------------------------------------------
void Calc::angleCorrelation(const CLocalMatrixStudiable& inputStudiable1, const CLocalMatrixStudiable& inputStudiable2, CLocalMatrixStudiable& resultStudiable)
{
    long nBins = inputStudiable1.Height();
    for (int i = 0; i < 180; i++)
    {
        for (int j = 0; j < 180; j++)
        {
            long sum = 0;
            for (int k = 0; k < nBins; k++)
            {
                sum += inputStudiable1(k, i) * inputStudiable2(k, j);
            }
            resultStudiable(i, j) = sum;
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Cost matrix computation
// 
// Computes the cost matrix of a matrix, here used on angle correlation matrix
// inputStudiable : An angle correlation matrix (size has to be 180-180)
// resultStudiable : The cost matrix corresponding to the input (size has to be 180-180)
//------------------------------------------------------------------------------------------------------------------------
void Calc::costMatrix(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable)
{
    for (int i = 0; i < 180; i++)
    {
        for (int j = 0; j < 180; j++)
        {
            if (inputStudiable(i, j) == 0)
            {
                resultStudiable(i, j) = 0;
            }
            else
            {
                resultStudiable(i, j) = 2000000 / inputStudiable(i, j);
            }
        }
    }
}


void Calc::normalizeCostMatrix(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable)
{
    long min = LONG_MAX;
    long indexRow = 0;
    long indexCol = 0;
    for (int i = 0; i < 180; i++)
    {
        for (int j = 0; j < 180; j++)
        {
            if (inputStudiable(i, j) < min)
            {
                min = inputStudiable(i, j);
                indexRow = i;
                indexCol = j;
            }
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Padding function
// 
// Turns a rectangle studiable to a square studiable using padding
// inputStudiable : A studiable of any shape
// resultStudiable : A square matrix representing the input studiable after padding
//------------------------------------------------------------------------------------------------------------------------
void Calc::padding(const CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable)
{
    long m = inputStudiable.Width();
    long n = inputStudiable.Height();

    long diff = abs(m - n);
    long diffSurDeux = diff / 2;

    long x = LONG_MIN;

    if (m > n)
    {
        for (int i = 0; i < diffSurDeux; i++)
        {
            for (int j = 0; j < m; j++)
            {
                resultStudiable(i, j) = LONG_MIN;
            }
        }
        for (int i = diffSurDeux; i < n+diffSurDeux; i++)
        {
            for (int j = 0; j < m; j++)
            {
                resultStudiable(i, j) = inputStudiable(i-diffSurDeux, j);
            }
        }
        for (int i = n + diffSurDeux; i < n + diff; i++)
        {
            for (int j = 0; j < m; j++)
            {
                resultStudiable(i, j) = LONG_MIN;
            }
        }
    }
    else
    {
        for (int j = 0; j < diffSurDeux; j++)
        {
            for (int i = 0; i < n; i++)
            {
                resultStudiable(i, j) = LONG_MIN;
            }
        }
        for (int j = diffSurDeux; j < m + diffSurDeux; j++)
        {
            for (int i = 0; i < n; i++)
            {
                resultStudiable(i, j) = inputStudiable(i, j - diffSurDeux);
            }
        }
        for (int j = m + diffSurDeux; j < m + diff; j++)
        {
            for (int i = 0; i < n; i++)
            {
                resultStudiable(i, j) = LONG_MIN;
            }
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------
// Inverse Radon transform
// 
// Inverse Radon transform transform 2D written in C++, turns a sinogramm into a 2D image
//------------------------------------------------------------------------------------------------------------------------
/*void Calc::iRadon(CLocalMatrixStudiable& inputStudiable, CLocalMatrixStudiable& resultStudiable)
{
    long M = resultStudiable.Width();
    long N = resultStudiable.Height();

    double xprev = 0;
    double yprev = 0;

    for (long theta = 0; theta < 180; theta++)
    {
        for (int x = 0; x < M; x++)
        {
            for (int y = 0; y < N; y++)
            {
                double sum = 0;
                long index = 0;
                double thetarad = theta * PI / 180;
                for (long rho = 0; rho < floor(N*sqrt(2)); rho++)
                {
                    xprev = x * cos(thetarad) - y * sin(thetarad) + 128;
                    yprev = x * sin(thetarad) + y * cos(thetarad) + 128;
                    if (xprev >= 0 && yprev >= 0 && xprev < M && yprev < N)
                    {
                        if (inputStudiable(rho, theta) != LONG_MIN)
                        {
                            sum += inputStudiable(rho, theta) * interpolation_bilineaire(xprev, yprev, inputStudiable);
                            index++;
                        }
                    }
                    resultStudiable(y, x) += sum/180;
                }
            }
        }
    }
}


double Calc::interpolation_bilineaire(long x, long y, CLocalMatrixStudiable& inputStudiable)
{
    long x1 = floor(x);
    long x2 = x1 + 1;
    long y1 = floor(y);
    long y2 = y1 + 1;

    if (x2 >= inputStudiable.Width())
    {
        x2 = x1;
    }

    if (y2 >= inputStudiable.Height())
    {
        y2 = y1;
    }

    double q11 = inputStudiable(y1, x1);
    double q12 = inputStudiable(y1, x2);
    double q21 = inputStudiable(y2, x1);
    double q22 = inputStudiable(y2, x2);

    double valeur_interpolee = (q11 * (x2 - x) * (y2 - y) +
        q21 * (x - x1) * (y2 - y) +
        q12 * (x2 - x) * (y - y1) +
        q22 * (x - x1) * (y - y1));

    return valeur_interpolee;
}*/


void Calc::RadonV4(const CLocalMatrixStudiable& data, CLocalMatrixStudiable& output)
{
    long nn = data.Width();
    long n2 = 2 * nn;
    CLocalMatrixStudiable a(n2, nn);
    CLocalMatrixStudiable b(n2, nn);
    CLocalMatrixStudiable c(n2, nn);
    CLocalMatrixStudiable d(n2, nn);
    //if (nn < 2 || nn & (nn - 1)) throw("nn must be power of 2 in Radon");
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) a(j + nn, i) = data(j, i);
    xformpiece(a);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) b(j + nn, i) = data(i, j);
    xformpiece(b);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) c(j + nn, i) = data(nn - 1 - i, j);
    xformpiece(c);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) d(j + nn, i) = data(nn - 1 - j, i);
    xformpiece(d);

    for (long j = 0; j < n2; j++) for (long i = 0; i < nn; i++) output(j, i) = c(n2-1-j, i);
    for (long j = 0; j < n2; j++) for (long i = 0; i < nn; i++) output(j, i+nn) = d(j, nn-1-i);
    for (long j = 0; j < n2; j++) for (long i = 0; i < nn; i++) output(j+nn, i + 2*nn) = a(n2-1-j, i);
    for (long j = 0; j < n2; j++) for (long i = 0; i < nn; i++) output(j+nn, i + 3*nn) = b(j, nn-1-i);
}


void Calc::xformpiece(CLocalMatrixStudiable& a)
{
    long mm, n, ng, i, j, k, ii, nr = a.Height(), nn = a.Width(), noff = nr - nn;
    std::vector<long> b(nn);
    for (j = 0; j < noff; j++) for (i = 0; i < nn; i++) a(j, i) = 0.;
    mm = 1;
    ng = nn >> 1;
    while (ng) {
        for (j = 0; j < nr; j++) { // loop over rows
            ii = k = 0; // ii will be 2*n*mm+i
            for (n = 0; n < ng; n++) {
                for (i = 0; i < mm; i++) {
                    b[k + 1] = a(j, ii);
                    b[k] = b[k + 1];
                    if ((j + i) + 1 < nr) b[k + 1] = b[k + 1] + a((j + i) + 1, ii + mm);
                    if (j + i < nr) b[k] = b[k] + a(j + i, ii + mm);
                    ii++;
                    k += 2;
                }
                ii += mm;
            }
            for (k = 0; k < nn; k++) a(j, k) = b[k];
        }
        mm += mm;
        ng >>= 1;
    }
}

/*void Calc::backproject(CLocalMatrixStudiable& ans)
{
    Int i, j;
    long N = ans.Height();
    long M = ans.Width();
    long nn = N;
    long n2 = 2 * nn;
    CLocalMatrixStudiable a(n2, nn);
    CLocalMatrixStudiable b(n2, nn);
    CLocalMatrixStudiable c(n2, nn);
    CLocalMatrixStudiable d(n2, nn);
    CLocalMatrixStudiable tmp(nn, nn);
    //if (ans.ncols() != nn) ans.resize(nn, nn);
    Doub val = 0.25 / (nn - 1.);
    backpiece(ans, a);
    backpiece(tmp, b);
    for (j = 0; j < nn; j++) for (i = 0; i < nn; i++) ans(j, i) += tmp(i, j);
    backpiece(tmp, c);
    for (j = 0; j < nn; j++) for (i = 0; i < nn; i++) ans(j, i) += tmp(i, nn - 1 - j);
    backpiece(tmp, d);
    for (j = 0; j < nn; j++) for (i = 0; i < nn; i++) ans(j, i) += tmp(nn - 1 - j, i);
    for (j = 0; j < nn; j++) for (i = 0; i < nn; i++) ans(j, i) *= val;
}


void Calc::backpiece(CLocalMatrixStudiable& ans, CLocalMatrixStudiable& a)
{
    Int mm, ng, n, i, ii, j, k, nn, n2;
    CLocalMatrixStudiable b(n2, nn); // 
    std::vector<long> bb(nn);
    mm = nn >> 1;
    ng = 1;
    for (j = 0; j < n2; j++) for (i = 0; i < nn; i++) b(j, i) = a(j, i); // copy a only to avoid overwriting it
    while (mm) {
        for (j = n2 - 1; j > nn - (mm + mm); j--) { // loop over rows
            for (i = 0; i < nn; i++) bb[i] = b(j, i);
            ii = 0; // ii will be 2*n*mm+i
            for (n = 0; n < ng; n++) {
                k = ii; // k will be 2*n*mm
                for (i = 0; i < mm; i++) {
                    b(j, ii) = bb[k + i + i] + bb[k + i + i + 1];
                    if (j + i < n2) b(j + i, mm + ii) = bb[k + i + i] + b(j - 1, k + i + i + 1);
                    ii++;
                }
                ii += mm;
            }
        }
        ng += ng;
        mm >>= 1;
    }
    for (j = 0; j < nn; j++) for (i = 0; i < nn; i++) ans(j, i) = b(j + nn, i);
}


void Calc::inverse_recursive(CLocalMatrixStudiable& img, CLocalMatrixStudiable& rad) {
    Int i, j, ii, jj, k, nn = img.Width(), n2 = 2*nn, nh = nn / 2;
    CLocalMatrixStudiable a(n2, nn);
    CLocalMatrixStudiable b(n2, nn);
    CLocalMatrixStudiable c(n2, nn);
    CLocalMatrixStudiable d(n2, nn);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) a(j + nn, i) = rad(j, i);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) b(j + nn, i) = rad(i, j);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) c(j + nn, i) = rad(nn - 1 - i, j);
    for (long j = 0; j < nn; j++) for (long i = 0; i < nn; i++) d(j + nn, i) = rad(nn - 1 - j, i);
    if (rad.Width() == 1) {
        img(0, 0) = a(1, 0);
    }
    else {
        { //scope
            CLocalMatrixStudiable radh(4*nh, 4*nh);
            CLocalMatrixStudiable imgh(nh, nh);
            inverse_recursive(imgh, radh); // recursive call!
            for (jj = 0, j = 0; j < nh; j++, jj += 2) {
                for (ii = 0, i = 0; i < nh; i++, ii += 2) {
                    img(jj, ii) = img(jj + 1, ii) = img(jj, ii + 1) = img(jj + 1, ii + 1) = imgh(j, i);
                }
            }
        } // endscope
        mprove(img, rad, 2);
    }
}


void Calc::mprove(CLocalMatrixStudiable& img, CLocalMatrixStudiable& rad, Int method) {
    Int i, j, nn = img.Height(), n2 = 2 * nn;
    CLocalMatrixStudiable img1(nn, nn);
    CLocalMatrixStudiable rad1(img.Height(), img.Width());
    for (long i = 0; i < img.Height(); i++) for (long j = 0; j < img.Width(); j++) rad1(i, j) = img(i, j);
    rad1 -= rad;
    if (method == 1) {
        inverse_recursive(img1, rad1);
    }
    else if (method == 2) {
        rad1.backproject(img1);
        //hipass(img1);
    }
    else throw ("no such method in mprove");
    for (j = 0; j < nn; j++) for (i = 0; i < nn; i++) img[j][i] -= img1[j][i];
}*/


//------------------------------------------------------------------------------------------------------------------------
// Alpha factor computation
// 
// Compute the scaling factor (Alpha) between the two inputs using Radon transforms
// inputStudiable1 : The matrix representing the first studiable
// inputStudiable2 : The matrix representing the second studiable
// alpha1 : The scale factor of the first studiable
// alpha2 : The scale factor of the second studiable
// thetaComputed : The theta difference computed between the two studiables
//------------------------------------------------------------------------------------------------------------------------
double Calc::Alpha(const CLocalMatrixStudiable& inputStudiable1, CLocalMatrixStudiable& inputStudiable2, double alpha1, double alpha2, long thetaComputed)
{
    double sum1, sum2, inter1, inter2, res, sumtest1, sumtest2, sumtest3, sumtest4, sumtest5, min;
    sum1 = 0;
    sum2 = 0;
    sumtest1 = 0;
    sumtest2 = 0;
    sumtest3 = 0;
    sumtest4 = 0;
    sumtest5 = 0;
    min = inputStudiable1.GetZMinLocal();
    long N = inputStudiable1.Height();
    long M = inputStudiable2.Height();
    long diff = abs(inputStudiable1.Height() - inputStudiable2.Height()) / 2;
    if (N > M)
    {
        for (int i = 0; i < M; i++)
        {
            sum1 += inputStudiable1(i + diff, 0) + min;
            sum2 += inputStudiable2(i, 0 + thetaComputed) + min;
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            sum1 += inputStudiable1(i, 0) + min;
            sum2 += inputStudiable2(i + diff, 0 + thetaComputed) + min;
        }
    }
    inter1 = sum1 * alpha1;
    inter2 = sum2 * alpha2;
    res = abs(inter2/inter1);
    return res;
}


//------------------------------------------------------------------------------------------------------------------------
// Alpha factor computation
// 
// Compute the scaling factor (Alpha) between the two inputs using Radon transforms
//------------------------------------------------------------------------------------------------------------------------
double Calc::AlphaPrecise(const CLocalMatrixStudiable& inputStudiable1, CLocalMatrixStudiable& inputStudiable2, double alpha1, double alpha2, long thetaComputed, long precision)
{
    double res, tempRes;
    res = 0;
    for (int i = 0; i < precision; i++)
    {
        tempRes = Alpha(inputStudiable1, inputStudiable2, alpha1, alpha2, (thetaComputed+i*10)%360);
        res += tempRes;
    }
    return res/precision;
}


//------------------------------------------------------------------------------------------------------------------------
// Xf factor computation
// 
// Compute the center of mass on the X axis (Xf) of the input
//------------------------------------------------------------------------------------------------------------------------
double Calc::TranslateXf(const CLocalMatrixStudiable& inputStudiable, long theta)
{
    double num1, num2, denom, res, rhomax, rhomin, thetarad, costheta, sintheta, min;
    thetarad = PI * (double)theta / 180;
    min = abs(inputStudiable.GetZMinLocal());
    costheta = cos(thetarad);
    sintheta = sin(thetarad);
    num1 = 0;
    num2 = 0;
    denom = 0;
    if (inputStudiable.Height() % 2 == 0)
    {
        rhomax = inputStudiable.Height()/2;
        rhomin = -rhomax+1;
    }
    else
    {
        rhomax = inputStudiable.Height()/2;
        rhomin = -rhomax;
    }
    for (int i = (long)rhomin; i < (long)rhomax+1; i++)
    {
        num1 += (double)(inputStudiable(i - (long)rhomin, theta) + min) * i;
        num2 += (double)(inputStudiable(i- (long)rhomin, (theta + 90) % 360)+min) * i;
        denom += (inputStudiable(i- (long)rhomin, theta)+min);
    }
    res = (num1/denom)*costheta-(num2/denom)*sintheta;
    double inter = res + (double)inputStudiable.Height()/2;
    return inter / sqrt(2);
}


//------------------------------------------------------------------------------------------------------------------------
// Yf factor computation
// 
// Compute the center of mass on the Y axis (Yf) of the input
//------------------------------------------------------------------------------------------------------------------------
double Calc::TranslateYf(const CLocalMatrixStudiable& inputStudiable, long theta)
{
    double num1, num2, denom, res, rhomax, rhomin, thetarad, costheta, sintheta, min;
    thetarad = PI * (double)theta / 180;
    min = abs(inputStudiable.GetZMinLocal());
    costheta = cos(thetarad);
    sintheta = sin(thetarad);
    num1 = 0;
    num2 = 0;
    denom = 0;
    if (inputStudiable.Height() % 2 == 0)
    {
        rhomax = inputStudiable.Height() / 2;
        rhomin = -rhomax + 1;
    }
    else
    {
        rhomax = inputStudiable.Height() / 2;
        rhomin = -rhomax;
    }
    for (int i = (long)rhomin; i < (long)rhomax + 1; i++)
    {
        num1 += (double)(inputStudiable(i - (long)rhomin, theta) + min) * i;
        num2 += (double)(inputStudiable(i - (long)rhomin, (theta + 90) % 360) + min) * i;
        denom += (inputStudiable(i - (long)rhomin, theta) + min);
    }
    res = (num1 / denom) * sintheta + (num2 / denom) * costheta;
    double inter = res + (double)inputStudiable.Height() / 2;
    return inter / sqrt(2);
}


//------------------------------------------------------------------------------------------------------------------------
// Theta factor computation
// 
// Compute the rotation factor (Theta) between the two inputs using Radon transforms
//------------------------------------------------------------------------------------------------------------------------
long Calc::Theta(const CLocalMatrixStudiable& inputStudiable1, const CLocalMatrixStudiable& inputStudiable2)
{
    long N1 = inputStudiable1.Width();
    long M1 = inputStudiable1.Height();
    long N2 = inputStudiable2.Width();
    long M2 = inputStudiable2.Height();
    double max = max(M1, M2);
    double min = min(M1, M2);
    double minVal = 0;
    long diff = (long)abs(min - max)/2;
    double sum;
    long id = 0;
    if (M1 > M2)
    {
        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < min; j++)
            {
                minVal += abs(inputStudiable1(j+diff, i) - inputStudiable2(j, i));
            }
        }
        for (int thetazero = 0; thetazero < 360; thetazero++)
        {
            sum = 0;
            for (int i = 0; i < N1; i++)
            {
                for (int j = 0; j < min; j++)
                {
                    sum += abs(inputStudiable1(j+diff, i) - inputStudiable2(j, (i+thetazero)%360));
                }
            }
            if (sum < minVal)
            {
                minVal = sum;
                id = thetazero;
            }
        }
    }
    else
    {
        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < min; j++)
            {
                minVal += abs(inputStudiable1(j, i) - inputStudiable2(j+diff, i));
            }
        }
        for (int thetazero = 0; thetazero < 360; thetazero++)
        {
            sum = 0;
            for (int i = 0; i < N1; i++)
            {
                for (int j = 0; j < min; j++)
                {
                    sum += abs(inputStudiable1(j, i) - inputStudiable2(j+diff, (i + thetazero) % 360));
                }
            }
            if (sum < minVal)
            {
                minVal = sum;
                id = thetazero;
            }
        }
    }
    return id;
}
