/*  Canny-Edge-Detector.c  (or canny.c) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

typedef struct Record{
    char a[2];
    int width, height, graylevel;
}Record;
Record record;

typedef struct Mask{
    double X[MAXMASK][MAXMASK];
    double Y[MAXMASK][MAXMASK];
}Mask;

typedef struct GradientVector{
    double X[PICSIZE][PICSIZE];
    double Y[PICSIZE][PICSIZE];
}GradientVector;

void canny(FILE*, FILE*, double, double);
void readImageToBuffer(int[PICSIZE][PICSIZE], FILE*);
Mask generateSmoothenerMask(int, double);
GradientVector calculateGradients(int[PICSIZE][PICSIZE], Mask, int);
double calculateMagnitudes(GradientVector, double, int, double[PICSIZE][PICSIZE]);
void scaleImageWithRespectTo(double maxVal, double[PICSIZE][PICSIZE]);
void writeImageToFile(double[PICSIZE][PICSIZE], FILE*);
double findMaxValue(double[PICSIZE][PICSIZE]);


int main(int argc, char** argv){
    // Command-Line Arguements
    FILE* inputImage=fopen(argv[1],"rb");
    FILE* outputImage=fopen(argv[2],"wb");
    double sigma = atof(argv[3]);
    double threshold = atof(argv[4]);
    
    fread(&record, sizeof(record), 1, inputImage);
    fwrite(&record, sizeof(record), 1, outputImage);

    canny(inputImage, outputImage, sigma, threshold);
}

void canny(FILE* inputImage, FILE* outputImage, double sigma, double threshold){
    Mask mask;
    int maskRadius = (int)(sigma * 3);
    int pic[PICSIZE][PICSIZE];
    double magnitudes[PICSIZE][PICSIZE];
    double maxMagnitude;
    GradientVector gradients;

    readImageToBuffer(pic, inputImage);
    mask = generateSmoothenerMask(maskRadius, sigma);
    gradients = calculateGradients(pic, mask, maskRadius);
    calculateMagnitudes(gradients, threshold, maskRadius, magnitudes);
    maxMagnitude = findMaxValue(magnitudes);
    scaleImageWithRespectTo(maxMagnitude, magnitudes);
    writeImageToFile(magnitudes, outputImage);
}

void readImageToBuffer(int pic[PICSIZE][PICSIZE], FILE* fileStream){
    for (int i = 0; i < 256; i++){
        for (int j = 0; j < 256; j++){
            pic[i][j] = getc (fileStream);
        }
    }   
}

Mask generateSmoothenerMask(int maskRadius, double sigma){
    int centerX = MAXMASK/2;
    int centerY = MAXMASK/2;
    Mask mask;
    for (int p = -maskRadius; p <= maskRadius; p++){
        for (int q = -maskRadius; q <= maskRadius; q++){
            mask.X[p+centerY][q+centerX] = q*(exp(-(q*q + p*p)/(2*sigma*sigma)));
            mask.Y[p+centerY][q+centerX] = p*(exp(-(q*q + p*p)/(2*sigma*sigma)));
        }
    }
    return mask;
}

GradientVector calculateGradients(int pic[PICSIZE][PICSIZE], Mask mask, int maskRadius){
    GradientVector gradients;
    int centerX = MAXMASK/2;
    int centerY = MAXMASK/2;
    for (int i = maskRadius; i <= 255-maskRadius; i++){
        for (int j = maskRadius; j<= 255-maskRadius; j++){
            int sumX = 0;
            int sumY = 0;
            for (int p = -maskRadius; p <= maskRadius; p++){
                for (int q = -maskRadius; q <= maskRadius; q++){
                    sumX += pic[i+p][j+q] * mask.X[p+centerY][q+centerX];
                    sumY += pic[i+p][j+q] * mask.Y[p+centerY][q+centerX];
                }
            }
            gradients.X[i][j] = sumX;
            gradients.Y[i][j] = sumY;
        }
    }   
    return gradients;
}

//Stores magnitudes and returns the maximum 
double calculateMagnitudes(GradientVector gradients, double threshold, int maskRadius, 
                            double magnitudes[PICSIZE][PICSIZE]){
    double maxVal = 0;
    for (int i = maskRadius; i < PICSIZE-maskRadius ; i++){
        for (int j = maskRadius; j < PICSIZE-maskRadius; j++){
            double magnitude = sqrt((double)((gradients.X[i][j]*gradients.X[i][j]) +
                        (gradients.Y[i][j]*gradients.Y[i][j])));
            if (magnitude > threshold){
                magnitudes[i][j] = magnitude;
            }else{
                magnitudes[i][j] = 0;
            }
        }
    }
    return maxVal;
}

void scaleImageWithRespectTo(double maxVal, double magnitudes[PICSIZE][PICSIZE]){
    for (int i = 0; i < PICSIZE; i++){
        for (int j = 0; j < PICSIZE; j++){
            magnitudes[i][j] = (magnitudes[i][j] / maxVal) * (PICSIZE-1);            
        }
    }
}

void writeImageToFile(double image[PICSIZE][PICSIZE], FILE* outputFile){
    for (int i = 0; i < PICSIZE; i++){
        for (int j = 0; j < PICSIZE; j++){
            fprintf(outputFile,"%c",(char)((int)(image[i][j])));
        }
    }
}

double findMaxValue(double array[PICSIZE][PICSIZE]){
    double max = 0;
    for (int i = 0; i < PICSIZE; i++){
        for (int j = 0; j < PICSIZE; j++){
            if (array[i][j] > max){
                max = array[i][j];
            }
        }
    }
    return max;
}
