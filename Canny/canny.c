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

typedef struct Mask{
    double X[MAXMASK][MAXMASK];
    double Y[MAXMASK][MAXMASK];
    int radius;
}Mask;

typedef struct GradientVector{
    double X[PICSIZE][PICSIZE];
    double Y[PICSIZE][PICSIZE];
}GradientVector;

void canny(FILE*, FILE*, double, double);
void readImageTo(int[PICSIZE][PICSIZE], FILE*);
Mask generateSmoothenerMask(double);
GradientVector calculateGradients(int[PICSIZE][PICSIZE], Mask);
double calculateMagnitudes(GradientVector, double, int, double[PICSIZE][PICSIZE]);
void scaleImageWithRespectTo(double maxVal, double[PICSIZE][PICSIZE]);
void applyPeakFilterTo(double[PICSIZE][PICSIZE], GradientVector, int);
void writeTo(FILE*, double[PICSIZE][PICSIZE]);
double findMaxValue(double[PICSIZE][PICSIZE]);

int main(int argc, char** argv){
    // Command-Line Arguements
    FILE* inputImage=fopen(argv[1],"rb");
    FILE* outputImage=fopen(argv[2],"wb");
    double sigma = atof(argv[3]);
    double threshold = atof(argv[4]);
    
    Record record;
    fread(&record, sizeof(record), 1, inputImage);
    fwrite(&record, sizeof(record), 1, outputImage);

    canny(inputImage, outputImage, sigma, threshold);
}

void canny(FILE* inputImageFile, FILE* outputImageFile, double sigma, double threshold){
    Mask mask;
    mask.radius = (int)(sigma * 3);
    int picBuffer[PICSIZE][PICSIZE];
    double itsMagnitudes[PICSIZE][PICSIZE];
    double maxMagnitude;
    GradientVector usingGradients;

    readImageTo(picBuffer, inputImageFile);
    mask = generateSmoothenerMask(sigma);
    usingGradients = calculateGradients(picBuffer, mask);
    calculateMagnitudes(usingGradients, threshold, mask.radius, itsMagnitudes);
    maxMagnitude = findMaxValue(itsMagnitudes);
    scaleImageWithRespectTo(maxMagnitude, itsMagnitudes);
    applyPeakFilterTo(itsMagnitudes, usingGradients, mask.radius);

    writeTo(outputImageFile, itsMagnitudes);
}

void readImageTo(int buffer[PICSIZE][PICSIZE], FILE* fileStream){
    for (int i = 0; i < 256; i++){
        for (int j = 0; j < 256; j++){
            buffer[i][j] = getc (fileStream);
        }
    }   
}

Mask generateSmoothenerMask(double sigma){
    int maskRadius = (int)(sigma*3);
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

GradientVector calculateGradients(int pic[PICSIZE][PICSIZE], Mask mask){
    GradientVector gradients;
    int centerX = MAXMASK/2;
    int centerY = MAXMASK/2;
    for (int i = mask.radius; i <= 255-mask.radius; i++){
        for (int j = mask.radius; j<= 255-mask.radius; j++){
            int sumX = 0;
            int sumY = 0;
            for (int p = -mask.radius; p <= mask.radius; p++){
                for (int q = -mask.radius; q <= mask.radius; q++){
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

void applyPeakFilterTo(double magnitudes[PICSIZE][PICSIZE], GradientVector gradients,
                        int maskRadius){
    for(int i = maskRadius; i < 256-maskRadius; i++){ 
        for(int j = maskRadius; j < 256-maskRadius; j++){
            if (gradients.X[i][j] == 0.0) { 
                gradients.X[i][j] = 0.00001;
            }
            double slope = gradients.Y[i][j]/gradients.X[i][j];
            if (slope <= 0.4142 && slope > -0.4142){
                if (magnitudes[i][j] > magnitudes[i][j-1] && magnitudes[i][j] > magnitudes[i][j+1]){ 
                    magnitudes[i][j] = 255;
                } 
                else{
                    magnitudes[i][j] = 0;
                }
            }
            else if (slope <= 2.4142 && slope > 0.4142){
                if (magnitudes[i][j] > magnitudes[i-1][j-1] && magnitudes[i][j] > magnitudes[i+1][j+1]){
                    magnitudes[i][j] = 255; 
                }
                else{
                    magnitudes[i][j] = 0;
                }
            }
            else if (slope <= -0.4142 && slope > -2.4142){
                if (magnitudes[i][j] > magnitudes[i+1][j-1] && magnitudes[i][j] > magnitudes[i-1][j+1]){
                    magnitudes[i][j] = 255;
                } 
                else{
                    magnitudes[i][j] = 0;
                }
            }else{
                if (magnitudes[i][j] > magnitudes[i-1][j] && magnitudes[i][j] > magnitudes[i+1][j]){
                    magnitudes[i][j] = 255;
                } 
                else{
                    magnitudes[i][j] = 0;
                }
            }
        }
    }
}

void writeTo(FILE* outputFile, double image[PICSIZE][PICSIZE]){
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

