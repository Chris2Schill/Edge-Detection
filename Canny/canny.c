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

typedef struct Threshold{
    double high;
    double low;
}Threshold;

void canny(FILE*, FILE*, double);
void readImageTo(int[PICSIZE][PICSIZE], FILE*);
Mask generateSmoothenerMask(double);
GradientVector calculateGradients(int[PICSIZE][PICSIZE], Mask);
double calculateMagnitudes(GradientVector, int, double[PICSIZE][PICSIZE]);
void scaleImageWithRespectTo(double maxVal, double[PICSIZE][PICSIZE]);
void findPeaks(double[PICSIZE][PICSIZE], GradientVector, int, double[PICSIZE][PICSIZE]);
Threshold getThreshold(double[PICSIZE][PICSIZE]);
void applyHysteresisThresholdTo(double[PICSIZE][PICSIZE], GradientVector, int,
                        double[PICSIZE][PICSIZE], Threshold, double[PICSIZE][PICSIZE]);
void writeTo(FILE*, double[PICSIZE][PICSIZE]);
double findMaxValue(double[PICSIZE][PICSIZE]);

int main(int argc, char** argv){
    FILE* inputImage=fopen(argv[1],"rb");
    FILE* outputImage=fopen(argv[2],"wb");
    double sigma = atof(argv[3]);
    
    Record record;
    fread(&record, sizeof(record), 1, inputImage);
    fwrite(&record, sizeof(record), 1, outputImage);

    canny(inputImage, outputImage, sigma);
}

void canny(FILE* inputImageFile, FILE* outputImageFile, double sigma){
    Mask mask;
    int picBuffer[PICSIZE][PICSIZE];
    double itsMagnitudes[PICSIZE][PICSIZE];
    double maxMagnitude;
    GradientVector usingGradients;
    double peaks[PICSIZE][PICSIZE];
    double final[PICSIZE][PICSIZE];

    readImageTo(picBuffer, inputImageFile);
    mask = generateSmoothenerMask(sigma);
    usingGradients = calculateGradients(picBuffer, mask);
    calculateMagnitudes(usingGradients, mask.radius, itsMagnitudes);
    maxMagnitude = findMaxValue(itsMagnitudes);
    scaleImageWithRespectTo(maxMagnitude, itsMagnitudes);

    findPeaks(itsMagnitudes, usingGradients, mask.radius, peaks);
    Threshold threshold = getThreshold(itsMagnitudes);
    applyHysteresisThresholdTo(itsMagnitudes, usingGradients, mask.radius, 
                                peaks, threshold, final);
    writeTo(outputImageFile, final);
}

void readImageTo(int buffer[PICSIZE][PICSIZE], FILE* fileStream){
    for (int i = 0; i < PICSIZE; i++){
        for (int j = 0; j < PICSIZE; j++){
            buffer[i][j] = getc (fileStream);
        }
    }   
}

Mask generateSmoothenerMask(double sigma){
    int centerX = MAXMASK/2;
    int centerY = MAXMASK/2;
    Mask mask;
    mask.radius = (int)(sigma*3);
    for (int p = -mask.radius; p <= mask.radius; p++){
        for (int q = -mask.radius; q <= mask.radius; q++){
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
double calculateMagnitudes(GradientVector gradients, int maskRadius, 
                            double magnitudes[PICSIZE][PICSIZE]){
    double maxVal = 0;
    for (int i = maskRadius; i < PICSIZE-maskRadius ; i++){
        for (int j = maskRadius; j < PICSIZE-maskRadius; j++){
            double magnitude = sqrt((double)((gradients.X[i][j]*gradients.X[i][j]) +
                        (gradients.Y[i][j]*gradients.Y[i][j])));
            magnitudes[i][j] = magnitude;
            /*
            if (magnitude > threshold){
                magnitudes[i][j] = magnitude;
            }else{
                magnitudes[i][j] = 0;
            }
            */
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

void applyHysteresisThresholdTo(double magnitudes[PICSIZE][PICSIZE], GradientVector gradients,
                        int maskRadius, double peaks[PICSIZE][PICSIZE], Threshold threshold,
                        double final[PICSIZE][PICSIZE]){
    int moreToDo = 1;
    for (int i = 0; i < PICSIZE; i++){
        for (int j = 0; j < PICSIZE; j++){
            if (peaks[i][j] == 255){
                if (magnitudes[i][j] > threshold.high){
                    peaks[i][j] = 0;
                    final[i][j] = 255;
                }
                else if (magnitudes[i][j] < threshold.low){
                    peaks[i][j] = 0;
                    final[i][j] = 0;
                }
                
            }
        }
    }
    while (moreToDo){
        moreToDo = 0;
        for (int i = 0; i < PICSIZE; i++){
            for (int j = 0; j < PICSIZE; j++){
                if (peaks[i][j] == 255){
                    for (int p = -1; p <= 1; p++){
                        for (int q = -1; q <= 1; q++){
                            if (final[i+p][j+q] == 255){
                                peaks[i][j] = 0;
                                final[i][j] = 255;
                                moreToDo = 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

void findPeaks(double magnitudes[PICSIZE][PICSIZE], GradientVector gradients,
               int maskRadius, double peaks[PICSIZE][PICSIZE]){

    for(int i = maskRadius; i < 256-maskRadius; i++){ 
        for(int j = maskRadius; j < 256-maskRadius; j++){
            if((gradients.X[i][j]) == 0.0) {
                gradients.Y[i][j] = 0.00001;
            }
            double slope = gradients.Y[i][j]/gradients.X[i][j];
            double magnitude = magnitudes[i][j];
            if (slope <= 0.4142 && slope > -0.4142){
                if (magnitude > magnitudes[i][j-1] && magnitude > magnitudes[i][j+1]){ 
                    peaks[i][j] = 255;
                } 
            }
            else if (slope <= 2.4142 && slope > 0.4142){
                if (magnitude > magnitudes[i-1][j-1] && magnitude > magnitudes[i+1][j+1]){
                    peaks[i][j] = 255; 
                }
                
            }
            else if (slope <= -0.4142 && slope > -2.4142){
                if (magnitude > magnitudes[i+1][j-1] && magnitude > magnitudes[i-1][j+1]){ 
                    peaks[i][j] = 255;
                } 
            }else{
                if (magnitude > magnitudes[i-1][j] && magnitude > magnitudes[i+1][j]){ 
                    peaks[i][j] = 255;
                } 
            }
        
        }  
    }
}

void writeTo(FILE* outputFile, double image[PICSIZE][PICSIZE]){
    for (int i = 0; i < PICSIZE; i++){
        for (int j = 0; j < PICSIZE; j++){
            fprintf(outputFile,"%c",(char)((double)(image[i][j])));
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

Threshold getThreshold(double magnitudes[PICSIZE][PICSIZE]){
    int frequencies[PICSIZE];
    int i, j;
    for (i = 0; i < PICSIZE; i++) {
       frequencies[i] = 0;
    }

    for (i = 0; i < PICSIZE; i++){
        for (j = 0; j < PICSIZE; j++){
            (frequencies[(int)magnitudes[i][j]])++;
        }
    }

    double percent = 0.10;
    double cutOff = percent * PICSIZE * PICSIZE;
    int AreaOfTops = 0;
    for (i = PICSIZE-1; i >= 0; i--){
        AreaOfTops += frequencies[i]; 
        if (AreaOfTops > cutOff){
            break;
        }
    }
    Threshold newThreshold;
    newThreshold.high = i;
    newThreshold.low = i * 0.35;   
    printf("High: %f, Low: %f\n", newThreshold.high, newThreshold.low);
    return newThreshold;
}
