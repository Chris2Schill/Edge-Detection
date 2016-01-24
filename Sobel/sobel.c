/* Sobel.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PIC_SIZE 256
#define MASK_SIZE 3

typedef struct Record{
    char a[2];
    int width, height, graylevel;
}Record;
Record record;


void sobel(FILE*, FILE*, double);
void readImageToPicBuffer(FILE*);
void calculateGradients();
void calculateMagnitudes(double*, double);
void adjustImageValuesWithRespectTo(double maxVal);
void writeToFile(FILE* outputFile);
void displayPicValues(double**);

int pic[PIC_SIZE][PIC_SIZE];
int outpicx[PIC_SIZE][PIC_SIZE];
int outpicy[PIC_SIZE][PIC_SIZE];
int maskx[MASK_SIZE][MASK_SIZE] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[MASK_SIZE][MASK_SIZE] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[PIC_SIZE][PIC_SIZE];


int main(int argc,char** argv){
    FILE *inputFile, *outputFile1, *outputFile2, *outputFile3;

    inputFile=fopen(argv[1],"rb");
    outputFile1=fopen(argv[2],"wb");
    outputFile2=fopen(argv[3],"wb");
    outputFile3=fopen(argv[4],"wb");
    double firstThreshold = atof(argv[5]);
    double secondThreshold = atof(argv[6]);
    double thirdThreshold = atof(argv[7]);

    fread(&record, sizeof(record), 1, inputFile);

    fwrite(&record, sizeof(record), 1, outputFile1);
    fwrite(&record, sizeof(record), 1, outputFile2);
    fwrite(&record, sizeof(record), 1, outputFile3);

    sobel(inputFile,outputFile1, firstThreshold);

    fclose(inputFile);
    inputFile = fopen(argv[1],"rb");
    fread(&record, sizeof(record), 1, inputFile);

    sobel(inputFile,outputFile2, secondThreshold);

    fclose(inputFile);
    inputFile = fopen(argv[1],"rb");
    fread(&record, sizeof(record), 1, inputFile);
    sobel(inputFile,outputFile3, thirdThreshold);
    return 0;
}

void sobel(FILE* inputFile, FILE* outputFile, double threshold){
    double maxVal = 0;
    readImageToPicBuffer(inputFile);
    calculateGradients();
    calculateMagnitudes(&maxVal, threshold);
    adjustImageValuesWithRespectTo(maxVal);
    writeToFile(outputFile);
}

void readImageToPicBuffer(FILE *fp){
    for (int i = 0; i < PIC_SIZE; i++){ 
        for (int j = 0; j < PIC_SIZE; j++){
            pic[i][j]  =  getc(fp);
            pic[i][j]  &= 0377;
        }
    }   
}

void calculateGradients(){
    int mr = 1;
    for (int i = mr; i < PIC_SIZE-mr; i++){
        for (int j = mr; j < PIC_SIZE-mr; j++){
            int sum1 = 0;
            int sum2 = 0;
            for (int p = -mr; p <= mr; p++){
                for (int q = -mr; q <= mr; q++){
                    sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                    sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
            }
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }
}

// Also stores the maximum value in maxVal
void calculateMagnitudes(double* maxVal, double threshold){
    int mr = 1;
    for (int i = mr; i < PIC_SIZE-mr ; i++){
        for (int j = mr; j < PIC_SIZE-mr; j++){
            double magnitude = sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                        (outpicy[i][j]*outpicy[i][j])));
            if (magnitude > threshold){
                ival[i][j] = magnitude;
            }else{
                ival[i][j] = 0;
            }

            if (ival[i][j] > *maxVal){ 
                *maxVal = ival[i][j];
            }
        }
    }
}

void adjustImageValuesWithRespectTo(double maxVal){
    for (int i = 0; i < PIC_SIZE; i++){
        for (int j = 0; j < PIC_SIZE; j++){
            ival[i][j] = (ival[i][j] / maxVal) * (PIC_SIZE-1);            
        }
    }
}

void writeToFile(FILE* outputFile){
    for (int i = 0; i < PIC_SIZE; i++){
        for (int j = 0; j < PIC_SIZE; j++){
            fprintf(outputFile,"%c",(char)((int)(ival[i][j])));
        }
    }
}

void displayPicValues(double **array){
    for (int i = 0; i < PIC_SIZE; i++){
        for (int j = 0; j < PIC_SIZE; j++){
            printf("%f ", array[i][j]);
        }
        printf("\n");
    }
}
