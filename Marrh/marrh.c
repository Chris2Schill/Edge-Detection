#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <stdlib.h>
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

typedef struct Record{
    char a[2];
    int width, height, graylevel;
}Record;
Record record;

void marrh(FILE*, FILE*, double);
void readImageToPicBuffer(FILE*);
void calculateSmoothenerMask(double*, int, double, int, int);
void convolution(int,int,int);
void getMaxAndMinVals(int, double*, double*, double*);
void writeToFile(FILE* outputFile1, int maxval);
void calculateGradients();
void calculateMagnitudes(double*, double);
void adjustImageValuesWithRespectTo(double maxVal);

int    pic[PICSIZE][PICSIZE];
double outpic1[PICSIZE][PICSIZE];
double outpic2[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double mask[MAXMASK][MAXMASK];
double conv[PICSIZE][PICSIZE];


int main(int argc, char** argv){

    int     i,j,p,q,s,t;
    double  maskval,sum,maxival,minival,maxval;

    // Command-Line Arguements
    FILE* inputFile=fopen(argv[1],"rb");
    FILE* outputFile1=fopen(argv[2],"wb");
    FILE* outputFile2=fopen(argv[3],"wb");
    double sig = atof(argv[4]);
    double ZEROTOL = atof(argv[5]);

    int mr = (int)(sig * 3);
    int centx = MAXMASK/2;
    int centy = MAXMASK/2;

    fread(&record, sizeof(record), 1, inputFile);
    fwrite(&record, sizeof(record), 1, outputFile1);
    fwrite(&record, sizeof(record), 1, outputFile2);

    readImageToPicBuffer(inputFile);
    calculateSmoothenerMask(&maskval, mr, sig, centx, centy);
    convolution(mr, centx, centy);

    maxval  = 0;
    maxival = 0;
    minival = 255;

    getMaxAndMinVals(mr, &maxival, &minival, &maxval);
    writeToFile(outputFile1, maxval);
    for (i = mr; i <= 255-mr; i++){
        for (j = mr; j <= 255-mr; j++){
            outpic2[i][j] = 0;
            if (conv[i][j] > ZEROTOL){
                for (p = -1; p <= 1; p++){
                    for (q = -1;q <= 1; q++){
                        if (conv[i+p][j+q] < -ZEROTOL){
                            outpic2[i][j] = 255;
                        }
                    }
                }
            }
            else if ((fabs)(conv[i][j]) < ZEROTOL){
                if (((conv[i+1][j] > ZEROTOL) && (conv[i-1][j] < -ZEROTOL)) ||
                        ((conv[i+1][j] < -ZEROTOL) && (conv[i-1][j] > ZEROTOL))){
                    outpic2[i][j] = 255;
                }
                else if (((conv[i][j+1] > ZEROTOL) && (conv[i][j-1] < -ZEROTOL))   ||
                        ((conv[i][j+1] < -ZEROTOL) && (conv[i][j-1] > ZEROTOL))){
                    outpic2[i][j] = 255;
                }
                else if (((conv[i+1][j+1] > ZEROTOL) && (conv[i-1][j-1] < -ZEROTOL))   ||
                        ((conv[i+1][j+1] < -ZEROTOL) && (conv[i-1][j-1] > ZEROTOL))){
                    outpic2[i][j] = 255;

                }
                else if (((conv[i+1][j-1] > ZEROTOL) && (conv[i-1][j+1] < -ZEROTOL))   ||
                        ((conv[i+1][j-1] < -ZEROTOL) && (conv[i-1][j+1] > ZEROTOL))){
                    outpic2[i][j] = 255;
                }
            }
        }
    }

    for (i = 0; i < 256; i++){
        for (j = 0; j < 256; j++){
            if (outpic2[i][j] == 255){
                outpic2[i][j]=0;
            }
            else {
                outpic2[i][j]=255;
            }
            fprintf(outputFile2,"%c",(char)((int)(outpic2[i][j])));
        }
    }
}

void readImageToPicBuffer(FILE* fileStream){
    for (int i = 0; i < 256; i++){
        for (int j = 0; j < 256; j++){
            pic[i][j]  =  getc (fileStream);
        }
    }   
}

void calculateSmoothenerMask(double* maskval, int mr, double sig, int centx, int centy){
    for (int p = -mr; p <= mr; p++){
        for (int q = -mr; q <= mr; q++){
            *maskval = ((2-(((p*p)+(q*q))/(sig*sig)))*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
            mask[p+centy][q+centx] = *maskval;
        }
    }   
}

void convolution(int mr, int centx, int centy){
    for (int i = mr; i <= 255-mr; i++){
        for (int j = mr; j<= 255-mr; j++){
            int sum = 0;
            for (int p = -mr; p <= mr; p++){
                for (int q = -mr; q <= mr; q++){
                    sum += pic[i+p][j+q] * mask[p+centy][q+centx];
                }
            }
            outpic1[i][j] = sum;
            conv[i][j] = sum;
        }
    }   
}

void getMaxAndMinVals(int mr, double* maxival, double* minival, double* maxval){
    for (int i = mr; i < 256-mr; i++){
        for (int j = mr; j < 256-mr; j++){
//            printf("i: %d, j: %d, mr: %d\n",i,j,mr);
            if (outpic1[i][j] > *maxival){
                *maxival = outpic1[i][j];
            }
            if (outpic1[i][j] < *minival){
                *minival = outpic1[i][j];
            }
        }
    }

    if (fabs(*maxival) > fabs(*minival)){
        *maxval = fabs(*maxival);
    }
    else{
        *maxval = fabs(*minival);
    }
}

void writeToFile(FILE* outputFile1, int maxval){
    for (int i = 0; i < 256; i++){
        for (int j = 0; j < 256; j++){
            outpic1[i][j] = ((((outpic1[i][j]) / maxval) + 1) * 127);
            fprintf(outputFile1,"%c",(char)((int)(outpic1[i][j])));
        }
    }       
}
