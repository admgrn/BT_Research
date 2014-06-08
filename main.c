#include <stdio.h>
#include <stdlib.h>

/* Possible TODO: Add these as input arugments */
const unsigned int BT_LEN = 625;
const unsigned int CTS_t = 60;
const unsigned int SIFS = 10;
const unsigned int DUMMY_t = 50;

const unsigned int MAX_SAMPLE_SIZE = 100;
const float BELIEF_RAND10_MAX_P = 0.01;

typedef struct
{
    unsigned int length;
    unsigned int* samples;
    
} Samples;

int ReadFoundTime(const char*, Samples**, unsigned int*);
Samples* AllocateSamples(const unsigned int);
Samples* GetCleanT(Samples*);
Samples* GetCleanF(Samples*);
void GetMaxValLoc(Samples*, unsigned int*, unsigned int*);
float Mean(Samples*);
int DetectBT(Samples*, unsigned int, unsigned int, unsigned int, float*, float*);


int main(int argc, const char* argv[])
{
    Samples* foundTime;
    Samples* cleanT;
    Samples* cleanF;
    unsigned int totalNumberSamples;
    unsigned int maxValue;
    unsigned int maxLocation;
    float sampleRatio;
    float cleanFAverage;
    
    if (argc == 2)
    {
        if (!ReadFoundTime(argv[1], &foundTime, &totalNumberSamples))
        {
            return 1;
        }
        
    }
    else
    {
        return 1;
    }
    
    cleanT = GetCleanT(foundTime);
    cleanF = GetCleanF(cleanT);
    GetMaxValLoc(cleanF, &maxValue, &maxLocation);
    
    if (DetectBT(cleanF, foundTime->length, maxValue, totalNumberSamples, &sampleRatio, &cleanFAverage))
    {
       printf("10 sample ratio %f, maxval %d, avg %f, We think there is BT!\n", sampleRatio,
              maxValue, cleanFAverage);
    }
    else
    {
        printf("10 sample ratio %f, maxval %d, avg %f, We think there is NOOOOOOOOOO BT!\n", sampleRatio,
               maxValue, cleanFAverage);

    }
    
    
    return 0;
}

int ReadFoundTime(const char* filename, Samples** foundTime, unsigned int* totalSamples)
{
    int i;
    int inputValue;

    FILE* inputStream = fopen(filename, "r+");
    inputValue = 0;
    
    if (inputStream == NULL)
    {
        return 0;
    }
    
    *foundTime = AllocateSamples(MAX_SAMPLE_SIZE);
    i = 0;
    
    if(fscanf(inputStream, "%d\n", totalSamples) == EOF)
    {
        return 0;
    }
    
    while(fscanf(inputStream, "%d\n", &inputValue) != EOF)
    {
        if (i >= MAX_SAMPLE_SIZE)
        {
            break;
        }
        
        (*foundTime)->samples[i++] = inputValue;
    }
    
    (*foundTime)->length = i;
    fclose(inputStream);
    return 1;
}

Samples* AllocateSamples(unsigned int size)
{
    Samples* sample;
    sample = malloc(sizeof(Samples));
    sample->samples = malloc(sizeof(unsigned int) * size);
    return sample;
}

Samples* GetCleanT(Samples* sample)
{
    unsigned int i;
    Samples* cleanT;
    
    cleanT = AllocateSamples(MAX_SAMPLE_SIZE);
    cleanT->length = sample->length;
    
    for (i = 0; i < sample->length; ++i)
    {
        /* TODO: Do we need to add 1 to all values, also negative 1 is subtracted to fix index */
        cleanT->samples[i] = (((sample->samples[i] - 1) - CTS_t - 2 * SIFS) % BT_LEN);
    }
    
    return cleanT;
}

Samples* GetCleanF(Samples* sample)
{
    Samples* cleanF;
    unsigned int D;
    unsigned int i, j;
    unsigned int T;
    unsigned int thisT;
    unsigned int thisBgn;
    unsigned int thisEnd;
    
    cleanF = AllocateSamples(BT_LEN);
    cleanF->length = BT_LEN;
    D = DUMMY_t + CTS_t + 2 * SIFS;
    T = BT_LEN - D;
    
    for (i = 0; i < cleanF->length; ++i)
    {
        cleanF->samples[i] = 0;
    }
    
    for (i = 0; i < sample->length; ++i)
    {
        thisT = sample->samples[i];
        
        if (thisT < T )
        {
            thisBgn = thisT;
            thisEnd = thisT + D - 1;    /* Do we need - 1 here */
            
            for (j = thisBgn; j < thisEnd; ++j)
            {
                cleanF->samples[j] = cleanF->samples[j] + 1;
            }
        }
        else
        {
            thisBgn = thisT;
            thisEnd = BT_LEN;
            
            for (j = thisBgn; j < thisEnd; ++ j)
            {
                cleanF->samples[j] = cleanF->samples[j] + 1;
            }
            
            thisBgn = 0;
            thisEnd = D - (BT_LEN - thisT);
            
           for (j = thisBgn; j < thisEnd; ++j)
            {
                cleanF->samples[j] = cleanF->samples[j] + 1;
            }
        }
    }
    
    return cleanF;
}

void GetMaxValLoc(Samples* sample, unsigned int* max, unsigned int* loc)
{
    unsigned int position = 0;
    unsigned int maximum = 0;
    unsigned int value;
    unsigned int i;
    
    for (i = 0; i < sample->length; ++i)
    {
        value = sample->samples[i];
        
        if (value > maximum)
        {
            maximum = value;
            position = i;
        }
    }
    
    *max = maximum;
    *loc = position;
}

float Mean(Samples* sample)
{
    unsigned int i;
    unsigned int value = 0;
    
    for (i = 0; i < sample->length; ++i)
    {
        value += sample->samples[i];
    }
    
    return (float)value / sample->length;
}

int DetectBT(Samples* sample, unsigned int found, unsigned int maxval, unsigned int total, float* ratio, float* average)
{
    float sampleRatio;
    float ave;
    
    sampleRatio = (float)found / total;
    ave = Mean(sample);
    
    
    if (sampleRatio > BELIEF_RAND10_MAX_P && maxval > ave * 2.5)
    {
        *ratio = sampleRatio;
        *average = ave;
        return 1;
    }
    else
    {
        *ratio = sampleRatio;
        *average = ave;
        return 0;
    }
    
}



