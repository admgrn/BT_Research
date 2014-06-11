#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

/* Possible TODO: Add these as input arugments */
const unsigned int BT_LEN = 625;		// Length of BT timeslot
const unsigned int CTS_t = 60;			// Time taken to transmit clear to send signal
const unsigned int SIFS = 10;			// Short Interframe Spacing
const unsigned int DUMMY_t = 50;		// Length of probe packet
const unsigned int MAX_SLOT_NUM = 6;		// Max allowable timeslots for a single BT transmission/ reponse
const unsigned int SCORE_SLOT = 2;

const unsigned int MAX_SAMPLE_SIZE = 100;
const float BELIEF_RAND10_MAX_P = 0.01;

typedef struct
{
    unsigned int length;
    unsigned int* samples;
} Samples;

typedef struct
{
    unsigned int min;
    unsigned int max;
} OutputBoundary;

int ReadFoundTime(const char*, Samples**, unsigned int*);

Samples* AllocateSamples(const unsigned int);

Samples* GetCleanT(Samples*);

Samples* GetCleanF(Samples*);

void GetMaxValLoc(Samples*, unsigned int*, unsigned int*);

float Mean(Samples*);

int DetectBT(Samples*, unsigned int, unsigned int, unsigned int);

OutputBoundary FindBoundary(Samples*, Samples*, Samples**, unsigned int*, unsigned int, unsigned int);

unsigned int FindSchedule(Samples*, Samples*, unsigned int[SCORE_SLOT][MAX_SLOT_NUM], unsigned int, unsigned int);

unsigned int Min(unsigned int[], unsigned int);

void DeallocateSamples(Samples*);

int main(int argc, const char* argv[])
{
    Samples* foundTime;
    Samples* cleanT;
    Samples* cleanF;
    Samples* shouldTake;
    OutputBoundary boundary;
    unsigned int scores[SCORE_SLOT][MAX_SLOT_NUM];
    unsigned int totalNumberSamples;
    unsigned int maxValue;
    unsigned int maxLocation;
    unsigned int totalTakeNum;
    unsigned int schedule;

    if (argc == 2)
    {
        if (!ReadFoundTime(argv[1], &foundTime, &totalNumberSamples))
        {
            return 1;
        }
    }
    else
    {
	#ifdef _DEBUG
        	printf("USAGE: %s <file>\n", argv[0]);
	#endif
        return 1;
    }

    cleanT = GetCleanT(foundTime);
    cleanF = GetCleanF(cleanT);
    GetMaxValLoc(cleanF, &maxValue, &maxLocation);

    if (!DetectBT(cleanF, foundTime->length, maxValue, totalNumberSamples))
    {
        return 0;
    }

    boundary = FindBoundary(foundTime, cleanF, &shouldTake, &totalTakeNum, maxLocation, maxValue);
    schedule = FindSchedule(foundTime, shouldTake, scores, maxLocation, totalTakeNum);

    #ifdef _DEBUG
      printf("outputrange =\n\n\t%d\t%d\n", boundary.min, boundary.max);
      printf("scores =\n\n");

      int i, j;

      for (i = 0; i < SCORE_SLOT; ++i)
      {
          for (j = 0; j < MAX_SLOT_NUM; ++j)
          {
              printf("\t%d",scores[i][j]);
          }
          printf("\n");
      }

      printf("we think the schedule is %d!!!\n", schedule);
    #endif

    DeallocateSamples(foundTime);
    DeallocateSamples(cleanT);
    DeallocateSamples(cleanF);
    DeallocateSamples(shouldTake);

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

void DeallocateSamples(Samples* sample)
{
    free(sample->samples);
    free(sample);
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

int DetectBT(Samples* sample, unsigned int found, unsigned int maxval, unsigned int total)
{
    float sampleRatio;
    float ave;

    sampleRatio = (float)found / total;
    ave = Mean(sample);

    if (sampleRatio > BELIEF_RAND10_MAX_P && maxval > ave * 2.5)
    {
        #ifdef _DEBUG
          printf("10 sample ratio %f, maxval %d, avg %f, We think there is BT!\n", sampleRatio, maxval, ave);
        #endif
        return 1;
    }
    else
    {
#ifdef _DEBUG
        printf("10 sample ratio %f, maxval %d, avg %f, We think there is NOOOOOOOOOO BT!\n", sampleRatio,
               maxval, ave);
#endif
        return 0;
    }
}

OutputBoundary FindBoundary(Samples* foundTime, Samples* cleanF, Samples** shouldTake,
                            unsigned int* totalTakeNum, unsigned int maxloc, unsigned int maxval)
{
    const unsigned int RIGHT_LEN = CTS_t + 2 * SIFS;
    const unsigned int RANGE_MARGIN = 10;
    const unsigned int EXPECTED_R_LEN = 6;
    const float TARGET_P = 0.005;
    const int incrementArray[2] = {-1,1};
    float measuredDirtyP;
    float oneP;
    unsigned int i;
    unsigned int D;
    unsigned int dummyTime;
    unsigned int stepdown;
    unsigned int currval;
    unsigned int currdist;
    int inc;
    int currpnt;
    uint8_t shouldTakeThisOne;
    OutputBoundary output;

    *shouldTake = AllocateSamples(foundTime->length);
    (*shouldTake)->length = foundTime->length;
    D = DUMMY_t + CTS_t + 2 * SIFS;

    for (i = 0; i < foundTime->length; ++i)
    {
        dummyTime = foundTime->samples[i] % BT_LEN;
        shouldTakeThisOne = 0;

        if (dummyTime < maxloc)
        {
            if (maxloc - dummyTime <= DUMMY_t + RANGE_MARGIN)
            {
                shouldTakeThisOne = 1;
            }
        }
        else
        {
            if (maxloc + BT_LEN - dummyTime <= DUMMY_t + RANGE_MARGIN)
            {
                shouldTakeThisOne = 1;
            }
        }
        
        if (dummyTime > maxloc)
        {
            if (dummyTime - maxloc <= RIGHT_LEN + RANGE_MARGIN)
            {
                shouldTakeThisOne = 1;
            }
        }
        else
        {
            if (dummyTime + BT_LEN - maxloc <= RIGHT_LEN + RANGE_MARGIN)
            {
                shouldTakeThisOne = 1;
            }
        }
        (*shouldTake)->samples[i] = shouldTakeThisOne;
    }
    
    *totalTakeNum = 0;
    
    for (i = 0; i < (*shouldTake)->length; ++i)
    {
        *totalTakeNum += (*shouldTake)->samples[i];
    }
    
    measuredDirtyP = (float)(foundTime->length - *totalTakeNum) / (BT_LEN - D);
    oneP = measuredDirtyP * EXPECTED_R_LEN;
    
    for (stepdown = 1; stepdown <= 4; ++stepdown)
    {
        if (powf(oneP, stepdown) < TARGET_P)
        {
            break;
        }
    }
    
    for (i = 0; i < 2; ++i)
    {
        inc = incrementArray[i];
        currpnt = maxloc;
        currval = maxval;
        currdist = 0;
        
        while (1)
        {
            currpnt = currpnt + inc;
            currdist = currdist + 1;
            
            if (currpnt < 0)
            {
                currpnt = BT_LEN - 1;
            }
            if (currpnt == BT_LEN)
            {
                currpnt = 0;
            }
            
            currval = cleanF->samples[currpnt];
            
            if (currval < (maxval - stepdown))
            {
#ifdef _DEBUG
                printf("broke from the loop when the val is %d, measureddirtyp %f, stepdown %d\n",
                       currval, measuredDirtyP, stepdown);
#endif
                break;
            }
        }
        
        if (i == 0)
            output.min = currpnt;
        else
            output.max = currpnt;
    }

    return output;
}

unsigned int FindSchedule(Samples* foundTime, Samples* shouldTake,
                  unsigned int scores[SCORE_SLOT][MAX_SLOT_NUM], unsigned int maxpos,
                  unsigned int totalTakeNum)
{
    const float MAGIC_C = 0.3;
    const unsigned int slenArray[] = {4,MAX_SLOT_NUM};
    int i, j, jj, jjj;
    unsigned int slen;
    unsigned int thisEstBgn;
    unsigned int slotIdx;
    
    unsigned int minEstFour;
    unsigned int minEstSix;
    
    
    /* Convert to an array literal? */
    for (i = 0; i < SCORE_SLOT; ++i)
    {
        for (j = 0; j < MAX_SLOT_NUM; ++j)
        {
            scores[i][j] = 0;
        }
    }
    
    for (j = 0; j < 2; ++j)
    {
        slen = slenArray[j];
        for (jj = 0; jj < slen; ++jj)
        {
            thisEstBgn = maxpos + jj * BT_LEN;       /* TODO: Need to subtract 1 */
            
            for (jjj = 0; jjj < foundTime->length; ++jjj)
            {
                if (shouldTake->samples[jjj] == 1)
                {
                    slotIdx = foundTime->samples[jjj] - thisEstBgn;
                    slotIdx = slotIdx % (slen * BT_LEN);
                    slotIdx = roundf(slotIdx / (float)BT_LEN);
                    
                    if (slotIdx != 0 && slotIdx != slen - 1 && slotIdx != slen)
                    {
                        scores[j][jj] = scores[j][jj] + 1;
                    }
                }
            }
        }
        /* TODO: Do we add 1 here */
        for (jj = slen; jj < MAX_SLOT_NUM; ++jj)
        {
            scores[j][jj] = 1000000;
        }
    }
    
    minEstFour = Min(scores[0], MAX_SLOT_NUM);
    minEstSix = Min(scores[1], MAX_SLOT_NUM);
    
    
    if (minEstFour > totalTakeNum * MAGIC_C && minEstSix > totalTakeNum * MAGIC_C)
    {
        return 2;
    }
    else
    {
        if ((float)minEstFour < 0.5 * minEstSix)
        {
            return 4;
        }
        else
        {
            return 6;
        }
    }
}

unsigned int Min(unsigned int* values, unsigned int size)
{
    unsigned int min = values[0];
    unsigned int i;
    
    for (i = 1; i < size; ++i)
    {
        if (values[i] < min)
        {
            
            min = values[i];
        }
    }
    
    return min;
}
