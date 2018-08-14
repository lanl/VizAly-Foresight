// Includes -------------------------------------------------------------------

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
    char oriFilePath[640], outputFilePath[640];
    int r5, r4, r3, r2, r1;
    int dataLength;
    int k;
    
    if (argc < 5)
    {
        printf ("Usage: %s [decimate 1/k data] [data type (-f or -d)] [data path] [dimension sizes...]\n", argv[0]);
        printf ("Example: %s 5 -f /home/dtao001/CESM-testdata/1800x3600/CLDLOW.dat 3600 1800\n", argv[0]);
        exit(0);
    }
    
    k = atoi(argv[1]);
    
    sprintf(oriFilePath, "%s", argv[3]);
    sprintf(outputFilePath, "%s.intrpl.out", oriFilePath);
    
    r1 = atoi(argv[4]);
    dataLength = r1;
    if (argc>=6)
    {
        r2 = atoi(argv[5]);
        dataLength *= r2;
    }
    if (argc>=7)
    {
        r3 = atoi(argv[6]);
        dataLength *= r3;
    }
    if (argc>=8)
    {
        r4 = atoi(argv[7]);
        dataLength *= r4;
    }
    if (argc>=9)
    {
        r5 = atoi(argv[8]);
        dataLength *= r5;
    }
    
    
    if (argv[2][0] == '-' && argv[2][1] == 'f')
    {
        
        float *data = (float*)malloc(dataLength*sizeof(float));
        float *decomp_data = (float*)malloc(dataLength*sizeof(float));
        
        FILE *pFile;
        pFile = fopen(oriFilePath, "rb");
        if (pFile == NULL)
        {
            printf ("Cannot open input file!\n");
            exit(1);
        }
        fread (data, sizeof(float), dataLength, pFile);
        fclose(pFile);
        
        int num_seg = (int)((dataLength-1)/k);
        
        int i, j;
        for (i = 0; i < num_seg; i++)
            for (j = 0; j < k; j++)
                decomp_data[i*k+j] = (data[(i+1)*k]-data[i*k])/k*j+data[i*k]; // linear interpolation
        
        for (j = 0; j < (dataLength-1)%k; j++)
            decomp_data[num_seg*k+j] = (data[dataLength-1]-data[num_seg*k])/(dataLength-1-num_seg*k)*j+data[num_seg*k];
        
        decomp_data[dataLength-1] = data[dataLength-1];
        
        pFile = fopen(outputFilePath, "wb");
        fwrite(decomp_data, sizeof(float), dataLength, pFile);
        fclose(pFile);
        
        free(data);
        free(decomp_data);
    }
    else
        if (argv[1][0] == '-' && argv[1][1] == 'd')
        {
            double *data = (double*)malloc(dataLength*sizeof(double));
            double *decomp_data = (double*)malloc(dataLength*sizeof(double));
            
            FILE *pFile;
            pFile = fopen(oriFilePath, "rb");
            if (pFile == NULL)
            {
                printf ("Cannot open input file!\n");
                exit(1);
            }
            fread (data, sizeof(double), dataLength, pFile);
            fclose(pFile);
            
            int num_seg = (int)((dataLength-1)/k);
            
            int i, j;
            for (i = 0; i < num_seg; i++)
                for (j = 0; j < k; j++)
                    decomp_data[i*k+j] = (data[(i+1)*k]-data[i*k])/k*j+data[i*k]; // linear interpolation
            
            for (j = 0; j < (dataLength-1)%k; j++)
                decomp_data[num_seg*k+j] = (data[dataLength-1]-data[num_seg*k])/(dataLength-1-num_seg*k)*j+data[num_seg*k];
            
            decomp_data[dataLength-1] = data[dataLength-1];
            
            
            pFile = fopen(outputFilePath, "wb");
            fwrite(decomp_data, sizeof(double), dataLength, pFile);
            fclose(pFile);
            
            free(data);
            free(decomp_data);
        }
        else
        {
            printf ("Usage: %s [decimate 1/k data] [data type (-f or -d)] [data path] [dimension sizes...]\n", argv[0]);
            printf ("Example: %s 5 -f /home/dtao001/CESM-testdata/1800x3600/CLDLOW.dat 3600 1800\n", argv[0]);
            exit(0);
        }
    
    return 0;
}
