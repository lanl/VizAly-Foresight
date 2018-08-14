//CHECK:
//What happens when ECQBits==1, or ECQBits==0 or ECQBits<0?
//Rounding? Scale eb by 0.99?

//FIX: convert all bit operations to unions (SHOULD BE FIXED)
//Possible improvement: Optimize bookkeeping bits
//Possible improvement: Guess the type (C/UC, Sparse/Not)
//Possible improvement: Get rid of writing/reading some of the indexes to in/out buffers
//Possible improvement: Get rid of all debug stuff, including Makefile debug flags
//Possible improvement: Get rid of "compressedBytes"
//Possible improvement: SparseCompressed, ECQBits=2: 1's and -1's can be represented by just 0 and 1, instead 10 and 11. 
//Possible improvement: SparseCompressed, ECQBits>2: Again: 1: 10, -1:11, Others: 0XX...XX 
//Possible improvement: WriteBitsFast: maybe remove some masks?
//Possible improvement: WriteBitsFast: Get rid of multiple calls!



#ifndef PASTRITOOLS_H
#define PASTRITOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h> //Just for debugging purposes!

#define MAX_PS_SIZE 100
#define MAX_ELEMENTS 10000
#define MAX_BUFSIZE 160000  //Should be a multiple of 8
#define D_W 0 //Debug switch: Write (input block)
#define D_R 0 //Debug switch: Read (compressed block)
#define D_G 0 //Debug switch: General
#define D_C 0 //Debug switch: C
//#define DEBUG 1 //Debug switch

//#define BOOKKEEPINGBITS 0 //Currently unused
//#define BOOKKEEPINGBITS 120 //Includes: mode, indexOffsets, compressedBytes, Pb_, ECQBits_ (8+64+32+8+8) 
//BOOKKEEPINGBITS is defined here, because if P & S is going to be used, they appear just after the bookkeeping part.
//This allows us to write P and S directly onto using outBuf.
  

// IMPORTANT NOTE:
//Read/Write up to 56 bits.
//More than that is not supported!


typedef struct parameters{
  double eb; //Error Bound
  
  uint8_t bf[4]; //Orbital types (basis function types). Typically in range [0,3]
  int idxRange[4];  //Ranges of indexes. idxRange[i]=(bf[i]+1)*(bf[i]+2)/2;
  
  int sbSize; //=idxRange[2]*idxRange[3];
  int sbNum;  //=idxRange[0]*idxRange[1];
  int bSize; //=sbSize*sbNum;
  
  uint16_t idxOffset[4]; //Index offset values
}parameters;

typedef union u_UI64I64D{
  uint64_t ui64;
  int64_t i64;
  double d;
} u_UI64I64D;

//Add or sub 0.5, depending on the sign:
static inline double ceil_FastD(double x){
  u_UI64I64D u1,half;
  u1.d=x;
  
  half.d=0.5;
  
  //printf("ceil_FastD:\nx=%lf  x=0x%lx\n",x,(*((uint64_t *)(&x))));
  //printf("sign(x):0x%lx\n", x);
  //printf("0.5:0x%lx\n", (*((uint64_t *)(&half))));
  half.ui64 |= (u1.ui64 & (int64_t)0x8000000000000000);
  //printf("sign(x)*0.5:0x%lx\n", (*((uint64_t *)(&half))));
  return x + half.d;
}



static inline double abs_FastD(double x){
  u_UI64I64D u1;
  u1.d=x;
  //(*((uint64_t *)(&x)))&=(int64_t)0x7FFFFFFFFFFFFFFF;
  u1.ui64&=(int64_t)0x7FFFFFFFFFFFFFFF;
  return u1.d;
}

static inline int64_t abs_FastI64(int64_t x){
  return (x^((x&(int64_t)0x8000000000000000)>>63))+((x&(int64_t)0x8000000000000000)!=0);
}


//Returns the min. bits needed to represent x.
//Same as: ceil(log2(abs(x))) 
//Actually to be completely safe, it correspond to: ceil(log2(abs(i)+1))+0.1
//+0.1 was for fixing rounding errors
//REMEMBER: To represent the whole range [-x:x], the number of bits required is bitsNeeded(x)+1
static inline int bitsNeededD(double x){
  u_UI64I64D u1;
  u1.d=x;
  return (((u1.ui64<<1)>>53)-1022) & (((x!=0)<<31)>>31);
  //uint64_t *ptr;
  //x=abs_FastD(x);
  //ptr=(uint64_t *)&x;
  //uint64_t retVal;
  //retVal=(((*ptr)<<1)>>53)-1022;  //This is good for all x except when x=0
  //return retVal&(((x!=0)<<31)>>31);  
}
static inline int bitsNeededUI64(uint64_t x){
  int shift;
  int res=0;
  
  //Get the absolute value of x:
  //x=(x^((x&(int64_t)0x8000000000000000)>>63))+((x&(int64_t)0x8000000000000000)!=0);
  //x=abs_FastI64(x);
  
  //printf("%d\n",(x&(uint64_t)0xFFFFFFFF00000000)!=0);
  shift=(((x&(uint64_t)0xFFFFFFFF00000000)!=0)*32);
  x>>=shift;
  res+=shift;
  
  //printf("%d\n",(x&(uint64_t)0x00000000FFFF0000)!=0);
  shift=(((x&(uint64_t)0x00000000FFFF0000)!=0)*16);
  x>>=shift;
  res+=shift;
  
  //printf("%d\n",(x&(uint64_t)0x000000000000FF00)!=0);
  shift=(((x&(uint64_t)0x000000000000FF00)!=0)*8);
  x>>=shift;
  res+=shift;
  
  //printf("%d\n",(x&(uint64_t)0x00000000000000F0)!=0);
  shift=(((x&(uint64_t)0x00000000000000F0)!=0)*4);
  x>>=shift;
  res+=shift;
  
  //printf("%d\n",(x&(uint64_t)0x000000000000000C)!=0);
  shift=(((x&(uint64_t)0x000000000000000C)!=0)*2);
  x>>=shift;
  res+=shift;
  
  //printf("%d\n",(x&(uint64_t)0x0000000000000002)!=0);
  shift=((x&(uint64_t)0x0000000000000002)!=0);
  x>>=shift;
  res+=shift;
  
  //printf("%d\n",(x&(uint64_t)0x0000000000000001)!=0);
  shift=((x&(uint64_t)0x0000000000000001)!=0);
  x>>=shift;
  res+=shift;
  
  //printf("BITS NEEDED: %d\n",res);
  return res;
}
static inline int bitsNeededI64(int64_t x){
  uint64_t ux;
  ux=abs_FastI64(x);
  return bitsNeededUI64(ux);
}

//Implementations(They are inline, so they should be in this header file)

static inline int myEndianType(){ //Should work for most cases. May not work at mixed endian systems.
  uint64_t n=1;
  if (*(uint8_t *)&n == 1){
    //cout<<"Little-Endian"<<endl;
    return 0;  //0 for little endian
  }
  else{
    //cout<<"Big-Endian"<<endl;
    return 1; //1 for big endian
  }
}

static inline void flipBytes64b(uint64_t *dataPtr){
  uint8_t *tempA;
  uint8_t temp8b;
  tempA=(uint8_t*)dataPtr;
  temp8b=tempA[7];
  tempA[7]=tempA[0];
  tempA[0]=temp8b;
  temp8b=tempA[6];
  tempA[6]=tempA[1];
  tempA[1]=temp8b;
  temp8b=tempA[5];
  tempA[5]=tempA[2];
  tempA[2]=temp8b;
  temp8b=tempA[4];
  tempA[4]=tempA[3];
  tempA[3]=temp8b;
  return;
}

//WARNING: readBits works properly only on Little Endian machines! (For Big Endians, some modifications are needed)

static inline uint64_t readBits_UI64(uint8_t* buffer,uint64_t *bitPosPtr,uint8_t numBits){ // numBits must be in range [0:56]
    uint64_t mask = ((uint64_t)0x0000000000000001<<numBits)-1;
    //cout<<"bitPos:"<<(*bitPosPtr)<<"\tbitPos>>3:"<<(*bitPosPtr>>3)<<endl;
    uint64_t temp64b = *(uint64_t*)(buffer + ( *bitPosPtr >> 3)); 
    //NOTE: bitPos>>3 is the same as bitPos/8
    temp64b >>= (*bitPosPtr) & (uint64_t)0x0000000000000007;
    
    //cout<<endl;
    //cout<<"bitpos>>3:"<<(bitPos>>3)<<" bitPos&0x7:"<<(bitPos & 0x00000007)<<" bitPos%8:"<<(bitPos%8)<<endl;
    //cout<<"Read:"<<(temp64b & mask)<<" temp64b:"<<temp64b<<" Mask:"<<mask<<" numBits:"<<numBits<<endl;
    
    (*bitPosPtr) += numBits;
    return (temp64b & mask);
};

static inline int64_t readBits_I64(uint8_t* buffer,uint64_t *bitPosPtr,uint8_t numBits){ // numBits must be in range [0:56]
  int64_t val;
  val=readBits_UI64(buffer,bitPosPtr,numBits);//Read value
  int64_t shiftAmount=64-numBits;
  val=(val<<shiftAmount)>>shiftAmount;//Sign correction
  return val;
}

//WARNING: readBits_EndianSafe is not tested on Big-Endian machines
static inline uint64_t readBits_EndianSafe(uint8_t* buffer,uint64_t *bitPosPtr,uint8_t numBits){ // numBits must be in range [0:56]
    uint64_t mask = ((uint64_t)0x0000000000000001<<numBits)-1;
    uint64_t temp64b = *(uint64_t*)(buffer + ((*bitPosPtr)>>3)); 
    //NOTE: (*bitPosPtr)>>3 is the same as (*bitPosPtr)/8
    if(myEndianType())
      flipBytes64b(&temp64b);
    temp64b >>= (*bitPosPtr) & (uint64_t)0x0000000000000007;
    (*bitPosPtr) += numBits;
    return temp64b & mask;
};

//WARNING: writeBits_Fast works properly only on Little Endian machines! (For Big Endians, some modifications are needed)
//The buffer should be initialized as 0's for this to work!
//Also, the range of data is not checked!(If data exceeds numBits, it may be cause problems)
static inline void writeBits_Fast(uint8_t* buffer,uint64_t *bitPosPtr,uint8_t numBits,int64_t data){
    //if(DEBUG){printf("writeBits_Fast: data:0x%lx %ld\n",data,data);} //DEBUG
    //if(DEBUG){printf("writeBits_Fast: numBits:0x%lx %ld\n",numBits,numBits);} //DEBUG
    uint64_t mask = ((uint64_t)0x0000000000000001<<numBits)-1;
    //if(DEBUG){printf("writeBits_Fast: mask:0x%lx %ld\n",mask,mask);} //DEBUG
    //if(DEBUG){printf("writeBits_Fast: data&mask:0x%lx %ld\n",((*(uint64_t*)&data)&mask),((*(uint64_t*)&data)&mask));} //DEBUG
    
    //if(DEBUG){printf("writeBits_Fast: buffer_O:0x%lx\n",*(uint64_t*)(buffer + ((*bitPosPtr)>>3)));} //DEBUG
    *(uint64_t*)(buffer + ((*bitPosPtr)>>3)) |= ((*(uint64_t*)&data)&mask) << ((*bitPosPtr) & (uint64_t)0x0000000000000007);
    //if(DEBUG){printf("writeBits_Fast: buffer_N:0x%lx\n",*(uint64_t*)(buffer + ((*bitPosPtr)>>3)));} //DEBUG

    
    (*bitPosPtr) += numBits;
};

//WARNING: writeBits_EndianSafe is not tested on Big-Endian machines
static inline void writeBits_EndianSafe(uint8_t* buffer,uint64_t *bitPosPtr,uint8_t numBits,uint64_t data){
    uint64_t mask = ((uint64_t)0x0000000000000001<<numBits)-1;
    data=data&mask;
    uint64_t temp64b_inBuffer=*(uint64_t*)(buffer + ((*bitPosPtr)>>3));
    uint64_t temp64b_outBuffer=data << ((*bitPosPtr) & (uint64_t)0x0000000000000007);
    if(myEndianType()){
      flipBytes64b(&temp64b_inBuffer);
    }
    temp64b_outBuffer |= temp64b_inBuffer;
    if(myEndianType()){
      flipBytes64b(&temp64b_outBuffer);
    }
    *(uint64_t*)(buffer + ((*bitPosPtr)>>3))=temp64b_outBuffer;  // "|=" may also work
    (*bitPosPtr) += numBits;
};

static inline int pastriCompress(uint8_t *inBuf,parameters *p,uint8_t *outBuf,int *numOutBytes){
  int ECQ1s, ECQOthers;
  int nonZeros;
  //uint64_t temp1; //DEBUG
  //Preprocess by calculating some parameters:
  //Calculate sbSize, sbNum, etc.:
  p->idxRange[0]=(p->bf[0]+1)*(p->bf[0]+2)/2;
  p->idxRange[1]=(p->bf[1]+1)*(p->bf[1]+2)/2;
  p->idxRange[2]=(p->bf[2]+1)*(p->bf[2]+2)/2;
  p->idxRange[3]=(p->bf[3]+1)*(p->bf[3]+2)/2;
  p->sbSize=p->idxRange[2]*p->idxRange[3];
  p->sbNum=p->idxRange[0]*p->idxRange[1];
  p->bSize=p->sbSize*p->sbNum;
  
  //if(DEBUG){printf("Parameters: bfs:%d %d %d %d eb:%.16lf\n",p->bf[0],p->bf[1],p->bf[2],p->bf[3],p->eb);}  //DEBUG
  //if(DEBUG){printf("Parameters: idxRanges:%d %d %d %d\n",p->idxRange[0],p->idxRange[1],p->idxRange[2],p->idxRange[3]);} //DEBUG
  //if(DEBUG){printf("Parameters: sbSize:%d sbNum:%d bSize:%d\n",p->sbSize,p->sbNum,p->bSize); }//DEBUG
  
  //Read some indexes and calculate index offsets:
  //A full read of the indexes would look like this:
  //i0Ptr=(uint16_t*)(inBuf);
  //i1Ptr=(uint16_t*)(inBuf+p->bSize*2);
  //i2Ptr=(uint16_t*)(inBuf+p->bSize*4);
  //i3Ptr=(uint16_t*)(inBuf+p->bSize*6);
  //But, we do not need all the indexes as they are just consecutive numbers.
  //We only need the offset information here.
  
  uint16_t *idx0,*idx1,*idx2,*idx3;
  idx0=(uint16_t*)(inBuf           );
  idx1=(uint16_t*)(inBuf+p->bSize*2);
  idx2=(uint16_t*)(inBuf+p->bSize*4);
  idx3=(uint16_t*)(inBuf+p->bSize*6);
  
 
  p->idxOffset[0]=(idx0[0]/p->idxRange[0])*p->idxRange[0];
  p->idxOffset[1]=(idx1[0]/p->idxRange[1])*p->idxRange[1];
  p->idxOffset[2]=(idx2[0]/p->idxRange[2])*p->idxRange[2];
  p->idxOffset[3]=(idx3[0]/p->idxRange[3])*p->idxRange[3];
  
  //if(DEBUG){printf("Parameters: idxOffsets:%d %d %d %d\n",p->idxOffset[0],p->idxOffset[1],p->idxOffset[2],p->idxOffset[3]);} //DEBUG
  
  int i,sb;
  double *data;
  data=(double*)(inBuf+p->bSize*8);
  //if(DEBUG){printf("data>eb :\n");for(i=0;i<p->bSize;i++){if(abs_FastD(data[i])>p->eb)printf("data[%d]=%.6e\n",i,data[i]);} }//DEBUG
  //if(DEBUG){printf("All Data:\n");for(i=0;i<p->bSize;i++){printf("data[%d]=%.6e\n",i,data[i]);}} //DEBUG
  
  
  //Start by finding the pattern.
  //Find the extremum point:
  double absExt=0; //Absolute value of Extremum
  int extIdx=-1; //Index of Extremum
  nonZeros=0;
  for(i=0;i<p->bSize;i++){
    //printf("data[%d] = %.16lf\n",i,data[i]);//DEBUG
    if(abs_FastD(data[i])>p->eb){
      nonZeros++;
      //if(DEBUG)printf("data[%d]:%.6e\n",i,data[i]); //DEBUG
    }
    if(abs_FastD(data[i])>absExt){
      absExt=abs_FastD(data[i]);
      extIdx=i;
    }
  }
  int patternIdx; //Starting Index of Pattern
  patternIdx=extIdx/p->sbSize*p->sbSize;
  double patternExt=data[extIdx];
  double binSize=2*p->eb;
  
  //if(DEBUG){printf("Extremum  : data[%d] = %.6e\n",extIdx,patternExt);} //DEBUG
  //if(DEBUG){printf("patternIdx: %d\n",patternIdx);} //DEBUG
  
  //if(DEBUG){for(i=0;i<p->sbSize;i++){printf("pattern[%d]=data[%d]=%.6e ???:%.6e Quantized:%d\n",i,patternIdx+i,data[patternIdx+i],ceil_FastD(data[patternIdx+i]/binSize),(int)ceil_FastD(data[patternIdx+i]/binSize)  );}   }//DEBUG
  
  
  
  //int64_t *patternQ=(int64_t*)(outBuf+15);  //Possible Improvement!
  int64_t patternQ[MAX_PS_SIZE];
  int64_t scalesQ[MAX_PS_SIZE];
  
  for(i=0;i<p->sbSize;i++){
    patternQ[i]=(int64_t)ceil_FastD(data[patternIdx+i]/binSize);
    if(D_W){printf("patternQ[%d]=%ld\n",i,patternQ[i]);}
  }
  
  
  int patternBits=bitsNeededD((abs_FastD(patternExt)/binSize)+1)+1;
  int scaleBits=patternBits;
  double scalesBinSize=1/(double)(((uint64_t)1<<(scaleBits-1))-1);
  //if(DEBUG){printf("(patternExt/binSize)+1: %.6e\n",(patternExt/binSize)+1);} //DEBUG
  //if(DEBUG){printf("scaleBits=patternBits: %d\n",scaleBits);} //DEBUG
  if(D_W){printf("scalesBinSize: %.6e\n",scalesBinSize);} //DEBUG
  
  
  //Calculate Scales.
  //The index part of the input buffer will be reused to hold Scale, Pattern, etc. values.
  int localExtIdx=extIdx%p->sbSize; //Local extremum index. This is not the actual extremum of the current sb, but rather the index that correspond to the global (block) extremum.
  //int64_t *scalesQ=(int64_t*)(outBuf+15+p->sbSize*8);  //Possible Improvement!
  int patternExtZero=(patternExt==0);
  //if(DEBUG){printf("patternExtZero: %d\n",patternExtZero);} //DEBUG
  for(sb=0;sb<p->sbNum;sb++){
    //scales[sb]=data[sb*p->sbSize+localExtIdx]/patternExt;
    //scales[sb]=patternExtZero ? 0 : data[sb*p->sbSize+localExtIdx]/patternExt;
    //assert(scales[sb]<=1);
    //if(DEBUG){printf("scales[%d]=[%d/%d]: %lf Quantized:%d \n",sb,sb*p->sbSize+localExtIdx,extIdx,patternExtZero ? 0 : data[sb*p->sbSize+localExtIdx]/patternExt,(int64_t)ceil_FastD((patternExtZero ? 0 : data[sb*p->sbSize+localExtIdx]/patternExt)/scalesBinSize));}//DEBUG
    scalesQ[sb]=(int64_t)ceil_FastD((patternExtZero ? 0 : data[sb*p->sbSize+localExtIdx]/patternExt)/scalesBinSize);
    if(D_W){printf("scalesQ[%d]=%ld\n",sb,scalesQ[sb]);}
  }
  //if(DEBUG){for(i=0;i<p->sbSize;i++){printf("scalesQ[%d]=%ld \n",i,scalesQ[i]);}} //DEBUG

  //int64_t *ECQ=(int64_t*)(outBuf+p->bSize*8); //ECQ is written into outBuf, just be careful when handling it.
  int64_t ECQ[MAX_ELEMENTS];
  //uint64_t wVal;
  
  uint64_t ECQExt=0; //The extremum of ECQ
  int _1DIdx;
  ECQ1s=0;
  ECQOthers=0;
  for(sb=0;sb<p->sbNum;sb++){
    for(i=0;i<p->sbSize;i++){
      _1DIdx=sb*p->sbSize+i;
      ECQ[_1DIdx]=(int64_t)ceil_FastD( (scalesQ[sb]*patternQ[i]*scalesBinSize*binSize-data[_1DIdx])/binSize );
      double absECQ=abs_FastD(ECQ[_1DIdx]);
      if(absECQ > ECQExt)
        ECQExt=absECQ;
      //if(DEBUG){printf("EC[%d]: %.6e Quantized:%ld \n",_1DIdx,(scalesQ[sb]*patternQ[i]*scalesBinSize*binSize-data[_1DIdx]),ECQ[_1DIdx]);} //DEBUG
      switch (ECQ[_1DIdx]){
        case 0:
          //ECQ0s++; //Currently not needed
          break;
        case 1:
          ECQ1s++;
          break;
        case -1:
          ECQ1s++;
          break;
        default:
          ECQOthers++;
          break;
      }
    }
  }
  
  //DEBUG: Self-check. Remove this later.
  for(sb=0;sb<p->sbNum;sb++){
    for(i=0;i<p->sbSize;i++){
      _1DIdx=sb*p->sbSize+i;
      double decompressed=scalesQ[sb]*patternQ[i]*scalesBinSize*binSize-ECQ[_1DIdx]*binSize;
      if(abs_FastD(decompressed-data[_1DIdx])>(p->eb)){
        printf("p->eb=%.6e\n",p->eb);
        printf("data[%d]=%.6e decompressed[%d]=%.6e diff=%.6e\n",_1DIdx,data[_1DIdx],_1DIdx,decompressed,abs_FastD(data[_1DIdx]-decompressed));
        assert(0);
      }
    }
  }
  
  //if(DEBUG){printf("|ECQExt|:%d bitsNeededUI64(ECQExt):%d\n",ECQExt,bitsNeededUI64(ECQExt)+1);} //DEBUG

  
  int ECQBits=bitsNeededUI64(ECQExt)+1;
  
  int _1DIdxBits=bitsNeededUI64(p->bSize);
  (*numOutBytes)=0;
  
  //Encode: 3 options:
  //Compressed, Sparse
  //Compressed, Non-Sparse
  //Uncompressed, Sparse : numZeros*
  
  //int UCSparseBits;  //Uncompressed, Sparse bits. Just like the original GAMESS data. Includes: mode, indexOffsets, nonZeros, indexes, data
  //int UCNonSparseBits;  //Uncompressed, NonSparse bits. Includes: mode, indexOffsets, data
  int CSparseBits;  //Includes: mode, indexOffsets, compressedBytes, patternBits, ECQBits,numOutliers,P, S, {Indexes(Sparse), ECQ}
  int CNonSparseBits;  //Includes: mode, indexOffsets, compressedBytes, patternBits, ECQBits,P, S, {ECQ}
  //int BOOKKEEPINGBITS=120; //Includes: mode, indexOffsets, compressedBytes, patternBits, ECQBits (8+64+32+8+8) //Moved to much earlier!
    
  //Consider: ECQ0s, ECQ1s, ECQOthers. Number of following values in ECQ: {0}, {1,-1}, { val<=-2, val>=2}
  //ECQ0s is actually not needed, but others are needed.

  //UCSparseBits = 8 + 16 +64+ nonZeros*128;
  //UCNonSparseBits = 8 + 64 + p->bSize*64;
  int numOutliers=ECQ1s+ECQOthers;
  if(ECQBits==2){
    CSparseBits = 8+64+32+8+8+16 + patternBits*p->sbSize + scaleBits*p->sbNum + ECQ1s*(_1DIdxBits+2);
    CNonSparseBits = 8+64+32+8+8+ patternBits*p->sbSize + scaleBits*p->sbNum + p->bSize + ECQ1s ;  //Or: ECQ0s+ECQ1s*2;
  }else{ //ECQBits>2
    //CSparseBits = 8+64+32+8+8+16 + patternBits*p->sbSize + scaleBits*p->sbNum + (ECQ1s+ECQOthers)*_1DIdxBits + ECQ1s*3+ ECQOthers*2;
    CSparseBits = 8+64+32+8+8+16 + patternBits*p->sbSize + scaleBits*p->sbNum + (numOutliers)*(_1DIdxBits+2) + ECQ1s;
    CNonSparseBits = 8+64+32+8+8+ patternBits*p->sbSize + scaleBits*p->sbNum + p->bSize + ECQ1s*2 + ECQOthers; //Or: ECQ0s+ECQ1s*3+ECQOthers*2
  }
  
  //int UCSparseBytes=(UCSparseBits+7)/8; 
  int UCSparseBytes=1 + 8 + 2 + nonZeros*16; 
  //int UCNonSparseBytes=(UCNonSparseBits+7)/8; 
  int UCNonSparseBytes=1 + 8 + p->bSize*8; 
  
  int CSparseBytes=(CSparseBits+7)/8; 
  int CNonSparseBytes=(CNonSparseBits+7)/8; 
  uint64_t bitPos=0;
  uint64_t bytePos=0;
  int i0,i1,i2,i3;
  
  *(uint16_t*)(&outBuf[1])=p->idxOffset[0];
  *(uint16_t*)(&outBuf[3])=p->idxOffset[1];
  *(uint16_t*)(&outBuf[5])=p->idxOffset[2];
  *(uint16_t*)(&outBuf[7])=p->idxOffset[3];
    
  if(D_W){printf("ECQ0s:%d ECQ1s:%d ECQOthers:%d Total:%d\n",p->bSize-ECQ1s-ECQOthers,ECQ1s,ECQOthers,p->bSize);} //DEBUG
  if(D_W){printf("numOutliers:%d\n",numOutliers);} //DEBUG
  
  //****************************************************************************************
  //if(0){ //DEBUG
  //W:UCSparse
  if((UCSparseBytes<UCNonSparseBytes) && (UCSparseBytes<CSparseBytes) && (UCSparseBytes<CNonSparseBytes) ){ 
    //Uncompressed, Sparse bits. Just like the original GAMESS data. Includes: mode, indexOffsets, nonZeros, indexes, data
    if(D_G){printf("UCSparse\n");} //DEBUG
    if(D_G)printf("ECQBits:%d\n",ECQBits); //DEBUG
    outBuf[0]=0; //mode
    
    *(uint16_t*)(&outBuf[9])=nonZeros;
    bytePos=11;//0:mode, 1-8:indexOffsets 9-10:NonZeros. So start from 11.
    
    for(i0=0;i0<p->idxRange[0];i0++)
      for(i1=0;i1<p->idxRange[1];i1++)
        for(i2=0;i2<p->idxRange[2];i2++)
          for(i3=0;i3<p->idxRange[3];i3++){
            _1DIdx=p->idxRange[3]*(i2+p->idxRange[2]*(i1+i0*p->idxRange[1]))+i3;
            if(abs_FastD(data[_1DIdx])>p->eb){
              *(uint16_t*)(&outBuf[bytePos])=i0+1+p->idxOffset[0];
              bytePos+=2;
              *(uint16_t*)(&outBuf[bytePos])=i1+1+p->idxOffset[1];
              bytePos+=2;
              *(uint16_t*)(&outBuf[bytePos])=i2+1+p->idxOffset[2];
              bytePos+=2;
              *(uint16_t*)(&outBuf[bytePos])=i3+1+p->idxOffset[3];
              bytePos+=2;
              
              *(double*)(&outBuf[bytePos])=data[_1DIdx];
              bytePos+=8;
            }
          }
    
    *numOutBytes=UCSparseBytes;
    if(D_G)printf("UCSparseBytes:%d \n",UCSparseBytes); //DEBUG
    
  //****************************************************************************************
  //}else if(0){ //DEBUG
  //W:UCNonSparse
  }else if((UCNonSparseBytes<UCSparseBytes) && (UCNonSparseBytes<CSparseBytes) && (UCNonSparseBytes<CNonSparseBytes) ){ 
    //Uncompressed, NonSparse bits. Includes: mode, indexOffsets, data
    if(D_G){printf("UCNonSparse\n");} //DEBUG
    if(D_G)printf("ECQBits:%d\n",ECQBits); //DEBUG
    outBuf[0]=1; //mode
    
    //void *memcpy(void *dest, const void *src, size_t n){...};
    memcpy(&outBuf[9], &inBuf[p->bSize*8], UCNonSparseBytes-9);
    *numOutBytes=UCNonSparseBytes;
    if(D_G)printf("UCNonSparseBytes:%d \n",UCNonSparseBytes); //DEBUG
    /*
    for(i=0;i<UCNonSparseBytes-17;i++){
      printf("%d ",inBuf[p->bSize*8+i]);
    }
    printf("\n");
    for(i=0;i<UCNonSparseBytes-17;i++){
      printf("%d ",outBuf[17+i]);
    }
    printf("\n");
    */
  //****************************************************************************************
  //}else if(1){ //DEBUG
  //W:CSparse
  }else if((CSparseBytes<UCNonSparseBytes) && (CSparseBytes<UCSparseBytes) && (CSparseBytes<CNonSparseBytes) ){ 
    //Includes: mode, indexOffsets, compressedBytes, patternBits, ECQBits,numOutliers,P, S, {Indexes(Sparse), ECQ}

    if(D_G){printf("CSparse\n");} //DEBUG
    if(D_G)printf("ECQBits:%d\n",ECQBits); //DEBUG
    //if(DEBUG){printf("patternBits:%d _1DIdxBits:%d\n",patternBits,_1DIdxBits);} //DEBUG
    outBuf[0]=2; //mode
    //outBuf bytes [1:8] are indexOffsets, which are already written
    //outBuf bytes [9:12] are reserved for compressedBytes.
    outBuf[13]=patternBits;
    outBuf[14]=ECQBits;
    
    //bitPos=15*8; //Currently, we are at the end of 15th byte.
    *(uint16_t*)(&outBuf[15])=numOutliers;
    bitPos=17*8; //Currently, we are at the end of 17th byte.
    
    //if(DEBUG){printf("bitPos_B:%ld\n",bitPos);} //DEBUG

    for(i=0;i<p->sbSize;i++){
      writeBits_Fast(outBuf,&bitPos,patternBits,patternQ[i]);//Pattern point
    }
    //if(DEBUG){printf("bitPos_P:%ld\n",bitPos);} //DEBUG
    for(i=0;i<p->sbNum;i++){
      writeBits_Fast(outBuf,&bitPos,scaleBits,scalesQ[i]);//Scale
    }
    //if(DEBUG){printf("bitPos_S:%ld\n",bitPos);} //DEBUG
    //if(DEBUG)printf("ECQBits:%d\n",ECQBits);
    switch(ECQBits){
      case 2:
        for(i=0;i<p->bSize;i++){
          switch(ECQ[i]){
            case 0:
              break;
            case 1:
              //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x0\n",i,ECQ[i]); //DEBUG
              writeBits_Fast(outBuf,&bitPos,_1DIdxBits,i);
              //writeBits_Fast(outBuf,&bitPos,2,0x10);
              //writeBits_Fast(outBuf,&bitPos,2,0);//0x00
              //writeBits_Fast(outBuf,&bitPos,2,0);//0x00
              writeBits_Fast(outBuf,&bitPos,1,0);//0x00
              break;
            case -1:
              //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x1\n",i,ECQ[i]); //DEBUG
              writeBits_Fast(outBuf,&bitPos,_1DIdxBits,i);
              //writeBits_Fast(outBuf,&bitPos,2,0x11);
              //writeBits_Fast(outBuf,&bitPos,2,1);//0x01
              //writeBits_Fast(outBuf,&bitPos,1,0);
              writeBits_Fast(outBuf,&bitPos,1,1);
              break;
            default:
              assert(0);
              break;
          }
        }
        break;
      default: //ECQBits>2
      for(i=0;i<p->bSize;i++){
        switch(ECQ[i]){
          case 0:
            break;
          case 1:
            //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x00\n",i,ECQ[i]); //DEBUG
            writeBits_Fast(outBuf,&bitPos,_1DIdxBits,i);
            //writeBits_Fast(outBuf,&bitPos,3,0);//0x000
            //writeBits_Fast(outBuf,&bitPos,1,0);
            writeBits_Fast(outBuf,&bitPos,1,0);
            writeBits_Fast(outBuf,&bitPos,1,0);
            break;
          case -1:
            //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x01\n",i,ECQ[i]); //DEBUG
            writeBits_Fast(outBuf,&bitPos,_1DIdxBits,i);
            //writeBits_Fast(outBuf,&bitPos,3,1);//0x001
            //writeBits_Fast(outBuf,&bitPos,1,0);
            writeBits_Fast(outBuf,&bitPos,1,0);
            writeBits_Fast(outBuf,&bitPos,1,1);
            break;
          default:
            //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x1 0x%lx\n",i,ECQ[i],ECQ[i]); //DEBUG
            writeBits_Fast(outBuf,&bitPos,_1DIdxBits,i);
            //writeBits_Fast(outBuf,&bitPos,2+ECQBits,((uint64_t)0x11<<ECQBits)|ECQ[i]);
            //writeBits_Fast(outBuf,&bitPos,2+ECQBits,(ECQ[i]&((uint64_t)0x00<<ECQBits))|((uint64_t)0x01<<ECQBits));
            //writeBits_Fast(outBuf,&bitPos,1,0);
            writeBits_Fast(outBuf,&bitPos,1,1);
            writeBits_Fast(outBuf,&bitPos,ECQBits,ECQ[i]);
            break;
        }
      }
      break;
    }
    
    //if(DEBUG){printf("bitPos_E:%ld\n",bitPos);} //DEBUG
    if(D_C){if(!((ECQBits>=2)||((ECQBits==1) && (numOutliers==0)))){printf("ERROR: ECQBits:%d numOutliers:%d This should not have happened!\n",ECQBits,numOutliers);assert(0);}} //DEBUG
          

    uint32_t bytePos=(bitPos+7)/8;
    *(uint32_t*)(&outBuf[9])=bytePos;
    
    if(D_G)printf("bitPos:%ld CSparseBits:%d bytePos:%d CSparseBytes:%d\n",bitPos,CSparseBits,bytePos,CSparseBytes); //DEBUG
    
    
  //****************************************************************************************
  //W:CNonSparse
  }else { 
    //Includes: mode, indexOffsets, compressedBytes, patternBits, ECQBits,P, S, {ECQ}
    if(D_G){printf("CNonSparse\n");} //DEBUG
    if(D_G)printf("ECQBits:%d\n",ECQBits); //DEBUG
    //if(DEBUG){printf("patternBits:%d _1DIdxBits:%d\n",patternBits,_1DIdxBits);} //DEBUG
    outBuf[0]=3; //mode
    //outBuf bytes [1:8] are indexOffsets, which are already written
    //outBuf bytes [9:12] are reserved for compressedBytes.
    outBuf[13]=patternBits;
    outBuf[14]=ECQBits;
    
    bitPos=15*8; //Currently, we are at the end of 15th byte.
    
    //if(DEBUG){printf("bitPos_B:%ld\n",bitPos);} //DEBUG

    for(i=0;i<p->sbSize;i++){
      writeBits_Fast(outBuf,&bitPos,patternBits,patternQ[i]);//Pattern point
    }
    //if(DEBUG){printf("bitPos_P:%ld\n",bitPos);} //DEBUG
    for(i=0;i<p->sbNum;i++){
      writeBits_Fast(outBuf,&bitPos,scaleBits,scalesQ[i]);//Scale
    }
    //if(DEBUG){printf("bitPos_S:%ld\n",bitPos);} //DEBUG
    //if(DEBUG)printf("ECQBits:%d\n",ECQBits);
    switch(ECQBits){
      case 2:
        for(i=0;i<p->bSize;i++){
          switch(ECQ[i]){
            case 0:
              //if(DEBUG)printf("Index:%d ECQ:%d Written:0x1\n",i,ECQ[i]); //DEBUG
              writeBits_Fast(outBuf,&bitPos,1,1);//0x1
              break;
            case 1:
              //if(DEBUG)printf("Index:%d ECQ:%d Written:0x00\n",i,ECQ[i]); //DEBUG
              //writeBits_Fast(outBuf,&bitPos,2,0);//0x00
              writeBits_Fast(outBuf,&bitPos,1,0);
              writeBits_Fast(outBuf,&bitPos,1,0);
              break;
            case -1:
              //if(DEBUG)printf("Index:%d ECQ:%d Written:0x01\n",i,ECQ[i]); //DEBUG
              //writeBits_Fast(outBuf,&bitPos,2,2); //0x01
              writeBits_Fast(outBuf,&bitPos,1,0);
              writeBits_Fast(outBuf,&bitPos,1,1);
              break;
            default:
              assert(0);
              break;
          }
        }
        break;
      default: //ECQBits>2
        //if(DEBUG) printf("AMG_W1:bitPos:%ld\n",bitPos); //DEBUG
        for(i=0;i<p->bSize;i++){
          //if(DEBUG){printf("AMG_W3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&outBuf[bitPos/8]));}; //DEBUG
          //if(DEBUG) printf("AMG_W2:bitPos:%ld\n",bitPos); //DEBUG
          //if(DEBUG) printf("ECQ[%d]:%ld\n",i,ECQ[i]); //DEBUG
          switch(ECQ[i]){
            case 0:
              //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x1\n",i,ECQ[i]); //DEBUG
              //if(DEBUG){printf("AMG_WB3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&outBuf[bitPos/8]));}; //DEBUG
              //temp1=bitPos;
              writeBits_Fast(outBuf,&bitPos,1,1);  //0x1
              //wVal=1; writeBits_Fast(outBuf,&bitPos,1,wVal); //0x1
              //if(DEBUG){printf("AMG_WA3:bitPos:%ld buffer[%ld]=0x%lx\n",temp1,temp1/8,*(uint64_t*)(&outBuf[temp1/8]));}; //DEBUG
              break;
            case 1:
              //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x000\n",i,ECQ[i]); //DEBUG
              //if(DEBUG){printf("AMG_WB3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&outBuf[bitPos/8]));}; //DEBUG
              //temp1=bitPos;
              //writeBits_Fast(outBuf,&bitPos,3,0); //0x000
              writeBits_Fast(outBuf,&bitPos,1,0);
              writeBits_Fast(outBuf,&bitPos,1,0);
              writeBits_Fast(outBuf,&bitPos,1,0);
              //wVal=0; writeBits_Fast(outBuf,&bitPos,3,wVal); //0x000
              //if(DEBUG){printf("AMG_WA3:bitPos:%ld buffer[%ld]=0x%lx\n",temp1,temp1/8,*(uint64_t*)(&outBuf[temp1/8]));}; //DEBUG
              break;
            case -1:
              //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x001\n",i,ECQ[i]); //DEBUG
              //if(DEBUG){printf("AMG_WB3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&outBuf[bitPos/8]));}; //DEBUG
              //temp1=bitPos;
              //writeBits_Fast(outBuf,&bitPos,3,8); //0x001
              writeBits_Fast(outBuf,&bitPos,1,0); 
              writeBits_Fast(outBuf,&bitPos,1,0); 
              writeBits_Fast(outBuf,&bitPos,1,1); 
              //wVal=8; writeBits_Fast(outBuf,&bitPos,3,wVal); //0x001
              //if(DEBUG){printf("AMG_WA3:bitPos:%ld buffer[%ld]=0x%lx\n",temp1,temp1/8,*(uint64_t*)(&outBuf[temp1/8]));}; //DEBUG
              break;
            default:
              //if(DEBUG)printf("Index:%d ECQ:%ld Written:0x01 0x%lx\n",i,ECQ[i]); //DEBUG
              //if(DEBUG){printf("AMG_WB3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&outBuf[bitPos/8]));}; //DEBUG
              //temp1=bitPos;
              //writeBits_Fast(outBuf,&bitPos,2,2); //0x01
              writeBits_Fast(outBuf,&bitPos,1,0); 
              writeBits_Fast(outBuf,&bitPos,1,1); 
              //wVal=2; writeBits_Fast(outBuf,&bitPos,2,wVal); //0x01
              writeBits_Fast(outBuf,&bitPos,ECQBits,ECQ[i]);
              //if(DEBUG){printf("AMG_WA3:bitPos:%ld buffer[%ld]=0x%lx\n",temp1,temp1/8,*(uint64_t*)(&outBuf[temp1/8]));}; //DEBUG
              break;
          }
        }
        break;
    }
    
    //if(DEBUG){printf("bitPos_E:%ld\n",bitPos);} //DEBUG
    if(D_C){if(!((ECQBits>=2)||((ECQBits==1) && (numOutliers==0)))){printf("ERROR: ECQBits:%d numOutliers:%d This should not have happened!\n",ECQBits,numOutliers);assert(0);}} //DEBUG
    
          

    uint32_t bytePos=(bitPos+7)/8;
    *(uint32_t*)(&outBuf[9])=bytePos;
    
    if(D_G)printf("bitPos:%ld CSparseBits:%d bytePos:%d CSparseBytes:%d\n",bitPos,CSparseBits,bytePos,CSparseBytes); //DEBUG

    
  }
  //for(i=213;i<233;i++)if(DEBUG)printf("AMG_WE:bitPos:%d buffer[%d]=0x%lx\n",i*8,i,*(uint64_t*)(&outBuf[i])); //DEBUG
  
  return 0;
}

static inline int pastriDecompress(uint8_t *inBuf,parameters *p,uint8_t *outBuf,int *numOutBytes){
  int j;
  //for(j=213;j<233;j++)if(DEBUG){printf("AMG_RB:bitPos:%d buffer[%d]=0x%lx\n",j*8,j,*(uint64_t*)(&inBuf[j]));}; //DEBUG
  
  
  //Preprocess by calculating some parameters:
  //Calculate sbSize, sbNum, etc.:
  p->idxRange[0]=(p->bf[0]+1)*(p->bf[0]+2)/2;
  p->idxRange[1]=(p->bf[1]+1)*(p->bf[1]+2)/2;
  p->idxRange[2]=(p->bf[2]+1)*(p->bf[2]+2)/2;
  p->idxRange[3]=(p->bf[3]+1)*(p->bf[3]+2)/2;
  p->sbSize=p->idxRange[2]*p->idxRange[3];
  p->sbNum=p->idxRange[0]*p->idxRange[1];
  p->bSize=p->sbSize*p->sbNum;
  
  int _1DIdxBits=bitsNeededUI64(p->bSize);
  
  double *data=(double*)(outBuf+p->bSize*8);
  int i0,i1,i2,i3;
  uint16_t *idx0,*idx1,*idx2,*idx3;
  int _1DIdx,nonZeros;
  
  int64_t patternQ[MAX_PS_SIZE]; 
  int64_t scalesQ[MAX_PS_SIZE];
  int64_t ECQTemp;
  uint64_t bytePos;
  uint64_t bitPos;
  uint64_t temp,temp2;
  //int sb,localIdx;
  double binSize,scalesBinSize;
  int patternBits,ECQBits,numOutliers;
  
  idx0=(uint16_t*)(outBuf           );
  idx1=(uint16_t*)(outBuf+p->bSize*2);
  idx2=(uint16_t*)(outBuf+p->bSize*4);
  idx3=(uint16_t*)(outBuf+p->bSize*6);
  p->idxOffset[0]=*(uint32_t*)(&inBuf[1]);
  p->idxOffset[1]=*(uint32_t*)(&inBuf[3]);
  p->idxOffset[2]=*(uint32_t*)(&inBuf[5]);
  p->idxOffset[3]=*(uint32_t*)(&inBuf[7]);
  
  for(i0=0;i0<p->idxRange[0];i0++)
    for(i1=0;i1<p->idxRange[1];i1++)
      for(i2=0;i2<p->idxRange[2];i2++)
        for(i3=0;i3<p->idxRange[3];i3++){
            //_1DIdx=i0*p->idxRange[1]*p->idxRange[2]*p->idxRange[3]+i1*p->idxRange[2]*p->idxRange[3]+i2*p->idxRange[3]+i3;
            _1DIdx=p->idxRange[3]*(i2+p->idxRange[2]*(i1+i0*p->idxRange[1]))+i3;
            idx0[_1DIdx]=i0+1+p->idxOffset[0];
            idx1[_1DIdx]=i1+1+p->idxOffset[1];
            idx2[_1DIdx]=i2+1+p->idxOffset[2];
            idx3[_1DIdx]=i3+1+p->idxOffset[3];
        }
  *numOutBytes=p->bSize*16;  
  
  //inBuf[0] is "mode"
  switch(inBuf[0]){
    //R:UCSparse
    case 0:
      if(D_G){printf("\nDC:UCSparse\n");} //DEBUG
      nonZeros=*(uint16_t*)(&inBuf[9]);
      bytePos=11;
      for(j=0;j<p->bSize;j++){
          data[j]=0;
      }
      for(j=0;j<nonZeros;j++){
        i0=*(uint16_t*)(&inBuf[bytePos])-1-p->idxOffset[0]; //i0
        bytePos+=2;
        i1=*(uint16_t*)(&inBuf[bytePos])-1-p->idxOffset[1]; //i1
        bytePos+=2;
        i2=*(uint16_t*)(&inBuf[bytePos])-1-p->idxOffset[2]; //i2
        bytePos+=2;
        i3=*(uint16_t*)(&inBuf[bytePos])-1-p->idxOffset[3]; //i3
        bytePos+=2;
        _1DIdx=p->idxRange[3]*(i2+p->idxRange[2]*(i1+i0*p->idxRange[1]))+i3;
        data[_1DIdx]=*(double*)(&inBuf[bytePos]);
        bytePos+=8; 
      }
      break;
    //R:UCNonSparse
    case 1:
      if(D_G){printf("\nDC:UCNonSparse\n");} //DEBUG
      memcpy(&outBuf[p->bSize*8], &inBuf[9], p->bSize*8);
      break;
    //R:CSparse
    case 2:
      if(D_G){printf("\nDC:CSparse\n");} //DEBUG
      //for(j=0;j<p->bSize;j++){
      //  data[j]=0;
      //}
      
      patternBits=inBuf[13];
      ECQBits=inBuf[14];
      
      if(D_R){printf("patternBits:%d ECQBits:%d _1DIdxBits:%d\n",patternBits,ECQBits,_1DIdxBits);} //DEBUG
      
      numOutliers=*(uint16_t*)(&inBuf[15]);
      if(D_R){printf("numOutliers:%d\n",numOutliers);} //DEBUG
      bitPos=17*8;

      scalesBinSize=1/(double)(((uint64_t)1<<(patternBits-1))-1);
  
      binSize=p->eb*2;
      
      if(D_R){printf("scalesBinSize:%.6e binSize:%.6e scalesBinSize*binSize:%.6e\n",scalesBinSize,binSize,scalesBinSize*binSize);} //DEBUG

      for(j=0;j<p->sbSize;j++){
        patternQ[j]=readBits_I64(inBuf,&bitPos,patternBits);//Pattern point
        if(D_R){printf("R:patternQ[%d]=%ld\n",j,patternQ[j]);}
      }
      for(j=0;j<p->sbSize;j++){
        scalesQ[j]=readBits_I64(inBuf,&bitPos,patternBits);//Scale
        if(D_R){printf("R:scalesQ[%d]=%ld\n",j,scalesQ[j]);}
      }
      
      for(j=0;j<p->bSize;j++){
        data[j]=scalesQ[j/p->sbSize]*patternQ[j%p->sbSize]*scalesBinSize*binSize;
      }
      
      
      switch(ECQBits){
        case 2:
          for(j=0;j<numOutliers;j++){
            //if(DEBUG){printf("readBits_UI64:%ld\n",readBits_UI64(inBuf,&bitPos,_1DIdxBits));} //DEBUG
            //if(DEBUG){printf("readBits_UI64:%ld\n",readBits_I64(inBuf,&bitPos,2));} //DEBUG
            
            _1DIdx=readBits_UI64(inBuf,&bitPos,_1DIdxBits);
            ECQTemp=readBits_I64(inBuf,&bitPos,1);
            ECQTemp= ((ECQTemp<<63)>>63)|(uint64_t)0x1;
            //if(D_R)printf("R:ECQ[%d]: %ld \n",_1DIdx,ECQTemp);
            //continue;
            //sb=_1DIdx/p->sbSize; 
            //localIdx=_1DIdx%p->sbSize;
            data[_1DIdx]-=ECQTemp*binSize;
            //if(DEBUG){printf("decompressed[%d]:%.6e\n",_1DIdx,data[_1DIdx]);} //DEBUG
          }
          break;
        default: //ECQBits>2
          if(D_C){if(!((ECQBits>=2)||((ECQBits==1) && (numOutliers==0)))){printf("ERROR: ECQBits:%d numOutliers:%d This should not have happened!\n",ECQBits,numOutliers);assert(0);}} //DEBUG
    
          for(j=0;j<numOutliers;j++){
            _1DIdx=readBits_UI64(inBuf,&bitPos,_1DIdxBits);
            //sb=_1DIdx/p->sbSize; 
            //localIdx=_1DIdx%p->sbSize;
            temp=readBits_UI64(inBuf,&bitPos,1);
            //if(DEBUG){printf("temp:%ld\n",temp);} //DEBUG
            switch(temp){
              case 0:  //+-1
                ECQTemp=readBits_I64(inBuf,&bitPos,1);
                ECQTemp= ((ECQTemp<<63)>>63)|(uint64_t)0x1;
                //if(DEBUG){printf("_1DIdx:%ld ECQTemp:0x%ld\n",_1DIdx,ECQTemp);} //DEBUG
                //if(D_R)printf("R:ECQ[%d]: %ld \n",_1DIdx,ECQTemp);
                break;
              case 1: //Others
                ECQTemp=readBits_I64(inBuf,&bitPos,ECQBits);
                //if(DEBUG){printf("_1DIdx:%ld ECQTemp:0x%ld\n",_1DIdx,ECQTemp);} //DEBUG
                //if(D_R)printf("R:ECQ[%d]: %ld \n",_1DIdx,ECQTemp);
                break;
              //default:
              //  printf("ERROR: Bad 2-bit value: 0x%lx",temp);
              // assert(0); //AMG
              //  break;
            }
            data[_1DIdx]-=ECQTemp*binSize;
            //if(DEBUG){printf("decompressed[%d]:%.6e\n",_1DIdx,data[_1DIdx]);} //DEBUG
          }
          break;
      }
      //static inline uint64_t readBits_UI64(uint8_t* buffer,uint64_t *bitPosPtr,uint64_t numBits){ // numBits must be in range [0:56]
      //patternQ=(int64_t*)(inBuf+15); 
      //scalesQ=(int64_t*)(inBuf+15+p->sbSize*8);
      break;
    //R:CNonSparse
    case 3:
      if(D_G){printf("\nDC:CNonSparse\n");} //DEBUG
      
      //for(j=0;j<p->bSize;j++){
      //  data[j]=0;
      //}
      patternBits=inBuf[13];
      ECQBits=inBuf[14];
      
      if(D_R){printf("patternBits:%d ECQBits:%d _1DIdxBits:%d\n",patternBits,ECQBits,_1DIdxBits);} //DEBUG
      
      bitPos=15*8;

      scalesBinSize=1/(double)(((uint64_t)1<<(patternBits-1))-1);
      binSize=p->eb*2;
      
      if(D_R){printf("scalesBinSize:%.6e binSize:%.6e scalesBinSize*binSize:%.6e\n",scalesBinSize,binSize,scalesBinSize*binSize);} //DEBUG

      for(j=0;j<p->sbSize;j++){
        patternQ[j]=readBits_I64(inBuf,&bitPos,patternBits);//Pattern point
        if(D_R){printf("R:patternQ[%d]=%ld\n",j,patternQ[j]);}
      }
      for(j=0;j<p->sbSize;j++){
        scalesQ[j]=readBits_I64(inBuf,&bitPos,patternBits);//Scale
        if(D_R){printf("R:scalesQ[%d]=%ld\n",j,scalesQ[j]);}
      }
      for(j=0;j<p->bSize;j++){
        data[j]=scalesQ[j/p->sbSize]*patternQ[j%p->sbSize]*scalesBinSize*binSize;
        //if(DEBUG){printf("DC:PS[%d]=%.6e\n",j,data[j]);}
      }

      switch(ECQBits){
        case 2:
          for(j=0;j<p->bSize;j++){
            //if(DEBUG){printf("readBits_UI64:%ld\n",readBits_UI64(inBuf,&bitPos,_1DIdxBits));} //DEBUG
            //if(DEBUG){printf("readBits_UI64:%ld\n",readBits_I64(inBuf,&bitPos,2));} //DEBUG
            //_1DIdx=readBits_UI64(inBuf,&bitPos,_1DIdxBits);
            temp=readBits_UI64(inBuf,&bitPos,1);
            switch(temp){
              case 0:
                ECQTemp=readBits_I64(inBuf,&bitPos,1);
                ECQTemp= ((ECQTemp<<63)>>63)|(uint64_t)0x1;
                break;
              case 1:
                ECQTemp=0;
                break;
              default:
                assert(0);
                break;
            }
            
            //if(DEBUG){printf("_1DIdx:%ld ECQTemp:0x%ld\n",_1DIdx,ECQTemp);} //DEBUG
            //continue;
            //sb=_1DIdx/p->sbSize; 
            //localIdx=_1DIdx%p->sbSize;
            data[j]-=ECQTemp*binSize;
            //if(DEBUG){printf("decompressed[%d]:%.6e\n",_1DIdx,data[_1DIdx]);} //DEBUG
          }
          break;
        default: //ECQBits>2
          //if(DEBUG)printf("AMG_R1:bitPos: %ld\n",bitPos);
          
          for(j=0;j<p->bSize;j++){
            //if(DEBUG){printf("AMG_R3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&inBuf[bitPos/8]));}; //DEBUG
            //if(DEBUG)printf("AMG_R2:bitPos: %ld\n",bitPos);

            //if(DEBUG){printf("readBits_UI64:%ld\n",readBits_UI64(inBuf,&bitPos,_1DIdxBits));} //DEBUG
            //if(DEBUG){printf("readBits_UI64:%ld\n",readBits_I64(inBuf,&bitPos,2));} //DEBUG
            //_1DIdx=readBits_UI64(inBuf,&bitPos,_1DIdxBits);
            temp=readBits_UI64(inBuf,&bitPos,1);
            //if(DEBUG){printf("AMG_R3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&inBuf[bitPos/8]));}; //DEBUG
            switch(temp){
              case 0:
                //if(DEBUG)printf("Read:0");
                temp2=readBits_UI64(inBuf,&bitPos,1);
                switch(temp2){
                  case 0:
                    //if(DEBUG)printf("0");
                    ECQTemp=readBits_I64(inBuf,&bitPos,1);
                    //if(DEBUG){printf("AMG_R3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&inBuf[bitPos/8]));}; //DEBUG
                    //if(DEBUG)printf("R:ECQTemp:%ld\n",ECQTemp);
                    ECQTemp= ((ECQTemp<<63)>>63)|(uint64_t)0x1;
                    //if(DEBUG)printf("R:ECQ[%d]: %ld\n",j,ECQTemp);
                    break;
                  case 1:
                    //if(DEBUG)printf("1\n");
                    ECQTemp=readBits_I64(inBuf,&bitPos,ECQBits);
                    //if(DEBUG){printf("AMG_R3:bitPos:%ld buffer[%ld]=0x%lx\n",bitPos,bitPos/8,*(uint64_t*)(&inBuf[bitPos/8]));}; //DEBUG
                    //if(DEBUG)printf("R:ECQ[%d]: %ld\n",j,ECQTemp);
                    break;
                  default:
                    assert(0);
                    break;
                }
                break;
              case 1:
                //if(DEBUG)printf("Read:1\n");
                ECQTemp=0;
                //if(DEBUG)printf("R:ECQ[%d]: %ld\n",j,ECQTemp);
                break;
              default:
                assert(0);
                break;
            }
            
            //if(DEBUG){printf("_1DIdx:%ld ECQTemp:0x%ld\n",_1DIdx,ECQTemp);} //DEBUG
            //continue;
            //sb=_1DIdx/p->sbSize; 
            //localIdx=_1DIdx%p->sbSize;
            data[j]-=ECQTemp*binSize;
            //if(DEBUG){printf("DC:data[%d]:%.6e\n",j,data[j]);} //DEBUG
          }
          break;
      }
      //static inline uint64_t readBits_UI64(uint8_t* buffer,uint64_t *bitPosPtr,uint64_t numBits){ // numBits must be in range [0:56]
      //patternQ=(int64_t*)(inBuf+15); 
      //scalesQ=(int64_t*)(inBuf+15+p->sbSize*8);
      break;
      
    default:
      assert(0);
      break;
  } 

  return 0;
}

//inBuf vs Decompressed
static inline int pastriCheck(uint8_t *inBuf,uint8_t *DC,parameters *p){
  int i;
  
  double *data=(double*)(inBuf+p->bSize*8);
  uint16_t *idx0=(uint16_t*)(inBuf           );
  uint16_t *idx1=(uint16_t*)(inBuf+p->bSize*2);
  uint16_t *idx2=(uint16_t*)(inBuf+p->bSize*4);
  uint16_t *idx3=(uint16_t*)(inBuf+p->bSize*6);

  double *data_dc=(double*)(DC+p->bSize*8);
  uint16_t *idx0_dc=(uint16_t*)(DC           );
  uint16_t *idx1_dc=(uint16_t*)(DC+p->bSize*2);
  uint16_t *idx2_dc=(uint16_t*)(DC+p->bSize*4);
  uint16_t *idx3_dc=(uint16_t*)(DC+p->bSize*6);
  
  //Comparing Indexes:
  for(i=0;i<p->bSize;i++){
    if(idx0[i]!=idx0_dc[i]){
      printf("idx0[%d]=%d  !=  %d=idx0_dc[%d]",i,idx0[i],idx0_dc[i],i);
      assert(0);
    }
    if(idx1[i]!=idx1_dc[i]){
      printf("idx1[%d]=%d  !=  %d=idx1_dc[%d]",i,idx1[i],idx1_dc[i],i);
      assert(0);
    }
    if(idx2[i]!=idx2_dc[i]){
      printf("idx2[%d]=%d  !=  %d=idx2_dc[%d]",i,idx2[i],idx2_dc[i],i);
      assert(0);
    }
    if(idx3[i]!=idx3_dc[i]){
      printf("idx3[%d]=%d  !=  %d=idx3_dc[%d]",i,idx3[i],idx3_dc[i],i);
      assert(0);
    }
  }
  
  //Comparing Data:
  for(i=0;i<p->bSize;i++){
    if(abs_FastD(data[i]-data_dc[i])>p->eb){
      printf("|data[%d]-data_dc[%d]|>eb : %.3e - %.3e = %.3e > %.3e",i,i,data[i],data_dc[i],abs_FastD(data[i]-data_dc[i]),p->eb);
      assert(0);
    }

  }
  return 0;
}

#endif





