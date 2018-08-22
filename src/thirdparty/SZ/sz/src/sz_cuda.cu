#include "cuda.h"
#include "cuda_runtime.h"
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include "sz_cuda.h"
#include <conf.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <printf.h>
#include <sqltypes.h>
#include <sz.h>
#include <sz_float.h>
#include <time.h>

// region CUDA macros
#ifdef __CUDACC__
#define GPU_DEVICE __device__
#define GPU_KERNEL __global__
#define GPU_HOST __host__
#define TIMER cudaEvent_t
//TODO (robertu#1) when using full cuda instead of thrust, this version of parallel launch is required {{{
//#define PARALLEL_LAUNCH(num_blocks,num_threads,function_name) function_name<<<(num_blocks), (num_threads)>>>
#define PARALLEL_LAUNCH(num_blocks,num_threads,function_name) function_name
//END TODO (robertu#1) }}}
inline void
start_timer(TIMER* start)
{
  cudaEventCreate(start);
  cudaEventRecord(*start, 0);
}
/**
 * @param start  -- timer to stop and deallocate
 * @return the time elapsed in ms
 */
inline double
stop_timer(TIMER* start)
{
  cudaEvent_t stop;
  cudaEventCreate(&stop);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  float elapsed;
  cudaEventElapsedTime(&elapsed, *start, stop);
  cudaEventDestroy(stop);
  cudaEventDestroy(*start);
  return elapsed;
}

#define NO_DEVICE -1

//TODO (robertu#2) use the openmp version until more cuda kernels are availible {{{
inline int
get_max_threads()
{
	return omp_get_max_threads();
}

inline void
set_max_threads(size_t threads)
{
  omp_set_num_threads(threads);
}

inline void
get_thread_id(int* idx)
{
	*idx = omp_get_thread_num();
}
/* code for cuda version when kernels are ready
inline int
get_max_threads()
{
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  for (int i = 0; i < num_devices; ++i) {
    struct cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, i);
    return properties.maxThreadsPerMultiProcessor;
  }
  return NO_DEVICE;
}
inline void
set_max_threads(size_t threads)
{
  (void)0; // NOOP on cuda silence compiler
}

void
get_thread_id(int* idx)
{
  *idx = blockIdx.x * blockDim.x + threadIdx.x;
}
*/
//END TODO (robertu#2) }}}

#else /*omp version*/
#define GPU_DEVICE
#define GPU_KERNEL
#define GPU_HOST
#define PARALLEL_LAUNCH(num_blocks,num_threads,function_name) function_name

#include <chrono>
#define TIMER std::chrono::high_resolution_clock::time_point

inline double stop_timer(TIMER* start)
{
	std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> duration = now - *start;
	return duration.count();
}

inline void start_timer(TIMER* start)
{
	*start = std::chrono::high_resolution_clock::now();
}

inline int
get_max_threads()
{
  return omp_get_max_threads();
}

inline void
set_max_threads(size_t threads)
{
  omp_set_num_threads(threads);
}

inline void
get_thread_id(int* idx)
{
  *idx = omp_get_thread_num();
}
#endif


// region Blockcount

struct BlockCount
{
  size_t split_index;
  size_t early_blockcount;
  size_t late_blockcount;
};
// endregion
// region MemoryBlocks
struct CompressionMemoryBlocks
{
  int* result_type;
  float* result_unpredictable_data;
  unsigned int* unpredictable_count;
  float* mean;
};

inline void
deleteCompressionMemoryBlocks(struct CompressionMemoryBlocks* memory)
{
  free((*memory).mean);
  free((*memory).result_unpredictable_data);
  free((*memory).unpredictable_count);
  free((*memory).result_type);
}

// endregion
// endregion
// region Prototypes


inline size_t get_offset(struct BlockCount const* const x, int const i);

inline struct CompressionMemoryBlocks
newCompressionMemoryBlocks(size_t num_blocks, size_t num_elements,
                           size_t unpred_data_max_size)
{
  struct CompressionMemoryBlocks memory;
  memory.result_type = (int*)malloc(num_elements * sizeof(int));
  memory.result_unpredictable_data =
    (float*)malloc(unpred_data_max_size * sizeof(float) * num_blocks);
  memory.unpredictable_count =
    (unsigned int*)malloc(num_blocks * sizeof(unsigned int));
  memory.mean = (float*)malloc(num_blocks * sizeof(float));
  return memory;
}

// endregion


unsigned char*
SZ_compress_float_1D_MDQ_cuda(float* oriData, size_t r1, double realPrecision,
                              size_t* comp_size)
{


  // determine quantization_intervals
  unsigned int quantization_intervals =
    quantization_intervals_1D(oriData, r1, realPrecision);

  // configure threads
  int thread_num = get_max_threads();
  size_t num_x = thread_num;
  set_max_threads(thread_num);

  // compute blockcounts
  struct BlockCount x;
  SZ_COMPUTE_BLOCKCOUNT(r1, num_x, x.split_index, x.early_blockcount,
                        x.late_blockcount);

  // compute work sizes
  size_t max_num_block_elements = x.early_blockcount;
  size_t num_blocks = num_x;
  size_t num_elements = r1;
  size_t unpred_data_max_size = max_num_block_elements;

  // allocate memory for compression
  struct CompressionMemoryBlocks memory =
    newCompressionMemoryBlocks(num_blocks, num_elements, unpred_data_max_size);

  // compress memory
  _sz_compress_float_1d_mdq_ra_block(oriData, r1, realPrecision, thread_num,
                                     unpred_data_max_size, &x, &memory);

  // build huffman encoding
  size_t nodeCount;
  unsigned char* treeBytes;
  unsigned int treeByteSize;
  buildHuffmanTree(thread_num, num_elements, &memory, &nodeCount, &treeBytes,
                   &treeByteSize);

  // compute number of unpredictable blocks
  size_t total_unpred = compute_total_unpred_gpu(num_blocks, &memory);

  // allocate output buffer
  unsigned char* result_pos;
  unsigned char* result = result_pos =
    (unsigned char*)malloc(compute_compressed_size(num_blocks, num_elements,
                                                   treeByteSize, total_unpred));

  // write compression metadata
  result_pos += initRandomAccessBytes(result_pos);
  result_pos = write_parallel_compresion_metadata(
    result_pos, thread_num, realPrecision, quantization_intervals, &memory,
    num_blocks, nodeCount, treeBytes, treeByteSize);
  free(treeBytes);

  // write unpredictable data
  size_t* unpred_offset =
    compute_unpred_offset(thread_num, num_blocks, &memory);
  copy_unpredictable(thread_num, unpred_data_max_size, &memory, result_pos,
                     unpred_offset);
  result_pos += total_unpred * sizeof(float);
  free(unpred_offset);

  // encode remaining data
  size_t* block_pos =
    (size_t*)result_pos; // block_pos exists to affect pointer math
  result_pos += num_blocks * sizeof(size_t);
  unsigned char* encoding_buffer =
    (unsigned char*)malloc(max_num_block_elements * sizeof(int) * num_blocks);
  encode_1D(oriData, thread_num, &x, max_num_block_elements, &memory, block_pos,
            encoding_buffer);

  // write encoded data
  size_t* block_offset =
    compute_block_offsets(thread_num, num_blocks, block_pos);
  copyEncodingBuffers(thread_num, max_num_block_elements, result_pos, block_pos,
                      encoding_buffer, block_offset);
  result_pos += block_offset[thread_num - 1] + block_pos[thread_num - 1];
  free(block_offset);

  // cleanup
  free(encoding_buffer);
  deleteCompressionMemoryBlocks(&memory);
  SZ_ReleaseHuffman();

  // return pointer to compressed data and size
  *comp_size = result_pos - result;
  return result;
}
unsigned char*
SZ_compress_float_2D_MDQ_cuda(float* oriData, size_t r1, size_t r2,
                              double realPrecision, size_t* comp_size)
{
  // determine quantization_intervals
  unsigned int quantization_intervals =
    quantization_intervals_2D(oriData, r1, r2, realPrecision);

  // configure threads
  int thread_num;
  size_t num_x;
  size_t num_y;
  config_threads_2D(&thread_num, &num_x, &num_y);

  // compute blockcounts
  struct BlockCount x, y;
  SZ_COMPUTE_BLOCKCOUNT(r1, num_x, x.split_index, x.early_blockcount,
                        x.late_blockcount);
  SZ_COMPUTE_BLOCKCOUNT(r2, num_y, y.split_index, y.early_blockcount,
                        y.late_blockcount);

  // compute work sizes
  size_t max_num_block_elements = x.early_blockcount * y.early_blockcount;
  size_t num_blocks = num_x * num_y;
  size_t num_elements = r1 * r2;
  size_t unpred_data_max_size = max_num_block_elements;
  size_t dim0_offset = r2;
  size_t buffer_size = y.early_blockcount * sizeof(float);

  // allocate memory for compression
  struct CompressionMemoryBlocks memory =
    newCompressionMemoryBlocks(num_blocks, num_elements, unpred_data_max_size);

  // compress memory
  PARALLEL_LAUNCH(1,thread_num,_sz_compress_float_2d_mdq_ra_block)(oriData, r1, r2, realPrecision, thread_num,
                                     num_y, unpred_data_max_size, dim0_offset,
                                     buffer_size, &x, &y, &memory);

  // build huffman encoding
  // TODO keep this one until I figure out why they pass arguments now
  // SZ_Reset(allNodes, stateNum);
  size_t nodeCount;
  unsigned char* treeBytes;
  unsigned int treeByteSize;
  buildHuffmanTree(thread_num, num_elements, &memory, &nodeCount, &treeBytes,
                   &treeByteSize);

  // compute number of unpredictable blocks
  size_t total_unpred = compute_total_unpred_gpu(num_blocks, &memory);

  // allocate output buffer
  unsigned char* result_pos;
  unsigned char* result = result_pos =
    (unsigned char*)malloc(compute_compressed_size(num_blocks, num_elements,
                                                   treeByteSize, total_unpred));

  // write compresion metadata
  result_pos += initRandomAccessBytes(result_pos);
  result_pos = write_parallel_compresion_metadata(
    result_pos, thread_num, realPrecision, quantization_intervals, &memory,
    num_blocks, nodeCount, treeBytes, treeByteSize);
  free(treeBytes);

  // write unpredictable data
  size_t* unpred_offset =
    compute_unpred_offset(thread_num, num_blocks, &memory);
  copy_unpredictable(thread_num, unpred_data_max_size, &memory, result_pos,
                     unpred_offset);
  result_pos += total_unpred * sizeof(float);
  free(unpred_offset);

  // encode remaining data
  size_t* block_pos = (size_t*)result_pos;
  result_pos += num_blocks * sizeof(size_t);
  unsigned char* encoding_buffer =
    (unsigned char*)malloc(max_num_block_elements * sizeof(int) * num_blocks);
  PARALLEL_LAUNCH(1,thread_num,encode_2D)(thread_num, num_y, &x, &y, max_num_block_elements, dim0_offset,
            &memory, block_pos, encoding_buffer);

  // write encoded data
  size_t* block_offset =
    compute_block_offsets(thread_num, num_blocks, block_pos);
  copyEncodingBuffers(thread_num, max_num_block_elements, result_pos, block_pos,
                      encoding_buffer, block_offset);
  result_pos += block_offset[thread_num - 1] + block_pos[thread_num - 1];
  free(block_offset);

  // cleanup
  free(encoding_buffer);
  deleteCompressionMemoryBlocks(&memory);
  SZ_ReleaseHuffman();

  *comp_size = result_pos - result;
  return result;
}
unsigned char*
SZ_compress_float_3D_MDQ_cuda(float* oriData, size_t r1, size_t r2, size_t r3,
                              double realPrecision, size_t* comp_size)
{


  // determine quantization_intervals
  unsigned int quantization_intervals =
    quantization_intervals_3D(oriData, r1, r2, r3, realPrecision);

  // configure threads
  int thread_num;
  size_t num_x, num_y, num_z;
  config_threads_3D(&thread_num, &num_x, &num_y, &num_z);

  // compute blockcounts
  struct BlockCount x, y, z;
  SZ_COMPUTE_BLOCKCOUNT(r1, num_x, x.split_index, x.early_blockcount,
                        x.late_blockcount);
  SZ_COMPUTE_BLOCKCOUNT(r2, num_y, y.split_index, y.early_blockcount,
                        y.late_blockcount);
  SZ_COMPUTE_BLOCKCOUNT(r3, num_z, z.split_index, z.early_blockcount,
                        z.late_blockcount);

  // compute work sizes
  size_t max_num_block_elements =
    x.early_blockcount * y.early_blockcount * z.early_blockcount;
  size_t num_blocks = num_x * num_y * num_z;
  size_t num_elements = r1 * r2 * r3;
  size_t unpred_data_max_size = max_num_block_elements;
  size_t dim0_offset = r2 * r3;
  size_t dim1_offset = r3;
  int num_yz = num_y * num_z;
  size_t buffer_size = y.early_blockcount * z.early_blockcount * sizeof(float);

  // allocate memory for compression
  struct CompressionMemoryBlocks memory =
    newCompressionMemoryBlocks(num_blocks, num_elements, unpred_data_max_size);

  // compress memory
  PARALLEL_LAUNCH(1,thread_num,_sz_compress_float_3d_mdq_ra_block)(
    oriData, r1, r2, r3, realPrecision, thread_num, num_z, unpred_data_max_size,
    dim0_offset, dim1_offset, num_yz, buffer_size, &x, &y, &z, &memory);

  // build huffman encoding
  size_t nodeCount;
  unsigned char* treeBytes;
  unsigned int treeByteSize;
  buildHuffmanTree(thread_num, num_elements, &memory, &nodeCount, &treeBytes,
                   &treeByteSize);

  // compute number of unpredictable blocks
  size_t total_unpred = compute_total_unpred_gpu(num_blocks, &memory);

  // allocate output buffer
  unsigned char* result_pos;
  unsigned char* result = result_pos =
    (unsigned char*)malloc(compute_compressed_size(num_blocks, num_elements,
                                                   treeByteSize, total_unpred));

  // write compression metadata
  result_pos += initRandomAccessBytes(result_pos);
  result_pos = write_parallel_compresion_metadata(
    result_pos, thread_num, realPrecision, quantization_intervals, &memory,
    num_blocks, nodeCount, treeBytes, treeByteSize);
  free(treeBytes);

  // write unpredictable data
  size_t* unpred_offset =
    compute_unpred_offset(thread_num, num_blocks, &memory);
  copy_unpredictable(thread_num, unpred_data_max_size, &memory, result_pos,
                     unpred_offset);
  result_pos += total_unpred * sizeof(float);
  free(unpred_offset);

  // encode remaining data
  size_t* block_pos = (size_t*)result_pos;
  result_pos += num_blocks * sizeof(size_t);
  unsigned char* encoding_buffer =
    (unsigned char*)malloc(max_num_block_elements * sizeof(int) * num_blocks);
  PARALLEL_LAUNCH(1,thread_num,encode_3D)(thread_num, num_z, &x, &y, &z, max_num_block_elements, dim0_offset,
            dim1_offset, num_yz, &memory, block_pos, encoding_buffer);

  // write encoded data
  size_t* block_offset =
    compute_block_offsets(thread_num, num_blocks, block_pos);
  copyEncodingBuffers(thread_num, max_num_block_elements, result_pos, block_pos,
                      encoding_buffer, block_offset);
  result_pos += block_offset[thread_num - 1] + block_pos[thread_num - 1];
  free(block_offset);

  // cleanup
  free(encoding_buffer);
  deleteCompressionMemoryBlocks(&memory);
  SZ_ReleaseHuffman();

  *comp_size = result_pos - result;
  return result;
}


void
decompressDataSeries_float_1D_cuda(float** data, size_t r1,
                                   unsigned char* comp_data)
{

  size_t num_elements = r1;

  *data = (float*)malloc(sizeof(float) * num_elements);

  unsigned char* comp_data_pos = comp_data;
  int thread_num = readIntBigEndian(&comp_data_pos);
  size_t num_x = thread_num;

  set_max_threads(thread_num);
  struct BlockCount x;
  SZ_COMPUTE_BLOCKCOUNT(r1, num_x, x.split_index, x.early_blockcount,
                        x.late_blockcount);

  size_t num_blocks = num_x;

  double realPrecision = bytesToDouble(comp_data_pos);
  comp_data_pos += 8;
  unsigned int intervals = readIntBigEndian(&comp_data_pos);

  updateQuantizationInfo(intervals);
  // intvRadius = (int)((tdps->intervals - 1)/ 2);

  struct CompressionMemoryBlocks memory;
  unsigned int tree_size = readIntBigEndian(&comp_data_pos);
  allNodes = readIntBigEndian(&comp_data_pos);
  stateNum = allNodes / 2;
  SZ_Reset();
  node root =
    reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos, allNodes);
  comp_data_pos += tree_size;

  unsigned int* unpred_count = (unsigned int*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(unsigned int);

  float* mean_pos = (float*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(float);

  memory.result_unpredictable_data = (float*)comp_data_pos;

  size_t total_unpred = 0;
  size_t* unpred_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  for (int i = 0; i < num_blocks; i++) {
    unpred_offset[i] = total_unpred;
    total_unpred += unpred_count[i];
  }

  comp_data_pos += total_unpred * sizeof(float);

  memory.result_type = (int*)malloc(num_elements * sizeof(int));
  // decode(comp_data_pos, num_elements, root, memory.result_type);
  size_t* block_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  size_t* block_pos = (size_t*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(size_t);
  block_offset[0] = 0;
  for (int t = 1; t < thread_num; t++) {
    block_offset[t] = block_pos[t - 1] + block_offset[t - 1];
  }
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t;
    size_t offset_x = get_offset(&x, i);
    size_t current_blockcount_x = get_current_blockcount(&x, i);
    size_t type_offset = offset_x;
    int* type = memory.result_type + type_offset;
    decode(comp_data_pos + block_offset[t], current_blockcount_x, root, type);
  }

#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t;
    size_t offset_x = get_offset(&x, i);
    float* data_pos = *data + offset_x;

    size_t current_blockcount_x = get_current_blockcount(&x, i);
    size_t type_offset = offset_x;
    int* type = memory.result_type + type_offset;

    float* unpredictable_data =
      memory.result_unpredictable_data + unpred_offset[t];
    float mean = mean_pos[t];
    int cur_unpred_count = decompressDataSeries_float_1D_RA_block(
      data_pos, mean, r1, current_blockcount_x, realPrecision, type,
      unpredictable_data);
  }

  free(memory.result_type);
  free(unpred_offset);
}

void
decompressDataSeries_float_2D_cuda(float** data, size_t r1, size_t r2,
                                   unsigned char* comp_data)
{
  // printf("num_block_elements %d num_blocks %d\n", max_num_block_elements,
  // num_blocks); fflush(stdout);
  TIMER timer;
  start_timer(&timer);

  size_t dim0_offset = r2;
  size_t num_elements = r1 * r2;

  *data = (float*)malloc(sizeof(float) * num_elements);

  unsigned char* comp_data_pos = comp_data;

  int thread_num = bytesToInt_bigEndian(comp_data_pos);
  comp_data_pos += 4;
  int thread_order = (int)log2(thread_num);
  size_t num_x, num_y;
  {
    int block_thread_order = thread_order / 2;
    switch (thread_order % 2) {
      case 0: {
        num_x = 1 << block_thread_order;
        num_y = 1 << block_thread_order;
        break;
      }
      case 1: {
        num_x = 1 << (block_thread_order + 1);
        num_y = 1 << block_thread_order;
        break;
      }
    }
  }
  printf("number of blocks: %zu %zu, thread_num %d\n", num_x, num_y,
         thread_num);
  set_max_threads(thread_num);
  struct BlockCount x, y;
  SZ_COMPUTE_BLOCKCOUNT(r1, num_x, x.split_index, x.early_blockcount,
                        x.late_blockcount);
  SZ_COMPUTE_BLOCKCOUNT(r2, num_y, y.split_index, y.early_blockcount,
                        y.late_blockcount);

  size_t num_blocks = num_x * num_y;

  double realPrecision = bytesToDouble(comp_data_pos);
  comp_data_pos += 8;
  unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
  comp_data_pos += 4;

  updateQuantizationInfo(intervals);
  // intvRadius = (int)((tdps->intervals - 1)/ 2);

  unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
  comp_data_pos += 4;
  allNodes = bytesToInt_bigEndian(comp_data_pos);
  stateNum = allNodes / 2;
  SZ_Reset();
  // printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
  // fflush(stdout);
  node root =
    reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos + 4, allNodes);

  struct CompressionMemoryBlocks memory;
  comp_data_pos += 4 + tree_size;
  unsigned int* unpred_count = (unsigned int*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(unsigned int);
  float* mean_pos = (float*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(float);
  memory.result_unpredictable_data = (float*)comp_data_pos;
  size_t total_unpred = 0;
  size_t* unpred_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  for (int i = 0; i < num_blocks; i++) {
    unpred_offset[i] = total_unpred;
    total_unpred += unpred_count[i];
  }
  comp_data_pos += total_unpred * sizeof(float);

  memory.result_type = (int*)malloc(num_elements * sizeof(int));
  // decode(comp_data_pos, num_elements, root, memory.result_type);
  size_t* block_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  size_t* block_pos = (size_t*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(size_t);
  block_offset[0] = 0;
  for (int t = 1; t < thread_num; t++) {
    block_offset[t] = block_pos[t - 1] + block_offset[t - 1];
  }
  printf("Read data info elapsed time: %.4f\n", stop_timer(&timer));
  start_timer(&timer);
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_y);
    int j = (t % num_y);
    size_t offset_x = get_offset(&x, i);
    size_t offset_y = get_offset(&y, j);

    size_t current_blockcount_x = get_current_blockcount(&x, i);
    size_t current_blockcount_y = get_current_blockcount(&y, j);

    size_t type_offset =
      offset_x * dim0_offset + offset_y * current_blockcount_x;
    int* type = memory.result_type + type_offset;
    decode(comp_data_pos + block_offset[t],
           current_blockcount_x * current_blockcount_y, root, type);
  }
  printf("Parallel Huffman decoding elapsed time: %.4f\n", stop_timer(&timer));
  start_timer(&timer);

#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / num_y;
    int j = (t % num_y);
    // printf("%d: %d %d %d\n", omp_get_thread_num(), i, j, k);
    size_t offset_x = get_offset(&x, i);
    size_t offset_y = get_offset(&y, j);
    float* data_pos = *data + offset_x * dim0_offset + offset_y;

    size_t current_blockcount_x = get_current_blockcount(&x, i);
    size_t current_blockcount_y = get_current_blockcount(&y, j);

    size_t type_offset =
      offset_x * dim0_offset + offset_y * current_blockcount_x;
    int* type = memory.result_type + type_offset;

    float* unpredictable_data =
      memory.result_unpredictable_data + unpred_offset[t];
    float mean = mean_pos[t];
    // printf("\n%d\ndata_offset: %ld\n", t, offset_x * dim0_offset + offset_y *
    // dim1_offset + offset_z); printf("memory.mean: %.2f\n", memory.mean);
    // for(int tmp=0; tmp<10; tmp++){
    // 	printf("%.2f ", unpredictable_data[tmp]);
    // }
    // printf("\n\n");
    int cur_unpred_count = decompressDataSeries_float_2D_RA_block(
      data_pos, mean, r1, r2, current_blockcount_x, current_blockcount_y,
      realPrecision, type, unpredictable_data);
  }
  printf("Parallel decompress elapsed time: %.4f\n", stop_timer(&timer));

  free(memory.result_type);
  free(unpred_offset);
}
void
decompressDataSeries_float_3D_cuda(float** data, size_t r1, size_t r2,
                                   size_t r3, unsigned char* comp_data)
{
  // printf("num_block_elements %d num_blocks %d\n", max_num_block_elements,
  // num_blocks); fflush(stdout);
  double elapsed_time = 0.0;
  elapsed_time = -omp_get_wtime();

  size_t dim0_offset = r2 * r3;
  size_t dim1_offset = r3;
  size_t num_elements = r1 * r2 * r3;

  *data = (float*)malloc(sizeof(float) * num_elements);

  unsigned char* comp_data_pos = comp_data;
  // int meta_data_offset = 3 + 1 + MetaDataByteLength;
  // comp_data_pos += meta_data_offset;

  int thread_num = bytesToInt_bigEndian(comp_data_pos);
  comp_data_pos += 4;
  int thread_order = (int)log2(thread_num);
  size_t num_x, num_y, num_z;
  {
    int block_thread_order = thread_order / 3;
    switch (thread_order % 3) {
      case 0: {
        num_x = 1 << block_thread_order;
        num_y = 1 << block_thread_order;
        num_z = 1 << block_thread_order;
        break;
      }
      case 1: {
        num_x = 1 << (block_thread_order + 1);
        num_y = 1 << block_thread_order;
        num_z = 1 << block_thread_order;
        break;
      }
      case 2: {
        num_x = 1 << (block_thread_order + 1);
        num_y = 1 << (block_thread_order + 1);
        num_z = 1 << block_thread_order;
        break;
      }
    }
  }
  printf("number of blocks: %zu %zu %zu, thread_num %d\n", num_x, num_y, num_z,
         thread_num);
  set_max_threads(thread_num);
  struct BlockCount x, y, z;
  SZ_COMPUTE_BLOCKCOUNT(r1, num_x, x.split_index, x.early_blockcount,
                        x.late_blockcount);
  SZ_COMPUTE_BLOCKCOUNT(r2, num_y, y.split_index, y.early_blockcount,
                        y.late_blockcount);
  SZ_COMPUTE_BLOCKCOUNT(r3, num_z, z.split_index, z.early_blockcount,
                        z.late_blockcount);

  size_t num_blocks = num_x * num_y * num_z;

  double realPrecision = bytesToDouble(comp_data_pos);
  comp_data_pos += 8;
  unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
  comp_data_pos += 4;

  updateQuantizationInfo(intervals);
  // intvRadius = (int)((tdps->intervals - 1)/ 2);

  unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
  comp_data_pos += 4;
  allNodes = bytesToInt_bigEndian(comp_data_pos);
  stateNum = allNodes / 2;
  SZ_Reset();
  // printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
  // fflush(stdout);
  node root =
    reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos + 4, allNodes);

  struct CompressionMemoryBlocks memory;
  comp_data_pos += 4 + tree_size;
  unsigned int* unpred_count = (unsigned int*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(unsigned int);
  float* mean_pos = (float*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(float);
  memory.result_unpredictable_data = (float*)comp_data_pos;
  size_t total_unpred = 0;
  size_t* unpred_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  for (int i = 0; i < num_blocks; i++) {
    unpred_offset[i] = total_unpred;
    total_unpred += unpred_count[i];
  }
  comp_data_pos += total_unpred * sizeof(float);

  memory.result_type = (int*)malloc(num_elements * sizeof(int));
  // decode(comp_data_pos, num_elements, root, memory.result_type);
  size_t* block_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  size_t* block_pos = (size_t*)comp_data_pos;
  comp_data_pos += num_blocks * sizeof(size_t);
  block_offset[0] = 0;
  for (int t = 1; t < thread_num; t++) {
    block_offset[t] = block_pos[t - 1] + block_offset[t - 1];
  }
  int num_yz = num_y * num_z;
  elapsed_time += omp_get_wtime();
  printf("Read data info elapsed time: %.4f\n", elapsed_time);
  elapsed_time = -omp_get_wtime();

#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_yz);
    int j = (t % num_yz) / num_z;
    int k = t % num_z;
    size_t offset_x = get_offset(&x, i);
    size_t offset_y = get_offset(&y, j);
    size_t offset_z = get_offset(&z, k);
    size_t current_blockcount_x = get_current_blockcount(&x, i);
    size_t current_blockcount_y = get_current_blockcount(&y, j);
    size_t current_blockcount_z = get_current_blockcount(&z, k);
    size_t type_offset = offset_x * dim0_offset +
                         offset_y * current_blockcount_x * dim1_offset +
                         offset_z * current_blockcount_x * current_blockcount_y;
    int* type = memory.result_type + type_offset;
    decode(comp_data_pos + block_offset[t],
           current_blockcount_x * current_blockcount_y * current_blockcount_z,
           root, type);
  }
  elapsed_time += omp_get_wtime();
  printf("Parallel Huffman decoding elapsed time: %.4f\n", elapsed_time);
  elapsed_time = -omp_get_wtime();

#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_yz);
    int j = (t % num_yz) / num_z;
    int k = t % num_z;
    // printf("%d: %d %d %d\n", omp_get_thread_num(), i, j, k);
    size_t offset_x = get_offset(&x, i);
    size_t offset_y = get_offset(&y, j);
    size_t offset_z = get_offset(&z, k);

    float* data_pos =
      *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

    size_t current_blockcount_x = get_current_blockcount(&x, i);
    size_t current_blockcount_y = get_current_blockcount(&y, j);
    size_t current_blockcount_z = get_current_blockcount(&z, k);

    size_t type_offset = offset_x * dim0_offset +
                         offset_y * current_blockcount_x * dim1_offset +
                         offset_z * current_blockcount_x * current_blockcount_y;
    int* type = memory.result_type + type_offset;

    float* unpredictable_data =
      memory.result_unpredictable_data + unpred_offset[t];
    float mean = mean_pos[t];
    int cur_unpred_count = decompressDataSeries_float_3D_RA_block(
      data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y,
      current_blockcount_z, realPrecision, type, unpredictable_data);
  }
  elapsed_time += omp_get_wtime();
  printf("Parallel decompress elapsed time: %.4f\n", elapsed_time);

  free(memory.result_type);
  free(unpred_offset);
}



inline void
_sz_compress_float_1d_mdq_ra_block(float* oriData, size_t r1,
                                   double realPrecision, int thread_num,
                                   size_t unpred_data_max_size,
                                   struct BlockCount* x,
                                   struct CompressionMemoryBlocks* memory)
{
  #pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {

	  int id;
	  get_thread_id(&id);
    size_t offset_x = get_offset(x, id);
    size_t current_blockcount_x = get_current_blockcount(x, id);

    float* data_pos = oriData + offset_x;

    size_t type_offset = offset_x;
    int* type = (*memory).result_type + type_offset;

    float* unpredictable_data =
      (*memory).result_unpredictable_data + id * unpred_data_max_size;
    (*memory).unpredictable_count[id] = SZ_compress_float_1D_MDQ_RA_block(
      data_pos, (*memory).mean + id, r1, current_blockcount_x, realPrecision,
      type, unpredictable_data);
  }
}

inline void
_sz_compress_float_2d_mdq_ra_block(float* oriData, size_t r1, size_t r2,
                                   double realPrecision, int thread_num,
                                   size_t num_y, size_t unpred_data_max_size,
                                   size_t dim0_offset, size_t buffer_size,
                                   struct BlockCount* x, struct BlockCount* y,
                                   struct CompressionMemoryBlocks* memory)
{
  float *P0, *P1; // buffer
  P0 = (float*)malloc(buffer_size * thread_num);
  P1 = (float*)malloc(buffer_size * thread_num);
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_y);
    int j = (t % num_y);
    size_t offset_x = get_offset(x, i);
    size_t offset_y = get_offset(x, j);
    float* data_pos = oriData + offset_x * dim0_offset + offset_y;

    size_t current_blockcount_x = get_current_blockcount(x, i);
    size_t current_blockcount_y = get_current_blockcount(y, j);
    size_t type_offset =
      offset_x * dim0_offset + offset_y * current_blockcount_x;
    int* type = (*memory).result_type + type_offset;

    float* unpredictable_data =
      (*memory).result_unpredictable_data + t * unpred_data_max_size;
    (*memory).unpredictable_count[t] = SZ_compress_float_2D_MDQ_RA_block(
      data_pos, (*memory).mean + t, r1, r2, current_blockcount_x,
      current_blockcount_y, realPrecision, P0 + (t * buffer_size),
      P1 + (t * buffer_size), type, unpredictable_data);
  }
  free(P0);
  free(P1);
}

inline void
_sz_compress_float_3d_mdq_ra_block(
  float* oriData, size_t r1, size_t r2, size_t r3, double realPrecision,
  int thread_num, size_t num_z, size_t unpred_data_max_size, size_t dim0_offset,
  size_t dim1_offset, int num_yz, size_t buffer_size, struct BlockCount* x,
  struct BlockCount* y, struct BlockCount* z,
  struct CompressionMemoryBlocks* memory)
{
  float* P0 = (float*)malloc(buffer_size * thread_num);
  float* P1 = (float*)malloc(buffer_size * thread_num);
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_yz);
    int j = (t % num_yz) / num_z;
    int k = t % num_z;
    size_t offset_x = get_offset(x, i);
    size_t offset_y = get_offset(y, j);
    size_t offset_z = get_offset(z, k);

    float* data_pos =
      oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

    size_t current_blockcount_x = get_current_blockcount(x, i);
    size_t current_blockcount_y = get_current_blockcount(y, j);
    size_t current_blockcount_z = get_current_blockcount(z, k);

    size_t type_offset = offset_x * dim0_offset +
                         offset_y * current_blockcount_x * dim1_offset +
                         offset_z * current_blockcount_x * current_blockcount_y;

    int* type = (*memory).result_type + type_offset;

    float* unpredictable_data =
      (*memory).result_unpredictable_data + t * unpred_data_max_size;
    (*memory).unpredictable_count[t] = SZ_compress_float_3D_MDQ_RA_block(
      data_pos, (*memory).mean + t, r1, r2, r3, current_blockcount_x,
      current_blockcount_y, current_blockcount_z, realPrecision,
      P0 + (t * buffer_size), P1 + (t * buffer_size), type,
      unpredictable_data);
  }
  free(P0);
  free(P1);
}
void
buildHuffmanTree(int thread_num, size_t num_elements,
                 struct CompressionMemoryBlocks* memory, size_t* nodeCount,
                 unsigned char** treeBytes, unsigned int* treeByteSize)
{
  (*nodeCount) = 0;
  SZ_Reset();
  Huffman_init_cuda((*memory).result_type, num_elements, thread_num);
  for (size_t i = 0; i < stateNum; i++)
    if (code[i])
      (*nodeCount)++;
  (*nodeCount) = (*nodeCount) * 2 - 1;
  *treeByteSize =
    convert_HuffTree_to_bytes_anyStates((int)*nodeCount, treeBytes);
}

inline void
copy_unpredictable(int thread_num, size_t unpred_data_max_size,
                   struct CompressionMemoryBlocks* memory,
                   unsigned char* result_pos, const size_t* unpred_offset)
{
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    float* unpredictable_data =
      (*memory).result_unpredictable_data + t * unpred_data_max_size;
    memcpy(result_pos + unpred_offset[t] * sizeof(float), unpredictable_data,
           (*memory).unpredictable_count[t] * sizeof(float));
  }
}

inline void
copyEncodingBuffers(int thread_num, size_t max_num_block_elements,
                    unsigned char* result_pos, size_t* block_pos,
                    unsigned char* encoding_buffer, const size_t* block_offset)
{
  #pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    memcpy(result_pos + block_offset[t],
           encoding_buffer + t * max_num_block_elements * sizeof(int),
           block_pos[t]);
  }
}
inline size_t
compute_compressed_size(size_t num_blocks, size_t num_elements,
                        unsigned int treeByteSize, size_t total_unpred)
{

  unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
  return meta_data_offset +                    // metadata
         sizeof(double) +                      // real precision
         sizeof(int) +                         // intervals
         sizeof(int) +                         // nodeCount
         treeByteSize +                        // huffman
         num_blocks * sizeof(unsigned short) + // block index
         num_blocks * sizeof(unsigned short) + // unpredictable count
         num_blocks * sizeof(float) +          // memory.mean
         total_unpred * sizeof(float) +        // unpred
         num_elements * sizeof(int);           // elements
}

inline size_t
compute_total_unpred(size_t num_blocks, struct CompressionMemoryBlocks* memory)
{
  size_t total_unpred = 0;
#pragma omp parallel for reduction(+:total_unpred)
  for (int i = 0; i < num_blocks; i++) {
    total_unpred += (*memory).unpredictable_count[i];
  }
  return total_unpred;
}

inline size_t
compute_total_unpred_gpu(size_t num_blocks, struct CompressionMemoryBlocks* memory)
{
	thrust::device_vector<unsigned int> dev(memory->unpredictable_count, memory->unpredictable_count + num_blocks);
	return thrust::reduce(dev.begin(), dev.end());
}

inline size_t*
compute_unpred_offset(int thread_num, size_t num_blocks,
                      struct CompressionMemoryBlocks* memory)
{
  size_t* unpred_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  unpred_offset[0] = 0;
  for (int t = 1; t < thread_num; t++) {
    unpred_offset[t] =
      (*memory).unpredictable_count[t - 1] + unpred_offset[t - 1];
  }
  return unpred_offset;
}

inline size_t*
compute_block_offsets(int thread_num, size_t num_blocks,
                      const size_t* block_pos)
{
  size_t* block_offset = (size_t*)malloc(num_blocks * sizeof(size_t));
  block_offset[0] = 0;
  for (int t = 1; t < thread_num; t++) {
    block_offset[t] = block_pos[t - 1] + block_offset[t - 1];
  }
  return block_offset;
}

inline void
config_threads_2D(int* thread_num, size_t* num_x, size_t* num_y)
{
  (*thread_num) = omp_get_max_threads();
  int thread_order = (int)log2((*thread_num));
  {
    int block_thread_order = thread_order / 2;
    switch (thread_order % 2) {
      case 0: {
        (*num_x) = 1 << block_thread_order;
        (*num_y) = 1 << block_thread_order;
        break;
      }
      case 1: {
        (*num_x) = 1 << (block_thread_order + 1);
        (*num_y) = 1 << block_thread_order;
        break;
      }
    }
    (*thread_num) = (*num_x) * (*num_y);
  }
  set_max_threads((*thread_num));
  // calculate block dims
}

inline void
config_threads_3D(int* thread_num, size_t* num_x, size_t* num_y, size_t* num_z)
{
  (*thread_num) = omp_get_max_threads();
  int thread_order = (int)log2((*thread_num));
  {
    int block_thread_order = thread_order / 3;
    switch (thread_order % 3) {
      case 0: {
        (*num_x) = 1 << block_thread_order;
        (*num_y) = 1 << block_thread_order;
        (*num_z) = 1 << block_thread_order;
        break;
      }
      case 1: {
        (*num_x) = 1 << (block_thread_order + 1);
        (*num_y) = 1 << block_thread_order;
        (*num_z) = 1 << block_thread_order;
        break;
      }
      case 2: {
        (*num_x) = 1 << (block_thread_order + 1);
        (*num_y) = 1 << (block_thread_order + 1);
        (*num_z) = 1 << block_thread_order;
        break;
      }
    }
    (*thread_num) = (*num_x) * (*num_y) * (*num_z);
  }
  set_max_threads((*thread_num));
}

inline size_t
get_current_blockcount(struct BlockCount const* const x, int const i)
{
  return ((i < x->split_index) ? x->early_blockcount : x->late_blockcount);
}
inline size_t
get_offset(struct BlockCount const* const x, int const i)
{
  return ((i < x->split_index) ? i * x->early_blockcount
                               : i * x->late_blockcount + x->split_index);
}

inline void
encode_1D(float* oriData, int thread_num, struct BlockCount* x,
          size_t max_num_block_elements, struct CompressionMemoryBlocks* memory,
          size_t* block_pos, unsigned char* encoding_buffer)
{
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t;
    unsigned char* encoding_buffer_pos =
      encoding_buffer + t * max_num_block_elements * sizeof(int);
    size_t enCodeSize = 0;
    size_t offset_x = get_offset(x, i);
    size_t current_blockcount_x = get_current_blockcount(x, i);
    size_t current_block_elements = current_blockcount_x;
    size_t type_offset = offset_x;
    int* type = (*memory).result_type + type_offset;
    encode(type, current_block_elements, encoding_buffer_pos, &enCodeSize);
    block_pos[t] = enCodeSize;
  }
}
inline void
encode_2D(int thread_num, size_t num_y, struct BlockCount* x,
          struct BlockCount* y, size_t max_num_block_elements,
          size_t dim0_offset, struct CompressionMemoryBlocks* memory,
          size_t* block_pos, unsigned char* encoding_buffer)
{
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_y);
    int j = (t % num_y);
    unsigned char* encoding_buffer_pos =
      encoding_buffer + t * max_num_block_elements * sizeof(int);
    size_t enCodeSize = 0;

    size_t offset_x = get_offset(x, i);
    size_t offset_y = get_offset(y, j);

    size_t current_blockcount_x = get_current_blockcount(x, i);
    size_t current_blockcount_y = get_current_blockcount(y, j);

    size_t current_block_elements = current_blockcount_x * current_blockcount_y;
    size_t type_offset =
      offset_x * dim0_offset + offset_y * current_blockcount_x;
    int* type = (*memory).result_type + type_offset;
    encode(type, current_block_elements, encoding_buffer_pos, &enCodeSize);
    block_pos[t] = enCodeSize;
  }
}
inline void
encode_3D(int thread_num, size_t num_z, struct BlockCount* x,
          struct BlockCount* y, struct BlockCount* z,
          size_t max_num_block_elements, size_t dim0_offset, size_t dim1_offset,
          int num_yz, struct CompressionMemoryBlocks* memory, size_t* block_pos,
          unsigned char* encoding_buffer)
{
#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int i = t / (num_yz);
    int j = (t % num_yz) / num_z;
    int k = t % num_z;
    unsigned char* encoding_buffer_pos =
      encoding_buffer + t * max_num_block_elements * sizeof(int);
    size_t enCodeSize = 0;
    size_t offset_x = get_offset(x, i);
    size_t offset_y = get_offset(y, j);
    size_t offset_z = get_offset(z, k);

    size_t current_blockcount_x = get_current_blockcount(x, i);
    size_t current_blockcount_y = get_current_blockcount(y, j);
    size_t current_blockcount_z = get_current_blockcount(z, k);

    size_t current_block_elements =
      current_blockcount_x * current_blockcount_y * current_blockcount_z;
    size_t type_offset = offset_x * dim0_offset +
                         offset_y * current_blockcount_x * dim1_offset +
                         offset_z * current_blockcount_x * current_blockcount_y;
    int* type = (*memory).result_type + type_offset;
    encode(type, current_block_elements, encoding_buffer_pos, &enCodeSize);
    block_pos[t] = enCodeSize;
  }
}

inline unsigned int
quantization_intervals_1D(float* oriData, size_t r1, double realPrecision)
{
  unsigned int quantization_intervals;
  if (optQuantMode == 1) {
    quantization_intervals =
      optimize_intervals_float_1D(oriData, r1, realPrecision);
    updateQuantizationInfo(quantization_intervals);
  } else {
    quantization_intervals = intvCapacity;
  }
  return quantization_intervals;
}
inline unsigned int
quantization_intervals_2D(float* oriData, size_t r1, size_t r2,
                          double realPrecision)
{
  unsigned int quantization_intervals;
  if (optQuantMode == 1) {
    quantization_intervals =
      optimize_intervals_float_2D_opt(oriData, r1, r2, realPrecision);
    printf("2D number of bins: %d\nerror bound %.20f\n", quantization_intervals,
           realPrecision);
    updateQuantizationInfo(quantization_intervals);
  } else {
    quantization_intervals = intvCapacity;
  }
  return quantization_intervals;
}
inline unsigned int
quantization_intervals_3D(float* oriData, size_t r1, size_t r2, size_t r3,
                          double realPrecision)
{
  unsigned int quantization_intervals;
  if (optQuantMode == 1) {
    // quantization_intervals = optimize_intervals_float_3D(oriData, r1,
    // realPrecision);
    quantization_intervals =
      optimize_intervals_float_3D_opt(oriData, r1, r2, r3, realPrecision);
    printf("3D number of bins: %d\nerror bound %.20f\n", quantization_intervals,
           realPrecision);
    // exit(0);
    updateQuantizationInfo(quantization_intervals);
  } else {
    quantization_intervals = intvCapacity;
  }
  return quantization_intervals;
}

inline unsigned char*
write_parallel_compresion_metadata(unsigned char* result_pos, int thread_num,
                                   double realPrecision,
                                   unsigned int quantization_intervals,
                                   struct CompressionMemoryBlocks* memory,
                                   size_t num_blocks, size_t nodeCount,
                                   const unsigned char* treeBytes,
                                   unsigned int treeByteSize)
{
  result_pos = writeIntBigEndian(result_pos, thread_num);
  result_pos = writeDoubleBigEndian(result_pos, realPrecision);
  result_pos = writeIntBigEndian(result_pos, quantization_intervals);
  result_pos = writeIntBigEndian(result_pos, treeByteSize);
  result_pos = writeIntBigEndian(result_pos, nodeCount);
  result_pos = writeBytes(result_pos, (unsigned char*)treeBytes, treeByteSize);
  result_pos =
    writeBytes(result_pos, (unsigned char*)(*memory).unpredictable_count,
               num_blocks * sizeof(unsigned int));
  result_pos = writeBytes(result_pos, (unsigned char*)(*memory).mean,
                          num_blocks * sizeof(float));
  return result_pos;
}

inline unsigned char*
writeBytes(unsigned char* output, const unsigned char* bytes, unsigned int size)
{
  memcpy(output, bytes, size);
  output += size;
  return output;
}
inline unsigned char*
writeDoubleBigEndian(unsigned char* output, double d)
{
  doubleToBytes(output, d);
  output += 8;
  return output;
}
inline unsigned char*
writeIntBigEndian(unsigned char* output, int i)
{
  intToBytes_bigEndian(output, i);
  output += 4;
  return output;
}

inline int
readIntBigEndian(unsigned char** data)
{
  int ret = bytesToInt_bigEndian(*data);
  (*data) += 4;
  return ret;
}

void
Huffman_init_cuda(int* s, size_t length, int thread_num)
{

  size_t i;
  size_t* freq = (size_t*)malloc(thread_num * allNodes * sizeof(size_t));
  memset(freq, 0, thread_num * allNodes * sizeof(size_t));
  size_t block_size = (length - 1) / thread_num + 1;
  size_t block_residue = length - (thread_num - 1) * block_size;

#pragma omp parallel for
  for (int t = 0; t < thread_num; t++) {
    int* s_pos = s + t * block_size;
    size_t* freq_pos = freq + t * allNodes;
    if (t < thread_num - 1) {
      for (size_t i = 0; i < block_size; i++) {
        freq_pos[s_pos[i]]++;
      }
    } else {
      for (size_t i = 0; i < block_residue; i++) {
        freq_pos[s_pos[i]]++;
      }
    }
  }

  size_t* freq_pos = freq + allNodes;
  for (int t = 1; t < thread_num; t++) {
    for (i = 0; i < allNodes; i++) {
      freq[i] += freq_pos[i];
    }
    freq_pos += allNodes;
  }

  for (i = 0; i < allNodes; i++)
    if (freq[i])
      qinsert(new_node(freq[i], i, 0, 0));

  while (qend > 2)
    qinsert(new_node(0, 0, qremove(), qremove()));

  build_code(qq[1], 0, 0, 0);
  free(freq);
}
