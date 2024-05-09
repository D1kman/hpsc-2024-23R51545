#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucket_sort(int *bucket, int *key,  int range, int n) {
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  if(id < range) {
  bucket[id] = 0;
  }
  __syncthreads();

  atomicAdd(&bucket[key[id]],1);
  __syncthreads();
  int offset = 0;
  for (int i = 0; i < range; i++) {
	if (id  < bucket[i]+offset && id >= offset) {
	  key[id] = i;
	}
	  offset += bucket[i];  
  }

}

int main() {
  int n = 50;
  int range = 5;
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
 
  
  bucket_sort<<<1, n>>>(bucket, key, range, n);
  
  cudaDeviceSynchronize();
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");

  cudaFree(bucket);
}
