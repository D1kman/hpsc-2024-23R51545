#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <x86intrin.h>
#include <vector>

void print_m512(__m512 vec) {
    float values[16];
    _mm512_store_ps(values, vec);

    for (int i = 0; i < 16; i++) {
        printf("%f ", values[i]);
    }
    printf("\n");
}

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], array[2*N], brray[N], irray[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  
  for (int i=0; i<N; i++){
	  array[i]=i;
  	  array[i+N]=(-1);
  }		  
  for(int i=0; i<N; i++) {
    
	__m512 avec = _mm512_load_ps(array);   
        __m512 bvec = _mm512_setzero_ps();
	__m512 ivec = _mm512_set1_ps(i);
        __mmask16 mask = _mm512_cmp_ps_mask(avec, ivec, _MM_CMPINT_NE);
	__mmask16 mask2 = _mm512_cmp_ps_mask(avec, bvec, _MM_CMPINT_GE);
        __m512 xjvec = _mm512_load_ps(x);
	__m512 xivec = _mm512_set1_ps(x[i]);
	__m512 yjvec = _mm512_load_ps(y);
	__m512 yivec = _mm512_set1_ps(y[i]);
	
	__m512 rxvec = _mm512_sub_ps(xivec, xjvec);
	__m512 ryvec = _mm512_sub_ps(yivec, yjvec);
	
	__m512 rx2vec = _mm512_mul_ps(rxvec, rxvec);
	__m512 ry2vec = _mm512_mul_ps(ryvec, ryvec);
	__m512 rvec = _mm512_add_ps(rx2vec, ry2vec);
	
	rvec = _mm512_rsqrt14_ps(rvec);
	ivec = _mm512_set1_ps(1);
	rvec = _mm512_div_ps(ivec, rvec);
	__m512 r2vec  = _mm512_mul_ps(rvec, rvec);
	rvec = _mm512_mul_ps(rvec, r2vec);

	__m512 mvec = _mm512_load_ps(m);
	rxvec = _mm512_mul_ps(rxvec, mvec);
	ryvec = _mm512_mul_ps(ryvec, mvec);
	
	rxvec = _mm512_div_ps(rxvec, rvec);
	ryvec = _mm512_div_ps(ryvec, rvec);
	rxvec = _mm512_mask_blend_ps(mask, bvec, rxvec);
	ryvec = _mm512_mask_blend_ps(mask, bvec, ryvec);
	rxvec = _mm512_mask_blend_ps(mask2, bvec, rxvec);
	ryvec = _mm512_mask_blend_ps(mask2, bvec, ryvec);
	
	fx[i] = _mm512_reduce_add_ps(rxvec);
	fy[i] = _mm512_reduce_add_ps(ryvec);
	fx[i] = -fx[i];
	fy[i] = -fy[i];
//        float rx = x[i] - x[j];
//        float ry = y[i] - y[j];
//        float r = std::sqrt(rx * rx + ry * ry);
//        fx[i] -= rx * m[j] / (r * r * r);
//        fy[i] -= ry * m[j] / (r * r * r);
     
    
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
