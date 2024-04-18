#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <x86intrin.h>
#include <vector>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], array[N], brray[N], irray[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  
  for (int i=0; i<N; i++)
	  array[i]=i;
  for(int i=0; i<N; i++) {
    
	__m512 avec = _mm512_load_ps(array);   
        __m512 bvec = _mm512_setzero_ps();
	__m512 ivec = _mm512_set1_ps(i);
        __mmask16 mask = _mm512_cmp_ps_mask(avec, ivec, _MM_CMPINT_NE);
	 
        __m512 xjvec = _mm512_load_ps(x);
	__m512 xivec = _mm512_set1_ps(x[i]);
	__m512 yjvec = _mm512_load_ps(y);
	__m512 yivec = _mm512_set1_ps(y[i]);

//	_mm512_mask_blend_ps(mask, bvec, xjvec);
	__m512 rxvec = _mm512_sub_ps(xivec, xjvec);
	__m512 ryvec = _mm512_sub_ps(yivec, yjvec);

	__m512 rx2vec = _mm512_mul_ps(rxvec, rxvec);
	__m512 ry2vec = _mm512_mul_ps(ryvec, ryvec);
	__m512 rvec = _mm512_add_ps(rx2vec, ry2vec);
	rvec = _mm512_rsqrt14_ps(rvec);



        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
     
    
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
