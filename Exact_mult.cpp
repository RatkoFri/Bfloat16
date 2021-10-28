#include <stdint.h>
#include <stdio.h>



void float2bfloat(const float src, float& dst) {
	const uint16_t* p = reinterpret_cast<const uint16_t*>(&src);
	uint16_t* q = reinterpret_cast<uint16_t*>(&dst);
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
  q[0] = p[0];
  q[1] = 0;
#else
	q[0] = 0;
	q[1] = p[1];
#endif
}


int main(){
	float a = 2.5;
	float b = 3.1245;
	float tmp_b,tmp_a;

	float2bfloat(b,tmp_b);
	float2bfloat(a,tmp_a);

	printf("True product %f\n", a*b);
	printf("Reduced precision product %f", tmp_a*tmp_b);

	return 0;
}