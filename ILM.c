#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


uint8_t LOD(uint8_t val){
    uint32_t n = 0, x;
    x = val;
    if (x <= 0x0000ffff) n += 16, x <<= 16;
    if (x <= 0x00ffffff) n += 8, x <<= 8;
    if (x <= 0x0fffffff) n += 4, x <<= 4;
    if (x <= 0x3fffffff) n += 2, x <<= 2;
    if (x <= 0x7fffffff) n++;
    return 31 - n;
}

uint16_t ILM(uint8_t a, uint8_t b, uint8_t iter){
    /*
        a, b -> input operands,
        iter -> number of iterations
        only two iterations supported
    */
    if (a == 0 || b == 0) return 0;

    uint8_t Ka, Kb; 
    Ka = LOD(a);
    Kb = LOD(b);

    uint8_t ResA, ResB, Res2B;
    ResA = a ^ (1 << Ka);
    ResB = b ^ (1 << Kb);

    uint16_t prod0, prod1;
    prod0 = a * (1<<Kb) + ResB * (1<<Ka);
    prod1 = 0;
    if(iter == 2){
        Ka = LOD(ResA);
        Kb = LOD(ResB);
        Res2B = ResB ^ (1 << Kb);
        prod1 = ResA * (1<<Kb) + Res2B * (1<<Ka);
    }

    return prod0 + prod1;
}

int main(){
    srand(time(0));
    uint8_t a = rand();
    uint8_t b = rand();
    printf("%d", ILM(a,b,2));
    uint32_t i, max;
    uint16_t p_exact, p_approx;
    FILE *log;

    log = fopen("log/ILM_log.txt", "a");
    char str1[80];

    max = 20;
    for(i=0;i<max;i++){
        a = rand();
        b = rand();
        p_exact =  a * b;
        p_approx = ILM(a,b,2);

        sprintf(str1, "A : %d, B : %d, TP : %d, AP : %d \n", a,b,p_exact,p_approx);
	    fputs(str1, log);
        

    }
    return 0;
}

