/*Author : Sohan Das, Indian Statistical Institute, KOlkata*/
#ifndef ECC_DOT_H
#define ECC_DOT_H

void add(long long int a[], long long int b[], long long int c[], int size);  // c = a + b
void mult(long long int a[], long long int b[], long long int m[], int size); // m = a * b
void sub(long long int a[], long long int b[], long long int s[], int size);  // s = a - b  with assumption a > b
void karat_1M1(long long int a[], long long int b[], long long int k_ab[], int size); // where size = sizeof(intput) = sizeof(a) = sizeof(b)
/* Actually we are considerinf A = a0 + a1*x, B = b0 + b1*x,
    therefore   AB = (a0*a1) + {(a0 + a1)*(b0 + b1) - a0*a1 - b0*b1}*x + (b0*b1)*x^2 
    but since we are considering 2^28 base so (b0*b1) may be greater than (2^28 - 1), then we have to shift those bits into next block, so we are taking size of AB = 2*size(a) = 2*2 = 4 */ 
void baseManagement(long long int input[], long long int outputbase[], int size, int base);  // This will prepare all the array elements on 2^28 base 
void karat_3M3(long long int a[], long long int b[], long long int k_ab[], int size);
// This karat_3M3 multiplies 3 degree polynomial with 3 degree polynomial using KARATSUBA multiplication
void karat_4M4(long long int a[], long long int b[], long long int k_ab[], int size);
//The karat_9M9 multiplies two 9 degree poly with 10 coefficients using KARATSUBA multiplication 
void karat_9M9(long long int a[], long long int b[], long long int k_ab[], int size);

void sub_mod_p(long long int a[], long long int b[], long long int s[], int size);  // subtraction in Z_p
int check(long long int x[], long long int y[], int size); // returns 1 when x>y, 2 when x<y , 3 when x=y

void Barrett(long long int x[], long long int x_Barrett_p[], int size);                                                                                // Barrett reduction with ( mod p)
void sqm_mod_p(long long int y[], long long int d[], long long int z[], int size);                                                                     // z = y^b (mod p), size = 10,here
void ECC_addition(long long int x1[10], long long int y1[10], long long int x2[10], long long int y2[10], long long int x3[10], long long int y3[10]); // (x1,y1)+(x2,y2) = (x3,y3) ECC Addition
void Doubling_Adding(long long int K[10], long long int p_x[10], long long int p_y[10], long long int k_Px[10], long long int k_Py[10] );   // K*(p_x , p_y) = (k_Px , k_Py)

#endif
