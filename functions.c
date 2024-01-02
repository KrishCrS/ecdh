/* ----------//////////----------- Function Definitions -----------///////////----------*/
/*Implemented by SOHAN DAS, an M.Tech(CrS) student of Indian Statistical Institute, KOLKATA*/

#include <stdio.h>
#include "ECC.h"

long long int p[10] = {268435455, 268435455, 268435455, 4095, 0, 0, 16777216, 0, 268435455, 15}; // 256-bit prime
// p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff (NIST Suggested)
long long int mu[11] = {327680, 3145728, 0, 0, 268435455, 268435439, 268435199, 268431359, 268435455, 1048575, 16777216}; // precomputed "mu" for Barrett reduction


void add(long long int a[], long long int b[], long long int c[], int size) // a and b will be added and the sum will be stored in c
{
    long long int carry = 0;
    for(int i =0; i<size; i++)
        c[i] = 0;
    for (int i = 0; i < size; i++)
    {
        c[i] = a[i] + b[i] + carry;
        carry = (c[i] >> 28) & 1;
        c[i] = c[i] & ((1 << 28) - 1); // (1<<28) - 1 = 2^28 - 1

        // printf(" %lld", c[i]);
    }
}


void mult(long long int a[], long long int b[], long long int m[], int size) // a and b will be multiplied and the product will be stored in m
{
    long long int restof28bit;
    for (int i = 0; i < 2 * size; i++)
    {
        m[i] = 0;
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            m[i + j] = m[i + j] + (a[i] * b[j]);
            restof28bit = m[i + j] >> 28;
            m[i + j] = m[i + j] & ((1 << 28) - 1);     // (1<<28) - 1 = 2^28 - 1
            m[i + j + 1] = m[i + j + 1] + restof28bit; // Carry forwarding
        }
    }

    // for (int k = 0; k < 2 * size; k++)
    //     printf(" %lld", m[k]);
}


void sub(long long int a_o[], long long int b_o[], long long int s[], int size) // s = a_o - b_o with the assumption that a_o > b_o , o stands for original
{   long long int a[size];
    long long int b[size];
    int i;
    for (i = 0; i < size; i++)  // initialization
    {
        s[i] = 0;
        a[i] = a_o[i];
        b[i] = b_o[i];
    }
        
    // for (i = 0; i < size; i++)
    // {
    //     if (a[i] < b[i])
    //     {
    //         a[i] = a[i] + (1 << 28);
    //         b[i + 1] = b[i + 1] + 1;
    //     }
    //     s[i] = a[i] - b[i];

    //     // printf(" %lld", s[i]);
    // }

    /* Using 2's complemnet method */
    long long int I[size], b_2s_com[size];
    for (i = 0; i < size; i++)
    {
        I[i] = 0;
        b[i] = (~b[i])&((1<<28)-1);  // 1's complement is done 
    }
    I[0] = 1;  // I = 1
    add(b, I, b_2s_com, size);  // 2's complement of b is done
    add(a, b_2s_com, s, size);  // s = a + b_2s_com
}





/*  ------------------  KARATSUBA MULTIPLICATIONS  -------------   */


// This karat_1M1 multiplies 1 degree polynomial with 1 degree polynomial using KARATSUBA multiplication
void karat_1M1(long long int a[], long long int b[], long long int k_ab[], int size)  // size = sizeof(a) = 2 = degree+1, k_ab = AB, where A = a0+a1*x , B = b0+b1*x
{
    int i;
    long long int coeffi[2*size]; // it will store coefficients of x^i 's for i = 0,1,2
    // coeffi[3] = 0 is taken for getting advantage in building for loop in base management part below
    for(i=0; i<2*size; i++)  // initialization
    {
        k_ab[i] = 0;
        coeffi[i] = 0;
    }
    coeffi[0] = a[0]*b[0];  //coefficient of x^0
    coeffi[2] = a[1]*b[1];  //coefficient of x^2
    coeffi[1] = (a[0] + a[1])*(b[0] + b[1]) - (coeffi[0] + coeffi[2]); //coefficient of x^1

    // base management part
    baseManagement(coeffi, k_ab, 2*size, 28);
}

// This will prepare all the array elements on 2^28 base 
void baseManagement(long long int input[], long long int outputbase[], int size, int base)  // size means sizeof(outputbase) = sizeof(input)
{
    int i;
    long long int inp_coeffi[size], carry = 0;
    for (i=0; i<size; i++)   // initialization 
    {
        inp_coeffi[i] = input[i];
        outputbase[i] = 0;
    }
    outputbase[0] = inp_coeffi[0]&((1<<base)-1);    // (1<<28) - 1 = 111.....11 (base number of 1's)
    for(i = 1; i<size; i++)
    {
        carry = inp_coeffi[i-1]>>base;  // carry from previous block
        inp_coeffi[i] = inp_coeffi[i]+carry;  // update the current block
        outputbase[i] = inp_coeffi[i]&((1<<base)-1);  // ((1<<size)-1) = 1111...11 , size no. of 1's
    }
}



/* This multiplies two 3 degree plonomials using   karat_1M1  :  
  (a0 + a1*x + a2*x^2 + a3*x^3)*(b0 + b1*x + b2*x^2 + b3*x^3)
= (A1 + A2*x^2)*(B1 + B2*x^2)  , where A1 = a0 + a1*x, A2 = a2 + a3*x, B1 = b0 + b1*x, B2 = b2 + b3*x
= A1*B1 + {(A1+A2)*(B1+B2)}*x^2 + (A2*B2)*x^4  */
void karat_3M3(long long int a[], long long int b[], long long int k_ab[], int size)  // here size will be 4 = 3+1 = degree + 1
{
    int i, j;
    long long int A1[2] = {a[0], a[1]}, A2[2] = {a[2],a[3]}, B1[2] = {b[0], b[1]}, B2[2] = {b[2], b[3]};  //initialization
    long long int coeffi[3][4] = {0};  //initalization 
    karat_1M1(A1,B1,coeffi[0],2);  // A1B1  is calculated & stored in coeffi[0]
    karat_1M1(A2, B2, coeffi[2], 2);  // A2B2  is calculated & stored in coeffi[2]

    long long int sum_A1_A2[2] = {(a[0]+a[2]), (a[1]+a[3])}, sum_B1_B2[2] = {(b[0]+b[2]), (b[1]+b[3])};
    long long int temp[4] = {0}, temp2[4] = {0};
    karat_1M1(sum_A1_A2, sum_B1_B2, temp, 2);  // temp = (A1 + A2)*(B1 + B2) is done
    sub(temp, coeffi[0], temp2, 4);   // temp2 = (A1 + A2)*(B1 + B2) - A1B1  is done
    sub(temp2, coeffi[2], coeffi[1], 4);   // coeffi[2] = (A1 + A2)*(B1 + B2) - A1B1 - A2B2  is done

    long long int pre_k_ab[8] = {0}; // it will manage postion of all the coefficients before base management, actually makes an 8 bit array from 3*4 = 12 bit array
    for(j = 0; j<4; j++)
        pre_k_ab[j] += coeffi[0][j];
    for(i=1; i<3; i++)
        for(j = 0; j<4; j++)
            pre_k_ab[((1<<i)+j)] = pre_k_ab[((1<<i)+j)] + coeffi[i][j];   // (1<<i) + j = 2^i +j

    baseManagement(pre_k_ab, k_ab, 8, 28);
}


/* This karat_4M4() multiplies two 4 degree polynomial using karat_3M3() as follows :
  (a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4)*(b0 + b1*x + b2*x^2 + b3*x^3 + b4*x^4)
= (c0 + C1*x)*(d0 + D1*x)  , where c0 = a0, C1 = a1 + a2*x + a3*x^2 + a4*x^3  &  d0 = b0, D1 = b1 + b2*x + b3*x^2 + b4*x^3
= c0d0 + {(c0 + C1)*(d0 + D1) - c0d0 - C1D1}*x + C1D1*x^2
 */
void karat_4M4(long long int a[], long long int b[], long long int k_ab[], int size)
{
    int i;
    long long int c0 = a[0], d0 = b[0];
    long long int C1[4] = {0}, D1[4] = {0};
    for(i=0; i<4; i++)
    {
        C1[i] = a[i+1];
        D1[i] = b[i+1];
    }
    long long int c0d0 = 0, C1D1[8] = {0}, sum_c0_C1[4] = {0}, sum_d0_D1[4] = {0}, pr_of_sums[8] = {0}, coeffi_of_x[8] = {0}, coef_x[8] = {0}; 
    c0d0 = a[0]*b[0];
    karat_3M3(C1, D1, C1D1, size-1);   //C1D1 = C1 * D1 is done
    
    C1[0] = C1[0] + c0;   // c0 + C1  is done
    D1[0] = D1[0] + d0;   // d0 + D1  is done
    baseManagement(C1, sum_c0_C1, 4, 28); 
    baseManagement(D1, sum_d0_D1, 4, 28);
    karat_3M3(sum_c0_C1, sum_d0_D1, pr_of_sums, size-1);   // pr_of_sums = (c0 + C1)*(d0 + D1)  is done
    sub(pr_of_sums, C1D1, coeffi_of_x, 8);   // coeffi_of_x = (c0 + C1)*(d0 + D1) - C1D1  is done
    long long int cd[8] = {0}, ext_c0d0[8] = {0};   // extended for doing subtraction in same size
    cd[0] = c0d0;   //to perform subtraction we have put down the same base, so we manage the base
    baseManagement(cd, ext_c0d0, 8, 28);  // same base number is gotten

    sub(coeffi_of_x, ext_c0d0, coef_x, 8);   // coef_x = (c0 + C1)*(d0 + D1) - C1D1 - c0d0  is done
    long long int arrange_coeffi_ab[10] = {0};  
    // sum up coefficients of corresponding degrees
    arrange_coeffi_ab[0] = ext_c0d0[0];
    arrange_coeffi_ab[1] = ext_c0d0[1]+coef_x[0];
    for(i=2; i<8; i++)
        arrange_coeffi_ab[i] = ext_c0d0[i] + coef_x[i-1] + C1D1[i-2];
    arrange_coeffi_ab[8] = coef_x[7] + C1D1[6];
    arrange_coeffi_ab[9] = C1D1[7];
    
    baseManagement(arrange_coeffi_ab, k_ab, 10, 28);

 }

/* This will multiply two 9 degree polynomial using karat_4M4 as follows :
  (a0+a1x+a2x^2+....+a9x^9)*(b0+b1x+b2x^2+....+b9x^9) 
= (E0 + E1*x^5)*(F0 + F1*x^5) , where E0 = a0+a1x+a2x^2+a3x^3+a4x^4, E1 = a5+a6x+a7x^2+a8x^3+a9x^4 
                                      F0 = b0+b1x+b2x^2+b3x^3+b4x^4, E1 = b5+b6x+b7x^2+b8x^3+b9x^4  
= (E0F0) + {(E0+E1)*(F0+F1) - E0F0 - E1F1}X^5 + E1F1*x^10 */

void karat_9M9(long long int a[], long long int b[], long long int k_ab[], int size)
{
    int i,j;
    long long int E[2][5] = {0}, F[2][5] = {0};
    for(i=0;i<2;i++)   // these are initializations
        for(j=0; j<5; j++)
        {   if(i==1)        // So,  
                i +=5;
            E[i][j] = a[i+j];
            F[i][j] = b[i+j];
        }
    long long int E0F0[10] = {0}, E1F1[10] = {0}, Coeff_x5[10] = {0};
    karat_4M4(E[0], F[0], E0F0, 5);   // E0F0 = E0 * F0  calculated
    karat_4M4(E[1], F[1], E1F1, 5);   // E1F1 = E1*F1  is calculated
    long long int temp1[10] = {0}, temp2[10] = {0};   // to manage the base
    for(i=0; i<10; i++)
        {
            temp1[i] = E0F0[i];
            E0F0[i] = 0;

            temp2[i] = E1F1[i];
            E1F1[i] = 0;
        }
    baseManagement(temp1, E0F0, 10, 28);   // E0F0 = E0 * F0 is calculated
    baseManagement(temp2, E1F1, 10, 28);   // E1F1 = E1 * F1 is calculated
    long long int sum_of_E0_E1[5] = {0}, sum_of_F0_F1[5] = {0};
    for(i = 0; i<5; i++)
    {
        sum_of_E0_E1[i] = E[0][i]+F[0][i];
        sum_of_F0_F1[i] = E[1][i]+F[1][i];
        temp1[i] = 0;    // reset 
        temp2[i] = 0 ;   //reset
    }
    long long int product[10] = {0};  // to store product of above two sums
    karat_4M4(sum_of_E0_E1, sum_of_F0_F1, product, 10);

    sub(product, E0F0, temp1,10);   // temp1 = (E0+E1)*(F0+F1) - E0F0
    sub(temp1, E1F1, temp2,10);  // temp2 = (E0+E1)*(F0+F1) - E0F0 - E1F1
    
    // coeffiection management using correct degree
    long long int final[20] = {0};
    final[0] = E0F0[0];
    final[1] = E0F0[1];
    final[2] = E0F0[2];
    final[3] = E0F0[3];
    final[4] = E0F0[4];
    final[5] = E0F0[5] + temp2[0];
    final[6] = E0F0[6] + temp2[1];
    final[7] = E0F0[7] + temp2[2];
    final[8] = E0F0[8] + temp2[3];
    final[9] = E0F0[9] + temp2[4];
    final[10] = temp2[5] + E1F1[0];
    final[11] = temp2[6] + E1F1[1];
    final[12] = temp2[7] + E1F1[2];
    final[13] = temp2[8] + E1F1[3];
    final[14] = temp2[9] + E1F1[4];
    final[15] = E1F1[5];
    final[16] = E1F1[6];
    final[17] = E1F1[7];
    final[18] = E1F1[8];
    final[19] = E1F1[9];

    baseManagement(final, k_ab, 20, 28);
}

/*   -----------------  KARATSUBA MULTIPLICATION ENDS ---------------*/



// check whether  x > y or x = y or x < y
// return  1 when x>y
//         2 when x<y
//         3 when x=y
int check(long long int x[], long long int y[], int size)
{
    int flag = 3, i;
    for (i = size - 1; i >= 0; i--)
    {
        if (x[i] > y[i])
        {
            flag = 1;
            return flag;
        }
        if (x[i] < y[i])
        {
            flag = 2;
            return flag;
        }
    }
    return flag;
}



void sub_mod_p(long long int a[], long long int b[], long long int s[], int size) // s = a - b  if a>b , s = p - (b-a) if a<b
{
    long long int pr[size];
    for (int i = 0; i < 10; i++)
        pr[i] = p[i];
    if(size>10)
    {
        for (int i = 10; i < size; i++)
        pr[i] = 0;
    }
    int flag = check(a,b,size);  // flag = 1 when a>b , flag = 2 when a<b and flag = 3 when a=b
    if ((flag == 1) || (flag == 3)) // when a >= b
        sub(a, b, s, size);
    if (flag == 2) // when a < b , then subtraction in mod p will be p - (b - a)
    {
        long long int temp[size]; // for storing b - a
        for (int i = 0; i < size; i++)
            temp[i] = 0;
        sub(b, a, temp, size);  // temp = b - a
        sub(pr, temp, s, size); // s = p - temp = p - (b - a) = p + (a - b)
        // printf("\n in sub_mod_p flag 2 occurs ");
    }
}


void Barrett(long long int x[], long long int x_Barret_p[], int size) // Barrett reduction with ( mod p), size means size of input x[]
{
    long long int q[11]; // since size of "mu" is 28*11 bits
    for (int i = 0; i < 11; i++)
    {
        q[i] = 0; // initialization
    }
    // q = floor(x/{B^(k-1)}) , k = number of 28-bit blocks used to store the value of p, so k = 10, here and B = 2^28, so to divide x by B^9 we will right shift x by 9*28 bits and store it in q
    int i = 0;
    while (i < (size - 9))
    {
        q[i] = x[i + 9];
        i++;
    }                      // q = floor(x/{B^(k-1)}) is done.
    long long int q22[22]; // to store  q*mu , both of them are of 11*28 bits so their product may be of 2*11*28 bit long
    for (int i = 0; i < 22; i++)
    {
        q22[i] = 0; // initialization
    }
    mult(q, mu, q22, 11);  // q22 = q*mu  is done
    long long int q11[11]; // since q11 = floor(q22/{B^(k+1)}) , q22 is of 22*28 bits and k = 10
    for (int i = 0; i < 11; i++)
    {
        q11[i] = 0; // initialization
    }
    int j = 0; // initialization
    while (j < 11)
    {
        q11[j] = q22[j + 11];
        j++;
    } // q11 = floor(q22/{B^(k+1)}) is done
    // r = x - q11p11 (mod B^(k+1)), B = 2^28, k= 10

    long long int q11p11[22] = {0}; // sizeof(q11) = 11 block, sizeof(p) = 10, so we have to extend it to 11 blocks and we say it p11 also have to extend sizeof(x) = 20 to 22 blocks
    long long int p11[11];          // extending block size of p from 10 to 11
    for (int i = 0; i < 10; i++)    // since sizeof(p) is 10 blocks means 10 long long int is used
    {
        p11[i] = p[i];
    }
    p11[10] = 0;
    mult(p11, q11, q11p11, 11);  // q11p11 = p11 * q11

    // r = x - q11p11 (mod B^(k+1)), first we will calculate x(mod B^(k+1)) and q11p11(mod B^(k+1)) then will subtract in  mod B^(k+1)
    long long int x11[11] = {0}; // for storing the value x (mod B^(k+1)), k=10 here
    long long int qp[11] = {0};  // for storing the value q11p11 (mod B^(k+1)), k=10 here
    for (int i = 0; i < 11; i++)
    {
        x11[i] = x[i];
        qp[i] = q11p11[i];
    }
    long long int r11[11] = {0};  // r11 = x11 - qp (mod B^(k+1))
    int flg = check(x11, qp, 11);  // flg = 1 if x11>qp, flg = 2 if x11<qp, flag = 3 if x11=qp
    if((flg == 1) || (flg == 3))  // if x11 >= qp
        sub(x11, qp, r11, 11);
    if(flg == 2)  // if x11 < qp
        {
            long long int tmp[11] = {0};
            sub(qp, x11, tmp, 11);  // temp = qp - x11  is done 
            //to do mod B^(k+1) we will use 2's complement method as follows
            for(int i = 0; i<11; i++)
            {
                tmp[i] = (~tmp[i])&((1<<28)-1);  // 1's complement in 28-bit is done
            }
            long long int I[11] = {0};
            I[0] = 1;  
            add(tmp, I, r11, 11);  // 2's complement is done
        }

    // Checking whether r >= p or not
    int flag = check(r11, p11, 11); // flag = 1 if r11>p11, flag = 2 if r11<p11, flag = 3 if r11=p11

    while ((flag == 1) || (flag == 3)) // do successive subtraction when r >= p
    {
        // printf("\n Entered in Barrett While loop   flag = %d\n", flag);
        long long int temp[11] = {0}; // to copy r11
        for (int i = 0; i < 11; i++)
        {
            temp[i] = r11[i];
            r11[i] = 0;
        }
        sub(temp, p11, r11, 11); // r11 = temp - p11 = r11(old) - p11

        // Checking whether r >= p or not
        flag = check(r11, p11, 11); // flag = 1 when r > p, flag = 2 when r < p, flag = 0 when r =p
    }
    static int count = 0;
    count++;
    // printf("\n You are checking Barrett  %d times", count);
    for (int i = 0; i < 10; i++)
    {
        x_Barret_p[i] = r11[i];
        // printf("\n x_Barrett_p[%d] = %lld", i, x_Barret_p[i]);
    }
}

// Square and multiply ( used for finding inverse in Z*_p)
void sqm_mod_p(long long int y[], long long int d[], long long int z[], int size)
{
    z[0] = 1;
    long long int z2[2 * size]; // for storing the value of z^2
    for (int i = 0; i < 2 * size; i++)
    {
        z2[i] = 0; // initialization
        // printf(" z2[%d] = %lld\n", i, z2[i]);
    }
    for (int i = size - 1; i >= 0; i--) // since sizeof(d) = 10 block = 10 *28 bit
    {
        for (int j = 27; j >= 0; j--) // since each block is of 28-bit long
        {
            mult(z, z, z2, size); // z2 = z*z, 10 is the size of inputs z's
            for (int j = 0; j < 10; j++)
                z[j] = 0;               // initialization of before storing values gotten from Barrett
            Barrett(z2, z, 2 * size);   // z = z^2 (mod p) is done
            long long int zy[2 * size]; // for storing z*y
            for (int i = 0; i < 2 * size; i++)
                zy[i] = 0; // initialization
            if (((d[i] >> j) & 1) == 1)
            {
                mult(z, y, zy, size); // zy = z*y
                for (int k = 0; k < 10; k++)
                    z[k] = 0;             // initialization before storing values gotten from Barrett
                Barrett(zy, z, 2 * size); // z = zy (mod p) is done
            }
        }
    }
}

void ECC_addition(long long int x_1[10], long long int y_1[10], long long int x_2[10], long long int y_2[10], long long int x3[10], long long int y3[10]) // (x1,y1)+(x2,y2) = (x3,y3) ECC Addition
{
    int i;
    long long int x1[10], y1[10], x2[10], y2[10];
    for (i = 0; i < 10; i++)
    {
        x3[i] = 0;
        y3[i] = 0;
        x1[i] = x_1[i];
        y1[i] = y_1[i];
        x2[i] = x_2[i];
        y2[i] = y_2[i];
    }
    long long int A[20] = {0}; // A = p-3, we take it's size 20 blocks since we will add it with 3x1^2, which is of 20 block
    for (i = 0; i < 10; i++)
        A[i] = p[i];
    A[0] = p[0] - 3; // A = (p - 3) is done
    // Now we will compare whether P1=(x1,y1) and P2=(x2,y2) are equal or not
    // lamda = (y2 - y1)/(x2 - x1) ,when P1 != P2
    // lamda = (3x1^2 + A)/(2*y1) ,when P1 = P2
    long long int lamda[10] = {0};
    // for comparing P1, P2 we will set flag_x, flag_y using check() function which
    // returns 1 when x>y
    //         2 when x<y
    //         3 when x=y  on comparing x and y in check
    int flag_x = check(x1, x2, 10);
    int flag_y = check(y1, y2, 10);
    // printf("flag_x =  %d ,  flag_y = %d", flag_x, flag_y);
    // if P1 == P2 calculate lamda = (3x1^2 + A)/(2*y1)
    if ((flag_x == 3) && (flag_y == 3))
    {
        long long int x1_sq[20] = {0}; // x_sq is for storing x1^2
        mult(x1, x1, x1_sq, 10);
        long long int x1_sq_2sums[20] = {0};   // will be used to store 2*x1^2
        add(x1_sq, x1_sq, x1_sq_2sums, 20);    // 2*x1^2 is stored in x1_sq_2sums
        long long int x1_sq_3sums[20] = {0};   // will be used to store 3*x1^2
        add(x1_sq, x1_sq_2sums, x1_sq_3sums, 20); // 3*x1^2 is done and stored in x1_sq_3sums
        // Now we have add 3*x1^2 with A
        long long int sum_3x1_sq_and_A[20] = {0};
        add(x1_sq_3sums, A, sum_3x1_sq_and_A, 20); // 3*x1^2 + A is done
        // use Barrett to reduce it into (mod p)
        long long int addedVal_mod_p[10] = {0};  // will be used to store reduced value of sum_3x1_sq_and_A in mod p
        Barrett(sum_3x1_sq_and_A, addedVal_mod_p, 20);
        //Calculate inverse of 2*y1 in Z*_p, we will use sqm_mod_p() to find inverse
        long long int y1_2sums[10] = {0};
        add(y1,y1,y1_2sums,10);
        long long int inverse_of_2y1[10] = {0};  // y^(p-2) will be the inverse of y in Z*_p
        long long int p_2[10] = {268435453, 268435455, 268435455, 4095, 0, 0, 16777216, 0, 268435455, 15};  //p_2 is actually  p-2
        sqm_mod_p(y1_2sums, p_2, inverse_of_2y1, 10);
        long long int ad_in[20] = {0};  // to store   addedVal_mod_p * inverse_of_2y1
        mult(addedVal_mod_p, inverse_of_2y1, ad_in, 10);
        // ad_in (mod p) will be stored in lamda
        Barrett(ad_in, lamda, 20);
    }

    else
    {
        long long int y2_y1[10] = {0};  // fro y2-y1
        sub_mod_p(y2,y1,y2_y1,10);  // (y2 - y1) (mod p) is done
        long long int x2_x1[10] = {0};  // for x2-x1
        sub_mod_p(x2,x1,x2_x1,10);  // (x2 - x1) (mod p) is done    
        
        // calculate inverse of x2-x1 in Z*_p using sqm_mod_p()
        long long int inv_x2_x1[10] = {0};  // (x2-x1)^(p-2) will be the inverse of (x2-x1) in Z*_p
        long long int p_2[10] = {268435453, 268435455, 268435455, 4095, 0, 0, 16777216, 0, 268435455, 15};  //p_2 is actually  p-2
        sqm_mod_p(x2_x1, p_2, inv_x2_x1, 10);   // calculating inverse of (x2-x1)
        long long int pr_in[20] = {0};  // pr_in  stores (y2_y1 * inv_x2_x1) i.e. (y2 - y1)/(x2 - x1)
        mult(y2_y1, inv_x2_x1, pr_in, 10);
        Barrett(pr_in, lamda, 20);

        // for(i = 0; i < 10; i++)
        //     printf("\n %lld",lamda[i]);
    }
    /*    VALUE OF LAMDA IS COMPUTED SUCCESSFULLY   */

    /* Calculation of x3 :  x3 = lamda^2 - x1 - x2 */
    long long int lamda_2[20] = {0};
    mult(lamda, lamda, lamda_2, 10);   // lamda^2  is done
    long long int add_x1_x2[10] = {0};  
    
    add(x1, x2, add_x1_x2, 10);   // x1+x2  is done

    // Extend the size of add_x1_x2 to 20 blocks
    long long int ext_add_x1_x2[20] = {0};
    for(i =0; i<10; i++)
        ext_add_x1_x2[i] = add_x1_x2[i];   // copied
    long long int ext_x3[20] = {0};
    sub(lamda_2, ext_add_x1_x2, ext_x3, 20); 
    Barrett(ext_x3, x3, 20);   // x3 is calculated successfully

    /* Calculation of y3 =  lamda(x1 - x3) - y1 */
    long long int x1_x3[10] = {0};  // for soring (x1 - x3)
    sub_mod_p(x1, x3, x1_x3, 10);   // (x1 - x3) is done
    long long int pr_lam[20] = {0};   // for storing lamda*(x1 - x3)
    mult(lamda, x1_x3, pr_lam, 10);   // lamda*(x1 - x3) is done
    long long int ext_y1[20] = {0};
    for(i=0; i<10; i++)
        ext_y1[i] = y1[i];    // y1 is extended to 20 block
    long long int ext_y3[20] = {0};   // for storing y3 in ext_y3 = pr_lam - ext_y1
    sub_mod_p(pr_lam, ext_y1, ext_y3, 20);  // y3 =  lamda(x1 - x3) - y1  is done, now reduce it
    Barrett(ext_y3, y3, 20);   // y3 is calculated

    // for(i = 0; i < 10; i++)
    //         printf("\n %lld     %lld",x3[i], y3[i]); 
}

// Doubling and Adding for computing K*(p_x , p_y) = (k_Px , k_Py)
void Doubling_Adding(long long int K[10], long long int p_x[10], long long int p_y[10], long long int k_Px[10], long long int k_Py[10] )   // K*(p_x , p_y) = (k_Px , k_Py)
{
    // initialization
    long long int k[10] = {0}, px[10] = {0}, py[10] = {0};
    long long int z1[10], z2[10];
    int i, r;
    for(i=0; i<10; i++)
    {
        k[i] = K[i];
        px[i] = p_x[i];
        py[i] = p_y[i];
        z1[i] = 0;
        z2[i] = 0;

        k_Px[i] = 0;
        k_Py[i] = 0;
    }
    int j, flag = 0, first_bit = 0, bit_value;
    long long int z3[10] = {0}, z4[10] = {0};  // will be used to store (temporary) sum of ECC_addition
    for(i = 9; i>=0; i--)
    {
        for(j = 27; j>=0; j--)
        {
            bit_value = ((k[i]>>j)&1);
            if((bit_value + flag)==0)
                continue;  // we will not perform any operation untill we get the first non-zero bit
            if((bit_value+first_bit)==1)
            {
                flag = 1;   // first non-zero bit is gotten, so flag is set to 1, also in later case if we get bit_vlue = 0 then this if(...) will not be executed but flag will remain set to 1
                first_bit = 3;  // since we want only once
            }
            if(((bit_value + flag)>0) && ((bit_value + flag)<3))  // if bit_value+flag = 1 or 2, that means first time getting non-zero bit
            {
                flag = 6;   // since we want to execute this only at the first time of getting non-zero bit
                for(r=0; r<10; r++)
                {
                    z1[r] = px[r];
                    z2[r] = py[r];
                }
                continue;  // since we want this only once, at the first timeof getting non-zero bit
            }
            if((bit_value + flag)>5)    // this means this bit is gotten after first non-zero bit
            {
                ECC_addition(z1,z2,z1,z2,z3,z4);  // (z3,z4) = 2*(z1,z2)  is done
                for(int l = 0; l<10; l++)
                {
                    z1[l] = z3[l];
                    z2[l] = z4[l];   // new(z1,z2) = 2 * old(z1,z2) 
                    z3[l] = 0;
                    z4[l] = 0;
                }
                if(bit_value)
                {
                    ECC_addition(z1,z2,px, py,z3,z4);  // if  bit_vlaue = 1 then  (z3,z4) = (z1,z2) + (px,py)
                    for(int l = 0; l<10; l++)
                    {
                        z1[l] = z3[l];
                        z2[l] = z4[l];   // new(z1,z2) =  old(z1,z2) + (px,py) 
                        z3[l] = 0;
                        z4[l] = 0;
                    }
                }
            }

        }
    }
    for(i=0;i<10;i++)
    {
        k_Px[i] = z1[i];
        k_Py[i] = z2[i];
    }
}