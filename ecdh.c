/* Elliptic Curve Key Exchange is implemented here. I used 2^28 base and NIST suggested 256-bit ECC prime p-256 = 2^256 - 2^224 + 2^192 + 2^96 - 1. All the required functions are written in functions.c file

Implemented by SOHAN DAS, an M.Tech(CrS) student of Indian Statistical Institute, KOLKATA
*/

#include <stdio.h>
#include "ECC.h"

int main()
{
    // P is a point on the EC, here we choose P as generator of the the EC group
    long long int px[10] = {144229014, 169055325, 187932916, 131601118, 15890179, 240532036, 133741798, 236110884, 186110450, 6};
    long long int py[10] = {129978869, 191104643, 22990539, 53835443, 236334030, 78102777, 193914859, 266446841, 266552034, 4};
    // these above two remain unchanged in this project

    // Let's consider Alice choose her private key  a ,a 256-bit random number from Z_p
    long long int a[10] = {20190170,
                           140767076,
                           152046517,
                           163557714,
                           89552916,
                           1848764,
                           78838121,
                           137173203,
                           3663413,
                           0};
    // Let's consider Bob choose her private key  b ,a 256-bit random number from Z_p
    long long int b[10] = {246941407,
                           92511354,
                           119248285,
                           186369831,
                           262066571,
                           191977025,
                           185213022,
                           100583089,
                           2104427,
                           0};
    // Alice's will make public a*(px, py), so let's compute that first
    long long int aPx[10] = {0};
    long long int aPy[10] = {0};
    Doubling_Adding(a, px, py, aPx, aPy); // a*(px, py) is computed
    // Bob will take this Alice's a*(px, py) and computes b*(a*P), P=(px,py)
    long long int baPx[10] = {0};
    long long int baPy[10] = {0};
    Doubling_Adding(b, aPx, aPy, baPx, baPy);

    // Now similarly Bob will make public b*(px,py), so let'c compute that
    long long int bPx[10] = {0}, bPy[10] = {10};
    Doubling_Adding(b, px, py, bPx, bPy); // b*(px, py) is copmputed
    // Alice will take this Bob's  b*(px, py) and computes a*(b*P)
    long long int abPx[10] = {0}, abPy[10] = {0};
    Doubling_Adding(a, bPx, bPy, abPx, abPy);

    // If (abPx, abPy) matches with the (baPx, baPy) then we get success means this product becomes a private key for both of them for further communication

    printf("\n Alice's computation and her key is :");
    for (int i = 0; i < 10; i++)
    {
        printf("\n abPx[%d] = %lld         abPy[%d] = %lld ", i, abPx[i], i, abPy[i]);
    }
    printf("\n\n Bob's computation and his key is :");
    for (int i = 0; i < 10; i++)
    {
        printf("\n baPx[%d] = %lld         baPy[%d] = %lld ", i, baPx[i], i, baPy[i]);
    }
    int count = 0;
    for (int i = 0; i < 10; i++)
        if (abPx[i] == baPx[i])
            count++;
    if (count == 10)
        printf("\n Shared key generatd \n\n");

    return 0;
}