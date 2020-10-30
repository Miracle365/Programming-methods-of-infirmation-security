#include <stdio.h>        // IO
#include <openssl/bn.h>   // BIGNUM lib from openssl
#include "Parameters.h"    // Lab parameters
#include "Point.h"     // Point functions
#include "Mont_Curve.h"  // Elliptic curve functions

int main ()
{
    printf("\n\e[36mOpenSSL (libcrypto) Montgomery elliptic curve lab\e[0m\n\n");

    BIGNUM *power = BN_new(),    						// power of point
    *test_x = BN_new(),                         // var to test algorythm
    *test_y = BN_new(),
    *max_rand = NULL,     						// range of power randint
    *one = NULL;          						// 1

    struct Point point  = {NULL, NULL, NULL};
    struct Point point1 = {NULL, NULL, NULL};
    struct Point point2 = {NULL, NULL, NULL};
    struct Point point3 = {NULL, NULL, NULL};

    struct Mont_curv curve = {NULL, NULL, NULL, NULL, NULL};

    BN_dec2bn(&one, "1");               				// make some "one" number
    Mont_Curv_Constr(&curve);   				// initialize elliptic curve with parameters
    BN_dec2bn(&max_rand, MAX_RAND);  				// initialize random range


    point_constr(&point);               				// initialize id-tc26-gost-3410-2012-256-paramSetA point
    point_constr(&point1);                              // Extra point for test 4(four Ps)
    point_constr(&point2);
    point_constr(&point3);
    printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


    BN_rand_range(power, max_rand);     				// make cryptorandom integer from 0 to max_rand -- power
    //BN_dec2bn(&power, "18791202894425416709884867221754140282839356181413286370677768873537813477604");
    printf("[INFO] Random power is %s\n", BN_bn2dec(power));


    Mont_Curv_Ladder(&curve, &point, power);
    printf("\e[32m[ANSWER] X ** %s =\n         = %s\e[0m\n", BN_bn2dec(power), BN_bn2dec(point.x));
    printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


    printf("\e[32m[CHECK 1] Point is ");
    if (Mont_Curv_Point_Check(&curve, &point) != 1)
        printf("NOT ");
    printf("on curve\e[0m\n");


    printf("\e[32m[CHECK 2] Point with Q power is ");
    Mont_Curv_Ladder(&curve, &point, curve.q);
    if (!BN_is_one(point.x))
        printf("NOT ");
    printf("neutral\e[0m\n");
    printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));


    BN_copy(test_x, point.x);
    BN_copy(test_y, point.y);
    BN_add(power, curve.q, one);
    Mont_Curv_Ladder(&curve, &point, power);
    printf("\e[32m[CHECK 3.1] Point with Q+1 power is ");
    if ((BN_cmp(point.x, test_x) != 0) )
        printf("NOT ");
    printf("the same!\e[0m\n");


    BIGNUM *two = NULL;
    BN_dec2bn(&two, "2");
    BN_sub(power, curve.q, two);

    Mont_Curv_Ladder(&curve, &point, power);

    printf("\e[32m[CHECK 3.2] Point with Q-1 power is ");
    if ((BN_cmp(point.x, test_x) != 0))
        printf("NOT ");
    printf("the same!\e[0m\n");

    printf("[INFO] X = %s,\n       Z = %s\n", BN_bn2dec(point.x), BN_bn2dec(point.z));



    //! Test 4
    printf("\e[32m[CHECK 4] ");
    BIGNUM *k1 = BN_new(),              // Initialize k1, k2
            *k2 = BN_new(),
            *res = BN_new();
    BN_rand_range(k1, max_rand);        //Generate k1 k2
    BN_rand_range(k2, max_rand);

    if(BN_cmp(k2, k1) == 1)
        if(BN_cmp(point2.x, point3.x))
            BN_swap(k1, k2);

    Mont_Curv_Ladder(&curve, &point1, k1);
    Mont_Curv_Ladder(&curve, &point2, k2);
    BN_add(res, k1, k2);
    Mont_Curv_Ladder(&curve, &point3, res);


    point_add(&point2, &point1, &point, curve.p);

    if(!BN_cmp(point2.x, point3.x))
        printf("Test passed!");
    else
        printf("Test failed!");





    //! Clear
    BN_free(k1); BN_free(k2); BN_free(res); BN_free(one); BN_free(two);
    point_clear(&point);
    point_clear(&point1);
    point_clear(&point2);
    point_clear(&point3);
    BN_free(test_x); BN_free(test_y);
    BN_free(power);
    Mont_Curv_Clear(&curve);


}
