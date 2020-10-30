#include <stdio.h>
#include <openssl/bn.h>
#include "Point.h"
#include "Parameters.h"

struct Mont_curv {
    BIGNUM *A, *B, *C, *p, *q;
};

/*
 * Initialize Montgomery curve
 * B*y^2=x^3+A*x^2+x mod p
 */
void Mont_Curv_Constr(struct Mont_curv *curve){

    BIGNUM *e     = NULL,
            *d     = BN_new(),
            *sum   = BN_new(),
            *diff  = BN_new(),
            *a     = NULL;

    BN_CTX *tmp = BN_CTX_new();

    BN_hex2bn(&curve->q, Q_PARAMETER);                      // init q parameter
    BN_hex2bn(&curve->p, P_PARAMETR);                        // init p parameter
    BN_hex2bn(&e, E_PARAMETER);                             // init e parameter
    BN_hex2bn(&d, D_PARAMETER);                             // init d parameter

    // A = 2 * (e + d) / (e - d)
    BN_add(sum, e, d);
    BN_sub(diff, e, d);
    BN_dec2bn(&curve->A, "2");
    BN_mul(curve->A, curve->A, sum, tmp);
    BN_mod_inverse(diff, diff, curve->p, tmp);
    BN_mod_mul(curve->A, curve->A, diff, curve->p, tmp);

    // B = 4 / (e - d)
    BN_dec2bn(&curve->B, "4");
    BN_mod_mul(curve->B, curve->B, diff, curve->p, tmp);

    // C = (A - 2) / 4
    BN_dec2bn(&curve->C, "2");
    BN_sub(curve->C, curve->A, curve->C);
    BN_dec2bn(&a, "4");
    BN_mod_inverse(a, a, curve->p, tmp);
    BN_mod_mul(curve->C, curve->C, a, curve->p, tmp);

    BN_free(a); BN_free(d); BN_free(e); BN_free(sum); BN_free(diff);
    BN_CTX_free(tmp);
}

void Mont_Curv_Ladder(const struct Mont_curv *curve, struct Point *point, const BIGNUM *power){

    BN_CTX *tmp = BN_CTX_new();                          // tmp element
    int size = BN_num_bits(power);                       // bit count of power

    struct Point r = {BN_new(), BN_new(), BN_new()},            // tmp points
    q = {NULL, NULL, NULL};
    point_neutral(&q);                                // init neutral point

    BN_copy(r.x, point->x);                              // copy point to temp point
    BN_copy(r.z, point->z);

    for (int i = size - 1; i >= 0; i--) {                // for every bit
        if (BN_is_bit_set(power, i)) {                   // if bit is 1
            point_add(&q, &r, point, curve->p);       // q = r (+) point mod p
            point_double(&r, curve->C, curve->p);     // r = r (*) 2 mod p
        } else {
            point_add(&r, &q, point, curve->p);       // r = q (+) point mod p
            point_double(&q, curve->C, curve->p);     // q = q (*) 2 mod p
        }
    }

    BN_copy(point->x, q.x);                              // copy q to point
    BN_copy(point->z, q.z);

    // Afin coords transformation
    // x_res = x * z^(-1) mod p
    BIGNUM *invert = BN_new();
    if (!BN_is_zero(point->z)){
        BN_mod_inverse(invert, point->z, curve->p, tmp);
        BN_mod_mul(point->x, point->x, invert, curve->p, tmp);
        BN_mod_mul(point->y, point->y, invert, curve->p, tmp);
        BN_mod_mul(point->z, point->z, invert, curve->p, tmp);
    } else {
        BN_mod_inverse(invert, point->x, curve->p, tmp);
        BN_mod_mul(point->x, point->x, invert, curve->p, tmp);
        BN_mod_mul(point->y, point->y, invert, curve->p, tmp);
    }

    BN_free(invert);
    point_clear(&r);
    point_clear(&q);
    BN_CTX_free(tmp);
}


/*
 * Clears memory
 */
void Mont_Curv_Clear(struct Mont_curv *curve){
    BN_free(curve->q);
    BN_free(curve->p);
    BN_free(curve->A);
    BN_free(curve->B);
    BN_free(curve->C);
}


/*
 * Checks if point is on curve
 */
int Mont_Curv_Point_Check(const struct Mont_curv *curve, const struct Point *point){

    BIGNUM *invert = BN_new(),
            *n = BN_new(),
            *pow_2 = NULL,
            *yY = BN_new(),
            *mul_res = BN_new();
    BN_CTX *tmp = BN_CTX_new();

    BN_dec2bn(&pow_2, "2");

    BN_exp(n, point->x, pow_2, tmp);                 // n = x^2
    BN_mul(mul_res, curve->A, point->x, tmp);        // mul_res = curve.a * x
    BN_add(n, n, mul_res);                           // n = n + curve.a * x
    BN_mul(n, n, point->x, tmp);                     // n = (n + curve.a * x) * x
    BN_add(n, n, point->x);                          // n = (n + curve.a * x) * x + x
    BN_mod_inverse(invert, curve->B, curve->p, tmp); // find invert B
    BN_mul(n, n, invert, tmp);                       // n = ((n + curve.a * x) * x + x) * B^-1

    BN_mul(yY, point->y, point->y, tmp);
    return BN_cmp(yY, n);
}
