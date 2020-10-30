#include <stdio.h>
#include <openssl/bn.h>

struct Mont_curv {
    BIGNUM *A, *B, *C, *p, *q;
};

void Mont_Curv_Constr(struct Mont_curv *curve);
void Mont_Curv_Ladder(const struct Mont_curv *curve, struct Point *point, const BIGNUM *power);
void Mont_Curv_Clear(struct Mont_curv *curve);
int Mont_Curv_Point_Check(const struct Mont_curv *curve, const struct Point *point);
