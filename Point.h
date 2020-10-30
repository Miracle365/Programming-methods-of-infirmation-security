#include <openssl/bn.h>

struct Point {
    BIGNUM *x, *z, *y;
};

void point_neutral(struct Point *point);
void point_constr(struct mi_point *point);
void point_add(struct Point *q, const struct Point *r, const struct Point *p1, const BIGNUM *p);
void point_double(struct Point *point, const BIGNUM *c, const BIGNUM *p);
void point_clear(struct Point *point);

