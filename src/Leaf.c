#include "beps.h"

Leaf leaf_mul(Leaf x, Leaf y) {
    Leaf ans;
    
    ans.o_sunlit = x.o_sunlit * y.o_sunlit;
    ans.o_shaded = x.o_shaded * y.o_shaded;
    ans.u_sunlit = x.u_sunlit * y.u_sunlit;
    ans.u_shaded = x.u_shaded * y.u_shaded;
    return ans;
}

double leaf_sum(Leaf x) {
    return (x.o_sunlit + x.o_shaded + x.u_sunlit + x.u_shaded);
}

double clamp(double x, double low, double high) {
    return max(min(x, high), low);
}

void init_leaf_struct(Leaf* x, Leaf replacement) {
    x->o_sunlit = replacement.o_sunlit;
    x->o_shaded = replacement.o_shaded;
    x->u_sunlit = replacement.u_sunlit;
    x->u_shaded = replacement.u_shaded;
}

void init_leaf_dbl(Leaf* x, double replacement) {
    x->o_sunlit = replacement;
    x->o_shaded = replacement;
    x->u_sunlit = replacement;
    x->u_shaded = replacement;
}

void init_leaf_dbl2(Leaf* x, double overstory, double understory) {
    x->o_sunlit = overstory;
    x->o_shaded = overstory;
    x->u_sunlit = understory;
    x->u_shaded = understory;
}
