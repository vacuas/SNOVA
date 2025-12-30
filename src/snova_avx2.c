#include "snova.h"

#if SNOVA_q == 16
#if SNOVA_l == 4
#include "snova_opt_16.c"
#else
#include "snova_avx2_16.c"
#endif
#elif defined(SNOVA_r) && (SNOVA_r != SNOVA_l)
#include "snova_avx2_rect.c"
#else
#if __GNUC__ >= 15
#include "snova_opt_q.c"
#else
#include "snova_avx2_q.c"
#endif
#endif
