#include "snova.h"

#if SNOVA_q == 16
#if SNOVA_r != SNOVA_l
#include "snova_rect_16.c"
#elif SNOVA_l == 4
#include "snova_opt_16.c"
#else
#include "snova_avx2_16.c"
#endif
#elif defined(SNOVA_r) && (SNOVA_r != SNOVA_l)
#include "snova_rect_q.c"
#elif SNOVA_l == 2
#include "snova_opt_q2.c"
#else
#include "snova_opt_q.c"
#endif
