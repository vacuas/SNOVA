#include "snova.h"

#if SNOVA_q == 16
#include "snova_opt_16.c"
#elif defined(SNOVA_r) && (SNOVA_r != SNOVA_l)
#include "snova_rect_q.c"
#else
#include "snova_opt_q.c"
#endif
