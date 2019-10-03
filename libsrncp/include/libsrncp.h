#include "helper.h"

#ifdef __cplusplus
extern "C" {
#endif

EXPORTED_FUNCTION void Unwrap3D(int h, int w, int d, double* phase, int* mask, int mask_largest_unwrapped_group);
EXPORTED_FUNCTION void Unwrap4D(int h, int w, int d, int s, double* phase, int* mask, int mask_largest_unwrapped_group);

#ifdef __cplusplus
};
#endif    

