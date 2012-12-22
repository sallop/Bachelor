     #include <fpu_control.h>
     static void __attribute__ ((constructor))
     trapfpe ()
     {
       __setfpucw (_FPU_DEFAULT &
                   ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM));
     }

