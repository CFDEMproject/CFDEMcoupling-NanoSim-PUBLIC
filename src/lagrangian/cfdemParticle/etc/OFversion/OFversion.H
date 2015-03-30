#define version23

// features of 2.1 work also in 2.3
#if defined(version23)
    #define version21
    #define version221
#endif

// features of 2.1 work also in 2.2
#if defined(version22)
    #define version21
    #define version221
#endif
