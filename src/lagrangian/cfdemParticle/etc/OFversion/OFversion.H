//#define version24Dev // use this for cfdemSolverTwoPhaseEuler
#define version24  // use this for all other solvers
//#define version23
//#define version22
//#define version21
//#define version16ext
//#define version15

// features of 2.4Dev work also in Dev
#if defined(versionDev)
    #define version24Dev
#endif

// basically use 24x settings + some dev features (e.g. new turbulence model structure)
#if defined(version24Dev)
    #define version24
#endif

// features of 2.4 work also in 2.3
#if defined(version24)
    #define version23
#endif

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
