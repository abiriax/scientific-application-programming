#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     1

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  USER_PAR_00             0

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1
#define  UNIT_LENGTH             1
#define  UNIT_VELOCITY           1

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    NO
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             DEFAULT
