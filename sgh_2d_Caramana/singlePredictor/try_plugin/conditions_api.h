#ifndef CONDITIONS_API_H
#define CONDITIONS_API_H

#include <stddef.h>

#define CFD_CASE_API_VERSION 1

typedef struct {
    double rho;
    double u;
    double v;
    double p;
} State;

typedef struct {
    double x;
    double y;
    double z;
} Point;

typedef enum {
    FACE_LEFT,
    FACE_RIGHT,
    FACE_BOTTOM,
    FACE_TOP
} BoundaryFace;

typedef struct {
    size_t api_version;

    const char *case_name;

	// function pointers
    State (*initial_condition)(Point x);

    State (*boundary_condition)(
        Point x,
        double time,
        State interior,
        BoundaryFace face
	);
} Conditions;


/* opaque loader object */
typedef struct ConditionModule ConditionModule;


/* Load a condition plugin.
   Returns NULL on failure. */
ConditionModule *condition_module_load(const char *filename);


/* Get its Conditions interface. */
const Conditions *
condition_module_conditions(const ConditionModule *module);


/* Unload it. */
void condition_module_free(ConditionModule *module);


/* Error from the most recent load operation */
const char *condition_module_error(void);


#endif
