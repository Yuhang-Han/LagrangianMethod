/* cavity.c */

#include "conditions_api.h"

static State cavity_initial_condition(Point x)
{
    return (State){
        .rho = 1.0,
        .u   = 0.0,
        .v   = 0.0,
        .p   = 1.0
    };
}

static State cavity_boundary_condition(
    Point x,
    double time,
    State interior,
    BoundaryFace face)
{
    State ghost = interior;

    if (face == FACE_TOP) {
        /* moving lid */
        ghost.u = 1.0;
        ghost.v = 0.0;
    } else {
        /* stationary wall */
        ghost.u = 0.0;
        ghost.v = 0.0;
    }

    return ghost;
}

static const Conditions conditions = {
	.api_version		= 1,
	.case_name          = "Lid-driven cavity",
    .initial_condition  = cavity_initial_condition,
    .boundary_condition = cavity_boundary_condition
};
