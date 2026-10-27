/* sod.c */

#include "conditions_api.h"

static State sod_initial_condition(Point x)
{
    if (x.x < 0.5) {
        return (State){
            .rho = 1.0,
            .u   = 0.0,
            .v   = 0.0,
            .p   = 1.0
        };
    }

    return (State){
        .rho = 0.125,
        .u   = 0.0,
        .v   = 0.0,
        .p   = 0.1
    };
}

static State sod_boundary_condition(
    Point x,
    double time,
    State interior,
    BoundaryFace face)
{
    return interior;
}

const Conditions conditions = {
	.api_version = 1,
	.case_name          = "Sod shock tube",
    .initial_condition  = sod_initial_condition,
    .boundary_condition = sod_boundary_condition
};
