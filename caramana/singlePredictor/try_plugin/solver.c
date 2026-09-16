#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

#include "conditions_api.h"

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s case.so\n", argv[0]);
        return 1;
    }

	ConditionModule *module = condition_module_load(argv[1]);
    if (module == NULL) {
        fprintf(stderr,
                "Cannot load conditions: %s\n",
                condition_module_error());

        return EXIT_FAILURE;
    }

    const Conditions *cond = condition_module_conditions(module);

    printf("Running case: %s\n", cond->case_name);
    printf("API version: %zu\n", cond->api_version);


    /* solver initialization */

    Point x = {0.25, 0.0, 0.0};
    State q = cond->initial_condition(x);

    /* main solver */
    /*
        ...
        q = cond->boundary_condition(...);
        ...
    */

    

    condition_module_free(module);
}
