#include "conditions_api.h"

#include <dlfcn.h>
#include <stdlib.h>
#include <stdio.h>


struct ConditionModule {
    void *handle;
    const Conditions *conditions;
};


static char error_message[256];


ConditionModule *
condition_module_load(const char *filename)
{
    ConditionModule *module = malloc(sizeof(*module));

    if (module == NULL) {
        snprintf(error_message,
                 sizeof(error_message),
                 "Cannot allocate ConditionModule");
        return NULL;
    }


    module->handle =
        dlopen(filename, RTLD_NOW | RTLD_LOCAL);

    if (module->handle == NULL) {
        snprintf(error_message,
                 sizeof(error_message),
                 "%s", dlerror());

        free(module);
        return NULL;
    }


    dlerror();  /* clear previous error */

    module->conditions =
        (const Conditions *)
        dlsym(module->handle, "conditions");

    const char *error = dlerror();

    if (error != NULL) {
        snprintf(error_message,
                 sizeof(error_message),
                 "%s", error);

        dlclose(module->handle);
        free(module);
        return NULL;
    }


    if (module->conditions->api_version
            != CFD_CASE_API_VERSION) {

        snprintf(error_message,
                 sizeof(error_message),
                 "Case API version mismatch: "
                 "solver=%d, case=%zu",
                 CFD_CASE_API_VERSION,
                 module->conditions->api_version);

        dlclose(module->handle);
        free(module);
        return NULL;
    }


    return module;
}


const Conditions *
condition_module_conditions(
    const ConditionModule *module)
{
    return module->conditions;
}


void
condition_module_free(ConditionModule *module)
{
    if (module == NULL)
        return;

    if (module->handle != NULL)
        dlclose(module->handle);

    free(module);
}


const char *
condition_module_error(void)
{
    return error_message;
}
