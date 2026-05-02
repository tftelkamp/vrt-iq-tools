#include <stdio.h>

/* Fallback values in case the defines are not provided */
#ifndef GIT_BRANCH
#define GIT_BRANCH "unknown"
#endif

#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

#ifndef GIT_DATE
#define GIT_DATE "unknown"
#endif

int main(void)
{
    printf("vrt-iq-tools build information\n");
    printf("-----------------\n");
    printf("Branch:      %s\n", GIT_BRANCH);
    printf("Last Commit: %s\n", GIT_COMMIT);
    printf("Commit Date: %s\n", GIT_DATE);

    return 0;
}