#include "utils.h"

/**
 * @brief Function to return the time in milliseconds.
 *
 * @return double time.
 */
double timestamp(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double)(tp.tv_sec * 1000.0 + tp.tv_usec / 1000.0));
}

void outputPrinter() {
}
