#include "pch.h"
#include "Logger.h"
void LOG(const char *fmt, ...) {
#if defined(DEBUG) || defined(INFOLOG)
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
	printf("\n");
    va_end(args);
#endif
	return;
}
