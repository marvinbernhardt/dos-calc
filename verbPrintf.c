#ifndef VERBPRINTF
#define VERBPRINTF

// print only if verbosity is true
void verbPrintf(bool verbosity, const char *format, ...)
{
    va_list args;

    if (!verbosity)
        return;

    va_start(args, format);
    vfprintf (stdout, format, args);
    va_end(args);
}

#endif
