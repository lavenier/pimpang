
#define _GNU_SOURCE
#include <sched.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/sysinfo.h>
#include <math.h>
#include <immintrin.h>
#include <omp.h>

#include <semaphore.h>

#include <dpu.h>
#include <assert.h>
#include <dpu_log.h>

#include "constant.h"

