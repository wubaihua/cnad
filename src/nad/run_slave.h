#include <complex.h>
#include <stdint.h>
#include "def.h"
#include "constant.h"
#include "gmath.h"
#include <stdio.h>
#include <math.h>
#include "cJSON.h"
#include "msmodel.h"
#include "msmodelio.h"
#include <stdbool.h>
#include "def_host.h"
#ifdef sunway
    #include <slave.h>
    #include <athread.h>
#endif

// extern SLAVE_FUN(init_slave)();


#ifdef sunway
extern SLAVE_FUN(dynamics_slave)(struct set_host *seth);
#elif defined(x86)
void dynamics_slave(struct set_host *seth);
#endif



// extern SLAVE_FUN(data_transport)(int id);



// extern SLAVE_FUN(free_slave)();
