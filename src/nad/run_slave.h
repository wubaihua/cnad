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

#include <slave.h>
#include <athread.h>


// extern SLAVE_FUN(init_slave)();



extern SLAVE_FUN(dynamics_slave)(struct set_host *seth);



// extern SLAVE_FUN(data_transport)(int id);



// extern SLAVE_FUN(free_slave)();
