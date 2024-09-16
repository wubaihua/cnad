# 在神威海洋之光平台上进行主核编译操作，使用课题组编译的主核LAPACK/BLAS
# 由于使用与USTC的公用账号，请勿从环境变量中设置链接库！

# 主核版本的静态库路径
# 提供netlib的LAPACK标准版本 (in F77)，有C和C++接口。

LAPACK_HOST_LIB=/home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/LAPACKE/include/
LAPACK_HOST_LIB2=/home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/
modpath=src/mod/
modelpath=src/msmodel/

# 编译。神威平台上的mpifort后端使用swgfortran。

# mpifort demo.f90 -L$LAPACK_HOST_LIB -llapack -lrefblas -lm -mieee -o demo

# swgcc -mslave -c athread_get_id.c 
# # swgfortran -mslave -c lapack.f -L/home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/host/ -llapack -lrefblas -lm -mieee
# swgfortran -mslave -c var.f90 

# swgfortran -mslave -c run.f90 
# mpifort -mhost -c demo.f90  
# # mpifort -O3 -mhybrid demo.o run.o athread_get_id.o var.o -L/usr/sw/yyzlib/xMath -lxMath 
# mpifort -O3 -mhybrid demo.o run.o athread_get_id.o var.o  -L$LAPACK_HOST_LIB  -lcblasd -llapacked -llapackd -lrefblasd -lm -lm_slave -mieee -o test


swgcc -mhost -c src/mod/cJSON.c -o cJSON.o
swgcc -mhost -c src/mod/constant.c -o constant.o
swgcc -mslave -c src/mod/gmath.c -o gmath.o -I$LAPACK_HOST_LIB -I$modpath
swgcc -mhost -c src/mod/iofun.c -o iofun.o -I$LAPACK_HOST_LIB -I$modpath
swgcc -mslave -c src/mod/iofun_slave.c -o iofun_slave.o -I$LAPACK_HOST_LIB -I$modpath

swgcc -mhost -c src/msmodel/msmodelio.c -o msmodelio.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath
swgcc -mslave -c src/msmodel/sbm.c -o sbm.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath

swgcc -mslave -c src/msmodel/msmodel.c -o msmodel.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath

swgcc -mslave -c src/nad/def.c -o def.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath -lm -lm_slave -lgfortran_slave
swgcc -mslave -c src/nad/run_slave.c -o run_slave.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath
swgcc -mhost -c src/nad/def_host.c -o def_host.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath
mpicc -mhost -c src/nad/nad.c -o nad.o -I$LAPACK_HOST_LIB -I$modpath -I$modelpath

mpicc -O3 -mhybrid  *.o  -L$LAPACK_HOST_LIB2  -llapacked -lcblasd -llapackd -lrefblasd -lm -lm_slave -lgfortran_slave -mieee -o cnad

