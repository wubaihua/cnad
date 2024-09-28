cmake_minimum_required(VERSION 3.10)
project(cnad)

# add_compile_options(-Wall -Wextra)

if(x86)
    add_compile_definitions(x86)
    set(CMAKE_C_COMPILER "mpicc")

    # add_compile_options(-g -O0)
    # 设置包含目录
    include_directories(
        # /home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/LAPACKE/include/
        src/mod/
        src/msmodel/
        src/nad/
    )

    # 使用 swgcc 编译的文件
    add_library(cJSON OBJECT src/mod/cJSON.c)
    set_target_properties(cJSON PROPERTIES COMPILE_FLAGS " ")

    add_library(constant OBJECT src/mod/constant.c)
    set_target_properties(constant PROPERTIES COMPILE_FLAGS " ")

    add_library(gmath OBJECT src/mod/gmath.c)
    set_target_properties(gmath PROPERTIES COMPILE_FLAGS " ")

    add_library(iofun OBJECT src/mod/iofun.c)
    set_target_properties(iofun PROPERTIES COMPILE_FLAGS " ")

    add_library(iofun_slave OBJECT src/mod/iofun_slave.c)
    set_target_properties(iofun_slave PROPERTIES COMPILE_FLAGS " ")

    add_library(msmodelio OBJECT src/msmodel/msmodelio.c)
    set_target_properties(msmodelio PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(sbm OBJECT src/msmodel/sbm.c)
    set_target_properties(sbm PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(morse3 OBJECT src/msmodel/morse3.c)
    set_target_properties(morse3 PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(sem OBJECT src/msmodel/sem.c)
    set_target_properties(sem PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(fmo OBJECT src/msmodel/fmo.c)
    set_target_properties(fmo PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(msmodel OBJECT src/msmodel/msmodel.c)
    set_target_properties(msmodel PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(def OBJECT src/nad/def.c)
    set_target_properties(def PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/ -lm -lm_slave -lgfortran_slave")

    add_library(run_slave OBJECT src/nad/run_slave.c)
    set_target_properties(run_slave PROPERTIES COMPILE_FLAGS "-I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(def_host OBJECT src/nad/def_host.c)
    set_target_properties(def_host PROPERTIES COMPILE_FLAGS "-I src/msmodel/ -I src/mod/")

    set(CMAKE_C_COMPILER "mpicc")
    # 使用 mpicc 编译的文件
    add_library(nad OBJECT src/nad/nad.c)
    set_target_properties(nad PROPERTIES COMPILE_FLAGS "-I src/msmodel/ -I src/mod/")

    # 添加可执行文件
    add_executable(cnad $<TARGET_OBJECTS:cJSON> $<TARGET_OBJECTS:constant> $<TARGET_OBJECTS:gmath> $<TARGET_OBJECTS:iofun> $<TARGET_OBJECTS:iofun_slave> 
                        $<TARGET_OBJECTS:msmodelio> $<TARGET_OBJECTS:sbm> $<TARGET_OBJECTS:morse3> $<TARGET_OBJECTS:sem> $<TARGET_OBJECTS:fmo> $<TARGET_OBJECTS:msmodel> 
                        $<TARGET_OBJECTS:def> $<TARGET_OBJECTS:run_slave> $<TARGET_OBJECTS:def_host> $<TARGET_OBJECTS:nad>)

    # 链接库
    target_link_libraries(cnad
            lapacke 
            m
            gfortran
    )


    
elseif(sunway)
    add_compile_definitions(sunway)
    set(CMAKE_C_COMPILER "swgcc")


    # 设置包含目录
    include_directories(
        /home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/LAPACKE/include/
        src/mod/
        src/msmodel/
        src/nad/
    )

    # 使用 swgcc 编译的文件
    add_library(cJSON OBJECT src/mod/cJSON.c)
    set_target_properties(cJSON PROPERTIES COMPILE_FLAGS "-mhost")

    add_library(constant OBJECT src/mod/constant.c)
    set_target_properties(constant PROPERTIES COMPILE_FLAGS "-mhost")

    add_library(gmath OBJECT src/mod/gmath.c)
    set_target_properties(gmath PROPERTIES COMPILE_FLAGS "-mslave")

    add_library(iofun OBJECT src/mod/iofun.c)
    set_target_properties(iofun PROPERTIES COMPILE_FLAGS "-mhost")

    add_library(iofun_slave OBJECT src/mod/iofun_slave.c)
    set_target_properties(iofun_slave PROPERTIES COMPILE_FLAGS "-mslave")

    add_library(msmodelio OBJECT src/msmodel/msmodelio.c)
    set_target_properties(msmodelio PROPERTIES COMPILE_FLAGS "-mhost -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(sbm OBJECT src/msmodel/sbm.c)
    set_target_properties(sbm PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(morse3 OBJECT src/msmodel/morse3.c)
    set_target_properties(morse3 PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(sem OBJECT src/msmodel/sem.c)
    set_target_properties(sem PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(fmo OBJECT src/msmodel/fmo.c)
    set_target_properties(fmo PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(msmodel OBJECT src/msmodel/msmodel.c)
    set_target_properties(msmodel PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(def OBJECT src/nad/def.c)
    set_target_properties(def PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/ -lm -lm_slave -lgfortran_slave")

    add_library(run_slave OBJECT src/nad/run_slave.c)
    set_target_properties(run_slave PROPERTIES COMPILE_FLAGS "-mslave -I src/nad/ -I src/msmodel/ -I src/mod/")

    add_library(def_host OBJECT src/nad/def_host.c)
    set_target_properties(def_host PROPERTIES COMPILE_FLAGS "-mhost -I src/msmodel/ -I src/mod/")

    set(CMAKE_C_COMPILER "mpicc")
    # 使用 mpicc 编译的文件
    add_library(nad OBJECT src/nad/nad.c)
    set_target_properties(nad PROPERTIES COMPILE_FLAGS "-mhost -I src/msmodel/ -I src/mod/")

    # 添加可执行文件
    add_executable(cnad $<TARGET_OBJECTS:cJSON> $<TARGET_OBJECTS:constant> $<TARGET_OBJECTS:gmath> $<TARGET_OBJECTS:iofun> $<TARGET_OBJECTS:iofun_slave> 
                        $<TARGET_OBJECTS:msmodelio> $<TARGET_OBJECTS:sbm> $<TARGET_OBJECTS:morse3> $<TARGET_OBJECTS:sem> $<TARGET_OBJECTS:fmo> $<TARGET_OBJECTS:msmodel> 
                        $<TARGET_OBJECTS:def> $<TARGET_OBJECTS:run_slave> $<TARGET_OBJECTS:def_host> $<TARGET_OBJECTS:nad>)

    # 链接库
    target_link_libraries(cnad -mhybrid
        /home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/liblapacked.a
        /home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/libcblasd.a
        /home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/liblapackd.a
        /home/export/online1/mdt00/shisuan/swustcfd/liugroup/bqli/generalTests/05-LAPACK/build/device/librefblasd.a
        m
        m_slave
        gfortran_slave
        ieee
    )




endif()
