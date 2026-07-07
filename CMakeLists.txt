cmake_minimum_required(VERSION 2.8.9)

################################################################################
#                               User Inputs
################################################################################

# Specify the compiler with the FC environment variable or
# -DCMAKE_Fortran_COMPILER=... when configuring CMake. Do not set
# CMAKE_Fortran_COMPILER_ID manually; CMake fills it after detecting the
# selected compiler, and the build logic below relies on the real value.

# Set the target host configuration.  This is a cache variable so collaborators
# can use `-Dhostname=...` without editing this file.
# Supported computers: derecho, rockfish, panoramix, njord, marcc,
# darwin-brew, cheyenne, stampede, default.
# "delta" is kept only as a backward-compatible alias for older commands that
# used the same NVHPC/CUDA/Cray-style GPU profile name.
set(hostname "derecho" CACHE STRING "Target host configuration")
set_property(CACHE hostname PROPERTY STRINGS rockfish panoramix njord marcc
    darwin-brew cheyenne stampede derecho delta default)

# Read and write binary files as big_endian or little_endian
# Options: DEFAULT, LITTLE, or BIG
set(WRITE_ENDIAN "DEFAULT" CACHE STRING "Output binary endian mode")
set(READ_ENDIAN "DEFAULT" CACHE STRING "Input binary endian mode")
set_property(CACHE WRITE_ENDIAN PROPERTY STRINGS DEFAULT LITTLE BIG)
set_property(CACHE READ_ENDIAN PROPERTY STRINGS DEFAULT LITTLE BIG)

# Build feature switches.
#
# These declarations define valid CMake cache variables and safe global
# defaults.  Case-specific physics should be selected in the case compile
# script, for example test-cases/<case>/compile_derecho.sh, not by editing
# these defaults.
option(USE_MPI "Enable parallelization support." ON)
option(USE_CPS "Enable concurrent precursor simulations." OFF)
option(USE_HIT "Enable homogeneous isotropic turbulence input." OFF)
option(USE_LVLSET "Enable level set." OFF)
option(USE_TURBINES "Enable actuator disk model without rotation." OFF)
option(USE_ATM "Enable actuator line/section model" OFF)
option(USE_LES_GPU "Enable optimized OpenACC explicit-residency LES core." OFF)
option(USE_CPU_BUILD "CPU baseline: nvfortran accepts the dormant cuf syntax but links no GPU libs (runs on GPU-less nodes)." OFF)
set(USE_GPU_AWARE_MPI "AUTO" CACHE STRING
    "GPU-aware MPI for USE_LES_GPU halos/tridag pipeline: AUTO = detect Cray GTL (Derecho/Perlmutter/HPE Cray EX); ON = force PPGPU_AWARE_MPI on machines whose MPI is CUDA-aware without an extra lib (e.g. OpenMPI+UCX); OFF = host-staged fallback paths")
set_property(CACHE USE_GPU_AWARE_MPI PROPERTY STRINGS AUTO ON OFF)
option(USE_DYN_TN "Dynamically update T, the timescale for Lagrangian averaging, using Taylor timescale." OFF)
option(USE_SAFETYMODE "Set extra safety measures." ON)
option(USE_CGNS "Enable CGNS data output." OFF)
option(USE_SCALARS "Enable scalar transport module." OFF)
option(USE_SCALARS_GPU "Enable optimized OpenACC explicit-residency scalar transport path. Requires USE_SCALARS and USE_LES_GPU." OFF)

if (USE_CPS AND NOT USE_MPI)
    message(FATAL_ERROR "USE_CPS requires USE_MPI=ON because concurrent precursor synchronization uses mpi_defs.")
endif ()
if (USE_CPU_BUILD AND USE_LES_GPU)
    message(FATAL_ERROR "USE_CPU_BUILD requires USE_LES_GPU=OFF; CPU baseline builds do not link GPU runtime libraries.")
endif ()
if (USE_SCALARS_GPU AND NOT USE_SCALARS)
    message(FATAL_ERROR "USE_SCALARS_GPU requires USE_SCALARS=ON.")
endif ()
if (USE_SCALARS_GPU AND NOT USE_LES_GPU)
    message(FATAL_ERROR "USE_SCALARS_GPU requires USE_LES_GPU=ON.")
endif ()
if (USE_GPU_AWARE_MPI STREQUAL "ON" AND NOT USE_MPI)
    message(FATAL_ERROR "USE_GPU_AWARE_MPI=ON requires USE_MPI=ON.")
endif ()
if (USE_GPU_AWARE_MPI STREQUAL "ON" AND NOT USE_LES_GPU)
    message(FATAL_ERROR "USE_GPU_AWARE_MPI=ON requires USE_LES_GPU=ON.")
endif ()

# Enable verbose mode
set(CMAKE_VERBOSE_MAKEFILE ON)

################################################################################
#                         Set up build (do not change)
################################################################################
# Set compiler versions for Cmake
#  Intel
#if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    #set(compiler_name "ifort")
#    set(compiler_name "ifx")
#  GNU
#elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
#    set(compiler_name "gfortran")
#endif ()
#find_program(WHICH NAMES which)
#exec_program(${WHICH} ARGS ${compiler_name} OUTPUT_VARIABLE compiler_path)
#set(CMAKE_Fortran_COMPILER ${compiler_path})

#set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

# Project details. LESGO is a Fortran codebase; keeping the project language
# Fortran-only avoids unrelated C/CXX compiler probes on HPC module stacks.
project(lesgo Fortran)
set(exec_name "lesgo")

# Keep compiler-generated Fortran module files out of the source tree and in a
# single build-local include directory.  This also makes parallel module
# dependencies less sensitive to stale `.mod` files from previous builds.
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
include_directories("${CMAKE_Fortran_MODULE_DIRECTORY}")

# Library links
if (${hostname} STREQUAL "marcc")
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} -mkl)
    else()
        link_libraries(fftw3)
    endif()
else()
    link_libraries(fftw3)
endif()

if (${hostname} STREQUAL "rockfish")
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        set(CMAKE_Fortran_COMPILER ${compiler_path})
        set(CMAKE_C_COMPILER gcc)
        set(CMAKE_CXX_COMPILER g++)
    else()
        link_libraries(fftw3)
    endif()
endif()
#============================================
# Compiler/library detection and FFTW setup.
if (${hostname} STREQUAL "derecho" OR ${hostname} STREQUAL "delta")
    message(STATUS "Configuring FFTW paths for Derecho/Cray-style GPU profile...")
    if(DEFINED ENV{FFTW_ROOT})
        include_directories($ENV{FFTW_ROOT}/include)
        link_directories($ENV{FFTW_ROOT}/lib)
        message(STATUS "Found FFTW at: $ENV{FFTW_ROOT}")
    elseif(DEFINED ENV{FFTW_DIR})
        include_directories($ENV{FFTW_DIR}/include)
        link_directories($ENV{FFTW_DIR}/lib)
        message(STATUS "Found FFTW at: $ENV{FFTW_DIR}")
    else()
        message(WARNING "FFTW environment variables not found; load cray-fftw before configuring.")
    endif()
    link_libraries(fftw3)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(STATUS "Using Intel Compiler, linking MKL...")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -xHost -ipo -qmkl")
    include_directories($ENV{MKLROOT}/include)
    include_directories($ENV{MKLROOT}/include/fftw)

# NVIDIA compiler fallback paths for non-Derecho local builds.
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    message(STATUS "Using NVIDIA HPC Compiler, linking FFTW3...")
    if(DEFINED ENV{FFTW_ROOT})
        include_directories($ENV{FFTW_ROOT}/include)
        link_directories($ENV{FFTW_ROOT}/lib)
        message(STATUS "Found FFTW at: $ENV{FFTW_ROOT}")
    elseif(DEFINED ENV{FFTW_DIR})
        include_directories($ENV{FFTW_DIR}/include)
        link_directories($ENV{FFTW_DIR}/lib)
        message(STATUS "Found FFTW at: $ENV{FFTW_DIR}")
    else()
        message(STATUS "FFTW_ROOT/FFTW_DIR not set; relying on default FFTW search paths.")
    endif()
    link_libraries(fftw3)

else()
    message(STATUS "Using GNU Compiler, linking FFTW3...")
    if(DEFINED ENV{FFTW_ROOT})
        include_directories($ENV{FFTW_ROOT}/include)
        link_directories($ENV{FFTW_ROOT}/lib)
        message(STATUS "Found FFTW at: $ENV{FFTW_ROOT}")
    elseif(DEFINED ENV{FFTW_DIR})
        include_directories($ENV{FFTW_DIR}/include)
        link_directories($ENV{FFTW_DIR}/lib)
        message(STATUS "Found FFTW at: $ENV{FFTW_DIR}")
    else()
        message(STATUS "FFTW_ROOT/FFTW_DIR not set; relying on default FFTW search paths.")
    endif()
    link_libraries(fftw3)
endif()

# Common include and link directories
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
if ( NOT(${hostname} STREQUAL "darwin-brew") )
    # include_directories($ENV{HOME}/include /usr/include)
    include_directories(/usr/include)
    link_directories($ENV{HOME}/lib /usr/lib64)
endif ()

# Include host-specific paths
#  panoramix
if (${hostname} STREQUAL "panoramix")
    if (USE_CGNS)
        include_directories(/opt/cgns-3.2.1-hdf5-mpi/include
        /opt/hdf5-1.8.12-mpi/include /usr/local/hdf5/include /usr/local/cgns/include
        /act/fftw3/mpich/gcc-4.7.2/include /act/hdf5-1.8.15/mpich/gcc-4.7.2/include
        /act/cgns-3.3.0/mpich/gcc-4.7.2/include)
        link_directories(/opt/hdf5-1.8.12-mpi/lib  /opt/cgns-3.2.1-hdf5-mpi/lib
        /usr/local/hdf5/lib  /usr/local/cgns/lib /act/fftw3/mpich/gcc-4.7.2/lib
        /act/hdf5-1.8.15/mpich/gcc-4.7.2/lib /act/cgns-3.3.0/mpich/gcc-4.7.2/lib)
    else (USE_CGNS)
        include_directories(/act/fftw3/mpich/gcc-8.3.0/include)
        link_directories(/act/fftw3/mpich/gcc-4.7.2/lib)
    endif ()
#  marcc
elseif (${hostname} STREQUAL "marcc")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/software/apps/compilers/intel/mkl/include/fftw/)
        link_directories(/software/apps/compilers/intel/mkl/lib/intel64/)
    # GNU -- site FFTW path; optional CGNS/HDF5 use user-local prefixes below.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        include_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/include/)
        link_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/lib/)
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()

#  rockfish
elseif (${hostname} STREQUAL "rockfish")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/fftw-3.3.8-wxxony3sfn4vctui2jomkwlzs6jfnvqr/include/)
        link_directories(/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/fftw-3.3.8-wxxony3sfn4vctui2jomkwlzs6jfnvqr/lib/)
    # GNU -- site FFTW path; optional CGNS/HDF5 use user-local prefixes below.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        include_directories(/data/apps/linux-centos8-cascadelake/gcc-9.3.0/fftw-3.3.8-uxkh4w5nrgeviv2ooaorhcuolbtbosyi/include/)
        link_directories(/data/apps/linux-centos8-cascadelake/gcc-9.3.0/fftw-3.3.8-uxkh4w5nrgeviv2ooaorhcuolbtbosyi/lib/)
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()

# Cheyenne
elseif (${hostname} STREQUAL "cheyenne")
    # Intel
    if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
        include_directories(/glade/u/apps/ch/opt/fftw3-mpi/3.3.6-pl2/mpt/2.15f/intel/18.0.1/include/)
        link_directories(/glade/u/apps/ch/opt/fftw3-mpi/3.3.6-pl2/mpt/2.15f/intel/18.0.1/lib/)
    # GNU -- site FFTW path; optional CGNS/HDF5 use user-local prefixes below.
    elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        include_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/include/)
        link_directories(/software/apps/fftw3/3.3.8/gcc/5.5.0/lib/)
        if (USE_CGNS)
            include_directories($ENV{HOME}/cgns/include/)
            link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
        endif ()
    endif ()
# stampede
elseif (${hostname} STREQUAL "stampede")
    # Intel
        include_directories(/opt/apps/intel18/impi18_0/fftw3/3.3.8/include /opt/intel/mkl/include/intel64/lp64)
        link_directories(/opt/apps/intel18/impi18_0/fftw3/3.3.8/lib /opt/intel/mkl/lib/intel64_lin)
#  njord
elseif (${hostname} STREQUAL "njord")
    # nothing yet, but can add if needed
#  all other systems
elseif (${hostname} STREQUAL "default")
    # include_directories($ENV{HOME}/codes/fftw3/include)
    # link_directories($ENV{HOME}/codes/fftw3/lib)
    # nothing yet, but can add if needed
    if (USE_CGNS)
        include_directories($ENV{HOME}/cgns/include/)
        link_directories($ENV{HOME}/hdf5/lib/ $ENV{HOME}/cgns/lib/)
    endif ()
endif ()

# Common sources are grouped by solver ownership.  Keep these groups aligned
# with docs/source_file_inventory.md so collaborators can map CMake membership
# to the source-tree index.
set(COMMON_DRIVER_SOURCES
    cfl_util.f90 clocks.f90 finalize.f90 functions.f90 grid.f90 initial.f90
    initialize.f90 input_util.f90 init_random_seed.f90 main.f90 messages.f90
    param.f90 param_output.f90 pid.f90 sim_param.f90 string_util.f90
    types.f90)

set(COMMON_LES_CORE_SOURCES
    iwmles.f90 convec.f90 derivatives.f90 divstress_uv.f90 divstress_w.f90 emul_complex.f90
    fft.f90 interpolag_Sdep.f90 interpolag_Ssim.f90 lagrange_Sdep.f90
    lagrange_Ssim.f90 press_stag_array.f90 scaledep_dynamic.f90 sgs_param.f90
    sgs_stag_util.f90 std_dynamic.f90 test_filtermodule.f90 tridag_array.f90
    wallstress.f90)

set(COMMON_FORCING_SOURCES
    coriolis.f90 forcing.f90 fringe.f90 inflow.f90 shifted_inflow.f90 sponge.f90)

set(COMMON_DIAGNOSTIC_SOURCES
    cuda_mpi_debug.f90 io.f90 rmsdiv.f90 stat_defs.f90 time_average.f90)

set(Sources
    ${COMMON_DRIVER_SOURCES}
    ${COMMON_LES_CORE_SOURCES}
    ${COMMON_FORCING_SOURCES}
    ${COMMON_DIAGNOSTIC_SOURCES})

set(MPI_SOURCES
    mpi_transpose_mod.f90 mpi_defs.f90)

set(CPS_SOURCES
    concurrent_precursor.f90)

set(HIT_SOURCES
    hit_inflow_gpu.f90 hit_inflow.f90)

set(LVLSET_SOURCES
    level_set_base.f90 level_set.f90 linear_simple.f90 trees_pre_ls.f90
    trees_base_ls.f90 trees_setup_ls.f90 trees_io_ls.f90
    trees_global_fmask_ls.f90)

set(TURBINE_SOURCES
    turbines_gpu.f90 turbines.f90 turbine_indicator.f90)

set(ATM_SOURCES
    atm_base.f90 atm_input_util.f90 actuator_turbine_model.f90
    atm_lesgo_interface.f90)

set(LES_GPU_SOURCES
    fft_gpu.f90 tridag_gpu.f90 convec_gpu.f90 derivatives_gpu.f90
    press_gpu.f90 lagrange_Sdep_gpu.f90 sgs_gpu.f90)

set(SCALAR_SOURCES
    scalars.f90 stability.f90)

if (USE_MPI)
    add_definitions(-DPPMPI)

    # Select the MPI Fortran wrapper for the active compiler/platform.
    if (${hostname} STREQUAL "derecho" OR ${hostname} STREQUAL "delta")
        set(CMAKE_Fortran_COMPILER ftn)
        message(STATUS "Setting MPI compiler to ftn (Cray wrapper for Derecho/Cray-style profile)")
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
        # Intel builds use the modern MPI wrapper.
        set(CMAKE_Fortran_COMPILER mpiifx)
        message(STATUS "Setting MPI compiler to mpiifx (Intel)")
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
        # NVHPC builds use an MPI wrapper around nvfortran.
        set(CMAKE_Fortran_COMPILER mpif90)
        message(STATUS "Setting MPI compiler to mpif90 (NVIDIA)")
    else()
        # Other compilers, including GNU, use the default mpif90 wrapper.
        set(CMAKE_Fortran_COMPILER mpif90)
        message(STATUS "Setting MPI compiler to mpif90 (Default/GNU)")
    endif()

    list(APPEND Sources ${MPI_SOURCES})
    set(exec_name "${exec_name}-mpi")
endif (USE_MPI)

if (USE_CPS)
    add_definitions(-DPPCPS)
    list(APPEND Sources ${CPS_SOURCES})
    set(exec_name "${exec_name}-cps")
endif (USE_CPS)

if (USE_HIT)
    add_definitions(-DPPHIT)
    list(APPEND Sources ${HIT_SOURCES})
    set(exec_name "${exec_name}-HIT")
endif (USE_HIT)

if (USE_LVLSET)
    add_definitions(-DPPLVLSET)
    list(APPEND Sources ${LVLSET_SOURCES})
    set(exec_name "${exec_name}-ls")
endif (USE_LVLSET)

if (USE_TURBINES)
    add_definitions(-DPPTURBINES)
    list(APPEND Sources ${TURBINE_SOURCES})
    set(exec_name "${exec_name}-turbines")
endif (USE_TURBINES)

if (USE_ATM)
    add_definitions(-DPPATM)
    list(APPEND Sources ${ATM_SOURCES})
    set(exec_name "${exec_name}-ATM")
endif (USE_ATM)

if (USE_LES_GPU)
    add_definitions(-DPPLES_GPU -DPPCONVEC_GPU -DPPDERIVS_GPU -DPPPRESS_GPU -DPPSGS_GPU)
    list(APPEND Sources ${LES_GPU_SOURCES})
    link_libraries(nvhpcwrapnvtx)
    set(exec_name "${exec_name}-lesgpu")
endif (USE_LES_GPU)

if (USE_DYN_TN)
    add_definitions(-DPPDYN_TN)
    set(exec_name "${exec_name}-dyntn")
endif (USE_DYN_TN)

if (USE_SAFETYMODE)
    add_definitions(-DPPSAFETYMODE)
else (USE_SAFETYMODE)
    set(exec_name "${exec_name}-safety_off")
endif (USE_SAFETYMODE)

if (USE_CGNS)
    add_definitions(-DPPCGNS)
    if (DEFINED CGNS_ROOT)
        include_directories(BEFORE ${CGNS_ROOT}/include)
        link_directories(${CGNS_ROOT}/lib ${CGNS_ROOT}/lib64)
    elseif (DEFINED ENV{CGNS_ROOT})
        include_directories(BEFORE $ENV{CGNS_ROOT}/include)
        link_directories($ENV{CGNS_ROOT}/lib $ENV{CGNS_ROOT}/lib64)
    endif ()
    if (DEFINED HDF5_ROOT)
        include_directories(BEFORE ${HDF5_ROOT}/include)
        link_directories(${HDF5_ROOT}/lib ${HDF5_ROOT}/lib64)
    elseif (DEFINED ENV{HDF5_ROOT})
        include_directories(BEFORE $ENV{HDF5_ROOT}/include)
        link_directories($ENV{HDF5_ROOT}/lib $ENV{HDF5_ROOT}/lib64)
    endif ()
    link_libraries(cgns hdf5_hl hdf5)
    set(exec_name "${exec_name}-cgns")
endif (USE_CGNS)

if (USE_SCALARS)
    add_definitions(-DPPSCALARS)
    if (USE_SCALARS_GPU)
        add_definitions(-DPPSCALARS_GPU)
        set(exec_name "${exec_name}-scalgpu")
    endif ()
    set(exec_name "${exec_name}-scalars")
    list(APPEND Sources ${SCALAR_SOURCES})
endif (USE_SCALARS)

if (WRITE_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPWRITE_LITTLE_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "LITTLE")

if (WRITE_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPWRITE_BIG_ENDIAN)
endif (WRITE_ENDIAN STREQUAL "BIG")

if (READ_ENDIAN STREQUAL "LITTLE")
    add_definitions(-DPPREAD_LITTLE_ENDIAN)
endif (READ_ENDIAN STREQUAL "LITTLE")

if (READ_ENDIAN STREQUAL "BIG")
    add_definitions(-DPPREAD_BIG_ENDIAN)
endif (READ_ENDIAN STREQUAL "BIG")

# ===================================================================
# Final compiler flag selection.
# ===================================================================
if (${hostname} STREQUAL "derecho" OR ${hostname} STREQUAL "delta")
  if (USE_CPU_BUILD)
    # CPU baseline build: -cuda so the dormant cuf/managed
    # declarations still compile, but no -gpu= target and no -cudalib so the
    # binary does not link libcuda/cufft and loads on GPU-less CPU nodes.
    message(STATUS "Applying CPU-baseline NVHPC flags (USE_CPU_BUILD)")
    set (CMAKE_Fortran_FLAGS "-Mpreprocess -cuda -O3 -fast")
  else ()
    message(STATUS "Applying Derecho NVHPC/CUDA flags for A100")
    set (CMAKE_Fortran_FLAGS "-Mpreprocess -cuda -O3 -gpu=mem:separate,cc80,lineinfo -cudalib=cufft,cusolver,cublas -fast")
    if (USE_LES_GPU)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc=gpu -Minfo=accel")
        # GPU-aware MPI: defines PPGPU_AWARE_MPI so per-step halos and the
        # chunk-pipelined tridag pass device pointers (host_data use_device).
        # Without it the host-staged fallbacks compile - correct but the
        # pressure pipeline is lost (~6x slower pressure solve).
        #   AUTO: enable iff Cray's CUDA GTL exists (cray-mpich needs it linked)
        #   ON:   force-enable; USE_MPI must be ON and MPI must be CUDA-aware
        #         (OpenMPI+UCX needs no extra lib; Cray GTL is still linked
        #         when present)
        #   OFF:  host-staged paths only
        if (USE_MPI AND NOT USE_GPU_AWARE_MPI STREQUAL "OFF" AND
            (USE_GPU_AWARE_MPI STREQUAL "ON" OR EXISTS "/opt/cray/pe/lib64/libmpi_gtl_cuda.so"))
            if (EXISTS "/opt/cray/pe/lib64/libmpi_gtl_cuda.so")
                link_directories("/opt/cray/pe/lib64")
                link_libraries(mpi_gtl_cuda)
                set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wl,-rpath=/opt/cray/pe/lib64")
                message(STATUS "USE_LES_GPU: GPU-aware MPI ON (Cray libmpi_gtl_cuda linked)")
            else ()
                message(STATUS "USE_LES_GPU: GPU-aware MPI ON (forced; MPI must be CUDA-aware, e.g. OpenMPI+UCX)")
            endif ()
            add_definitions(-DPPGPU_AWARE_MPI)
        elseif (NOT USE_MPI)
            message(STATUS "GPU-aware MPI OFF (USE_MPI=OFF)")
        else ()
            message(WARNING "GPU-aware MPI OFF (USE_GPU_AWARE_MPI=${USE_GPU_AWARE_MPI}, no Cray GTL found) - host-staged halos + monolithic tridag; pressure ~6x slower")
        endif ()
    endif ()
  endif ()
elseif (CMAKE_Fortran_COMPILER MATCHES "nvfortran" OR CMAKE_Fortran_COMPILER MATCHES "pgf90" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    message(STATUS "Detected NVIDIA Compiler for Flags")
    set (CMAKE_Fortran_FLAGS "-Mpreprocess -cuda -O3 -gpu=mem:separate -cudalib=cufft,cusolver,cublas -fast")
    if (USE_LES_GPU)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc=gpu -Minfo=accel")
    endif ()
    #set (CMAKE_Fortran_FLAGS "-Mpreprocess  -O3 -gpu=mem:managed  -fast")

elseif (CMAKE_Fortran_COMPILER MATCHES "ifx" OR CMAKE_Fortran_COMPILER MATCHES "ifort" OR CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(STATUS "Detected Intel Compiler for Flags")
    set (CMAKE_Fortran_FLAGS "-cpp -O3 -xHost -ipo -qmkl")

else ()
    message(STATUS "Defaulting to GNU Compiler Flags")
    set (CMAKE_Fortran_FLAGS "-Mpreprocess -cuda -O3 -fast")
endif ()
# ===================================================================

# main lesgo executable
add_executable(${exec_name} ${Sources})
