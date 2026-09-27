cmake_minimum_required(VERSION 3.18)
project(updc_jwst_c VERSION 1.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(UPDC_NATIVE_ARCH "Enable native-CPU optimization" ON)

find_package(Threads REQUIRED)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(FFTW3F QUIET IMPORTED_TARGET fftw3f)
  pkg_check_modules(FFTW3F_THREADS QUIET IMPORTED_TARGET fftw3f_threads)
endif()

if(NOT FFTW3F_FOUND)
  find_path(FFTW3F_INCLUDE_DIR fftw3.h)
  find_library(FFTW3F_LIBRARY NAMES fftw3f)
  if(NOT FFTW3F_INCLUDE_DIR OR NOT FFTW3F_LIBRARY)
    message(FATAL_ERROR "FFTW3 single-precision development files (fftw3f) were not found")
  endif()
endif()

if(NOT FFTW3F_THREADS_FOUND)
  find_library(FFTW3F_THREADS_LIBRARY NAMES fftw3f_threads)
  if(NOT FFTW3F_THREADS_LIBRARY)
    message(FATAL_ERROR "FFTW3 single-precision threads library (fftw3f_threads) was not found")
  endif()
endif()


add_executable(updc_psf_projector psf_projector.c overlap.c)
target_compile_options(updc_psf_projector PRIVATE -O3 -Wall -Wextra)
target_link_libraries(updc_psf_projector PRIVATE Threads::Threads m)
add_executable(updc_jwst_c updc_core.c overlap.c)
target_compile_definitions(updc_jwst_c PRIVATE
  UPDC_FFTW_THREADS
  UPDC_PARALLEL
)
target_compile_options(updc_jwst_c PRIVATE -O3 -Wall -Wextra)
target_link_libraries(updc_jwst_c PRIVATE Threads::Threads m)

if(FFTW3F_FOUND)
  target_link_libraries(updc_jwst_c PRIVATE PkgConfig::FFTW3F)
else()
  target_include_directories(updc_jwst_c PRIVATE ${FFTW3F_INCLUDE_DIR})
  target_link_libraries(updc_jwst_c PRIVATE ${FFTW3F_LIBRARY})
endif()

if(FFTW3F_THREADS_FOUND)
  target_link_libraries(updc_jwst_c PRIVATE PkgConfig::FFTW3F_THREADS)
else()
  target_link_libraries(updc_jwst_c PRIVATE ${FFTW3F_THREADS_LIBRARY})
endif()

if(UPDC_NATIVE_ARCH)
  include(CheckCCompilerFlag)
  if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    check_c_compiler_flag("-mcpu=native" UPDC_HAVE_NATIVE_FLAG)
    if(UPDC_HAVE_NATIVE_FLAG)
      target_compile_options(updc_jwst_c PRIVATE -mcpu=native)
      target_compile_options(updc_psf_projector PRIVATE -mcpu=native)
    endif()
  else()
    check_c_compiler_flag("-march=native" UPDC_HAVE_NATIVE_FLAG)
    if(UPDC_HAVE_NATIVE_FLAG)
      target_compile_options(updc_jwst_c PRIVATE -march=native)
      target_compile_options(updc_psf_projector PRIVATE -march=native)
    endif()
  endif()
endif()

install(TARGETS updc_jwst_c updc_psf_projector RUNTIME DESTINATION bin)

enable_testing()
add_test(NAME core_self_test COMMAND updc_jwst_c --self-test)
add_test(NAME psf_projector_self_test COMMAND updc_psf_projector --self-test)
