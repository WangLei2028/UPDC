CC ?= cc
PKG_CONFIG ?= pkg-config
CFLAGS ?= -O3 -std=c11 -DNDEBUG -Wall -Wextra
CPPFLAGS += -DUPDC_FFTW_THREADS -DUPDC_PARALLEL
FFTW_CFLAGS := $(shell $(PKG_CONFIG) --cflags fftw3f 2>/dev/null)
FFTW_LIBS := $(shell $(PKG_CONFIG) --libs fftw3f 2>/dev/null)
FFTW_THREADS_LIBS := $(shell $(PKG_CONFIG) --libs fftw3f_threads 2>/dev/null)
ifeq ($(strip $(FFTW_THREADS_LIBS)),)
FFTW_THREADS_LIBS := -lfftw3f_threads
endif

all: updc_jwst_c updc_psf_projector

updc_jwst_c: updc_core.c overlap.c overlap_geometry.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(FFTW_CFLAGS) -pthread -o $@ updc_core.c overlap.c $(FFTW_THREADS_LIBS) $(FFTW_LIBS) -lm

updc_psf_projector: psf_projector.c overlap.c overlap_geometry.h
	$(CC) $(CFLAGS) -pthread -o $@ psf_projector.c overlap.c -lm

check: updc_jwst_c updc_psf_projector
	./updc_jwst_c --self-test
	./updc_psf_projector --self-test

clean:
	rm -f updc_jwst_c updc_psf_projector

.PHONY: all check clean
