target := pptb.x

ifndef MPI_ENABLE
MPI_ENABLE := 1
endif

ifeq (${MPI_ENABLE},1)
cc_h := $(shell which mpicxx)
else
cc_h := $(shell which g++)
endif

compflags :=
compflags += -DMPI_ENABLE=${MPI_ENABLE}

hflags :=
hflags += -fconcepts-diagnostics-depth=3 -fmax-errors=1
ifeq (${sanny},1)
hflags += -fsanitize=undefined,address -fstack-protector-all
endif
ifeq (${gee},1)
hflags += -g
endif

hflags += -std=c++20 -O3
hflags += ${compflags}
flags = ${hflags}

dflags :=
dflags += -x cu
dflags += -std=c++20 -O3 #-prec-div=false
dflags += -ccbin=${cc_h}
dflags += --extended-lambda
dflags += -w

ifeq (${gee},1)
dflags += -lineinfo
endif

#dflags += --expt-relaxed-constexpr
#dflags += --verbose --keep

cc_d = $(shell which nvcc)

cc := ${cc_h}
flags := ${hflags}

main:
	${cc} ${flags} -I${PPTB}/src -I${SPADE}/src -I${SCIDF}/src main.cc -o ${target}

run: main
	./${target}

clean:
	rm -f ${target}
