
OBJ = obj/param.o obj/matrix.o obj/protein.o obj/syncbuf.o obj/query.o obj/hit.o obj/database.o obj/index.o obj/align.o obj/stat.o obj/save.o obj/dpu_xfer.o

SRC = param.c matrix.c protein.c hit.c align.c query.c syncbuf.c index.c database.c stat.c save.c dpu_xfer.c

FLAG = -O3 -Wall

NB_TASKLET = 12

pang : host_pang.c $(OBJ) constant.h dpu_pang
	gcc --std=c99 $(FLAG) -march=native  -fopenmp host_pang.c -o pang $(OBJ) -lm `dpu-pkg-config --cflags --libs dpu`

dpu_pang : dpu_pang.c constant.h
	dpu-upmem-dpurte-clang -DNR_TASKLETS=$(NB_TASKLET) -DSTACK_SIZE_DEFAULT=2048 -o dpu_pang dpu_pang.c

struct.h : $(SRC) constant.h include.h
	python3 gen_struct.py $(SRC) > struct.h

obj/param.o : param.c struct.h
	gcc --std=c99 $(FLAG) -march=native -c param.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv param.o obj/param.o

obj/matrix.o : matrix.c struct.h
	gcc --std=c99 $(FLAG) -march=native -c matrix.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv matrix.o obj/matrix.o

obj/protein.o : protein.c struct.h
	gcc --std=c99 $(FLAG) -march=native -c protein.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv protein.o obj/protein.o

obj/syncbuf.o : syncbuf.c struct.h
	gcc --std=c99 $(FLAG) -march=native -c syncbuf.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv syncbuf.o obj/syncbuf.o

obj/query.o : query.c struct.h
	gcc --std=c99 $(FLAG) -march=native -c query.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv query.o obj/query.o

obj/hit.o : hit.c struct.h
	gcc --std=c99 $(FLAG) -march=native -fopenmp -c hit.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv hit.o obj/hit.o

obj/database.o : database.c struct.h
	gcc --std=c99 $(FLAG) -march=native -c database.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv database.o obj/database.o

obj/index.o : index.c struct.h
	gcc --std=c99 $(FLAG) -march=native -fopenmp -c index.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv index.o obj/index.o

obj/align.o : align.c struct.h
	gcc --std=c99 $(FLAG) -march=native -fopenmp -c align.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv align.o obj/align.o

obj/stat.o : stat.c struct.h
	gcc --std=c99 $(FLAG) -march=native  -c stat.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv stat.o obj/stat.o

obj/save.o : save.c struct.h
	gcc --std=c99 $(FLAG) -march=native  -c save.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv save.o obj/save.o

obj/dpu_xfer.o : dpu_xfer.c struct.h
	gcc --std=c99 $(FLAG) -march=native  -c dpu_xfer.c -lm `dpu-pkg-config --cflags --libs dpu`
	mv dpu_xfer.o obj/dpu_xfer.o

