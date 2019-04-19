BOOL,true = true
BOOL,1 = true
BOOL,false = false
BOOL,0 = false
BOOL, = true
OR,true,true  = true
OR,true,false = true
OR,false,true = true
OR,false,false = false

have.sse2.test, = false
have.sse2 = $(BOOL,$(have.sse2.test,$(SSE2_FLAG)))

have.avx.test, = false
have.avx = $(BOOL,$(have.avx.test,$(AVX_FLAG)))

INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR) $(PTHREAD_FLAGS)

simd.source.sse2.true = linearAlgebra_sse2.c moments_sse2.c
simd.source.sse2.false = 
simd.source.sse2 = $(simd.source.sse2.$(have.sse2))

simd.source.avx.true = linearAlgebra_avx.c
simd.source.avx.false = 

simd.sources = $(simd.source.sse2) $(simd.source.avx.$(have.avx))
simd.objects = $(simd.sources:.c=.o)

LOCAL_SOURCES=linearAlgebra.c blockingThreadManager.c hierarchicalThreadManager.c linearAlgebra.c moments.c simd.c thread.c $(simd.sources)
LOCAL_OBJECTS=linearAlgebra.o blockingThreadManager.o hierarchicalThreadManager.o linearAlgebra.o moments.o simd.o thread.o $(simd.objects)

all : library

library : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) ../misc.a $(LOCAL_OBJECTS)

clean :
	-rm -f config.h
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/misc/alloca.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/intrinsic.h :
$(INCLUDE_DIR)/misc/linearAlgebra.h : $(INCLUDE_DIR)/misc/stddef.h
$(INCLUDE_DIR)/misc/simd.h :
$(INCLUDE_DIR)/misc/stats.h : $(INCLUDE_DIR)/misc/stddef.h $(INCLUDE_DIR)/misc/thread.h
$(INCLUDE_DIR)/misc/stddef.h :
$(INCLUDE_DIR)/misc/thread.h : $(INCLUDE_DIR)/misc/stddef.h

blockingThreadManager.o : blockingThreadManager.c $(INCLUDE_DIR)/misc/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c blockingThreadManager.c -o blockingThreadManager.o

hierarchicalThreadManager.o : hierarchicalThreadManager.c $(INCLUDE_DIR)/misc/thread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c hierarchicalThreadManager.c -o hierarchicalThreadManager.o

linearAlgebra.o : linearAlgebra.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c linearAlgebra.c -o linearAlgebra.o

linearAlgebra.sse2.true = $(CC) $(ALL_CPPFLAGS) $(SSE2_FLAG) $(CFLAGS) -c linearAlgebra_sse2.c -o linearAlgebra_sse2.o
linearAlgebra.sse2.false = touch linearAlgebra_sse2.o
compile.linearAlgebra.sse2 = $(linearAlgebra.sse2.$(have.sse2))

linearAlgebra_sse2.o : linearAlgebra_sse2.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	@$(compile.linearAlgebra.sse2)

linearAlgebra.avx.true = $(CC) $(ALL_CPPFLAGS) $(AVX_FLAG) $(CFLAGS) -c linearAlgebra_avx.c -o linearAlgebra_avx.o
linearAlgebra.avx.false = touch linearAlgebra_avx.o
compile.linearAlgebra.avx = $(linearAlgebra.avx.$(have.avx))

linearAlgebra_avx.o : linearAlgebra_avx.c $(INCLUDE_DIR)/misc/linearAlgebra.h
	@$(compile.linearAlgebra.avx)

moments.o : moments.c $(INCLUDE_DIR)/misc/stats.h $(INCLUDE_DIR)/misc/simd.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c moments.c -o moments.o

moments.sse2.true = $(CC) $(ALL_CPPFLAGS) $(SSE2_FLAG) $(CFLAGS) -c moments_sse2.c -o moments_sse2.o
moments.sse2.false = touch moments_sse2.o
compile.moments.sse2 = $(moments.sse2.$(have.sse2))

moments_sse2.o : moments_sse2.c $(INCLUDE_DIR)/misc/stats.h
	@$(compile.moments.sse2)

simd.o : simd.c $(INCLUDE_DIR)/misc/simd.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c simd.c -o simd.o

thread.o : thread.c $(INCLUDE_DIR)/misc/thread.h pthread.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c thread.c -o thread.o

