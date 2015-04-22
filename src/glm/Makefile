INCLUDE_DIR=../include

PKG_CPPFLAGS=-I$(INCLUDE_DIR)
ALL_CPPFLAGS=$(R_XTRA_CPPFLAGS) $(PKG_CPPFLAGS) $(CPPFLAGS)

LOCAL_SOURCES = glm.c
LOCAL_OBJECTS = glm.o

all : libglm.a

libglm.a : $(LOCAL_OBJECTS)
	$(AR) $(ARFLAGS) libglm.a $(LOCAL_OBJECTS)

clean :
	-rm -f *.o
	-rm -f *.a

rebuild : clean all

$(INCLUDE_DIR)/glm.h :

glm.o : glm.c $(INCLUDE_DIR)/glm/glm.h
	$(CC) $(ALL_CPPFLAGS) $(CFLAGS) -c $< -o $@
