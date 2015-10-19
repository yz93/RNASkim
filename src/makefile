OBJS = rs_estimate.o rs_estimate_lib.o rs_selected_key.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

p1 : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o p1

rs_estimate.o : rs_estimate_lib.h rs_estimate.cc
	$(CC) $(CFLAGS) rs_estimate.cc

rs_estimate_lib.o : rs_estimate_lib.h rs_estimate_lib.cc
	$(CC) $(CFLAGS) rs_estimate_lib.cc

rs_selected_key.o : rs_selected_key.h rs_selected_key.cc
	$(CC) $(CFLAGS) rs_selected_key.cc

clean:
	\rm *.o *~ p1
