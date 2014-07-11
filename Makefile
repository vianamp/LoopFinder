CC = g++

DFLAG = -DDEBUG

CFLAGS = -O1 -O2 -O3 -Os

debug: LoopFinder.cpp
	@$(CC) $(DFLAG) $(CFLAGS) LoopFinder.cpp -o LoopFinder

dynamic: LoopFinder.cpp
	@$(CC) $(CFLAGS) LoopFinder.cpp -o LoopFinder

static:
	@$(CC) $(CFLAGS) LoopFinder.cpp -o LoopFinder