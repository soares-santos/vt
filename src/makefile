# Targets:
#   all           
#   vt		
#   voronoi.o        
#   triangle.o
#   nr.o
#   distclean	# delete all executables, *.o and *~ files.

# SRC is the dir containing the source files, BIN will contain the executables.

#SRC =$(VT_DIR)/src/
SRC =./
BIN =$(VT_DIR)/bin/

# CC should be set to the name of your favorite C compiler.

CC = gcc

# CSWITCHES is a list of all switches passed to the C compiler.  

CSWITCHES = -O

# libraries and include files

LIB = -L$(CFITSIO_DIR)/lib -L. -L$(HOME)/lib -L./lib -lcfitsio -lm -L$(VT_DIR)/lib
INC = -I$(CFITSIO_DIR)/include -I$(HOME)/include -I. -I./include -I$(VT_DIR)/include

# TRILIBDEFS is a list of definitions used to compile an object code version
#   of Triangle (triangle.o) to be called by another program.  The file
#   "triangle.h" contains detailed information on how to call triangle.o.
# The -DTRILIBRARY should always be used when compiling Triangle into an
#   object file. 
# The -DCDT_ONLY switch gets rid of all meshing algorithms above and beyond 
#   constrained Delaunay triangulation.  Specifically, defining the CDT_ONLY 
#   symbol eliminates the -r, -q, -a, -S, and -s switches.  The REDUCED and 
#   CDT_ONLY symbols may be particularly attractive when Triangle is called 
#   by another program that does not need all of Triangle's features.

TRILIBDEFS = -DTRILIBRARY -DREDUCED -DCDT_ONLY

# RM should be set to the name of your favorite rm (file deletion program).

RM = /bin/rm -f

# The action starts here.

all: $(BIN)vt
	$(CC) -o $(BIN)fitscopy $(SRC)fitscopy.c $(INC) $(LIB)
	$(CC) -o $(BIN)dumpfits $(SRC)dumpfits.c $(INC) $(LIB)
	$(CC) -o $(BIN)tabmerge $(SRC)tabmerge.c $(INC) $(LIB)
	$(CC) -o $(BIN)liststruc $(SRC)liststruc.c $(INC) $(LIB)
	$(CC) -o $(BIN)match_clusters $(SRC)match_clusters.c $(INC) $(LIB)
	$(CC) -o $(BIN)zbins $(SRC)zbins.cpp -lstdc++
	cp $(SRC)findclusters.sh $(BIN)findclusters.sh
	cp $(SRC)cleancat.sh $(BIN)cleancat.sh
	cp $(SRC)get_boundaries.sh $(BIN)get_boundaries.sh
	cp $(SRC)matchscript.sh $(BIN)matchscript.sh
	cp $(SRC)vtfind $(BIN)vtfind
	cp $(SRC)rankclusters.sh $(BIN)rankclusters.sh
	$(RM) $(BIN)*.o
	cp $(SRC)*.h $(VT_DIR)/include

$(BIN)vt: $(BIN)nr.o $(BIN)triangle.o $(BIN)voronoi.o $(SRC)vt.c
	$(CC) $(CSWITHCES) -o $(BIN)vt $(BIN)nr.o $(BIN)nrutil.o $(BIN)triangle.o $(BIN)voronoi.o $(SRC)vt.c $(INC) $(LIB)

$(BIN)voronoi.o: $(SRC)voronoi.c $(SRC)voronoi.h 
	$(CC) $(CSWITCHES) -c -o $(BIN)voronoi.o $(SRC)voronoi.c $(INC)

$(BIN)triangle.o: $(SRC)triangle.c $(SRC)triangle.h
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c -o $(BIN)triangle.o $(SRC)triangle.c

$(BIN)nr.o: $(SRC)nr.c $(SRC)nr.h $(SRC)nrutil.c $(SRC)nrutil.h
	$(CC) $(CSWITCHES) -c $(SRC)nr.c $(SRC)nrutil.c
	mv $(SRC)nr.o $(SRC)nrutil.o $(BIN)

match: $(SRC)match_clusters.c $(SRC)voronoi.h 
	$(CC) -o $(BIN)match_clusters $(SRC)match_clusters.c $(INC) $(LIB)

zbins: $(SRC)zbins.cpp  
	$(CC) -o $(BIN)zbins $(SRC)zbins.cpp -lstdc++

scripts: $(SRC)findclusters.sh $(SRC)cleancat.sh $(SRC)get_boundaries.sh $(SRC)matchscript.sh $(SRC)vtfind
	cp $(SRC)findclusters.sh $(BIN)findclusters.sh
	cp $(SRC)cleancat.sh $(BIN)cleancat.sh
	cp $(SRC)get_boundaries.sh $(BIN)get_boundaries.sh
	cp $(SRC)matchscript.sh $(BIN)matchscript.sh
	cp $(SRC)vtfind $(BIN)vtfind
	cp $(SRC)rankclusters.sh $(BIN)rankclusters.sh

distclean:
	$(RM) $(BIN)*.o $(BIN)vt $(SRC)*~ $(BIN)fitscopy $(BIN)dumpfits $(BIN)tabmerge $(BIN)liststruc
