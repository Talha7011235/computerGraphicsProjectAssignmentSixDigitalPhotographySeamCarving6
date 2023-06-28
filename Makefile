CC          = c++

#-----------------------------------------
#Optimization ----------------------------
OPT   = -O3 -Wno-deprecated

#-----------------------------------------

TARGETS = seamcarving

OBJECTS = carver.o

#-----------------------------------------

LIBS = 
INCS = -I/usr/local/include/eigen3 -I/usr/include/eigen3 -I'/mnt/d/C++/eigen-3.4.0'

CCOPTS = $(OPT) $(INCS)
LDOPTS = $(OPT) $(INCS) $(LIBS)

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)

test: seamcarving
	./seamcarving River.png River-remove-both.png 640 480
	./seamcarving River.png River-add-height.png 1124 1424
	./seamcarving Balloons.jpg Balloons-remove-height.jpg 640 480
	./seamcarving cabin-in.jpg cabin-remove-width.jpg 312 384
	./seamcarving cabin-in.jpg cabin-add-width.jpg 712 384
	./seamcarving cabin-in.jpg cabin-remove-width-add-height.jpg 312 584
	./seamcarving webPageImages/cabin-in_good_case_before_image.jpg webPageImages/cabin-in_good_case_after_image.jpg 512 584
	./seamcarving man.jpg man-threeQuarter-width.jpg 375 333
	./seamcarving man.jpg man-half-width.jpg 250 333
	./seamcarving man.jpg man-quarter-width.jpg 125 333
	./seamcarving man.jpg man-oneFifth-width.jpg 100 333
	./seamcarving webPageImages/man_bad_case_before_image.jpg webPageImages/man_bad_case_after_image.jpg 100 333

clean:
	/bin/rm -f *.o $(TARGETS)

seamcarving: $(OBJECTS) main6.o
	$(CC) $(OBJECTS) main6.o $(LDOPTS) -o seamcarving

main6.o: main6.cpp
	$(CC)  $(CCOPTS) -c main6.cpp

carver.o: carver.cpp
	$(CC)  $(CCOPTS) -c carver.cpp

.cpp.o:
	$(CC) $(CCOPTS) -c $<
