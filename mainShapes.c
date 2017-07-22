// Isaac Haseley - Nathan Mannes - 2017

/*
Compile and run with:
clang mainShapes.c 000pixel.o -lglfw -framework OpenGL && ./a.out
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "000pixel.h"
#include "linAl.c"
#include "spline.c"

int main(void) {
	if (pixInitialize(512, 512, "Final Project") != 0) {
		return 1;
    }
    else {
        int i = 0;
        int xValues[] = {60, 160, 260, 360, 460};
        int numPoints[] = {4, 10, 16, 40, 400};
        for (int i = 0; i < 5; i++) {
            splineEllipse(xValues[i], 430, 45, 45, numPoints[i], 0);
            splineEllipse(xValues[i], 300, 55, 45, numPoints[i], 0);
            splineEllipse(xValues[i], 150, 45, 75, numPoints[i], 0);
        }
        
        pixRun();
        return 0;
    } 
}