// Isaac Haseley - Nathan Mannes - 2017

/*
Compile and run with:
clang mainInteractive.c 000pixel.o -lglfw -framework OpenGL && ./a.out
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "000pixel.h"
#include "linAl.c"
#include "spline.c"

double xPos, yPos; // mousePosition
splineSpline spline; // spline

// Returns the desired index of an added point
int locateIndex(splineSpline *spline, int x) {
    for (int i = 0; i < spline->numPoints; i++) {
        if (x == spline->xValues[i]) {
            return -1;
        }
        if (x < spline->xValues[i]) {
            return i;
        }
    }
    return spline->numPoints;
}

// Inserts a point into the spline in an insertion sort kind of way
void addPoint(splineSpline *spline, int x, int y) {
    if (locateIndex(spline, x) == -1) {
        return;
    }
    int i;
    for (i = spline->numPoints; i > locateIndex(spline, x); i--) {
        spline->xValues[i] = spline->xValues[i - 1];
        spline->yValues[i] = spline->yValues[i - 1];
    }
    spline->xValues[i] = x;
    spline->yValues[i] = y;
    spline->numPoints += 1;
}

// Tracks mouse position
void handleMouseMove(double x, double y) {
    xPos = x;
    yPos = y;
}

// Adds new spline coordinates on click
void handleMouseDown(int button, int shiftIsDown, int controlIsDown,
    int altOptionIsDown, int superCommandIsDown) {
    if (button ==0){
        addPoint(&spline, (int)xPos, (int)yPos);
        pixClearRGB(0, 0, 0);
        splineRender(&spline, 1, 1);
    }
}

int main(void) {
	if (pixInitialize(512, 512, "Final Project") != 0) {
		return 1;
    }
    else {
		pixSetMouseDownHandler(handleMouseDown);
		pixSetMouseMoveHandler(handleMouseMove);
		      
        int xValues[MAX_SPLINE_LENGTH];
        int yValues[MAX_SPLINE_LENGTH];
        initializeSpline(&spline, 0, xValues, yValues);
                        
        pixRun();
        splineDestroy(&spline);
        return 0;
    } 
}
