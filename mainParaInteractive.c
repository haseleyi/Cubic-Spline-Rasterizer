// Isaac Haseley - Nathan Mannes - 2017

#include <stdio.h>
#include <math.h>
#include "000pixel.h"
#include "linAl.c"

#define MAX_SPLINE_LENGTH 500

typedef struct splineSpline splineSpline; 
struct splineSpline {
  int *xValues;
  int *yValues;
  int numPoints;
  int softBoundIndex; // For seeing which index we would theoretically use if we were to do stuff
  int hardBoundIndex; // For actually using an index
double **bcdX, **bcdY; // For storing the polynomial info
};

// GLOBAL VARIABLES
int xPos, yPos; // mousePosition
splineSpline spline; // spline

double distance(double x0, double y0, double x1, double y1) {
  return sqrt(((x1 - x0) * (x1 - x0)) + ((y1 - y0) * (y1 - y0)));
}

// Calculates y-value at a given x-value on a given cubic polynomial
double polynomial(double t, int a, double b, double c, double d) {
  return  a + 
    b * t +
    c * t * t +
    d * t * t * t;
}

// Sets spline's struct values
void initializeSpline(splineSpline *spline, int numPoints, int xValues[], 
                      int yValues[]) {
  spline->numPoints = numPoints;
  if (numPoints == 0) {
    int x[MAX_SPLINE_LENGTH];
    int y[MAX_SPLINE_LENGTH];
    spline->xValues = x;
    spline->yValues = y;
  }
  else {
    spline->softBoundIndex = -1;
    spline->hardBoundIndex = -1;
    spline->xValues = xValues;
    spline->yValues = yValues;
  }
}

// A point is softbound if it is within 25 pixels of a point on a spline,
// A point is hardbound if a softbound point is clicked
int splineIsBoundSoft(splineSpline *spline) {
  return spline->softBoundIndex;
}

int splineIsBoundHard(splineSpline *spline) {
  return spline->hardBoundIndex;
}

void splineBindHard(splineSpline *spline) {
  spline->hardBoundIndex = spline->softBoundIndex;
}

void splineUnbind(splineSpline *spline) {
  spline->hardBoundIndex = -1;
  spline->softBoundIndex = -1;
}

// Draws 9 pixels instead of 1
void thickPoint(double x, double y, double r, double g, double b) {
  pixSetRGB(x + 1, y, r, g, b);
  pixSetRGB(x,     y, r, g, b);
  pixSetRGB(x - 1, y, r, g, b);

  pixSetRGB(x + 1, y + 1, r, g, b);
  pixSetRGB(x, y + 1, r, g, b);
  pixSetRGB(x - 1, y + 1, r, g, b);

  pixSetRGB(x + 1, y - 1, r, g, b);
  pixSetRGB(x, y - 1, r, g, b);
  pixSetRGB(x - 1, y - 1, r, g, b);
}

// Draws a slightly thicker point
void thiccerPoint(double x, double y, double r, double g, double b) {

  thickPoint(x + 1, y + 1, r, g, b);
  thickPoint(x - 1, y + 1, r, g, b);
  thickPoint(x + 1, y - 1, r, g, b);
  thickPoint(x - 1, y - 1, r, g, b);
}

// DOES NOT assume spline's x-values are in increasing order!
void splineRender(splineSpline *spline) {
  // Calculates ht-values: the offset between the t values and the t value at the next point
  if (spline->numPoints == 0) {
    return;
  }

  // Write up the tValues
  double tValues[MAX_SPLINE_LENGTH];
  tValues[0] = 0;
  for (int i = 1; i < MAX_SPLINE_LENGTH; i++) {
    tValues[i] = i;
  }
    
  double ht[spline->numPoints - 1];
  for (int i = 0; i < spline->numPoints - 1; i++) {
    ht[i] = tValues[i + 1] - tValues[i];
  }
  // Builds tridiagonal matrix represented by a 1D array
  int tridiagonal[spline->numPoints * spline->numPoints];
  tridiagonal[0] = 1;
  tridiagonal[spline->numPoints * spline->numPoints - 1] = 1;
  for (int i = 1; i < spline->numPoints * spline->numPoints - 1; i++) {
    tridiagonal[i] = 0;
  }
  int hIndex = 0;
  for (int i = spline->numPoints; i < spline->numPoints * spline->numPoints 
         - spline->numPoints - 2; i += spline->numPoints + 1) {
    tridiagonal[i] = ht[hIndex];
    tridiagonal[i + 1] = 2 * (ht[hIndex] + ht[hIndex + 1]);
    tridiagonal[i + 2] = ht[hIndex + 1];
    hIndex++;
  }
    
  // Builds constant vector for system of equations,
  double constantsX[spline->numPoints], cX[spline->numPoints];
  double constantsY[spline->numPoints], cY[spline->numPoints];
  constantsX[0] = 0;
  constantsY[0] = 0;
  constantsX[spline->numPoints - 1] = 0;
  constantsY[spline->numPoints - 1] = 0;
  for (int i = 1; i < spline->numPoints - 1; i++) {
    constantsY[i] = (3 / ht[i]) * (spline->yValues[i + 1] - spline->yValues[i])
      - (3 / ht[i - 1]) * (spline->yValues[i] - spline->yValues[i - 1]);

    constantsX[i] = (3 / ht[i]) * (spline->xValues[i + 1] - spline->xValues[i])
      - (3 / ht[i - 1]) * (spline->xValues[i] - spline->xValues[i - 1]);
  }
    
  // Solves system of equations for polynomials' c-values
  // for t(x) and t(y)
  double inverse[spline->numPoints * spline->numPoints];
  invertTridiagonal(spline->numPoints, tridiagonal, inverse);
  matVecMult(spline->numPoints, inverse, constantsY, cY);
  matVecMult(spline->numPoints, inverse, constantsX, cX);
    
  // Uses c-values to solve for remaining polynomial constants
  // For t(x) and t(y)
  double bY[spline->numPoints - 1], dY[spline->numPoints - 1];
  double bX[spline->numPoints - 1], dX[spline->numPoints - 1];
  for (int i = 0; i < spline->numPoints - 1; i++) {
    bY[i] = (1 / ht[i]) * (spline->yValues[i + 1] - spline->yValues[i]) -
      (ht[i] / 3) * (2 * cY[i] + cY[i + 1]);
    dY[i] = (cY[i + 1] - cY[i]) / (3 * ht[i]);
    
    bX[i] = (1 / ht[i]) * (spline->xValues[i + 1] - spline->xValues[i]) -
      (ht[i] / 3) * (2 * cX[i] + cX[i + 1]);
    dX[i] = (cX[i + 1] - cX[i]) / (3 * ht[i]);
  }
    
  // For each x, uses corresponding polynomial to find y-value
  // Calls pixSetRGB on xy-coordinate
  int *aY = spline->yValues;
  int *aX = spline->xValues;

  double *bcdX[3], *bcdY[3];
  bcdX[0] = bX;
  bcdX[1] = cX;
  bcdX[2] = dX;
  bcdY[0] = bY;
  bcdY[1] = cY;
  bcdY[2] = dY;

  spline->bcdX = bcdX;
  spline->bcdY = bcdY;
  double x, y; // For sanity
  // Now we have pointers to the A+B(t) + C(t^2) + D(t^3)
  // And we have pointers to the A'+B'(t) + C'(t^2) + D'(t^3)
  double t = 0;
  for (int i = 0; i < spline->numPoints - 1; i += 1) {
    for (t = 0; t < 1; t += .001) {
      x = polynomial(t, aX[i], bcdX[0][i], bcdX[1][i], bcdX[2][i]);
      y = polynomial(t, aY[i], bcdY[0][i], bcdY[1][i], bcdY[2][i]);
      pixSetRGB(round(x), round(y), 0, t, 1 - t);
    }
  }
}

// For filling, we implement a crossing number algorithm
int crossingNumber(int x, int y, int numVals, double xVals[], double yVals[]) {
  int x1, x2, y1, y2, intersections;
  x1 = xVals[0];
  x2 = xVals[numVals - 2];
  y1 = yVals[0];
  y2 = yVals[numVals - 2];
  // Hard code the first/last one
  if ((y > y1 && y < y2) || (y < y1 && y > y2)) { 
    intersections++;
  }
  printf("scaption\n");
  fflush(stdout); 
  for (int i = 0; i < numVals; i++) {
    x1 = xVals[i];
    x2 = xVals[i + 1];
    y1 = yVals[i];
    y2 = yVals[i + 1];
    if ((y > y1 && y < y2) || (y < y1 && y > y2)) {
      intersections++;
    }
  }
  return intersections;
}

//to fill a spline, the first and last point MUST be the same!
void fill(splineSpline *spline, double pointsPerSpline) {
  if ((spline->xValues[0] == spline->xValues[spline->numPoints - 1] ) &&
     (spline->yValues[0] == spline->yValues[spline->numPoints - 1] )) {
    int numVals = pointsPerSpline * spline->numPoints;
    int index = 0;;
    double xVals[numVals], yVals[numVals], x, y, t;
    printf("s\n");
    fflush(stdout);
    for (int i = 0; i < spline->numPoints - 1; i += 1) {
      printf("%d\n", i);
      fflush(stdout);
      for (t = 0; t < 1; t+= (1/pointsPerSpline)) {
        printf("%f %d %f %f %f\n", t, spline->xValues[i], spline->bcdX[0][i],
                       spline->bcdX[1][i], spline->bcdX[2][i]);
        fflush(stdout);
        
        x = polynomial(t, spline->xValues[i], spline->bcdX[0][i],
                       spline->bcdX[1][i], spline->bcdX[2][i]);
        printf("%f\n", x);
        fflush(stdout);
        y = polynomial(t, spline->yValues[i], spline->bcdY[0][i],
                       spline->bcdY[1][i], spline->bcdY[2][i]);
        xVals[index] = x;
        yVals[index] = y;
        index++;
      }
    }
    for (int i = 0; i < 512; i++) {
      for (int j = 0; j < 512; j++) {
        if (crossingNumber(x, y, numVals, xVals, yVals) % 2 == 1) {
          pixSetRGB(x, y, 1, 1, 1);
        }
      }
    }
  }
}

// Returns the desired index of an added point
int locateIndex(splineSpline *spline, int x, int y) {
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

int closestPointIndex(splineSpline *spline, int x, int y) {
  double minDistance = 25; 
  double test;
  int index;
  for (int i = 0; i < spline->numPoints; i++) {
    test = distance(x, y, spline->xValues[i], spline->yValues[i]);
    if (test < minDistance) {
      minDistance = test;
      index = i;
    }
  }
  if (minDistance < 25) {
    return index;
  }
  else {
    return -1;
  }
}

void bindIndex(splineSpline *spline, int x, int y) {
  spline->softBoundIndex = closestPointIndex(spline, x, y);
  if (spline->softBoundIndex > -1) {
    pixClearRGB(0,0,0);
    splineRender(spline);
    thiccerPoint(spline->xValues[spline->softBoundIndex],
               spline->yValues[spline->softBoundIndex], 1, 0, 1);
  }
  else {
    pixClearRGB(0,0,0);
    splineRender(spline);
  }
}

// Inserts a point to the end of the spline
void addEndPoint(splineSpline *spline, int x, int y) {
  spline->xValues[spline->numPoints] = x;
  spline->yValues[spline->numPoints] = y;
  spline->numPoints += 1;
}

// Adds a point at a specific index
void insertionAddPoint(splineSpline *spline, int index, int x, int y) {
  for (int i = spline->numPoints; i > index; i--) {
    spline->xValues[i] = spline->xValues[i - 1];
    spline->yValues[i] = spline->yValues[i - 1];
  }
  spline->xValues[index] = x;
  spline->yValues[index] = y;
  spline->numPoints++;
}

void removeIndex(splineSpline *spline, int index) {
  for (int i = index; i < spline->numPoints - 1; i++) {
    spline->xValues[i] = spline->xValues[i + 1];
    spline->yValues[i] = spline->yValues[i + 1];
  }
  spline->numPoints--;
}

// KEYBOARD HANDLERS: SCREEN DOES NOT UPDATE UNTIL MOUSE POSITION IS UPDATED

void handleKeyDown(int key, int shiftIsDown, int controlIsDown,
		int altOptionIsDown, int superCommandIsDown) {
  if (key == 70){
    printf("This button would fill if we had gotten it working\n");
    fflush(stdout);
    //fill(&spline, 100);
  }
}

void handleMouseMove(double x, double y) {
  xPos = round(x);
  yPos = round(y);
  if (splineIsBoundHard(&spline) != -1) {
    pixClearRGB(0, 0, 0);
    removeIndex(&spline, spline.hardBoundIndex);
    insertionAddPoint(&spline, spline.hardBoundIndex, xPos, yPos);
    splineRender(&spline);
  }
  else { //"Binds" point index if it is within 25 pixels
    bindIndex(&spline, x, y);
  }
}

void handleMouseDown(int button, int shiftIsDown, int controlIsDown,
                     int altOptionIsDown, int superCommandIsDown) {
  if (button == 0) { //Adds point if no points are soft-bound, moves point if point is soft bound
    int bind = splineIsBoundSoft(&spline);
    if (bind != -1) {
        // Checks for the existence of a hardbound point
      if (splineIsBoundHard(&spline)> -1) {
        splineUnbind(&spline);
        splineRender(&spline);
      }
        // You click and there is a softbound point,
        // so it hard-binds a point
      else {
        splineBindHard(&spline);
        splineRender(&spline);
      }
    }
    // You click and there are no soft-bound points
    else {
      pixClearRGB(0, 0, 0);
      addEndPoint(&spline, xPos, yPos);
      splineRender(&spline);
      printf("Added a point\n");
    }
  }
  if (button == 1) { //Removes softbound point if right clicked
    if (splineIsBoundSoft(&spline) != -1) {
      pixClearRGB(0, 0, 0);
      removeIndex(&spline, spline.softBoundIndex);
      printf("Point has been removed\n");
      splineUnbind(&spline);
      splineRender(&spline);
    }
  }
  fflush(stdout);
}

int main(void) {
  if (pixInitialize(512, 512, "Cubic Splines") != 0) {
    return 1;
  }
  else {
    /*Key Handlers*/
    pixSetKeyDownHandler(handleKeyDown);
    pixSetMouseDownHandler(handleMouseDown);
    pixSetMouseMoveHandler(handleMouseMove);

    int xValues[MAX_SPLINE_LENGTH] = {200, 300};
    int yValues[MAX_SPLINE_LENGTH] = {256, 256};
    initializeSpline(&spline, 2, xValues, yValues);
    splineRender(&spline);    
    pixRun();
    return 0;
  } 
}
