// Isaac Haseley - Nathan Mannes - 2017


#define MAX_SPLINE_LENGTH 512

// Thanks for providing the inspiration for this line of code
typedef struct splineSpline splineSpline; 
struct splineSpline {
	int numPoints;
    double xTrans, yTrans;
    int *xValues;
    int *yValues;
    double *b, *c, *d;
    double *rot;
};

// Sets spline's struct values with defaults for rotation and translation
void initializeSpline(splineSpline *spline, int numPoints, int xValues[], 
    int yValues[]) {
    spline->numPoints = numPoints;
    spline->xTrans = 0;
    spline->yTrans = 0;
    spline->xValues = (int *)malloc(numPoints * sizeof(int));
    spline->yValues = (int *)malloc(numPoints * sizeof(int));
    spline->b = (double *)malloc((numPoints - 1) * sizeof(double));
    spline->c = (double *)malloc((numPoints - 1) * sizeof(double));
    spline->d = (double *)malloc((numPoints - 1) * sizeof(double));
    spline->rot = (double *)malloc(4 * sizeof(double));
    for (int i = 0; i < numPoints; i++) {
        spline->xValues[i] = xValues[i];
        spline->yValues[i] = yValues[i];
    }
    spline->rot[0] = 1;
    spline->rot[1] = 0;
    spline->rot[2] = 0;
    spline->rot[3] = 1;
}

// Sets spline's struct values, including custom rotation and translation
void initializeSplinePlus(splineSpline *spline, int numPoints, int xValues[], 
    int yValues[], double rot[], double xTrans, double yTrans) {
    spline->numPoints = numPoints;
    spline->xValues = xValues;
    spline->yValues = yValues;
    spline->rot = rot;
    spline->xTrans = xTrans;
    spline->yTrans = yTrans;
}

// Deallocates the resources backing the given spline
void splineDestroy(splineSpline *spline) {
    free(spline->xValues);
    free(spline->yValues);
    free(spline->b);
    free(spline->c);
    free(spline->d);
    free(spline->rot);
}

// Draws 9 pixels instead of 1
void drawThickPoint(double x, double y, double r, double g, double b) {
    pixSetRGB(x + 1, y, r, g, b);
    pixSetRGB(x, y, r, g, b);
    pixSetRGB(x - 1, y, r, g, b);

    pixSetRGB(x + 1, y + 1, r, g, b);
    pixSetRGB(x, y + 1, r, g, b);
    pixSetRGB(x - 1, y + 1, r, g, b);

    pixSetRGB(x + 1, y - 1, r, g, b);
    pixSetRGB(x, y - 1, r, g, b);
    pixSetRGB(x - 1, y - 1, r, g, b);
}

// Calculates y-value at a given x-value on a given cubic polynomial
// polyIndex is the spline's polynomial, indexed from 0
// xLeft is the x-value at which said polynomial begins
double polynomial(splineSpline *spline, int polyIndex, int xLeft, int x) {
    return  spline->yValues[polyIndex] + 
            spline->b[polyIndex] * (x - xLeft) +
            spline->c[polyIndex] * (x - xLeft) * (x - xLeft) +
            spline->d[polyIndex] * (x - xLeft) * (x - xLeft) * (x - xLeft);
}

// Draws vertical segment at x from yBottom to yTop
// Colors segment if fill == 1
// Draws thick curves if thick == 1
void drawVertical(int x, int yBottom, int yTop, int fill, int thick) {
    for (int y = yBottom + 1; y < yTop; y++) {
        if (fill == 1) {
            if (thick == 1) {
                drawThickPoint(x, y, 0, 200.0 / (double)x, (double)x / 500.0);
            }
            else {
                pixSetRGB(x, y, 0, 200.0 / (double)x, (double)x / 500.0);
            }
        }
        else {
            if (thick == 1) {
                drawThickPoint(x, y, 0, 0, 0);
            }
            else {
                pixSetRGB(x, y, 0, 0, 0);
            }
        }
    }
}

// Assumes spline's x-values are in increasing order
// Draws thick curves if thick == 1
// Colors spline coordinates red if red == 1
void splineRender(splineSpline *spline, int thick, int red) {
    
    // Calculates h-values
    double h[spline->numPoints - 1];
    for (int i = 0; i < spline->numPoints - 1; i++) {
        h[i] = spline->xValues[i + 1] - spline->xValues[i];
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
        tridiagonal[i] = h[hIndex];
        tridiagonal[i + 1] = 2 * (h[hIndex] + h[hIndex + 1]);
        tridiagonal[i + 2] = h[hIndex + 1];
        hIndex++;
    }
    
    // Builds constant vector for system of equations
    double constants[spline->numPoints], c[spline->numPoints];
    constants[0] = 0;
    constants[spline->numPoints - 1] = 0;
    for (int i = 1; i < spline->numPoints - 1; i++) {
        constants[i] = (3 / h[i]) * (spline->yValues[i + 1] - spline->yValues[i])
            - (3 / h[i - 1]) * (spline->yValues[i] - spline->yValues[i - 1]);
    }
    
    // Solves system of equations for polynomials' c-values
    double inverse[spline->numPoints * spline->numPoints];
    invertTridiagonal(spline->numPoints, tridiagonal, inverse);
    matVecMult(spline->numPoints, inverse, constants, c);
    spline->c = c;
    
    // Uses c-values to solve for remaining polynomial constants
    double b[spline->numPoints - 1], d[spline->numPoints - 1];
    for (int i = 0; i < spline->numPoints - 1; i++) {
        b[i] = (1 / h[i]) * (spline->yValues[i + 1] - spline->yValues[i]) -
            (h[i] / 3) * (2 * c[i] + c[i + 1]);
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    spline->b = b;
    spline->d = d;
    
    // For each x, uses corresponding polynomial to find y-value
    // Calls pixSetRGB on xy-coordinate
    int *a = spline->yValues;
    int previousX = spline->xValues[0];
    int previousY = spline->yValues[0];
    for (int i = 0; i < spline->numPoints - 1; i++) {
        for (double x = spline->xValues[i]; x <= spline->xValues[i + 1]; x++) {
            // Rasterize pixel on spline
            double currentY = round(polynomial(spline, i, spline->xValues[i], x));
            double currentX = x;
            
                // Starts on rotation code
                currentX -= spline->xTrans;
                currentY -= spline->yTrans;
                currentX = spline->rot[0] * currentX + spline->rot[1] * currentY;
                currentY = spline->rot[2] * currentX + spline->rot[3] * currentY;
                currentX = round(currentX + spline->xTrans);
                currentY = round(currentY + spline->yTrans);

            if (thick == 1) {
                drawThickPoint(currentX, currentY, 0, 200.0 / (double)currentX, (double)currentX / 500.0);
            }
            else {
                pixSetRGB(currentX, currentY, 0, 200.0 / (double)currentX, (double)currentX / 500.0);
            }
            // Rasterize pixel gap
            if (currentY > previousY) {
                drawVertical(currentX, previousY, currentY, 1, thick);
            }
            else if (currentY < previousY) {
                drawVertical(currentX, currentY, previousY, 1, thick);
            }
            previousX = currentX;
            previousY = currentY;
        }
    }
    if (red == 1) {
        for (int i = 0; i < spline->numPoints; i++) {
            drawThickPoint(spline->xValues[i], spline->yValues[i], 1, 0, 0);
        }
    }
}

// Work in progress
// We're stumped on the bug here - run mainFillingDebug and see prints below
// Rasterizes area between two splines
void fillBetween(splineSpline *spline1, splineSpline *spline2, int fill) {
    // Just after calling fillBetween, our topSpline constants are:
    printf("topSpline: b[0] = %f, c[0] = %f, d[0] = %f\n", spline1->b[0], 
        spline1->c[0], spline1->d[0]);
    printf("topSpline: b[1] = %f, c[1] = %f, d[1] = %f\n", spline1->b[1], 
        spline1->c[1], spline1->d[1]);
    fflush(stdout);
    /* 
    In our terminal, we see the following: 
    
    topSpline: b[0] = -1.500000, c[0] = 0.000000, d[0] = 0.000247
    topSpline: b[1] = 0.000000, c[1] = 0.000000, d[1] = 0.000000
    topSpline: b[0] = 0.000000, c[0] = 0.000000, d[0] = 0.000000
    topSpline: b[1] = 0.000000, c[1] = 0.000000, d[1] = 0.000000
    
    This indicates that between the end of splineRender and the 
        start of fillBetween, topSpline's b[0] and d[0] values
        change to 0. 
    */
    
    int xStart = spline1->xValues[0];
    if (spline2->xValues[0] > xStart) {
        xStart = spline2->xValues[0];
    }
    int xEnd = spline1->xValues[spline1->numPoints - 1];
    if (spline2->xValues[spline2->numPoints - 1] < xEnd) {
        xEnd = spline2->xValues[spline2->numPoints - 1];
    }
    for (int x = xStart; x <= xEnd; x++) {
        // Determine which polynomial x is on for each spline
        int polyIndex1 = 0;
        for (int i = 0; i < spline1->numPoints - 1; i++) {
            if (spline1->xValues[i] <= x && x <= spline1->xValues[i + 1]) {
                polyIndex1 = i;
                break;
            }
        }
        int polyIndex2 = 0;
        for (int i = 0; i < spline2->numPoints - 1; i++) {
            if (spline2->xValues[i] <= x && x <= spline2->xValues[i + 1]) {
                polyIndex2 = i;
                break;
            }
        }
                
        // Determine which x-value is on top, draw vertical
        int y1 = round(polynomial(spline1, polyIndex1, 
            spline1->xValues[polyIndex1], x));
        int y2 = polynomial(spline2, polyIndex2, 
            spline2->xValues[polyIndex2], x);        
        if (y1 < y2) {
            drawVertical(x, y1, y2, fill, 1);
        }
        else {
            drawVertical(x, y2, y1, fill, 1);
        }
    }
}

// Renders an ellipse centered at x, y with given radii
// Ellipse is composed of two semicircle splines, one on top of the other
// numPoints is the number of spline coordinates on the ellipse's border
// fill is used in mainFilling
void splineEllipse(double x, double y, double xRadius, double yRadius, 
    int numPoints, int fill) {
    
	// Calculates points for splines
    int pointsInHalfCircle = (numPoints / 2) + 1;
	double deltaT = M_PI / (pointsInHalfCircle - 1);
    int xValues[pointsInHalfCircle], yTop[pointsInHalfCircle], 
        yBottom[pointsInHalfCircle];
    xValues[0] = (int) round(x - xRadius);
    yTop[0] = (int) round(y);
    yBottom[0] = (int) round(y);  
    int index = 0;
	for (double t = M_PI - deltaT; t >= -.001; t -= deltaT) {    
        int newX = (int) round(xRadius * cos(t) + x);
        if (t > M_PI / 2) {
            if (newX != xValues[index]) {
                index++;
                xValues[index] = newX;
                yTop[index] = (int) round(yRadius * sin(t) + y);
                yBottom[index] = (int) round(yRadius * sin(-1 * t) + y);
            }
        }
        else {
            if (newX != xValues[index]) {
                index++;
            }
            xValues[index] = newX;
            yTop[index] = (int) round(yRadius * sin(t) + y);
            yBottom[index] = (int) round(yRadius * sin(-1 * t) + y);
        }
	}

	// Initializes and draws splines
	splineSpline topSpline, bottomSpline;
	initializeSpline(&topSpline, index + 1, xValues, yTop);
	splineRender(&topSpline, 0, 0);
	initializeSpline(&bottomSpline, index + 1, xValues, yBottom);
	splineRender(&bottomSpline, 0, 0);
    
    // Displays bug we've isolated in filling code
    if (fill == 1) {
        // Just before calling fillBetween, our topSpline constants are:
        printf("topSpline: b[0] = %f, c[0] = %f, d[0] = %f\n", topSpline.b[0], 
            topSpline.c[0], topSpline.d[0]);
        printf("topSpline: b[1] = %f, c[1] = %f, d[1] = %f\n", topSpline.b[1], 
            topSpline.c[1], topSpline.d[1]);
        fflush(stdout);
        fillBetween(&topSpline, &bottomSpline, fill);
    }
}