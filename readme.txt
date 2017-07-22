Isaac Haseley - Nathan Mannes - 2017

Our project investigates the math behind data-fitting cubic splines.
Our files are as follows:

0. linAl.c

    linAl contains linear algebra functions used to render splines.
    
1. spline.c

    This file defines a spline struct and several accompanying functions used by 
    our mains. 
    This first version of a cubic spline is "horizontal": it does not rotate back
    on itself.
    
    Works in progress:
        Our unfinished code for rotating splines is present but does not impede
        other functions.
        Our code for filling the area between two horizontal splines contains a 
        puzzling bug that we're tracking down via print statements.
        
2. mainInteractive.c

    Compile and run with:
    clang mainInteractive.c 000pixel.o -lglfw -framework OpenGL && ./a.out
    
    Directions:
        Click anywhere to add a point to the horizontal spline.
        Points added directly above or below existing points are ignored.
        
3. mainShapes.c

    Compile and run with:
    clang mainShapes.c 000pixel.o -lglfw -framework OpenGL && ./a.out

    mainShapes explores the use of horizontal splines to create ellipses.
    Each ellipse is composed of two semicircular splines, one above the other.
    The number of points on each spline increases across each of the five columns.

4. mainParaInteractive.c

    Compile and run with:
    clang mainParaInteractive.c 000pixel.o -lglfw -framework OpenGL && ./a.out
    
    This is an interactive spline generator that lets the user play
    with a spline. Blue pixels dominate the beginning of each polynomial, so
    the user can easily visualize each spline component.
    
    Directions:
        Left click far from any points to add points to the end of the spline. 	
        Left click close to a spline to move it around.
        Left click again to drop it in place.
        Right click close to a linchpin point to remove it.

Sources cited:

    We first implemented the math by staring at this document for like a week:
    http://www.physics.arizona.edu/~restrepo/475A/Notes/sourcea-/node35.html

    How we invert a tridiagonal matrix: http://www.mat.uc.pt/preprints/ps/p0516.pdf

    Derived inspiration for parametrizing from this code base:
    https://www.codeproject.com/articles/560163/csharp-cubic-spline-interpolation