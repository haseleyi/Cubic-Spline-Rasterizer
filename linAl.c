// Isaac Haseley - Nathan Mannes - 2017


// Multiplies a square matrix of arbitrary size by a vector of corresponding length
// mockMatrix is a matrix unpacked into an array
void matVecMult(int rowLength, double mockMatrix[], double inVector[],  
    double outVector[]) {
    double sum;
    for (int r = 0; r < rowLength; r++) {
        sum = 0;
        for (int c = 0; c < rowLength; c++) {
            sum += mockMatrix[r * rowLength + c] * inVector[c];
        }
        outVector[r] = sum;
    }
}


// Inverts a tridiagonal matrix of arbitrary size
// Second and third parameters are matrices unpacked into arrays
void invertTridiagonal(int rowLength, int tridiagonal[], double inverse[]) {
    
    // Strides through matrix and builds three diagonals 
    double a[rowLength], b[rowLength - 1], c[rowLength - 1];
    for (int i = 0; i < rowLength - 1; i++) {
        a[i] = tridiagonal[i * (rowLength + 1)];
        b[i] = tridiagonal[1 + i * (rowLength + 1)];
        c[i] = tridiagonal[rowLength + i * (rowLength + 1)];
    }
    a[rowLength - 1] = tridiagonal[rowLength * rowLength - 1];
    
    // Theta is indexed from 0 as per the definition of an inverse tridiagonal
    // Meanwhile, a, b, and c are indexed from 0 in contrast to the math's 
    //      starting index of 1, so we subtract 1 from their formula indices
    double theta[rowLength + 1];
    theta[0] = 1;
    theta[1] = a[0];
    for (int i = 2; i < rowLength + 1; i++) {
        theta[i] = a[i - 1] * theta[i - 1] - b[i - 2] * c[i - 2] * theta[i - 2];
    }
    
    // Phi is indexed from 1 as per the definition of an inverse tridiagonal
    double phi[rowLength + 2];
    phi[rowLength + 1] = 1;
    phi[rowLength] = a[rowLength - 1];
    for (int i = rowLength - 1; i > 0; i--) {
        phi[i] = a[i - 1] * phi[i + 1] - b[i - 1] * c[i - 1] * phi[i + 2];
    }
    
    // Calculate each value in the final inverse matrix
    for (int i = 0; i < rowLength; i++) {
        for (int j = 0; j < rowLength; j++) {
            
            // Calculates values below the diagonal
            if (i < j) {
                inverse[i * rowLength + j] = pow(-1, i + j);
                for (int x = i; x < j; x++) {
                    inverse[i * rowLength + j] *= b[x];
                }
                inverse[i * rowLength + j] *= theta[i] * phi[j + 2] / theta[rowLength];
            }
            
            // Calculates values on the diagonal
            else if (i == j) {
                inverse[i * rowLength + j] = theta[i] * phi[j + 2] / theta[rowLength];
            }
            
            // Calculates values above the diagonal
            else if (i > j) {
                inverse[i * rowLength + j] = pow(-1, i + j);
                for (int x = j; x < i; x++) {
                    inverse[i * rowLength + j] *= c[x];
                }
                inverse[i * rowLength + j] *= theta[j] * phi[i + 2] / theta[rowLength];
            }
        }
    }
}