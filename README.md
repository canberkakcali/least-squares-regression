# least-squares-regression
Linear Least Squares Regression with  Eigen Dense Solver using QR decomposition and LU decomposition.

### Usage:

    regressor [-d num] [--qr] <filename>

**-d**

    Specifies the degree of the polynomial. The default value is 2.
  
**--qr**

    Uses QR decomposition instead of LU decomposition. LU decomposition is the default behaviour.
  
**filename**

    File to read the points from.
  
  
  
### Examples:

  To find a 4th degree polynomial using the QR decomposition execute:
  
      `regressor -d 4 --qr file.txt`
  
  
  To find a 2nd degree polynomial using the LU decomposition execute:
  
      `regressor -d 2 file.txt`
    
  or
  
      `regressor file.txt`
    
