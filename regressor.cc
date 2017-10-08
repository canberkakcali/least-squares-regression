#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd; // Dynamic double matrix

#define MAX_DEGREE 10 // Don't change unless you know what you are doing.

// A simple structure to represent points
typedef struct { 
	double x;
	double y;
} point;


// Returns the leftmost point's x value
double findMinXPoint(vector<point> v) {
	double min = v.at(0).x;
	for(int i=1;i<v.size();i++) {
		if(v.at(i).x < min) min = v.at(i).x;
	}
	return min;
}

// Returns the rightmost point's x value
double findMaxXPoint(vector<point> v) {
	double max = v.at(0).x;
	for(int i=1;i<v.size();i++) {
		if(v.at(i).x > max) max = v.at(i).x;
	}
	return max;
}

void printUsage(char* cmd) {
	cout << "Usage: " << cmd << " [-d num] [--qr] <filename>" << endl <<
	"\t-d" << endl << "\t\tSpecifies the degree of the polynomial. The default value is 2." << endl <<
					  "\t\tNote that the maximum degree supported by this program is " << MAX_DEGREE << "." << endl <<
	"\t--qr" << endl << "\t\tUses the QR decomposition instead of the LU decomposition which is the default." << endl << endl <<
	"Examples: " << endl <<
	"\tTo find a 4th degree polynomial using the QR decomposition execute:" << endl <<
	"\t\t" << cmd << " -d 4 --qr input.txt" << endl <<
	"\tTo find a 2nd degree polynomial using the LU decomposition execute:" << endl <<
	"\t\t" << cmd << " input.txt" << endl;
}

int main(int argc, char** argv) {
		
		bool flagQR = false; // Flag to remember --qr argument
		int degree = 2; // Degree of the polynomial
		char* fileName; // File to read from
		
		/* 
		* 		PARAMETER CHECK
		*/
		
		
		if(argc == 1) { // If there are no arguments other than the program name, the input file is missing.
			cout << "Input file is required." << endl;
			printUsage(argv[0]);
			return EXIT_FAILURE;
		} else {
			for(int i=1;i<argc;i++) { // Read all the arguments passed.
				if(argv[i][0] == '-') { // Check if the argument is an option
					if(strcmp(argv[i], "--qr") == 0) { // Check if the option is "--qr"
						flagQR = true;

					} else if(argv[i][1] == 'd') { // Check if the option is "-d"
						char* p;
						degree = strtol(argv[i+1], &p, 10); // Convert the next argument to number
						if(*p) { // Check if the argument is numeric
							cout << "-d option must be followed by a number. You entered: " << argv[i+1] << endl;
							printUsage(argv[0]);
							return EXIT_FAILURE;
						} 
						if(degree <  0) {
							cout << "Degree cannot be a negative number." << endl;
							printUsage(argv[0]);
							return EXIT_FAILURE;
						}
						if(degree > MAX_DEGREE) {
							cout << "Maximum degree supported by this program is " << MAX_DEGREE << "." << endl <<
									"Degree is set to " << MAX_DEGREE << "." << endl;
							degree = MAX_DEGREE;
						}
						i++; // Skip the next argument (num)

					} else {
						cout << "Unrecognized option " << argv[i] << endl;
						printUsage(argv[0]);
						return EXIT_FAILURE;
					}
				} else { // If it is not "--qr", then it must be the filename
					fileName = argv[i];	
				}
			}
		}
		
		/* 
		* 		READ FILE
		*/
		
		ifstream infile(fileName); // open file stream for the input file 
		if(infile.fail()) { // Ensure there are no failures
			cout << "Error reading file: " << fileName << endl;
			return EXIT_FAILURE;
		}
		
		vector<point> points; // A dynamic list to store read points
		point buffer; // We will read every line and parse the numbers into a point storage
		
		while(infile >> buffer.x >> buffer.y) { // As long as there are input, parse the inputs into point buffer
			points.push_back(buffer); // Push the point into the dynamic list
		}
		infile.close(); // Reading is over, close the file.
		
		
		/*
		*		INITIALIZE THE NECESSARY MATRICES
		*/
		
		
		if(points.size() < 2) {
			cout << fileName << " has not enough points." << endl
				 << "Please input at least 2 points." << endl;
			return EXIT_FAILURE;
		}
		
		MatrixXd A(points.size(), degree+1); // for a 2nd degree polynomial: x^2  x^1  x^0  matrix => 3 columns, N rows => (N,3) matrix
		
		for(int i=0;i<points.size();i++) { // parse every point into the matrix

			for(int j=0; j<=degree; j++) {
				A(i,j) = pow( points.at(i).x , degree-j );
			}

		}
		cout << "X Matrix: " << endl << A << endl << endl; // Print the matrix
		
		
		MatrixXd v(points.size(), 1); // parse the y values of the points into the Result Vector => (N,1) matrix
        	for(int i=0;i<points.size();i++) {
        		v(i,0) = points.at(i).y;
		}
        	cout << "Result Vector: " << endl << v << endl << endl; // Print the vector
        
        
		/*
		*		CALCULATION
		*/


		MatrixXd Solution(degree+1,1); // A vector matrix to store a,b,c... solution

		if(flagQR) { // If the "--qr" parameter is passed, use QR decomposition
			
			Solution = A.fullPivHouseholderQr().solve(v);
			cout << "Linear Least Square System Solution (using QR decomposition): " << endl;
			
		} else { // Otherwise, use the LU decomposition
		
			Solution = A.fullPivLu().solve(v);
			cout << "Linear Least Square System Solution (using LU decomposition): " << endl;
			
		}
		
		
		/*
		*		GENERATE THE OUTPUT
		*/
		
		
		cout << "Number of points = " << points.size() << endl;
		
		for(int i=0; i<=degree; i++) {
			cout << (char)('a'+i) << " = " << Solution(i,0) << endl;
		}
		
		// Final function:
		cout << "y = (" << Solution(0,0) << ")x^" << degree;
		for(int i=1; i<=degree; i++) {
			cout << " + (" << Solution(i,0) << ")x^" << (degree-i);
		}

		// Error
		cout << endl << "The error is: " << (A*Solution - v).norm() << endl << endl; // SUM OF (yi - f(xi))^2 FOR ALL i: 1 to N   where yi = given y and  f(x) = generated final function 

			/* << "a = " << Solution(0,0) << endl
        	 << "b = " << Solution(1,0) << endl
        	 << "c = " << Solution(2,0) << endl
			 << "y = (" << Solution(0,0) << ")x^2 + (" << Solution(1,0) << ")x + (" << Solution(2,0) << ")" << endl // final function
			 << "The error is: " << (A*Solution - v).norm() << endl << endl; // SUM OF (yi - f(xi))^2 FOR ALL i: 1 to N   where yi = given y and  f(x) = generated final function
		*/
		
		ofstream outfile("graph.gp"); // Create a GNUPlot File
		
		if(outfile.fail()) { // Ensure there are no failures
			cout << "Error writing file graph.gp" << endl;
			
		} else {
			
			/* 
			*  Set the graph range, based on leftmost x value and rightmost x value of the points in the given input file.
			*  Example: (4,5) , (9,15) , (6,8) , (5,7)
			*  Leftmost point = (4,5)
			*  Rightmost point = (9,15)
			*  Range will be = [4-1 , 9+1] = [3,10]
			*/
			outfile << "set xrange [" << (int)(findMinXPoint(points) - 1) << ":" << (int)(findMaxXPoint(points) + 1) << "]" << endl;
			
			
			// Create a marker for all of the given points.
			for(int i=0; i<points.size(); i++) {
				outfile << "set label " << i+1 << " \"(" << points.at(i).x << ", " << points.at(i).y << ")\"  at " << points.at(i).x << "," << points.at(i).y << " point pointtype 1" << endl;
			}
			
			// Plot the function
			//outfile	<< "plot (" << Solution(0,0) << ")*x*x + (" << Solution(1,0) << ")*x + (" << Solution(2,0) << ")" << endl;
			outfile << "plot";
			for(int i=0; i<=degree; i++) {
				outfile << " +(" << Solution(i,0) << ")*x**" << (degree-i);
			}
			outfile << endl;


			// Inform the user
			cout 	<< "*** Graph file with the name graph.gp is created at current directory." << endl << endl;
			
			// Writing is done, close the file.
			outfile.close();
			
			// Draw the plot.
			system("gnuplot -p graph.gp");
		}
		
		
        return EXIT_SUCCESS;
}
