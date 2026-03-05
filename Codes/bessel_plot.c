/*
 * File: bessel_mult.c
 * Purpose: Computes and graphs Bessel functions from J_0(x) to J_6(x) using gnuplot.
 */

#include <stdio.h>  /* Includes standard I/O library to use the FILE type and printf/fprintf. */
#include <stdlib.h> /* Includes standard library to use the popen and pclose process functions. */
#include <math.h>   /* Includes the math library to access the j0, j1, and jn functions. */

int main(void) {    /* Defines the main entry point of the program, taking no arguments. */
    
    FILE *gnuplot_pipe; /* Declares a pointer variable to a FILE structure for the process pipe. */
    
    /* Opens a persistent pipe to the gnuplot process in write-only mode. */
    gnuplot_pipe = popen("gnuplot -persistent", "w"); 
    
    if (gnuplot_pipe == NULL) { /* Evaluates if the operating system failed to create the pipe. */
        printf("Error: Could not open pipe.\n"); /* Prints an error to standard output. */
        return 1; /* Returns an error code to the operating system, terminating the program. */
    } /* Closes the conditional block. */

    /* Sends a string to gnuplot to set the title of the graph. */
    fprintf(gnuplot_pipe, "set title 'Bessel Functions J_0(x) to J_6(x)'\n"); 
    /* Sends a string to gnuplot to label the x-axis. */
    fprintf(gnuplot_pipe, "set xlabel 'x'\n"); 
    /* Sends a string to gnuplot to label the y-axis. */
    fprintf(gnuplot_pipe, "set ylabel 'J_n(x)'\n"); 
    /* Instructs gnuplot to render a background grid for easier reading. */
    fprintf(gnuplot_pipe, "set grid\n"); 

    /* Constructs a single, multi-line gnuplot command to expect 7 distinct datasets inline. */
    fprintf(gnuplot_pipe, "plot '-' with lines title 'J_0(x)', " /* First dataset definition. */
                          "'-' with lines title 'J_1(x)', " /* Second dataset definition. */
                          "'-' with lines title 'J_2(x)', " /* Third dataset definition. */
                          "'-' with lines title 'J_3(x)', " /* Fourth dataset definition. */
                          "'-' with lines title 'J_4(x)', " /* Fifth dataset definition. */
                          "'-' with lines title 'J_5(x)', " /* Sixth dataset definition. */
                          "'-' with lines title 'J_6(x)'\n"); /* Seventh dataset and newline execution. */

    int n; /* Allocates 4 bytes on the stack for an integer variable to track the Bessel order. */
    double x; /* Allocates 8 bytes on the stack for a double variable to track the x-coordinate. */

    /* Initializes n to 0 and loops until n is 6, incrementing by 1 each iteration. */
    for (n = 0; n <= 6; n++) { 
        /* Initializes x to 0.0 and loops until x reaches 20.0, incrementing by 0.1. */
        for (x = 0.0; x <= 20.0; x += 0.1) { 
            /* Checks if the current order is 0 to use the specific j0 function. */
            if (n == 0) { 
                /* Computes J_0(x) and streams the x and y values to the gnuplot pipe. */
                fprintf(gnuplot_pipe, "%f %f\n", x, j0(x)); 
            } else if (n == 1) { /* Checks if the current order is 1 to use the specific j1 function. */
                /* Computes J_1(x) and streams the x and y values to the gnuplot pipe. */
                fprintf(gnuplot_pipe, "%f %f\n", x, j1(x)); 
            } else { /* Executes for all other values of n (2 through 6). */
                /* Computes J_n(x) using the general function and streams the values. */
                fprintf(gnuplot_pipe, "%f %f\n", x, jn(n, x)); 
            } /* Closes the conditional execution block. */
        } /* Closes the inner loop managing the x-coordinates. */
        
        /* Sends the 'e' character to notify gnuplot that the current dataset stream has ended. */
        fprintf(gnuplot_pipe, "e\n"); 
    } /* Closes the outer loop managing the Bessel function orders. */

    /* Instructs the operating system to close the pipe and free the associated resources. */
    pclose(gnuplot_pipe); 

    return 0; /* Returns a success code of 0 to the operating system. */
} /* Closes the main function block. */