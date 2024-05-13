#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 1000 // Maximum number of points on the curve

double g(double y) {
    return -4*y*y*y*y + 3*y*y + 0.2*y + 0.1;
}

double dg(double y) {
    return -16*y*y*y + 6*y + 0.2;
}

// void corrector(double h, int nc, double err, int nn) {
//     double x[MAX_POINTS];
//     int i;

//     // Initialization
//     x[0] = 0.0;

//     // Predictor-Corrector Method
//     for (i = 1; i < nc && i < MAX_POINTS; i++) {
//         x[i] = x[i - 1] + h; // Predictor step
//         int j;
//         for (j = 0; j < nn; j++) {
//             double f_xi = f(x[i]);
//             double df_xi = df(x[i]);
//             double dx = -f_xi / df_xi; // Corrector step
//             x[i] += dx;
//             if (fabs(dx) < err) break; // Stop iteration if error is within tolerance
//         }
//         if (i == nc - 1) break; // Reached maximum number of points
//     }

//     // Output the points
//     printf("Points on the curve:\n");
//     for (int j = 0; j < i; j++) {
//         printf("%f\n", x[j]);
//     }
// }

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Usage: %s <h> <nc> <err> <nn>\n", argv[0]);
        return 1;
    }

    double h, err;
    int nc, nn;

    // Convert command line arguments to appropriate data types
    h = atof(argv[1]);
    nc = atoi(argv[2]);
    err = atof(argv[3]);
    nn = atoi(argv[4]);

    int i;

    // Initialization
    double y = 0.0;

    // Find the starting point with the Newton method
    while (fabs(g(y)) > err) {
        y += -(g(y) / dg(y));
    }

    printf("Starting point: (0, %e).\n", y);

    return 0;
}
