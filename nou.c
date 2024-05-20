#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-10
#define MAX_ITER 1000

// Define the function f(x, y)
double f(double x, double y)
{
    return -4 * pow(y, 4) - 6 * pow(x, 2) * pow(y, 2) + 11 * x * pow(y, 2) + 3 * pow(y, 2) - 6 * pow(x, 4) + 21 * pow(x, 3) - 19 * pow(x, 2) + 0.2 * y + 0.1;
}

// Partial derivative of f with respect to x
double fx(double x, double y)
{
    return -12 * x * pow(y, 2) + 11 * pow(y, 2) - 24 * pow(x, 3) + 63 * pow(x, 2) - 38 * x;
}

// Partial derivative of f with respect to y
double fy(double x, double y)
{
    return -16 * pow(y, 3) - 12 * pow(x, 2) * y + 22 * x * y + 6 * y + 0.2;
}

// Newton method to find y for a given x such that f(x, y) = 0
double newton_method(double (*f)(double, double), double (*df)(double, double), double x0, double y_start)
{
    double y = y_start;
    while (fabs(f(x0, y)) > TOL)
    {
        y -= f(x0, y) / df(x0, y);
    }
    return y;
}

// Two dimensional Newton Method
double two_d_newton_method(double (*f)(double, double), double (*fx)(double, double),
                           double (*fy)(double, double), double *pred, double *prev, double h)
{
    // F is (f, (x-x0)^2 + (y-x1)^2 - h^2)
    // Create DF(x0, x1) to calc its inverse.
    double DF[4];
    DF[0] = fx(pred[0], pred[1]);
    DF[1] = fy(pred[0], pred[1]);
    DF[2] = 2 * (pred[0] - prev[0]);
    DF[3] = 2 * (pred[1] - prev[1]);
    
    double det = determinant(DF[0], DF[1], DF[2], DF[3]);
    // If the determinant is zero (or close enough), assume the matrix is not invertible
    if (fabs(det) < TOL) {
        printf("The determinant too close to zero!");
    }

    double inverse[4];
    inverse[0] = DF[3] / det;
    inverse[1] = -DF[1] / det;
    inverse[2] = -DF[2] / det;
    inverse[3] = DF[0] / det;

    // Update the current point
    double result[2];
    result[0] = f(pred[0], pred[1]);
    result[1] = pow((pred[0] - prev[0]), 2) + pow((pred[1] - prev[1]), 2) - h * h;
    // Here we do pred = pred - inverse(pred) * F(pred)
    pred[0] = pred[0] - (inverse[0] * result[0] + inverse[1] * result[1]);
    pred[1] = pred[1] - (inverse[2] * result[0] + inverse[3] * result[1]);

    // TODO: Make it iterable and stop when F(x0, x1) is close enough to 0
}

// Predictor-corrector method
void predictor_corrector(double (*func)(double, double), double h, int nc, double err, int nn, double *points)
{
    double x = 0.0, y = newton_method(func, fy, 0.0, 0.0); // Start at (0, y0)
    printf("Starting y: %e", y);
    points[0] = x;
    points[1] = y;

    for (int i = 1; i < nc; ++i)
    {
        double x_pred = x + h;
        double y_pred = y; // Initial prediction

        // Corrector step
        for (int j = 0; j < nn; ++j)
        {
            double f_val = func(x_pred, y_pred);
            double fy_val = fy(x_pred, y_pred);
            if (fabs(f_val) < err)
            {
                break;
            }
            y_pred = y_pred - f_val / fy_val;
        }

        if (fabs(func(x_pred, y_pred)) >= err)
        {
            printf("Failed to converge at step %d\n", i);
            break;
        }

        x = x_pred;
        y = y_pred;
        points[2 * i] = x;
        points[2 * i + 1] = y;

        // Check if we are close to the starting point (closed curve)
        if (i > 1 && fabs(x - points[0]) < h && fabs(y - points[1]) < h)
        {
            printf("Closed curve completed\n");
            break;
        }
    }
}

int main()
{
    double h, err;
    int nc, nn;

    // Example input values, these should be taken from user input
    h = 0.01;      // Step size
    nc = 160;      // Max number of points
    err = 1e-6;    // Maximum allowed error
    nn = 10000000; // Number of iterations in corrector method

    double *points = (double *)malloc(2 * nc * sizeof(double));
    if (points == NULL)
    {
        printf("Memory allocation failed\n");
        return 1;
    }

    predictor_corrector(f, h, nc, err, nn, points);

    // Write points to file
    FILE *file = fopen("curve_points.txt", "w");
    if (file == NULL)
    {
        printf("File opening failed\n");
        return 1;
    }

    for (int i = 0; i < nc; ++i)
    {
        fprintf(file, "%lf %lf\n", points[2 * i], points[2 * i + 1]);
    }

    fclose(file);
    free(points);

    return 0;
}
