#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define TOL 1e-10

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

// Function to calculate the determinant of a 2x2 matrix
double determinant_2x2(double matrix[4])
{
    return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

// Newton method to find y for a given x such that f(x, y) = 0
double newton_method(double (*f)(double, double), double (*df)(double, double), double x0, double y_start, double err)
{
    double y = y_start;
    while (fabs(f(x0, y)) > err)
    {
        y -= f(x0, y) / df(x0, y);
    }
    return y;
}

// Two dimensional Newton Method
double two_d_newton_method(double (*f)(double, double), double (*fx)(double, double),
                           double (*fy)(double, double), double *pred, double *prev, double h, double err, int nn)
{
    int iter = 0;
    double DF[4], inverse[4], result[2], det;

    while (iter < nn)
    {
        // Compute the Jacobian matrix DF
        DF[0] = fx(pred[0], pred[1]);
        DF[1] = fy(pred[0], pred[1]);
        DF[2] = 2 * (pred[0] - prev[0]);
        DF[3] = 2 * (pred[1] - prev[1]);

        // Compute the determinant of DF
        det = determinant_2x2(DF);
        if (fabs(det) < TOL)
        {
            printf("The determinant is too close to zero!\n");
            return -1;
        }

        // Compute the inverse of the Jacobian matrix DF
        inverse[0] = DF[3] / det;
        inverse[1] = -DF[1] / det;
        inverse[2] = -DF[2] / det;
        inverse[3] = DF[0] / det;

        // Compute the function F at the current prediction
        result[0] = f(pred[0], pred[1]);
        result[1] = pow((pred[0] - prev[0]), 2) + pow((pred[1] - prev[1]), 2) - h * h;

        // Check for convergence
        if (sqrt(result[0] * result[0] + result[1] * result[1]) < err)
        {
            return 0; // Converged successfully
        }

        // Update the current point using the Newton-Raphson step
        pred[0] = pred[0] - (inverse[0] * result[0] + inverse[1] * result[1]);
        pred[1] = pred[1] - (inverse[2] * result[0] + inverse[3] * result[1]);

        iter++;
    }

    printf("Maximum iterations (nn) reached without convergence.\n");
    return -1;
}

void tangent_vector(double (*fx)(double, double), double (*fy)(double, double), double x, double y, double *tangent_x, double *tangent_y)
{
    // Calculate the gradient
    *tangent_x = fx(x, y);
    *tangent_y = fy(x, y);

    // Calculate an orthogonal vector to the gradient
    double temp = *tangent_x;
    *tangent_x = -*tangent_y;
    *tangent_y = temp;

    // Compute the norm of the orthogonal vector
    double norm = sqrt(*tangent_x * *tangent_x + *tangent_y * *tangent_y);

    // Normalize the orthogonal vector to get the unit tangent vector
    *tangent_x = *tangent_x / norm;
    *tangent_y = *tangent_y / norm;
}

// Predictor-corrector method
void predictor_corrector(double (*f)(double, double),
                         double (*fx)(double, double),
                         double (*fy)(double, double),
                         double h, int nc, double err, int nn, double *points, bool verbose)
{
    // Use prev to store the previous point, pred to store the prediction
    double prev[2], pred[2];

    // The first point is found taking x = 0 and using Newton's method to find y
    prev[0] = 0.0;
    prev[1] = newton_method(f, fy, 0.0, 0.0, err);
    if (verbose)
    {
        printf("Starting point: (%e, %e)\n", prev[0], prev[1]);
    }

    // Save the first point to the points array
    points[0] = prev[0];
    points[1] = prev[1];

    double tangent_x, tangent_y;

    for (int i = 1; i < nc; ++i)
    {
        // Find the tangent vector at (x, y)
        tangent_vector(fx, fy, prev[0], prev[1], &tangent_x, &tangent_y);
        if (verbose)
        {
            printf("The number %d unitary tangent vector is: (%e, %e)\n", i, tangent_x, tangent_y);
        }

        // Use the tangent vector to predict the next point
        pred[0] = prev[0] + h * tangent_x;
        pred[1] = prev[1] + h * tangent_y;
        if (verbose)
        {
            printf("Initial prediction: (%e, %e)\n", pred[0], pred[1]);
        }

        // Use the two dimensional Newton method to correct the prediction
        two_d_newton_method(f, fx, fy, pred, prev, h, err, nn);
        if (verbose)
        {
            printf("Adjusted prediction: (%e, %e)\n", pred[0], pred[1]);
        }

        // Update the previous prediction for the next iteration
        prev[0] = pred[0];
        prev[1] = pred[1];

        // Save the predicted point to the points array
        points[2 * i] = pred[0];
        points[2 * i + 1] = pred[1];

        // Check if we are close to the starting point (closed curve)
        if (i > 1 && sqrt(pow(pred[0] - points[0], 2) + pow(pred[1] - points[1], 2)) < h)
        {
            printf("Closed curve completed\n");
            break;
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("Usage: %s <h> <nc> <err> <nn>\n", argv[0]);
        return 1;
    }

    double h, err;
    int nc, nn;

    // Convert command line arguments to appropriate data types
    h = atof(argv[1]);   // Step size
    nc = atoi(argv[2]);  // Max number of points
    err = atof(argv[3]); // Maximum allowed error
    nn = atoi(argv[4]);  // Number of iterations in corrector method

    if (h <= 0 || nc <= 0 || err <= 0 || nn <= 0)
    {
        printf("All input parameters must be positive.\n");
        return 1;
    }

    double *points = (double *)malloc(2 * nc * sizeof(double));
    if (points == NULL)
    {
        printf("Memory allocation failed\n");
        return 1;
    }

    predictor_corrector(f, fx, fy, h, nc, err, nn, points, false);

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
