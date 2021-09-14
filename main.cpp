#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>

using namespace::std;

double k(int i, int xy_grid, int xy_interval);
double potential(double x, double y);
double test_potential(double x, double y);
double Psi0(double x, double y);
double test_Psi0(double x, double y);
int index(int i, int j, int N);
void time_evolution(char* filename, int xy_interval, double dxy, int xy_grid, int time_grid, double dt,
                    double kx, double ky, double norm, int x01, int y00, double V, double*** Psi);
double test_potential_harm(double x, double y, int i, int j, int xy_grid, int xy_interval);

int main() {
    // Initialize variables etc.
    char filename[9] = "SSSS.csv";
    const int xy_interval = 20;
    // dxy must be smaller or equal to 0.1 and there must exist an integer z so that z * dxy = 0.1 for indices to be
    // able to resolve x = 0.1 and y = 0.0. It is also assumed that dx = dy!
    const double dxy = 0.1;
    // xy_grid must be of size (xy_interval / dxy + 1) to be symmetric and
    // include endpoints of the interval -10 <= x,y <= 10
    const int xy_grid = xy_interval / dxy + 1;
    // dt = 0.1 seems to be a good time step
    const double dt = 0.1;
    const int time_grid = 100 / dt;
    double kx = {};
    double ky = {};
    double norm =  xy_grid * xy_grid;
    // Indices for x = 0.1 and y = 0.0 when dxy = 0.1
    int x01 = {xy_grid / 2 + int(0.1 / dxy)};
    int y00 = {xy_grid / 2};
    // Initialize potential and wave function array
    double V = {};
    // Allocate memory for the wave function Psi[xy_grid][xy_grid][2]
    auto*** Psi = (double***)malloc(xy_grid * sizeof(double**));
    for (int i = 0; i < xy_grid; i++) {
        Psi[i] = (double**)malloc(xy_grid * sizeof(double*));
            for (int j = 0; j < xy_grid; j++) {
                Psi[i][j] = (double *) malloc(2 * sizeof(double));
            }
    }
    cout << "Simulate time evolution of Psi(x = " << -10 + dxy * x01 << ", y = " << -10 + dxy * y00 << ")" << endl;
    // Run the simulation
    time_evolution(filename, xy_interval, dxy, xy_grid, time_grid, dt, kx, ky, norm, x01, y00, V, Psi);
    return 0;
}

// Function to create wave vector k
double k(int i, int xy_grid, int xy_interval) {
    double k;
    if (i < (xy_grid / 2 + 1)) {
        k = i * (2 * M_PI) / xy_interval;
    }
    else  {
        k = 2 * M_PI  / xy_interval * (i -  xy_grid);
    }
    return k;
}

// Function to create the potential V(x,y)
double potential(double x, double y) {
    double V;
    V = - 5 * pow( (1 + pow((x / 5), 2) + pow((y / 4), 2) ), - 4);
    return V;
}

// Function to create a test_potential V(x,y) = 0
double test_potential(double x, double y) {
    double V;
    V = 0.0;
    return V;
}

// Function to create a test_potential V(x,y) = 0
double test_potential_harm(double x, double y, int i, int j, int xy_grid, int xy_interval) {
    double V;
    V = k(i, xy_grid, xy_interval) * k(j, xy_grid, xy_interval) * (pow(x, 2) + pow(y, 2)) / 2;
    return V;
}

// Function to convert matrix indices to array indices (row)
int index(int i, int j, int N) {
    return i + N * j;
}

// Function to create Psi(t=0)
double Psi0(double x, double y) {
    double Psi;
    Psi = (1 / M_PI) * exp(- ( pow((x - 1), 2) + pow((y - 1), 2) ) );
    return Psi;
}

// Function to create test wave function
double test_Psi0(double x, double y) {
    double Psi;
    Psi = sin(x);
    return Psi;
}

// Function to perform the time evolution of the spectral split step solver
void time_evolution(char* filename, int xy_interval, double dxy, int xy_grid, int time_grid, double dt,
                    double kx, double ky, double norm, int x01, int y00, double V, double*** Psi) {
    // Loop to create wave function (and k according to 0, ... n, -n, ... -1)
    for (int i = 0; i < xy_grid; i++) {
        for (int j = 0; j < xy_grid; j++) {
            // Change Psi0 to test_Psi0 for a simple sine function
            Psi[i][j][0] = test_Psi0(- 10 + dxy * i, -10 + dxy * j);
            Psi[i][j][1] = 0;
        }
    }

    // Create fftw plans and arrays
    auto *forward = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xy_grid * xy_grid));
    auto *backward = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xy_grid * xy_grid));
    fftw_plan forward_plan = fftw_plan_dft_2d(xy_grid, xy_grid, forward, forward, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_2d(xy_grid, xy_grid, backward, backward, FFTW_BACKWARD, FFTW_ESTIMATE);
    double ecos;
    double esin;
    // Prepare output file with header and wave function at t = 0
    ofstream file(filename);
    file << "time, Real part, Imaginary part, time_grid = " << time_grid << endl;
    file << "0.0" << "," << Psi[x01][y00][0] << "," << Psi[x01][y00][1] << endl;

    // Evolve Psi in time
    for (int T = 1; T < time_grid; T++) {
        // Create exp(-i * t/2 * V(0)) * Psi(0) (use trigonometric functions)
        for (int i = 0; i < xy_grid; i++) {
            for (int j = 0; j < xy_grid; j ++) {
                // Change potential to test_potential for V = 0. This needs to be changed in the last step as well
                V = test_potential_harm(- 10 + dxy * i, -10 + dxy * j, i, j, xy_grid, xy_interval);
                ecos = cos(dt / 2 * V);
                esin = -sin(dt / 2 * V);
                forward[index(i, j, xy_grid)][0] = ecos * Psi[i][j][0] - esin * Psi[i][j][1];
                forward[index(i, j, xy_grid)][1] = ecos * Psi[i][j][1] + esin * Psi[i][j][0];
            }
        }
        // Forward Fourier transform
        fftw_execute(forward_plan);
        // Create exp(-i * k^2 * t) * output from above (use trigonometric functions)
        for (int i = 0; i < xy_grid; i++) {
            for (int j = 0; j < xy_grid; j ++) {
                kx = k(i, xy_grid, xy_interval);
                ky = k(j, xy_grid, xy_interval);
                ecos = cos(dt * (pow(kx, 2) + pow(ky, 2)));
                esin = -sin(dt * (pow(kx, 2) + pow(ky, 2)));
                backward[index(i, j, xy_grid)][0] =
                        ecos * forward[index(i, j, xy_grid)][0] - esin * forward[index(i, j, xy_grid)][1];
                backward[index(i, j, xy_grid)][1] =
                        ecos * forward[index(i, j, xy_grid)][1] + esin * forward[index(i, j, xy_grid)][0];
            }
        }
        // Backward Fourier transform
        fftw_execute(backward_plan);
        // Create exp(-i * t / 2 * V) * output from backward Fourier transform (use trigonometric functions)
        for (int i = 0; i < xy_grid; i++) {
            for (int j = 0; j < xy_grid; j ++) {
                // Change potential to test_potential for V = 0. This needs to be changed in the first step as well
                V = test_potential_harm(- 10 + dxy * i, -10 + dxy * j, i, j, xy_grid, xy_interval);
                ecos = cos(dt / 2 * V);
                esin = -sin(dt / 2 * V);
                Psi[i][j][0] = (ecos * backward[index(i, j, xy_grid)][0]- esin * backward[index(i, j, xy_grid)][1])
                               / norm;
                Psi[i][j][1] = (ecos * backward[index(i, j, xy_grid)][1] + esin * backward[index(i, j, xy_grid)][0])
                               / norm;
            }
        }
        // Write the temporal spectrum of Psi(0.1,0) to file
        file << dt * T << "," << Psi[x01][y00][0] << "," << Psi[x01][y00][1] << endl;
    }
    // Clear variables
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(forward);
    fftw_free(backward);
    free(Psi);
}
