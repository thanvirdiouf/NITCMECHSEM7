#include <iostream>
#include <fstream>
#include <cmath>

// Function to set the initial condition
double initial_condition(double x) {
    return 1 - x + std::sin(2 * M_PI * x);
}

int main() {
    // Parameters
    double alpha = 0.5;
    double L = 1.0;
    int Nx = 100;  // Number of spatial points
    int Nt = 1000;  // Number of time steps
    double dx = L / (Nx - 1);
    double dt = alpha * dx * dx / 2;  // Stability requirement

    // Boundary conditions
    double T_left = 1.0;
    double T_right = 0.0;

    // Initialization
    double x_values[Nx];
    double T[Nt][Nx];

    // Set initial condition
    for (int i = 0; i < Nx; ++i) {
        x_values[i] = i * dx;
        T[0][i] = initial_condition(x_values[i]);
    }

    // FTCS scheme
    for (int n = 0; n < Nt - 1; ++n) {
        for (int i = 1; i < Nx - 1; ++i) {
            T[n + 1][i] = T[n][i] + alpha * (T[n][i + 1] - 2 * T[n][i] + T[n][i - 1]);
        }

        // Boundary conditions
        T[n + 1][0] = T_left;
        T[n + 1][Nx - 1] = T_right;
    }

    // Write results to a CSV file
    std::ofstream outputFile("heat_conduction_results.csv");
    if (outputFile.is_open()) {
        for (int n = 0; n < Nt; ++n) {
            for (int i = 0; i < Nx; ++i) {
                outputFile << x_values[i] << "," << T[n][i];
                if (i < Nx - 1) {
                    outputFile << ",";
                }
            }
            outputFile << "\n";
        }
        outputFile.close();
        std::cout << "Results written to heat_conduction_results.csv" << std::endl;
    } else {
        std::cerr << "Unable to open the output file." << std::endl;
        return 1;
    }

    return 0;
}
