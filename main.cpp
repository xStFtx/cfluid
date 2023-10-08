#include <iostream>
#include <vector>
#include <cmath>

const int Lx = 100;          // Grid size in the x-direction
const int Ly = 50;           // Grid size in the y-direction
const int Q = 9;             // Number of velocity directions
const double tau = 0.6;      // Relaxation time
const double density_in = 1.0;
const double density_out = 0.0;
const double initial_velocity = 0.1;

class FluidSimulation {
public:
    FluidSimulation() {
        initialize();
    }

    // Initialize the simulation
    void initialize() {
        // Initialize density, velocity, and distribution functions
        density.assign(Lx, std::vector<double>(Ly, density_in));
        ux.assign(Lx, std::vector<double>(Ly, 0.0));
        uy.assign(Lx, std::vector<double>(Ly, 0.0));
        f.assign(Lx, std::vector<std::vector<double>>(Ly, std::vector<double>(Q, 0.0)));

        // Set initial velocity in the x-direction
        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                ux[x][y] = initial_velocity;
            }
        }
    }

    // Perform a single time step of the simulation
    void step() {
        // Calculate macroscopic variables (density and velocity)
        calculateMacroscopicVariables();

        // Collision step
        collisionStep();

        // Streaming step
        streamingStep();
    }

    // Visualize the simulation results (simplified, does not provide actual visualization)
    void visualize() {
        // Print or visualize the simulation results
    }

private:
    // Grid data (density, velocity components, and distribution functions)
    std::vector<std::vector<double>> density;
    std::vector<std::vector<double>> ux;
    std::vector<std::vector<double>> uy;
    std::vector<std::vector<std::vector<double>>> f;

    // Helper functions for collision and streaming steps
    void calculateMacroscopicVariables() {
        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                density[x][y] = 0.0;
                ux[x][y] = 0.0;
                uy[x][y] = 0.0;
                for (int k = 0; k < Q; k++) {
                    density[x][y] += f[x][y][k];
                    ux[x][y] += f[x][y][k] * c[k][0];
                    uy[x][y] += f[x][y][k] * c[k][1];
                }
                ux[x][y] /= density[x][y];
                uy[x][y] /= density[x][y];
            }
        }
    }

    void collisionStep() {
        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                double local_density = density[x][y];
                double usqr = ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y];
                for (int k = 0; k < Q; k++) {
                    double cu = c[k][0] * ux[x][y] + c[k][1] * uy[x][y];
                    double feq_val = feq(k, local_density, ux[x][y], uy[x][y], usqr);
                    f[x][y][k] = f[x][y][k] - (f[x][y][k] - feq_val) / tau;
                }
            }
        }
    }

    void streamingStep() {
        std::vector<std::vector<std::vector<double>>> f_new(Lx, std::vector<std::vector<double>>(Ly, std::vector<double>(Q, 0.0)));
        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                for (int k = 0; k < Q; k++) {
                    int new_x = (x + c[k][0] + Lx) % Lx;
                    int new_y = (y + c[k][1] + Ly) % Ly;
                    f_new[new_x][new_y][k] = f[x][y][k];
                }
            }
        }
        f = f_new;
    }

    double feq(int k, double rho, double u_x, double u_y, double usqr) {
        double cu = c[k][0] * u_x + c[k][1] * u_y;
        return w[k] * rho * (1.0 + 3.0 * cu + 9.0 / 2.0 * cu * cu - 3.0 / 2.0 * usqr);
    }

    // Lattice weights and directions (D2Q9 model)
    double w[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    int c[Q][2] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};
};

int main() {
    FluidSimulation simulation;

    // Main time-stepping loop
    for (int t = 0; t < 1000; t++) {
        simulation.step();
    }

    // Visualize the results
    simulation.visualize();

    return 0;
}
