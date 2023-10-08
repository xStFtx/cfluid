#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <omp.h>
#include <SFML/Graphics.hpp>

const int Lx = 100;
const int Ly = 50;
const int Q = 9;
const double tau = 0.6;
const double density_in = 1.0;
const double density_out = 0.0;
const double initial_velocity = 0.1;

class FluidSimulation {
public:
    FluidSimulation() {
        initialize();
    }

    void initialize() {
        density.assign(Lx, std::vector<double>(Ly, density_in));
        ux.assign(Lx, std::vector<double>(Ly, 0.0));
        uy.assign(Lx, std::vector<double>(Ly, 0.0));
        f.assign(Lx, std::vector<std::vector<double>>(Ly, std::vector<double>(Q, 0.0)));

        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                ux[x][y] = initial_velocity;
            }
        }
    }

    void step() {
        calculateMacroscopicVariables();
        applyBoundaryConditions();
        collisionStep();
        streamingStep();
    }

    void visualize(sf::RenderWindow& window) {
        // Render the density field and boundary conditions
        for (int x = 0; x < Lx; x++) {
            for (int y = 0; y < Ly; y++) {
                sf::RectangleShape pixel(sf::Vector2f(4.0f, 4.0f));
                pixel.setPosition(x * 4.0f, y * 4.0f);

                // Set color based on density
                int colorValue = static_cast<int>(density[x][y] * 255);
                sf::Color color(colorValue, colorValue, colorValue);

                // Check if it's a boundary cell and set a distinct color
                if (isBoundaryCell(x, y)) {
                    color = sf::Color::Red; // Change to any color you like for boundary cells
                }

                pixel.setFillColor(color);
                window.draw(pixel);
            }
        }
    }

private:
    std::vector<std::vector<double>> density;
    std::vector<std::vector<double>> ux;
    std::vector<std::vector<double>> uy;
    std::vector<std::vector<std::vector<double>>> f;

    // Constants for the lattice Boltzmann method
    const double w[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    const int c[Q][2] = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}, {1, 1}, {-1, 1}, {-1, -1}, {1, -1}};

    void calculateMacroscopicVariables() {
        #pragma omp parallel for collapse(2)
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

    void applyBoundaryConditions() {
        // Implement boundary conditions here
        // For simplicity, set no-slip boundary conditions at walls
        for (int x = 0; x < Lx; x++) {
            ux[x][0] = 0.0;
            uy[x][0] = 0.0;
            ux[x][Ly - 1] = 0.0;
            uy[x][Ly - 1] = 0.0;
        }
        for (int y = 0; y < Ly; y++) {
            ux[0][y] = 0.0;
            uy[0][y] = 0.0;
            ux[Lx - 1][y] = 0.0;
            uy[Lx - 1][y] = 0.0;
        }
    }

    void collisionStep() {
        #pragma omp parallel for collapse(2)
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
        #pragma omp parallel for collapse(3)
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

    bool isBoundaryCell(int x, int y) {
        // Check if (x, y) is on the boundary, e.g., if it's on the edges of the domain
        return x == 0 || x == Lx - 1 || y == 0 || y == Ly - 1;
    }
};

int main() {
    std::cout << "Running..." << std::endl;
    FluidSimulation simulation;

    sf::RenderWindow window(sf::VideoMode(Lx * 4, Ly * 4), "Fluid Simulation");
    window.setFramerateLimit(60);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        simulation.step();

        window.clear();
        simulation.visualize(window);
        window.display();
    }

    return 0;
}
