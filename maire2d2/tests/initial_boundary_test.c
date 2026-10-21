#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "../initial_bound_condition.h"

static void expect_near(double actual, double expected)
{
    assert(isfinite(actual));
    assert(fabs(actual - expected) <= 2e-12 * fmax(1.0, fabs(expected)));
}

static void test_coordinates(void)
{
    /* Nonzero origins must also work below the first physical node. */
    const struct {
        size_t node_x, node_y;
        double x, y;
    } samples[] = {
        {0, 0, -0.4, -0.9},
        {2, 3,  0.0,  0.0},
        {3, 4,  0.2,  0.3},
        {1, 3, -0.2,  0.0},
        {2, 2,  0.0, -0.3},
        {5, 7,  0.6,  1.2}
    };

    for (size_t i = 0; i < sizeof samples / sizeof samples[0]; ++i) {
        double coor[2] = {NAN, NAN};
        initializeNodeCoor(samples[i].node_x, samples[i].node_y,
                           0.2, 0.3, 2, 3, coor);
        expect_near(coor[0], samples[i].x);
        expect_near(coor[1], samples[i].y);

        coor[0] = coor[1] = NAN;
        boundarilizeCoor(samples[i].node_x, samples[i].node_y,
                        0.2, 0.3, 2, 3, coor);
        expect_near(coor[0], samples[i].x);
        expect_near(coor[1], samples[i].y);
    }
}

static void check_initial_energy(size_t nx, size_t ny,
                                 size_t start_x, size_t start_y)
{
    const size_t nodes_x = nx + start_x + 2;
    const size_t nodes_y = ny + start_y + 2;
    const double dx = 1.2 / (double)nx;
    const double dy = 1.2 / (double)ny;
    const double area = dx * dy;
    double coor[nodes_y][nodes_x][2];

    /* The initializer must not depend on uninitialized ghost vertices. */
    for (size_t y = 0; y < nodes_y; ++y)
        for (size_t x = 0; x < nodes_x; ++x)
            coor[y][x][0] = coor[y][x][1] = NAN;
    for (size_t y = start_y; y <= start_y + ny; ++y)
        for (size_t x = start_x; x <= start_x + nx; ++x)
            initializeNodeCoor(x, y, dx, dy, start_x, start_y, coor[y][x]);

    size_t heated_cells = 0;
    double total_energy = 0.0;
    for (size_t y = start_y; y < start_y + ny; ++y) {
        for (size_t x = start_x; x < start_x + nx; ++x) {
            const double volume = initializeCellVolume(dx, dy);
            const double density = initializeCellDensity();
            const double mass = initializeCellMass(volume, density);
            double momentum[2] = {NAN, NAN};
            initializeCellSpMomentum(momentum);
            const double pressure = initializeCellPresure(
                x, y, dx, dy, nodes_x, nodes_y, coor, 1.4);
            const double energy = initializeCellSpEnergy(
                density, momentum, pressure, 1.4);

            expect_near(volume, area);
            expect_near(density, 1.0);
            expect_near(mass, area);
            expect_near(momentum[0], 0.0);
            expect_near(momentum[1], 0.0);
            expect_near(pressure, 0.4 * density * energy);
            if (x == start_x && y == start_y) {
                /* Maire p. 109: the full stated energy goes in this cell. */
                expect_near(mass * energy, 0.244816);
                expect_near(pressure, 0.4 * 0.244816 / area);
                if (nx == 30 && ny == 30) {
                    expect_near(pressure, 61.204);
                    expect_near(energy, 153.01);
                }
            } else {
                expect_near(pressure, 1e-6);
                expect_near(energy, 2.5e-6);
            }
            heated_cells += pressure > 1e-6;
            total_energy += mass * energy;
        }
    }

    assert(heated_cells == 1);
    expect_near(total_energy,
                0.244816 + ((double)(nx * ny) - 1.0) * area * 2.5e-6);
    if (nx == 30 && ny == 30)
        expect_near(total_energy, 0.244819596);
}

static void test_initial_state(void)
{
    assert(SEDOV_NUM_CELLS == 30);
    expect_near(SEDOV_GAMMA, 1.4);
    expect_near(SEDOV_FINAL_TIME, 1.0);
    expect_near(SEDOV_DOMAIN_LENGTH, 1.2);
    expect_near(SEDOV_AMBIENT_DENSITY, 1.0);
    expect_near(SEDOV_AMBIENT_PRESSURE, 1e-6);
    expect_near(SEDOV_EXPLOSION_ENERGY, 0.244816);

    check_initial_energy(30, 30, 1, 1);
    check_initial_energy(60, 60, 1, 1);
    check_initial_energy(7, 11, 3, 2);
    check_initial_energy(1, 1, 1, 1);

    /* Specific total energy includes kinetic energy for nonzero velocity. */
    double momentum[2] = {4.0, -3.0};
    expect_near(initializeCellSpEnergy(2.0, momentum, 3.0, 1.4), 16.25);
    expect_near(boundarilizeDensity(), 1.0);
    expect_near(boundarilizePressure(), 1e-6);
    expect_near(boundarilizeVolume(0.2, 0.3), 0.06);
    expect_near(boundarilizeMass(0.2, 0.3), 0.06);
    expect_near(boundarilizeSpEnergy(1.4), 2.5e-6);
    boundarilizeSpMomentum(momentum);
    expect_near(momentum[0], 0.0);
    expect_near(momentum[1], 0.0);
}

static void test_cell_boundaries(void)
{
    enum { NX = 4, NY = 3 };
    double mass[NY + 2][NX + 2];
    double density[NY + 2][NX + 2];
    double pressure[NY + 2][NX + 2];
    double volume[NY + 2][NX + 2];
    double energy[NY + 2][NX + 2];
    double momentum[NY + 2][NX + 2][2];
    struct state {
        double mass, density, pressure, volume, energy, ux, uy;
    } original[NY + 2][NX + 2];

    for (size_t y = 0; y < NY + 2; ++y) {
        for (size_t x = 0; x < NX + 2; ++x) {
            mass[y][x] = density[y][x] = pressure[y][x] = NAN;
            volume[y][x] = energy[y][x] = NAN;
            momentum[y][x][0] = momentum[y][x][1] = NAN;
            if (x == 0 || x == NX + 1 || y == 0 || y == NY + 1)
                continue;
            /* A changing interior state, including both velocity components. */
            const double rho = 1.0 + 0.1 * (double)x + 0.2 * (double)y;
            const double v = 0.01 * (10.0 + (double)x + 2.0 * (double)y);
            const double p = 0.02 + 0.5 * (double)x + 0.7 * (double)y;
            const double ux = 0.3 * (double)x + 0.2 * (double)y;
            const double uy = -0.4 * (double)x + 0.1 * (double)y;
            const double e = p / (0.4 * rho) + 0.5 * (ux * ux + uy * uy);
            original[y][x] = (struct state){rho * v, rho, p, v, e, ux, uy};
            mass[y][x] = rho * v;
            density[y][x] = rho;
            pressure[y][x] = p;
            volume[y][x] = v;
            energy[y][x] = e;
            momentum[y][x][0] = ux;
            momentum[y][x][1] = uy;
        }
    }

    applySedovCellBoundaries(NX, NY, 0.2, 0.3, 1.4,
                            mass, density, pressure, volume, energy, momentum);

    for (size_t y = 0; y < NY + 2; ++y) {
        for (size_t x = 0; x < NX + 2; ++x) {
            struct state expected;
            if (x == NX + 1 || y == NY + 1) {
                expected = (struct state){0.06, 1.0, 1e-6, 0.06, 2.5e-6, 0.0, 0.0};
            } else {
                expected = original[y == 0 ? 1 : y][x == 0 ? 1 : x];
                if (x == 0)
                    expected.ux = -expected.ux;
                if (y == 0)
                    expected.uy = -expected.uy;
            }
            expect_near(mass[y][x], expected.mass);
            expect_near(density[y][x], expected.density);
            expect_near(pressure[y][x], expected.pressure);
            expect_near(volume[y][x], expected.volume);
            expect_near(energy[y][x], expected.energy);
            expect_near(momentum[y][x][0], expected.ux);
            expect_near(momentum[y][x][1], expected.uy);
            expect_near(density[y][x], mass[y][x] / volume[y][x]);
            const double kinetic = 0.5 * (momentum[y][x][0] * momentum[y][x][0] +
                                          momentum[y][x][1] * momentum[y][x][1]);
            expect_near(pressure[y][x], 0.4 * density[y][x] * (energy[y][x] - kinetic));
        }
    }

    /* Symmetry must follow the current solution on every application. */
    momentum[2][1][0] = 7.0;
    momentum[2][1][1] = -3.0;
    pressure[2][1] = 9.0;
    energy[2][1] = 9.0 / (0.4 * density[2][1]) + 29.0;
    applySedovCellBoundaries(NX, NY, 0.2, 0.3, 1.4,
                            mass, density, pressure, volume, energy, momentum);
    expect_near(momentum[2][0][0], -7.0);
    expect_near(momentum[2][0][1], -3.0);
    expect_near(pressure[2][0], 9.0);
    expect_near(energy[2][0], energy[2][1]);
}

static void test_node_boundaries(void)
{
    enum { NX = 4, NY = 3 };
    double coor[NY + 3][NX + 3][2];
    double original[NY + 3][NX + 3][2];
    for (size_t y = 0; y < NY + 3; ++y)
        for (size_t x = 0; x < NX + 3; ++x)
            coor[y][x][0] = coor[y][x][1] = NAN;

    /* A deformed mesh whose physical symmetry axes are still x=0 and y=0. */
    for (size_t y = 1; y <= NY + 1; ++y) {
        for (size_t x = 1; x <= NX + 1; ++x) {
            const double i = (double)x - 1.0;
            const double j = (double)y - 1.0;
            coor[y][x][0] = i * (0.2 + 0.02 * j) + 0.005 * i * i;
            coor[y][x][1] = j * (0.3 + 0.015 * i) + 0.003 * j * j;
        }
    }

    for (size_t pass = 0; pass < 2; ++pass) {
        for (size_t y = 1; y <= NY + 1; ++y) {
            for (size_t x = 1; x <= NX + 1; ++x) {
                coor[y][x][0] *= 1.1;
                coor[y][x][1] *= 0.9;
                original[y][x][0] = coor[y][x][0];
                original[y][x][1] = coor[y][x][1];
            }
        }

        applySedovNodeBoundaries(NX, NY, coor);
        for (size_t y = 1; y <= NY + 1; ++y)
            for (size_t x = 1; x <= NX + 1; ++x)
                for (size_t component = 0; component < 2; ++component)
                    expect_near(coor[y][x][component], original[y][x][component]);

        /* Check the reflection/extrapolation identities, including all corners. */
        for (size_t y = 0; y < NY + 3; ++y) {
            expect_near(coor[y][0][0] + coor[y][2][0], 0.0);
            expect_near(coor[y][0][1], coor[y][2][1]);
            for (size_t component = 0; component < 2; ++component)
                expect_near(coor[y][NX + 2][component] -
                            2.0 * coor[y][NX + 1][component] + coor[y][NX][component], 0.0);
        }
        for (size_t x = 0; x < NX + 3; ++x) {
            expect_near(coor[0][x][0], coor[2][x][0]);
            expect_near(coor[0][x][1] + coor[2][x][1], 0.0);
            for (size_t component = 0; component < 2; ++component)
                expect_near(coor[NY + 2][x][component] -
                            2.0 * coor[NY + 1][x][component] + coor[NY][x][component], 0.0);
        }
    }
}

static void test_nodal_symmetry(void)
{
    const struct {
        size_t x, y;
        double ux, uy;
    } samples[] = {
        {3, 2, 0.0,   0.0},
        {3, 5, 0.0, -11.0},
        {6, 2, 7.0,   0.0},
        {6, 5, 7.0, -11.0}
    };
    for (size_t i = 0; i < sizeof samples / sizeof samples[0]; ++i) {
        double velocity[2] = {7.0, -11.0};
        applySedovVelocityBoundary(samples[i].x, samples[i].y, 3, 2, velocity);
        expect_near(velocity[0], samples[i].ux);
        expect_near(velocity[1], samples[i].uy);
    }
}

int main(void)
{
    test_coordinates();
    test_initial_state();
    test_cell_boundaries();
    test_node_boundaries();
    test_nodal_symmetry();
    puts("initial_boundary_test: PASS (Sedov energy, EOS, coordinates, and symmetry)");
    return 0;
}
