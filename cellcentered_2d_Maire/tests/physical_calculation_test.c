#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "../physical_calculation.h"

static void expectNear(double actual, double expected)
{
    assert(isfinite(actual));
    assert(fabs(actual - expected) <= 3e-12 * fmax(1.0, fabs(expected)));
}

static double quadraticForm(double matrix[2][2], double vector[2])
{
    double product[2];
    matMulVec2(matrix, vector, product);
    return vecMulVec2(vector, product);
}

static void testRectangularCornersAndDispatch(void)
{
    const double width = 3.0;
    const double height = 0.8;
    const double impedance = 2.5;
    double q[4][2] = {{0.0, 0.0}, {width, 0.0}, {width, height}, {0.0, height}};

    for (size_t corner = 0; corner < 4; ++corner) {
        const size_t next = (corner + 1) % 4;
        const size_t prev = (corner + 3) % 4;
        double matrix[2][2], selected[2][2], old_glace[2][2], normal[2];
        const double length = calCornerLength(q[corner], q[next], q[prev]);
        calCornerNormal(q[corner], q[next], q[prev], normal);

        calMatCorner_eucclhyd(q[corner], q[next], q[prev], impedance, matrix);
        /* Each face contributes half its length to its normal direction. */
        expectNear(matrix[0][0], 0.5 * impedance * height);
        expectNear(matrix[1][1], 0.5 * impedance * width);
        expectNear(matrix[0][1], 0.0);
        expectNear(matrix[1][0], 0.0);

        calMatCorner(CORNER_SCHEME_EUCCLHYD, q[corner], q[next], q[prev], impedance, selected);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 2; ++j)
                expectNear(selected[i][j], matrix[i][j]);

        calMatCorner_glace(length, normal, impedance, old_glace);
        calMatCorner(CORNER_SCHEME_GLACE, q[corner], q[next], q[prev], impedance, selected);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 2; ++j)
                expectNear(selected[i][j], old_glace[i][j]);
    }
}

static void testSkewCornerAndRotation(void)
{
    double p[2] = {0.0, 0.0};
    double next[2] = {3.0, 0.0};
    double prev[2] = {1.0, 2.0};
    const double impedance = 1.7;
    double matrix[2][2];
    calMatCorner_eucclhyd(p, next, prev, impedance, matrix);

    expectNear(matrix[0][1], matrix[1][0]);
    expectNear(matrix[0][0] + matrix[1][1], 0.5 * impedance * (3.0 + sqrt(5.0)));
    const double determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    /* det(M)=z^2*(3/2)*(sqrt(5)/2)*sin(theta)^2, sin(theta)^2=4/5. */
    expectNear(determinant, 3.0 * impedance * impedance / sqrt(5.0));
    assert(matrix[0][0] > 0.0 && determinant > 0.0);
    double directions[3][2] = {{1.0, 0.0}, {0.0, 1.0}, {2.0, -3.0}};
    for (size_t i = 0; i < 3; ++i)
        assert(quadraticForm(matrix, directions[i]) > 0.0);

    /* A rigid rotation and translation must rotate the tensor, not change it. */
    const double rotation[2][2] = {{0.6, -0.8}, {0.8, 0.6}};
    const double translation[2] = {7.25, -3.5};
    double rotated[3][2], rotated_matrix[2][2];
    double *vertices[3] = {p, next, prev};
    for (size_t vertex = 0; vertex < 3; ++vertex)
        for (size_t i = 0; i < 2; ++i)
            rotated[vertex][i] = rotation[i][0] * vertices[vertex][0]
                                 + rotation[i][1] * vertices[vertex][1] + translation[i];
    calMatCorner_eucclhyd(rotated[0], rotated[1], rotated[2], impedance, rotated_matrix);
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            double expected = 0.0;
            for (size_t k = 0; k < 2; ++k)
                for (size_t l = 0; l < 2; ++l)
                    expected += rotation[i][k] * matrix[k][l] * rotation[j][l];
            expectNear(rotated_matrix[i][j], expected);
        }
    }
}

static void testDirectionalDissipation(void)
{
    double p[2] = {0.0, 0.0}, next[2] = {2.0, 0.0}, prev[2] = {0.0, 1.0};
    double eucclhyd[2][2], glace[2][2];
    calMatCorner(CORNER_SCHEME_EUCCLHYD, p, next, prev, 3.0, eucclhyd);
    calMatCorner(CORNER_SCHEME_GLACE, p, next, prev, 3.0, glace);

    double x_direction[2] = {1.0, 0.0};
    double y_direction[2] = {0.0, 1.0};
    double corner_tangent[2] = {2.0, -1.0};
    expectNear(quadraticForm(eucclhyd, x_direction), 1.5);
    expectNear(quadraticForm(eucclhyd, y_direction), 3.0);
    /* GLACE cannot damp a jump perpendicular to its combined corner normal. */
    expectNear(quadraticForm(glace, corner_tangent), 0.0);
    expectNear(quadraticForm(eucclhyd, corner_tangent), 9.0);
}

static void testCornerForce(void)
{
    double p[2] = {0.0, 0.0}, next[2] = {2.0, 0.0}, prev[2] = {0.0, 1.0};
    double matrix[2][2], normal[2], force[2];
    const double length = calCornerLength(p, next, prev);
    calCornerNormal(p, next, prev, normal);
    calMatCorner_eucclhyd(p, next, prev, 3.0, matrix);

    double cell_velocity[2] = {0.4, -0.8};
    double node_velocity[2] = {1.4, 1.2};
    calCornerForce(force, matrix, cell_velocity, 7.0, length, normal, node_velocity);
    /* At this SW corner: p*l*n=(-3.5,-7), M*delta_u=(1.5,6). */
    expectNear(force[0], -5.0);
    expectNear(force[1], -13.0);

    cell_velocity[0] += 11.0;
    node_velocity[0] += 11.0;
    cell_velocity[1] -= 9.0;
    node_velocity[1] -= 9.0;
    calCornerForce(force, matrix, cell_velocity, 7.0, length, normal, node_velocity);
    expectNear(force[0], -5.0);
    expectNear(force[1], -13.0);

    calCornerForce(force, matrix, cell_velocity, 7.0, length, normal, cell_velocity);
    expectNear(force[0], -3.5);
    expectNear(force[1], -7.0);

    calCornerForce(force, matrix, cell_velocity, 0.0, length, normal, node_velocity);
    double relative_velocity[2] = {1.0, 2.0};
    expectNear(-vecMulVec2(relative_velocity, force), 13.5);
}

static void checkAcousticRiemann(double left, double right, double bottom, double top,
                                size_t axis)
{
    double p[2] = {0.0, 0.0};
    /* Incident cells SW, SE, NW, NE; next/prev follow each cell's CCW order. */
    double next[4][2] = {{-left, 0.0}, {0.0, -bottom}, {0.0, top}, {right, 0.0}};
    double prev[4][2] = {{0.0, -bottom}, {right, 0.0}, {-left, 0.0}, {0.0, top}};
    const size_t x_sides[4] = {0, 1, 0, 1};
    const size_t y_sides[4] = {0, 0, 1, 1};
    const double impedances[2] = {2.3, 0.7};
    const double pressures[2] = {5.0, 1.0};
    const double velocities[2] = {0.4, -0.6};
    double corner_matrices[4][2][2], normals[4][2], lengths[4], cell_velocities[4][2];
    double nodal_matrix[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double rhs[2] = {0.0, 0.0};

    for (size_t cell = 0; cell < 4; ++cell) {
        const size_t side = axis == 0 ? x_sides[cell] : y_sides[cell];
        cell_velocities[cell][axis] = velocities[side];
        cell_velocities[cell][1 - axis] = -0.25;
        lengths[cell] = calCornerLength(p, next[cell], prev[cell]);
        calCornerNormal(p, next[cell], prev[cell], normals[cell]);
        calMatCorner(CORNER_SCHEME_EUCCLHYD, p, next[cell], prev[cell],
                     impedances[side], corner_matrices[cell]);
        attMat2(nodal_matrix, corner_matrices[cell]);
        double velocity_contribution[2];
        matMulVec2(corner_matrices[cell], cell_velocities[cell], velocity_contribution);
        for (size_t component = 0; component < 2; ++component)
            rhs[component] += lengths[cell] * pressures[side] * normals[cell][component]
                              + velocity_contribution[component];
    }

    double inverse[2][2], node_velocity[2];
    inverseMatrix2(nodal_matrix, inverse);
    matMulVec2(inverse, rhs, node_velocity);
    /* Acoustic Riemann velocity: (z_L*u_L+z_R*u_R+p_L-p_R)/(z_L+z_R). */
    expectNear(node_velocity[axis], 1.5);
    expectNear(node_velocity[1 - axis], -0.25);

    double total_force[2] = {0.0, 0.0};
    for (size_t cell = 0; cell < 4; ++cell) {
        const size_t side = axis == 0 ? x_sides[cell] : y_sides[cell];
        double force[2];
        calCornerForce(force, corner_matrices[cell], cell_velocities[cell], pressures[side],
                       lengths[cell], normals[cell], node_velocity);
        attVec2(total_force, force);
    }
    expectNear(total_force[0], 0.0);
    expectNear(total_force[1], 0.0);
}

static void testAcousticRiemannReduction(void)
{
    const double geometries[][4] = {
        {1.0, 1.0, 1.0, 1.0},
        {0.07, 1.8, 2.5, 0.12},
        {3.0, 0.2, 0.1, 2.0}
    };
    for (size_t geometry = 0; geometry < sizeof geometries / sizeof geometries[0]; ++geometry)
        for (size_t axis = 0; axis < 2; ++axis)
            checkAcousticRiemann(geometries[geometry][0], geometries[geometry][1],
                                geometries[geometry][2], geometries[geometry][3], axis);
}

int main(void)
{
    testRectangularCornersAndDispatch();
    testSkewCornerAndRotation();
    testDirectionalDissipation();
    testCornerForce();
    testAcousticRiemannReduction();
    puts("physical_calculation_test: PASS (EUCCLHYD tensors, dissipation, Riemann reduction, and forces)");
    return 0;
}
