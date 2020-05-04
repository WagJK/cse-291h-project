#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include <vector>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>
#include "particle.h"

using namespace glm;
using namespace std;

class Tetrahedron {

public:
    // R is rest state
    // P is current position
    // 0, 1, 2 needs to be anti-clockwise
    Tetrahedron(Particle* P0, Particle* P1, Particle* P2, Particle* P3, float mass, float E, float v) :
        P0(P0), P1(P1), P2(P2), P3(P3), E(E), v(v)
    {
        // add mass to every particle
        P0->SetMass(P0->GetMass() + mass / 4.0f);
        P1->SetMass(P1->GetMass() + mass / 4.0f);
        P2->SetMass(P2->GetMass() + mass / 4.0f);
        P3->SetMass(P3->GetMass() + mass / 4.0f);

        // store rest state
        r0 = P0->GetPosition();
        r1 = P1->GetPosition();
        r2 = P2->GetPosition();
        r3 = P3->GetPosition();

        // compute lame consts
        mu = E * v / ((1 + v) * (1 - 2 * v));
        lambda = E / (2 * (1 + v));

        // compute face vec
        // 1, 3, 2
        n0 = (1.0 / 2.0) * cross(P3->GetPosition() - P1->GetPosition(), P2->GetPosition() - P1->GetPosition());
        // 0, 2, 3
        n1 = (1.0 / 2.0) * cross(P2->GetPosition() - P0->GetPosition(), P3->GetPosition() - P0->GetPosition());
        // 0, 3, 1
        n2 = (1.0 / 2.0) * cross(P3->GetPosition() - P0->GetPosition(), P1->GetPosition() - P0->GetPosition());
        // 0, 1, 2
        n3 = (1.0 / 2.0) * cross(P1->GetPosition() - P0->GetPosition(), P2->GetPosition() - P0->GetPosition());

        // compute R & R_inv
        R[0] = r0 - r3;
        R[1] = r1 - r3;
        R[2] = r2 - r3;
        R_inv = inverse(R);
    }

    void print_vec(vec3 x) {
        for (int j = 0; j < 3; j++) printf("%f ", x[j]);
        printf("\n");
    }

    void print_mat(mat3 x) {
        for (int i = 0; i < 3; i++) print_vec(x[i]);
        printf("\n");
    }

    mat3 ComputeT() {
        mat3 T;
        T[0] = P0->GetPosition() - P3->GetPosition();
        T[1] = P1->GetPosition() - P3->GetPosition();
        T[2] = P2->GetPosition() - P3->GetPosition();
        return T;
    }

    mat3 ComputeF() {
        mat3 T = ComputeT();
        return T * R_inv;
    }

    mat3 ComputeStrainGreen() {
        mat3 F = ComputeF(), I(1.0f);
        mat3 eps = (1.0 / 2.0) * (transpose(F) * F - I);
        return eps;
    }

    mat3 ComputeStressIso() {
        mat3 eps = ComputeStrainGreen(), I(1.0f);
        mat3 omg = 2.0 * mu * eps + lambda * (eps[0][0] + eps[1][1] + eps[2][2]) * I;
        return omg;
    }

    void ApplyGForces() {
        vec3 g(0.0f, -9.8f, 0.0f);
        P0->ApplyPermForce(P0->GetMass() * g);
        P1->ApplyPermForce(P0->GetMass() * g);
        P2->ApplyPermForce(P0->GetMass() * g);
        P3->ApplyPermForce(P0->GetMass() * g);
    }

    void ApplyTempForces(vec3& f) {
        P0->ApplyTempForce(f);
        P1->ApplyTempForce(f);
        P2->ApplyTempForce(f);
        P3->ApplyTempForce(f);
    }

    void ClearTempForces() {
        P0->ClearTempForce();
        P1->ClearTempForce();
        P2->ClearTempForce();
        P3->ClearTempForce();
    }

    void ComputeForces() {
        mat3 F = ComputeF();
        mat3 omg = ComputeStressIso();
        // print_mat(omg);

        vec3 f0 = -F * omg * n0;
        vec3 f1 = -F * omg * n1;
        vec3 f2 = -F * omg * n2;
        vec3 f3 = -F * omg * n3;
        // print_vec(f0);

        P0->ApplyTempForce(f0);
        P1->ApplyTempForce(f1);
        P2->ApplyTempForce(f2);
        P3->ApplyTempForce(f3);
    }


private:
    Particle* P0, * P1, * P2, * P3;
    float E, v, mu, lambda;
    vec3 r0, r1, r2, r3;
    vec3 n0, n1, n2, n3;
    mat3 R, R_inv;
};

#endif