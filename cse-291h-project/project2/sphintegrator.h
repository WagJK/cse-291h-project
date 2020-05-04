#ifndef SPHINTEGRATOR_H
#define SPHINTEGRATOR_H

#include "sph.h"

// Derived integrators
class SPHIntegrator {
private:
    float Penalty;
public:
    SPHIntegrator(float penalty) : Penalty(penalty) {}

    void Integrate(SPHSystem& system, float timestep) {
        // Get positions & velocities
        int numDofs = system.getSize();
        vector<vec3> pos(numDofs), vel(numDofs), acc(numDofs);
        vector<Particle*> ps(numDofs);

        system.getPositions(pos);
        system.getVelocities(vel);
        system.getParticles(ps);
        for (int i = 0; i < numDofs; i++) {
            if (pos[i][1] <= 0) {
                vec3 penalty_force(0, -Penalty * pos[i][1], 0);
                ps[i]->applyTempForce(penalty_force);
            }
        }
        system.getAccelations(acc);
        // Forward Euler step
        for (int i = 0; i < numDofs; i++) {
            vel[i] += acc[i] * timestep;
            pos[i] += vel[i] * timestep;
            vel[i] = vel[i] * 0.999;
        } // Store results
        system.setPositions(pos);
        system.setVelocities(vel);
    }
};

#endif