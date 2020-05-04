#ifndef CONCRETEBLOCKITG_H
#define CONCRETEBLOCKITG_H

#include "concreteblock.h"

// Derived integrators
class ConcreteBlockIntegrator {
public:
    ConcreteBlockIntegrator(float penalty) : Penalty(penalty){}

    void Integrate(ConcreteBlockSystem& system, float timestep) {
        // Get positions & velocities
        int numDofs = system.GetNumDofs();
        vector<vec3> pos(numDofs), vel(numDofs), acc(numDofs);
        vector<Particle*> ps(numDofs);

        system.GetPositions(pos);
        system.GetVelocities(vel);
        system.GetParticles(ps);
        for (int i = 0; i < numDofs; i++) {
            if (pos[i][1] <= 0) {
                vec3 penalty_force(0, -Penalty * pos[i][1], 0);
                ps[i]->ApplyTempForce(penalty_force);
            }
        }
        system.ComputeAccelerations(acc);
        // Forward Euler step
        for (int i = 0; i < numDofs; i++) {
            vel[i] += acc[i] * timestep;
            pos[i] += vel[i] * timestep;
            vel[i] = vel[i] * 0.999;
        } // Store results
        system.SetPositions(pos);
        system.SetVelocities(vel);
    }
private:
    float Penalty;
};

#endif