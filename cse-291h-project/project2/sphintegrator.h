#ifndef SPHINTEGRATOR_H
#define SPHINTEGRATOR_H

#include "sph.h"

// Derived integrators
class SPHIntegrator {
private:
    float penalty;
public:
    SPHIntegrator(float penalty) : penalty(penalty) {}

    void integrate(SPHSystem& system, float timestep) {
        vector<Particle*>* ps = system.getParticles();
        for (Particle* p : *ps)
            p->integrate(timestep, penalty);
    }
};

#endif