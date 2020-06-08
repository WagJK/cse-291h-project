#ifndef SPHINTEGRATOR_H
#define SPHINTEGRATOR_H

#include "sph.h"
#include "omp.h"

// Derived integrators
class SPHIntegrator {
private:
    float penalty;
public:
    SPHIntegrator(float penalty) : penalty(penalty) {}

    void integrate(SPHSystem& system, float timestep) {
        vector<Particle*>* ps = system.getParticles();
#pragma omp parallel for
        for (int i = 0 ; i < ps->size() ; i++)
            ps->at(i)->integrate(timestep, penalty);
    }
};

#endif