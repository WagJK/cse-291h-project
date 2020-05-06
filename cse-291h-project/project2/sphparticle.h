#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H

#include <vector>
#include <map>
#include <cmath> 
#include <algorithm>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>

using namespace glm;
using namespace std;

class Particle {
private:
    vec3 pos;
    vec3 vel;
    vec3 permF;
    vec3 tempF;
    float mass;
    float dens;
    float pres;
    vec3 Fpres;
    vec3 Fvisc;
    bool built_basics;
    bool built_forces;

public:
    Particle(vec3 pos, float mass) : pos(pos), mass(mass) {
        vel *= 0;
        tempF *= 0;
        permF *= 0;
        built_basics = false;
        dens = 0;
        pres = 0;
        built_forces = false;
        Fpres = { 0, 0, 0 };
        Fvisc = { 0, 0, 0 };
    }

    float getMass() { return mass; }

    float getDensity() { 
        if (!built_basics) throw "basics not built!";
        return dens; 
    }

    float getPressure() {
        if (!built_basics) throw "basics not built!"; 
        return pres; 
    }

    vec3 getFpres() {
        if (!built_forces) throw "forces not built!";
        return Fpres;
    }

    vec3 getFvisc() {
        if (!built_forces) throw "forces not built!";
        return Fvisc;
    }

    vec3 getPosition() { return pos; }

    vec3 getVelocity() { return vel; }

    void setMass(float mass) { this->mass = mass; }

    void setDensity(float dens) { this->dens = dens;  }

    void setPressure(float pres) { this->pres = pres; }

    void setFpres(vec3 Fpres) { this->Fpres = Fpres; }

    void setFvisc(vec3 Fvisc) { this->Fvisc = Fvisc; }

    void setPosition(vec3 pos) { this->pos = pos; }

    void setVelocity(vec3 vel) { this->vel = vel; }

    void set_built_basics(bool flag) { this->built_basics = flag; }

    void set_built_forces(bool flag) { this->built_forces = flag; }

    void clearTempForce() { tempF *= 0; }

    void clearPermForce() { permF *= 0; }

    void applyPermForce(vec3 f) { permF += f;  }

    void applyTempForce(vec3 f) { tempF += f; }

    vec3 computeTempAcceleration() { return (1 / mass) * (tempF); }

    vec3 computeAcceleration() { return (1 / mass) * (tempF + permF); }

    void integrate(float deltaTime) {
        vec3 accel = computeAcceleration();
        vel += accel * deltaTime;
        pos += vel * deltaTime;
    }
};

#endif