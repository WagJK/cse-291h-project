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
    vector<Particle*>* neighbors;
    bool flagBuiltBasics;
    bool flagBuiltForces;
    bool flagBuiltNeighbors;

public:
    Particle(vec3 pos, float mass) : pos(pos), mass(mass) {
        vel *= 0;
        tempF *= 0;
        permF *= 0;
        flagBuiltBasics = false;
        dens = 0;
        pres = 0;
        flagBuiltForces = false;
        Fpres = { 0, 0, 0 };
        Fvisc = { 0, 0, 0 };
        flagBuiltNeighbors = false;
        neighbors = new vector<Particle*>();
    }

    float getMass() { return mass; }
    void setMass(float mass) { this->mass = mass; }

    float getDensity() { 
        if (!flagBuiltBasics) throw "basics not built!";
        return dens; 
    }
    void setDensity(float dens) { this->dens = dens; }

    float getPressure() {
        if (!flagBuiltBasics) throw "basics not built!"; 
        return pres; 
    }
    void setPressure(float pres) { this->pres = pres; }

    vec3 getFpres() {
        if (!flagBuiltForces) throw "forces not built!";
        return Fpres;
    }
    void setFpres(vec3 Fpres) { this->Fpres = Fpres; }

    vec3 getFvisc() {
        if (!flagBuiltForces) throw "forces not built!";
        return Fvisc;
    }
    void setFvisc(vec3 Fvisc) { this->Fvisc = Fvisc; }

    vec3 getPosition() { return pos; }
    void setPosition(vec3 pos) { this->pos = pos; }

    vec3 getVelocity() { return vel; }
    void setVelocity(vec3 vel) { this->vel = vel; }

    void setNeighbors(vector<Particle*>* ps) { neighbors = ps; }
    vector<Particle*>* getNeighbors(bool check) {
        if (check && !flagBuiltNeighbors)
            throw "neighbors not built!";
        return neighbors;
    }

    void setBuiltBasicsFlag(bool flag) { this->flagBuiltBasics = flag; }
    void setBuiltForcesFlag(bool flag) { this->flagBuiltForces = flag; }
    void setBuiltNeighborsFlag(bool flag) { this->flagBuiltNeighbors = flag; }

    void clearTempForce() { tempF *= 0; }
    void clearPermForce() { permF *= 0; }
    void applyPermForce(vec3 f) { permF += f;  }
    void applyTempForce(vec3 f) { tempF += f; }

    vec3 computeTempAcceleration() { return (1 / mass) * (tempF); }
    vec3 computeAcceleration() { return (1 / mass) * (tempF + permF); }

    void integrate(float deltaTime, float penalty) {
        vec3 acc = computeAcceleration();
        vel += acc * deltaTime;
        pos += vel * deltaTime;
        vel = vel * penalty;
    }
};

#endif