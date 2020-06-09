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

public:
    Particle(vec3 pos, float mass) : pos(pos), mass(mass) {
        vel *= 0;
        tempF *= 0;
        permF *= 0;
    }

    float getMass() { return mass; }
    void setMass(float mass) { this->mass = mass; }

    vec3 getPosition() { return pos; }
    void setPosition(vec3 pos) { this->pos = pos; }

    vec3 getVelocity() { return vel; }
    void setVelocity(vec3 vel) { this->vel = vel; }

    void clearTempForce() { tempF *= 0; }
    void clearPermForce() { permF *= 0; }
    void applyPermForce(vec3 f) { permF += f; }
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

class FluidParticle : public Particle {
private:
    float dens;
    vec3 norm;
    float pres;
    vec3 Fpres;
    vec3 Fvisc;

    vector<FluidParticle*>* neighbors;
    vector<vector<FluidParticle*>*>* neighborBlocks;

    float v_diff, kv, Ek;
    bool flagBuiltBasics;
    bool flagBuiltForces;
    bool flagBuiltNeighbors;

public:
    FluidParticle(vec3 pos, float mass) : Particle(pos, mass) {
        flagBuiltBasics = false;
        dens = 0;
        norm = { 0, 0, 0 };
        pres = 0;
        
        flagBuiltForces = false;
        Fpres = { 0, 0, 0 };
        Fvisc = { 0, 0, 0 };
        
        flagBuiltNeighbors = false;
        neighbors = new vector<FluidParticle*>();
        neighborBlocks = new vector<vector<FluidParticle*>*>();
        
        v_diff = 0; kv = 0; Ek = 0;
    }

    float getDensity() {
        if (!flagBuiltBasics) throw "basics not built!";
        return dens;
    }
    void setDensity(float dens) { this->dens = dens; }

    vec3 getNormal() {
        if (!flagBuiltBasics) throw "basics not built";
        return norm;
    }
    void setNormal(vec3 norm) { this->norm = norm; }

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

    vec3 getDiffusePotential() { return vec3(v_diff, kv, Ek); }
    void setDiffusePotential(float vi_diff, float ki, float Eki) {
        this->v_diff = vi_diff;
        this->kv = ki;
        this->Ek = Eki;
    }

    vector<FluidParticle*>* getNeighbors(bool check) {
        if (check && !flagBuiltNeighbors)
            throw "neighbors not built!";
        return neighbors;
    }
    vector<vector<FluidParticle*>*>* getNeighborBlocks() { return neighborBlocks; }

    void setBuiltBasicsFlag(bool flag) { this->flagBuiltBasics = flag; }
    void setBuiltForcesFlag(bool flag) { this->flagBuiltForces = flag; }
    void setBuiltNeighborsFlag(bool flag) { this->flagBuiltNeighbors = flag; }

};

class DiffuseParticle : public FluidParticle {

private:
    int lifetime;
    int diffuseType; // 0 - Fluid ; 1 - Spray ; 2 - Foam ; 3 - Bubble

public:
    const static int TYPE_SPRY = 1;
    const static int TYPE_FOAM = 2;
    const static int TYPE_BUBB = 3;

    DiffuseParticle(vec3 pos, float mass, int type, int lifetime) : FluidParticle(pos, mass) {
        this->diffuseType = type;
        this->lifetime = lifetime;
    }

    void setDiffuseType(int type) { diffuseType = type; }
    int getDiffuseType() { return diffuseType; }

    bool reduceLifetime() {
        lifetime--;
        return (lifetime <= 0);
    }
};
#endif