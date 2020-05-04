#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H

#include <vector>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>

using namespace glm;
using namespace std;

class Particle {
private:
    vec3 Position;
    vec3 Velocity;
    vec3 PermForce;
    vec3 TempForce;
    float Mass;

public:
    Particle(vec3 Position, float Mass) : Position(Position), Mass(Mass) {
        Velocity *= 0;
        TempForce *= 0;
        PermForce *= 0;
    }

    float getMass() { return Mass; }

    vec3 getPosition() { return Position; }

    vec3 getVelocity() { return Velocity; }

    void setMass(float Mass) { this->Mass = Mass; }

    void setPosition(vec3 Position) { this->Position = Position; }

    void setVelocity(vec3 Velocity) { this->Velocity = Velocity; }

    void clearTempForce() { TempForce *= 0; }

    void clearPermForce() { PermForce *= 0; }

    void applyPermForce(vec3 f) { PermForce += f;  }

    void applyTempForce(vec3 f) { TempForce += f; }

    vec3 computeTempAcceleration() { return (1 / Mass) * (TempForce); }

    vec3 computeAcceleration() { return (1 / Mass) * (TempForce + PermForce); }

    void integrate(float deltaTime) {
        vec3 accel = computeAcceleration();
        Velocity += accel * deltaTime;
        Position += Velocity * deltaTime;
    }
};

#endif