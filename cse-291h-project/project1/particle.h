#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>

using namespace glm;
using namespace std;

class Particle {
public:
    Particle(vec3 Position, float Mass) : Position(Position), Mass(Mass) {
        Velocity *= 0;
        TempForce *= 0;
        PermForce *= 0;
    }

    void ClearTempForce() { TempForce *= 0; }

    void ClearPermForce() { PermForce *= 0; }

    void ApplyPermForce(vec3 f) { PermForce += f;  }

    void ApplyTempForce(vec3 f) { TempForce += f; }

    void Integrate(float deltaTime) {
        vec3 accel = ComputeAcceleration();
        Velocity += accel * deltaTime;
        Position += Velocity * deltaTime;
    }

    float GetMass() { return Mass; }

    vec3 GetPosition() { return Position; }

    vec3 GetVelocity() { return Velocity; }

    void SetMass(float Mass) { this->Mass = Mass; }

    void SetPosition(vec3 Position) { this->Position = Position; }

    void SetVelocity(vec3 Velocity) { this->Velocity = Velocity; }

    vec3 ComputeTempAcceleration() { return (1 / Mass) * (TempForce); }

    vec3 ComputeAcceleration() { return (1 / Mass) * (TempForce + PermForce); }

private:
    vec3 Position;
    vec3 Velocity;
    vec3 PermForce;
    vec3 TempForce;
    float Mass;
};

#endif