#include <vector>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>

using namespace std;
using namespace glm;

// Base system
class PhysicsSystem {
public:
    virtual int GetNumDofs() { return -1; };
    virtual void GetPositions(float* pos) {};
    virtual void GetPositions(vector<vec3>& pos) {};
    virtual void GetVelocities(vector<vec3>& vel) {};
    virtual void SetPositions(const vector<vec3>& pos) {};
    virtual void SetVelocities(const vector<vec3>& vel) {};
    virtual void ApplyForces(vec3& f) {};
    virtual void ComputeForces() {};
    virtual void ComputeAccelerations(vector<vec3>& acc) {};
};

// Base integrator
class Integrator {
public:
    virtual void Integrate(PhysicsSystem& system, float timestep) {};
};


class Particle {
public:
    Particle(vec3 Position, float Mass) : 
        Position(Position), Mass(Mass) {
        Velocity *= 0;
        Force *= 0;
    }

    void ClearForce() { Force *= 0; }

    void ApplyForce(vec3& f) { Force += f; }
    
    void Integrate(float deltaTime) {
        vec3 accel = (1 / Mass) * Force;
        Velocity += accel * deltaTime;
        Position += Velocity * deltaTime;
        Force = vec3(0);
    }
    
    vec3 GetPosition() { return Position; }
    
    vec3 GetVelocity() { return Velocity; }
    
    void SetPosition(vec3 Position) { this->Position = Position; }
    
    void SetVelocity(vec3 Velocity) { this->Velocity = Velocity; }
    
    vec3 ComputeAcceleration() { return (1 / Mass) * Force; }

private:
    vec3 Position;
    vec3 Velocity;
    vec3 Force;
    float Mass;
};

class Tetrahedron {

public:
    // R is rest state
    // P is current position
    // 0, 1, 2 needs to be anti-clockwise
    Tetrahedron(
        Particle* R0, Particle* R1, Particle* R2, Particle* R3, 
        Particle* P0, Particle* P1, Particle* P2, Particle* P3, float E, float v) :
        R0(R0), R1(R1), R2(R2), R3(R3), P0(P0), P1(P1), P2(P2), P3(P3), E(E), v(v) 
    {
        // compute lame consts
        mu = E * v / ((1 + v) * (1 - 2 * v));
        lambda = E / (2 * (1 + v));
        // compute face vec
        n0 = (1.0 / 2.0) * cross(P2->GetPosition() - P1->GetPosition(), P3->GetPosition() - P1->GetPosition());
        n1 = (1.0 / 2.0) * cross(P2->GetPosition() - P0->GetPosition(), P3->GetPosition() - P0->GetPosition());
        n2 = (1.0 / 2.0) * cross(P1->GetPosition() - P0->GetPosition(), P3->GetPosition() - P0->GetPosition());
        n3 = (1.0 / 2.0) * cross(P1->GetPosition() - P0->GetPosition(), P2->GetPosition() - P0->GetPosition());
        // compute R & R_inv
        R[0] = R0->GetPosition() - R3->GetPosition();
        R[1] = R1->GetPosition() - R3->GetPosition();
        R[2] = R2->GetPosition() - R3->GetPosition();
        R_inv = inverse(R);
    }

    Tetrahedron(Particle* R0, Particle* R1, Particle* R2, Particle* R3, float E, float v) {
        Tetrahedron(R0, R1, R2, R3, R0, R1, R2, R3, E, v);
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
        mat3 F = ComputeF(), I{};
        mat3 eps = (1.0 / 2.0) * (transpose(F) * F - I);
        return eps;
    }

    mat3 ComputeStressIso() {
        mat3 eps = ComputeStrainGreen(), I{};
        mat3 omg = 2.0 * mu * eps + lambda * (eps[0][0] + eps[1][1] + eps[2][2]) * I;
        return omg;
    }

    void ApplyForces(vec3 &f) {
        P0->ApplyForce(f);
        P1->ApplyForce(f);
        P2->ApplyForce(f);
        P3->ApplyForce(f);
    }

    void ComputeForces() {
        mat3 F = ComputeF();
        mat3 omg = ComputeStressIso();
        vec3 f0 = F * omg * n0;
        vec3 f1 = F * omg * n1;
        vec3 f2 = F * omg * n2;
        vec3 f3 = F * omg * n3;
        P0->ApplyForce(f0);
        P1->ApplyForce(f1);
        P2->ApplyForce(f2);
        P3->ApplyForce(f3);
    }


private:
    Particle* R0, * R1, * R2, * R3;
    Particle* P0, * P1, * P2, * P3;
    float E, v, mu, lambda;
    vec3 n0, n1, n2, n3;
    mat3 R, R_inv;
};



class ConcreteBlockSystem : PhysicsSystem {
public:
    ConcreteBlockSystem(vec3 pos, vec3 size, vec3 gap, float mass, float E, float v) :
        pos(pos), size(size), gap(gap), mass(mass), E(E), v(v)
    {

        for (int ix = 0; ix < size[0]; ix++) {
            Ps.push_back(vector<vector<Particle>>());
            for (int iy = 0; iy < size[1]; iy++) {
                Ps[ix].push_back(vector<Particle>());
                for (int iz = 0; iz < size[2]; iz++) {
                    vec3 pos(pos[0] + gap[0] * ix, pos[1] + gap[1] * iy, pos[2] + gap[2] * iz);
                    Particle P(pos, mass);
                    Ps[ix][iy].push_back(P);
                }
            }
        }
        for (int ix = 0; ix < size[0] - 1; ix++) {
            for (int iy = 0; iy < size[1] - 1; iy++) {
                for (int iz = 0; iz < size[2] - 1; iz++) {
                    // 0,0,0 - 1,0,0 - 0,1,0 - 0,0,1
                    Tetrahedron top1(&Ps[ix][iy][iz], &Ps[ix + 1][iy][iz], &Ps[ix][iy + 1][iz], &Ps[ix][iy][iz + 1], E, v);
                    // 0,1,1 - 0,1,0 - 1,1,1 - 0,0,1
                    Tetrahedron top2(&Ps[ix][iy + 1][iz + 1], &Ps[ix][iy + 1][iz], &Ps[ix + 1][iy + 1][iz + 1], &Ps[ix][iy][iz + 1], E, v);
                    // 1,1,0 - 1,1,1 - 0,1,0 - 1,0,0
                    Tetrahedron bot1(&Ps[ix + 1][iy + 1][iz], &Ps[ix + 1][iy + 1][iz + 1], &Ps[ix][iy + 1][iz], &Ps[ix + 1][iy][iz], E, v);
                    // 1,0,1 - 0,0,1 - 1,1,1 - 1,0,0
                    Tetrahedron bot2(&Ps[ix + 1][iy][iz + 1], &Ps[ix][iy][iz + 1], &Ps[ix + 1][iy + 1][iz + 1], &Ps[ix + 1][iy][iz], E, v);
                    // center: 1,0,0 - 0,1,0 - 0,0,1 - 1,1,1
                    Tetrahedron cntr(&Ps[ix + 1][iy][iz], &Ps[ix][iy + 1][iz], &Ps[ix][iy][iz + 1], &Ps[ix + 1][iy + 1][iz + 1], E, v);
                    Ts.push_back(top1);
                    Ts.push_back(top2);
                    Ts.push_back(bot1);
                    Ts.push_back(bot2);
                    Ts.push_back(cntr);
                }
            }
        }
    }

    vec3 GetSize() { return size; }

    int GetNumDofs() {
        return size[0] * size[1] * size[2];
    }

    int GetNumLines() {
        int cnt_lines = 0;
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    if (ix != size[0] - 1) cnt_lines++;
                    if (iy != size[1] - 1) cnt_lines++;
                    if (iz != size[2] - 1) cnt_lines++;

                    if (ix != size[0] - 1 && iy != size[1] - 1) cnt_lines += 2;
                    if (iy != size[1] - 1 && iz != size[2] - 1) cnt_lines += 2;
                    if (ix != size[0] - 1 && iz != size[2] - 1) cnt_lines += 2;
                }
            }
        }
        return cnt_lines;
    }

    void GetLinesIndices(unsigned int* indices) {
        int cnt = 0;
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    int idx = ix * size[1] * size[2] + iy * size[2] + iz;
                    if (ix != size[0] - 1) {
                        indices[cnt++] = idx;
                        indices[cnt++] = idx + size[1] * size[2];
                    }
                    if (iy != size[1] - 1) {
                        indices[cnt++] = idx;
                        indices[cnt++] = idx + size[2];
                    }
                    if (iz != size[2] - 1) {
                        indices[cnt++] = idx;
                        indices[cnt++] = idx + 1;
                    }
                    if (ix != size[0] - 1 && iy != size[1] - 1) {
                        indices[cnt++] = idx;
                        indices[cnt++] = idx + size[1] * size[2] + size[2];
                        indices[cnt++] = idx + size[1] * size[2];
                        indices[cnt++] = idx + size[2];
                    }
                    if (iy != size[1] - 1 && iz != size[2] - 1) {
                        indices[cnt++] = idx;
                        indices[cnt++] = idx + size[2] + 1;
                        indices[cnt++] = idx + size[2];
                        indices[cnt++] = idx + 1;
                    }
                    if (ix != size[0] - 1 && iz != size[2] - 1) {
                        indices[cnt++] = idx;
                        indices[cnt++] = idx + size[1] * size[2] + 1;
                        indices[cnt++] = idx + size[1] * size[2];
                        indices[cnt++] = idx + 1;
                    }
                }
            }
        }
    }

    void GetPositions(float* pos) {
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    vec3 ppos = Ps[ix][iy][iz].GetPosition();
                    int idx = 3 * (ix * size[1] * size[2] + iy * size[2] + iz);
                    pos[idx + 0] = ppos[0];
                    pos[idx + 1] = ppos[1];
                    pos[idx + 2] = ppos[2];
                    // printf("%f %f %f\n", ppos[0], ppos[1], ppos[2]);
                }
            }
        }
    }

    void GetPositions(vector<vec3>& pos) {
        pos.clear();
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    pos.push_back(Ps[ix][iy][iz].GetPosition());
    }
    
    void GetVelocities(vector<vec3>& vel) {
        vel.clear();
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    vel.push_back(Ps[ix][iy][iz].GetVelocity());
    }

    void SetPositions(const vector<vec3>& pos) {
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].SetPosition(pos[ix * size[1] * size[2] + iy * size[2] + iz]);
    }
    
    void SetVelocities(vec3& vel) {
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].SetVelocity(vel);
    }

    void SetVelocities(const vector<vec3>& vel) {
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].SetVelocity(vel[ix * size[1] * size[2] + iy * size[2] + iz]);
    }

    void ApplyForces(vec3& f) {
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].ApplyForce(f);
    }

    void ComputeForces() {
        for (auto tthd : Ts) tthd.ComputeForces();
    }

    void ComputeAccelerations(vector<vec3>& acc) {
        acc.clear();
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    acc.push_back(Ps[ix][iy][iz].ComputeAcceleration());
                }
            }
        }
    }

    void PrintPositions() {
        printf("=========================================\n");
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    vec3 pos = Ps[ix][iy][iz].GetPosition();
                    printf("%.3f-%.3f-%.3f  ", pos[0], pos[1], pos[2]);
                } printf("\n");
            } printf("\n\n");
        }
    }

private:
    vector<vector<vector<Particle>>> Ps;
    vector<Tetrahedron> Ts;
    vec3 pos, size, gap;
    float E, v, mass;
};

// Derived integrators
class ConcreteBlockIntegrator {
public:
    void Integrate(ConcreteBlockSystem& system, float timestep) {
        // Get positions & velocities
        int numDofs = system.GetNumDofs();
        vector<vec3> pos(numDofs), vel(numDofs);
        system.GetPositions(pos);
        system.GetVelocities(vel);
        // Compute accelerations
        vector<vec3> acc(numDofs);
        system.ComputeAccelerations(acc);
        // Forward Euler step
        for (int i = 0; i < numDofs; i++) {
            vel[i] += acc[i] * timestep;
            pos[i] += vel[i] * timestep;
            if (pos[i][2] <= 0) {
                vel[i][2] = -vel[i][2];
                pos[i][2] = -pos[i][2];
            }
        }
        // Store results
        system.SetPositions(pos);
        system.SetVelocities(vel);
    }
};
