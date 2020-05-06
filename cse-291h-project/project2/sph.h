#ifndef SPH_H
#define SPH_H

#include <vector>
#include <cmath> 
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/scalar_multiplication.hpp>
#include "sphparticle.h"
#include "spatialhashtable.h"

#define PI 3.14159265358979323846264338327950288

using namespace std;
using namespace glm;

class SPHSystem {

private:
    SpatialHashTable* table;
    vector<Particle> Ps;
    vec3 spos, size, gap, m_d;
    float k, density0, supportRadius, smoothingRadius, penalty;

    float f(float q) {
        float ans;
        if (0 <= q && q < 1)
            ans = 2.0 / 3.0 - pow(q, 2) + 0.5 * pow(q, 3);
        else if (1 <= q && q < 2)
            ans = 1.0 / 6.0 * pow(2 - q, 3);
        else
            ans = 0.0;
        return ans * 3.0 / 2.0 / PI;
    }

    float df(float q) {
        float ans;
        if (0 <= q && q < 1)
            ans = - 2 * q + 3.0 / 2.0 * pow(q, 2);
        else if (1 <= q && q < 2)
            ans = - 0.5 * pow(2 - q, 2);
        else
            ans = 0.0;
        return ans * 3.0 / 2.0 / PI;
    }

    float W(Particle x_i, Particle x_j) {
        float q = distance(x_i.getPosition(), x_j.getPosition()) / smoothingRadius;
        return 1.0 / pow(smoothingRadius, 3) * f(q);
    }

    vec3 dW(Particle x_i, Particle x_j) {
        float q = distance(x_i.getPosition(), x_j.getPosition()) / smoothingRadius;
        vec3 norm = x_i.getPosition() - x_j.getPosition();
        norm = norm / sqrt(dot(norm, norm));
        return 1.0 / pow(smoothingRadius, 4) * df(q) * norm;
    }

    vector<Particle> getNeighbors(Particle p) {
        vector<Particle> ans;
        for (int i = 0; i < Ps.size(); i++) {
            if (distance(p.getPosition(), Ps[i].getPosition()) <= supportRadius)
                ans.push_back(Ps[i]);
        }
    }

    float computeDensity(Particle P) {
        float ans = 0;
        vector<Particle> neighbors = getNeighbors(P);
        for (int i = 0; i < neighbors.size(); i++)
            ans += neighbors[i].getMass() * W(P, neighbors[i]);
        return ans;
    }

    float computePressure(Particle P) {
        float density = computeDensity(P);
        return k * (pow(density / density0, 7) - 1);
    }

    vec3 computeFpres(Particle P) {
        vec3 dPres;
        float mass = P.getMass();
        float dens = computeDensity(P);
        float pres = computePressure(P);
        vector<Particle> neighbors = getNeighbors(P);
        for (int j = 0; j < neighbors.size(); j++) {
            float dens_j = computeDensity(neighbors[j]);
            float pres_j = computePressure(neighbors[j]);
            float mass_j = neighbors[j].getMass();
            dPres += mass_j * (pres / pow(dens, 2) + pres_j / pow(dens_j, 2)) * dW(P, neighbors[j]);
        } dPres *= dens;
        return - mass / dens * dPres;
    }

    vec3 computeFvisc(Particle P)
    {
        vec3 ddVisc;
        float mass = P.getMass();
        vec3 vel = P.getVelocity();
        vector<Particle> neighbors = getNeighbors(P);
        for (int j = 0; j < neighbors.size(); j++) {
            float dens_j = computeDensity(neighbors[j]);
            float mass_j = neighbors[j].getMass();
            vec3 x_ij = P.getPosition() - neighbors[j].getPosition();
            vec3 v_ij = P.getVelocity() - neighbors[j].getVelocity();
            vec3 dW_ij = dW(P, neighbors[j]);
            ddVisc += mass_j / dens_j * v_ij * (x_ij * dW_ij) / (dot(x_ij, x_ij) + 0.01 * pow(smoothingRadius, 2));
        } ddVisc *= 2;
        return mass * vel * ddVisc;
    }

public:
    SPHSystem(vec3 pos, vec3 size, vec3 gap, vec3 m_d, float penalty, 
        float k, float density0, float supportRadius, float smoothingRadius) :
        spos(pos), size(size), gap(gap), m_d(m_d), penalty(penalty), 
        k(k), density0(density0), supportRadius(supportRadius), smoothingRadius(smoothingRadius)
    {
        this->table = new SpatialHashTable(m_d);
        float mass = pow(smoothingRadius, 3) * density0;
        // float mass = pow(2.0 / 3.0 * smoothingRadius, 3) * density0;
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    // mat4 rotate_matrix(1.0f);
                    // rotate_matrix = rotate(rotate_matrix, radians(0.0f), vec3(0, 0, 1));
                    // rotate_matrix = rotate(rotate_matrix, radians(0.0f), vec3(1, 0, 0));
                    vec3 ppos(gap[0] * ix, gap[1] * iy, gap[2] * iz);
                    // ppos = rotate_matrix * vec4(ppos, 1);
                    ppos = ppos + pos;
                    Particle P(ppos, 0); 
                    Ps.push_back(P);
                }
            }
        }
    }

    void applyGForces() {
        vec3 g(0.0f, -9.8f, 0.0f);
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i].applyPermForce(Ps[i].getMass() * g);
    }

    void applySPHForces() {
        for (int i = 0 ; i < Ps.size() ; i++) {
            Ps[i].applyTempForce(computeFpres(Ps[i]));
            Ps[i].applyTempForce(computeFvisc(Ps[i]));
        }
    }

    void applyPenaltyForces() {
        for (int i = 0; i < Ps.size(); i++) {
            vec3 pos = Ps[i].getPosition();
            if (pos.y <= 0) {
                vec3 penalty_force(0, -penalty * pos.y, 0);
                Ps[i].applyTempForce(penalty_force);
            }
            if (pos.x <= 0) {
                vec3 penalty_force(-penalty * pos.x, 0, 0);
                Ps[i].applyTempForce(penalty_force);
            }
            if (pos.x > spos.x + gap.x * (size.x - 1)) {
                vec3 penalty_force(-penalty * (pos.x - (spos.x + gap.x * (size.x - 1))), 0, 0);
                Ps[i].applyTempForce(penalty_force);
            }
            if (pos.z <= 0) {
                vec3 penalty_force(0, 0, -penalty * pos.z);
                Ps[i].applyTempForce(penalty_force);
            }
            if (pos.z > spos.z + gap.z * (size.z - 1)) {
                vec3 penalty_force(0, 0, -penalty * (pos.z - (spos.z + gap.z * (size.z - 1))));
                Ps[i].applyTempForce(penalty_force);
            }
        }
    }

    void buildTable() {
        table->clear();
        table->build(Ps);
    }

    int getSize() { return Ps.size(); }

    int getNumDofs() { return getSize(); }

    void getParticles(vector<Particle*>& ps) {
        ps.clear();
        for (int i = 0 ; i < ps.size() ; i++) 
            ps.push_back(&Ps[i]);
    }

    void getPositions(float* pos) {
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i].getPosition();
            for (int j = 0 ; j < 3 ; j++) pos[3 * i + j] = ppos[j];
        }
    }

    void getPositions(vector<vec3>& pos) {
        pos.clear();
        for (int i = 0 ; i < Ps.size() ; i++) 
            pos.push_back(Ps[i].getPosition());
    }

    void getVelocities(float* vel) {
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 pvel = Ps[i].getVelocity();
            for (int j = 0 ; j < 3 ; j++) vel[3 * i + j] = pvel[j];
        }
    }

    void getVelocities(vector<vec3>& vel) {
        vel.clear();
        for (int i = 0 ; i < Ps.size() ; i++) 
            vel.push_back(Ps[i].getVelocity());
    }

    void setPositions(const vector<vec3>& pos) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i].setPosition(pos[i]);
    }

    void setVelocities(vec3& vel) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i].setVelocity(vel);
    }

    void setVelocities(const vector<vec3>& vel) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i].setVelocity(vel[i]);
    }

    void applyTempForces(vec3& f) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i].applyTempForce(f);
    }

    void clearTempForces() {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i].clearTempForce();
    }

    void getAccelations(vector<vec3>& accs) 
    {
        accs.clear();
        for (int i = 0 ; i < Ps.size() ; i++)
            accs.push_back(Ps[i].computeAcceleration());
    }

    void getPosAcc(float* posacc) {
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i].getPosition();
            vec3 pacc = ppos + Ps[i].computeTempAcceleration();
            for (int j = 0 ; j < 3 ; j++) posacc[6 * i + 0 + j] = ppos[j];
            for (int j = 0 ; j < 3 ; j++) posacc[6 * i + 3 + j] = pacc[j];
        }
    }

    void getPosAcc(vector<vec3>& posacc) {
        posacc.clear();
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i].getPosition();
            vec3 pacc = ppos + Ps[i].computeTempAcceleration();
            posacc.push_back(ppos);
            posacc.push_back(ppos + pacc);
        }
    }

    void printPositions() 
    {
        printf("=========================================\n");
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 pos = Ps[i].getPosition();
            printf("%.3f-%.3f-%.3f ", pos[0], pos[1], pos[2]);
        }
    }
};

#endif