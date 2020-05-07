#ifndef SPH_H
#define SPH_H

#include <vector>
#include <cstdlib>
#include <cmath> 
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/scalar_multiplication.hpp>
#include "sphparticle.h"
#include "sphmap.h"

#define LOG 0
#define SHOW_FPRES 0
#define SHOW_FVISC 0
#define PI 3.14159265358979323846264338327950288

using namespace std;
using namespace glm;

class SPHSystem {

private:
    SpatialHashTable* table;
    vector<Particle*> Ps;
    vec3 spos;
    vec3 size;
    vec3 gap;
    vec3 m_d;
    vec3 container_lb;
    vec3 container_ub;
    float k;
    float density0;
    float supportRadius;
    float smoothingRadius;
    float penalty;
    float viscosity;

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
            ans = -2 * q + 3.0 / 2.0 * pow(q, 2);
        else if (1 <= q && q < 2)
            ans = -0.5 * pow(2 - q, 2);
        else
            ans = 0.0;
        return ans * 3.0 / 2.0 / PI;
    }

    float W(Particle* x_i, Particle* x_j, float smoothingRadius) {
        float q = distance(x_i->getPosition(), x_j->getPosition()) / smoothingRadius;
        return 1.0 / pow(smoothingRadius, 3) * f(q);
    }

    vec3 dW(Particle* x_i, Particle* x_j, float smoothingRadius) {
        float q = distance(x_i->getPosition(), x_j->getPosition()) / smoothingRadius;
        vec3 norm = x_i->getPosition() - x_j->getPosition();
        norm = norm / sqrt(dot(norm, norm));
        return 1.0 / pow(smoothingRadius, 4) * df(q) * norm;
    }

public:
    SPHSystem(vec3 pos, vec3 size, vec3 gap, float m_d, vec3 container_lb, vec3 container_ub, 
        float penalty, float k, float density0, float supportRadius, float smoothingRadius, float viscosity) :
        spos(pos), size(size), gap(gap), m_d(m_d), container_lb(container_lb), container_ub(container_ub),
        penalty(penalty), k(k), density0(density0), supportRadius(supportRadius), smoothingRadius(smoothingRadius), viscosity(viscosity)
    {
        this->table = new SpatialHashTable(m_d);
        float mass = pow(supportRadius / 2.0, 3) * density0; // pow(2.0 / 3.0 * smoothingRadius, 3) * density0;

        vec3 center((size.x-1) * gap.x / 2, (size.y - 1) * gap.y / 2, (size.z - 1) * gap.z / 2);
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    vec3 ppos(gap[0] * ix, gap[1] * iy, gap[2] * iz);

                    if (distance(ppos, center) > (size.x - 1) * gap.x / 2 + 0.1) continue;

                    vec3 rand(0.01 * rand() / RAND_MAX, 0.01 * rand() / RAND_MAX, 0.01 * rand() / RAND_MAX);
                    ppos = ppos + pos + rand;
                    Ps.push_back(new Particle(ppos, mass));
                }
            }
        }
    }

    void build() {
        table->clear();
        table->build(Ps);
        for (int i = 0; i < Ps.size(); i++) computeNeighbors(Ps[i]);
        for (int i = 0; i < Ps.size(); i++) computeBasics(Ps[i]);
        for (int i = 0; i < Ps.size(); i++) computeForces(Ps[i]);
    }

    void computeNeighbors(Particle* p) {
        table->neighbors(p, p->getNeighbors(false));
        p->setBuiltNeighborsFlag(true);
    }

    void computeBasics(Particle* p) {
        float mass = p->getMass();
        vec3 pos = p->getPosition();
        vector<Particle*>* neighbors = p->getNeighbors(true);
        // compute density & pressure
        float dens = 0;
        for (int j = 0; j < neighbors->size(); j++) {
            Particle* p_j = neighbors->at(j);
            dens += p_j->getMass() * W(p, p_j, smoothingRadius);
        }
        float pres = k * (pow(dens / density0, 7) - 1);
        p->setDensity(dens);
        p->setPressure(pres);
        p->setBuiltBasicsFlag(true);
    }

    void computeForces(Particle* p) {
        vec3 pos = p->getPosition();
        vec3 vel = p->getVelocity();
        float mass = p->getMass();
        float dens = p->getDensity();
        float pres = p->getPressure();
        vector<Particle*>* neighbors = p->getNeighbors(true);
        // -----------------------
        // compute F_pressure
        // -----------------------
        vec3 dPres(0, 0, 0);
        if (SHOW_FPRES) {
            vec3 tpos = p->getPosition();
            printf("stats of %.0f-%.0f-%.0f:\tmass: %.2f\tdensity: %.2f\tpressure: %.2f\n",
                tpos.x, tpos.y, tpos.z, mass, dens, pres);
        }
        for (int j = 0; j < neighbors->size(); j++) {
            Particle* p_j = neighbors->at(j);
            float mass_j = p_j->getMass();
            float dens_j = p_j->getDensity();
            float pres_j = p_j->getPressure();
            vec3 dW_ij = dW(p, p_j, smoothingRadius);
            vec3 dPres_j = mass_j * (pres / pow(dens, 2) + pres_j / pow(dens_j, 2)) * dW_ij;
            dPres += dPres_j;
            if (SHOW_FPRES) {
                vec3 tpos = p_j->getPosition();
                printf("\tstats of %.0f-%.0f-%.0f:\tmass: %.2f  density: %.2f  pressure: %.2f  param: %.2f  dW: %.2f-%.2f-%.2f  dPres: %.2f-%.2f-%.2f\n",
                    tpos.x, tpos.y, tpos.z,
                    mass_j, dens_j, pres_j,
                    mass_j * (pres / pow(dens, 2) + pres_j / pow(dens_j, 2)),
                    dW_ij.x, dW_ij.y, dW_ij.z,
                    dPres_j.x, dPres_j.y, dPres_j.z
                );
            }
        }
        vec3 Fpres = - mass * dPres;
        if (SHOW_FPRES) {
            printf("dpres: %.4f %.4f %.4f\n", dPres.x, dPres.y, dPres.z);
            printf("Fpres: %.4f %.4f %.4f\n", Fpres.x, Fpres.y, Fpres.z);
        }
        // -----------------------
        // compute F_viscosity
        // -----------------------
        vec3 ddVisc(0, 0, 0);
        if (SHOW_FVISC) {
            vec3 tpos = pos;
            printf("stats of %.0f-%.0f-%.0f:\tmass: %.2f\tvel: %.2f-%.2f-%.2f\n",
                tpos.x, tpos.y, tpos.z, mass, vel.x, vel.y, vel.z);
        }
        for (int j = 0; j < neighbors->size(); j++) {
            Particle* p_j = neighbors->at(j);
            float mass_j = p_j->getMass();
            float dens_j = p_j->getDensity();
            vec3 x_ij = pos - p_j->getPosition();
            vec3 v_ij = vel - p_j->getVelocity();
            vec3 dW_ij = dW(p, p_j, smoothingRadius);
            vec3 ddVisc_j = mass_j / dens_j * v_ij * (dot(x_ij, dW_ij) / (dot(x_ij, x_ij) + 0.01 * pow(smoothingRadius, 2)));
            ddVisc += ddVisc_j;

            if (SHOW_FVISC) {
                vec3 tpos = p_j->getPosition();
                printf("\tstats of %.0f-%.0f-%.0f:\tmass: %.2f  density: %.2f  dW: %.2f-%.2f-%.2f  ddVisc: %.2f-%.2f-%.2f\n",
                    tpos.x, tpos.y, tpos.z,
                    mass_j, dens_j,
                    dW_ij.x, dW_ij.y, dW_ij.z,
                    ddVisc_j.x, ddVisc_j.y, ddVisc_j.z
                );
            }
        }
        vec3 Fvisc = 2 * mass * viscosity * ddVisc;
        if (SHOW_FVISC) {
            printf("ddVisc: %.4f %.4f %.4f\n", ddVisc.x, ddVisc.y, ddVisc.z);
            printf("Fvisc: %.4f %.4f %.4f\n", Fvisc.x, Fvisc.y, Fvisc.z);
        }
        p->setFpres(Fpres);
        p->setFvisc(Fvisc);
        p->setBuiltForcesFlag(true);
    }

    void applyGForces() {
        vec3 g(0.0f, -98.0f, 0.0f);
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i]->applyPermForce(Ps[i]->getMass() * g);
    }

    void applyTempForces(vec3& f) {
        for (int i = 0; i < Ps.size(); i++)
            Ps[i]->applyTempForce(f);
    }   

    void clearTempForces() {
        for (int i = 0; i < Ps.size(); i++)
            Ps[i]->clearTempForce();
    }

    void applySPHForces() {
        for (int i = 0 ; i < Ps.size() ; i++) {
            Ps[i]->applyTempForce(Ps[i]->getFpres());
            Ps[i]->applyTempForce(Ps[i]->getFvisc());
        }
    }

    void setContainer(vec3 container_lb, vec3 container_ub) {
        this->container_lb = container_lb;
        this->container_ub = container_ub;
    }

    void applyPenaltyForces() {
        for (int i = 0; i < Ps.size(); i++) {
            vec3 pos = Ps[i]->getPosition();
            if (pos.y <= container_lb.y) {
                vec3 penalty_force(0, -penalty * (pos.y - container_lb.y), 0);
                Ps[i]->applyTempForce(penalty_force);
            }
            if (pos.y > container_ub.y) {
                vec3 penalty_force(0, -penalty * (pos.y - container_ub.y), 0);
                Ps[i]->applyTempForce(penalty_force);
            }
            if (pos.x <= container_lb.x) {
                vec3 penalty_force(-penalty * (pos.x - container_lb.x), 0, 0);
                Ps[i]->applyTempForce(penalty_force);
            }
            if (pos.x > container_ub.x) {
                vec3 penalty_force(-penalty * (pos.x - container_ub.x), 0, 0);
                Ps[i]->applyTempForce(penalty_force);
            }
            if (pos.z <= container_lb.z) {
                vec3 penalty_force(0, 0, -penalty * (pos.z - container_lb.z));
                Ps[i]->applyTempForce(penalty_force);
            }
            if (pos.z > container_ub.z) {
                vec3 penalty_force(0, 0, -penalty * (pos.z - container_ub.z));
                Ps[i]->applyTempForce(penalty_force);
            }
        }
    }

    void getAccelations(vector<vec3>& accs) {
        accs.clear();
        for (int i = 0; i < Ps.size(); i++)
            accs.push_back(Ps[i]->computeAcceleration());
    }

    int getSize() { return Ps.size(); }

    int getNumDofs() { return getSize(); }

    vector<Particle*>* getParticles() { return &Ps; }

    void getPositions(float* pos) {
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i]->getPosition();
            for (int j = 0 ; j < 3 ; j++) pos[3 * i + j] = ppos[j];
        }
    }

    void getPositions(vector<vec3>& pos) {
        pos.clear();
        for (int i = 0 ; i < Ps.size() ; i++) 
            pos.push_back(Ps[i]->getPosition());
    }

    void getVelocities(float* vel) {
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 pvel = Ps[i]->getVelocity();
            for (int j = 0 ; j < 3 ; j++) vel[3 * i + j] = pvel[j];
        }
    }

    void getVelocities(vector<vec3>& vel) {
        vel.clear();
        for (int i = 0 ; i < Ps.size() ; i++) 
            vel.push_back(Ps[i]->getVelocity());
    }

    void setPositions(const vector<vec3>& pos) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i]->setPosition(pos[i]);
    }

    void setVelocities(vec3& vel) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i]->setVelocity(vel);
    }

    void setVelocities(const vector<vec3>& vel) {
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i]->setVelocity(vel[i]);
    }

    void getPosAcc(float* posacc) {
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i]->getPosition();
            vec3 pacc = ppos + Ps[i]->computeTempAcceleration();
            for (int j = 0 ; j < 3 ; j++) posacc[6 * i + 0 + j] = ppos[j];
            for (int j = 0 ; j < 3 ; j++) posacc[6 * i + 3 + j] = pacc[j];
        }
    }

    void getPosAcc(vector<vec3>& posacc) {
        posacc.clear();
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i]->getPosition();
            vec3 pacc = ppos + Ps[i]->computeTempAcceleration();
            posacc.push_back(ppos);
            posacc.push_back(ppos + pacc);
        }
    }

    void printPositions() 
    {
        printf("\n");
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 pos = Ps[i]->getPosition();
            printf("%.3f-%.3f-%.3f ", pos[0], pos[1], pos[2]);
        }
    }
};

#endif