#ifndef SPH_H
#define SPH_H

#include <omp.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/scalar_multiplication.hpp>
#include "sphparticle.h"
#include "sphmap.h"
#include "sphgrid.h"

#define LOG 0
#define SHOW_FPRES 0
#define SHOW_FVISC 0
#define PI 3.14159265358979323846264338327950288f

using namespace std;
using namespace glm;

class SPHSystem {

private:
    SpatialHashTable* table;
    SpatialGrid* grid;

    vector<FluidParticle*> Ps;
    vector<DiffuseParticle*> DPs;

    vec3 spos;
    vec3 size;
    vec3 gap;
    vec3 m_d;
    vec3 g_d;
    vec3 container_lb;
    vec3 container_ub;
    vec3 Imax;

    float k;
    float g;
    float dt;
    float density0;
    float sptRadius;
    float smtRadius;
    float penalty;
    float itgPenalty;
    float viscosity;
    float surfaceTension;
    float norm_h;
    float kdiffuse, kb, kd;

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

    float W(FluidParticle* x_i, FluidParticle* x_j) {
        float q = distance(x_i->getPosition(), x_j->getPosition()) / smtRadius;
        return 1.0 / pow(smtRadius, 3) * f(q);
    }

    vec3 dW(FluidParticle* x_i, FluidParticle* x_j) {
        float q = distance(x_i->getPosition(), x_j->getPosition()) / smtRadius;
        vec3 norm = x_i->getPosition() - x_j->getPosition();
        norm = norm / sqrt(dot(norm, norm));
        return 1.0 / pow(smtRadius, 4) * df(q) * norm;
    }

    float C(float r) {
        if (2 * r > smtRadius && r <= smtRadius)
            return 32.0f / PI / pow(smtRadius, 9) * (pow(smtRadius - r, 3) * pow(r, 3));
        else if (r > 0 && 2 * r <= smtRadius)
            return 32.0f / PI / pow(smtRadius, 9) * (2 * pow(smtRadius - r, 3) * pow(r, 3) - pow(smtRadius, 6) / 64.0f);
        else
            return 0.0f;
    }

    float Wdiffuse(float length_x_ij, float influenceRadius) {
        if (length_x_ij <= influenceRadius) {
            return (1.0f - length_x_ij / influenceRadius) * 3.0f / (PI * pow(influenceRadius, 3));
        } else return 0;
    }



public:
    SPHSystem(
        vec3 pos, vec3 vel, vec3 size, vec3 gap,
        float dt, float m_d, float g_d,
        vec3 container_lb, vec3 container_ub,
        float g, float penalty, float itgPenalty, float k, float density0,
        float sptRadius, float smtRadius, float viscosity, float surfaceTension, float norm_h,
        float kdiffuse, float kb, float kd, vec3 Imax
    ) :
        spos(pos), size(size), gap(gap), dt(dt), m_d(m_d), g_d(g_d),
        container_lb(container_lb), container_ub(container_ub),
        g(g), penalty(penalty), itgPenalty(itgPenalty), k(k), density0(density0), 
        sptRadius(sptRadius), smtRadius(smtRadius), viscosity(viscosity), surfaceTension(surfaceTension), norm_h(norm_h),
        kdiffuse(kdiffuse), kb(kb), kd(kd), Imax(Imax)
    {
        this->table = new SpatialHashTable(m_d);
        this->grid = new SpatialGrid(container_lb, container_ub, g_d, smtRadius);

        float mass = pow(sptRadius / 2.0, 3) * density0; 
                    // pow(2.0 / 3.0 * smtRadius, 3) * density0;

        vec3 center((size.x-1) * gap.x / 2, (size.y - 1) * gap.y / 2, (size.z - 1) * gap.z / 2);
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    vec3 ppos(gap[0] * ix, gap[1] * iy, gap[2] * iz);
                    if (distance(ppos, center) > (size.x - 1) * gap.x / 2 + 0.1) continue;
                    vec3 rand(0.02 * rand() / RAND_MAX, 0.02 * rand() / RAND_MAX, 0.02 * rand() / RAND_MAX);
                    ppos = ppos + pos + rand;
                    FluidParticle* p = new FluidParticle(ppos, mass);
                    p->setVelocity(vel);
                    Ps.push_back(p);
                }
            }
        }
    }

    void integrate() {
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++)
            Ps[i]->integrate(dt, itgPenalty);
#pragma omp parallel for
        for (int i = 0; i < DPs.size(); i++)
            DPs[i]->integrate(dt, itgPenalty);
    }

    void build() {
        table->clear();

        table->build(Ps);

        maintainFluidParticles();

        computeDiffusePotentials();

        maintainDiffuseParticles();
        
        grid->computeDensityTable(table);

        printf("Particles:%d\tDiffuse:%d\n", Ps.size(), DPs.size());
    }



    void maintainFluidParticles() {
        for (int i = 0; i < Ps.size(); i++) computeNeighborBlocks(Ps[i]);
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++) computeNeighbors(Ps[i]);
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++) computeBasics(Ps[i]);
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++) computeNormals(Ps[i]);
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++) computeForces(Ps[i]);
    }

    void computeNeighborBlocks(FluidParticle* p) {
        table->computeNeighborBlocks(p);
    }

    void computeNeighbors(FluidParticle* p) {
        table->computeNeighbors(p);
        p->setBuiltNeighborsFlag(true);
    }

    void computeBasics(FluidParticle* p) {
        float mass = p->getMass();
        vec3 pos = p->getPosition();
        vector<FluidParticle*>* neighbors = p->getNeighbors(true);
        // compute density & pressure
        float dens = 0; 
        for (int j = 0; j < neighbors->size(); j++) {
            FluidParticle* p_j = neighbors->at(j);
            dens += p_j->getMass() * W(p, p_j);
        }
        p->setDensity(dens);
        float pres = k * (pow(dens / density0, 7) - 1);
        p->setPressure(pres);
        p->setBuiltBasicsFlag(true);
    }

    void computeNormals(FluidParticle* p) {
        vector<FluidParticle*>* neighbors = p->getNeighbors(true);
        vec3 norm = { 0, 0, 0 };
        for (int j = 0; j < neighbors->size(); j++) {
            FluidParticle* p_j = neighbors->at(j);
            float mass_j = p_j->getMass();
            float dens_j = p_j->getDensity();
            norm += norm_h * mass_j / dens_j * dW(p, p_j);
            // norm += p_j->getMass() * dW(p, p_j);
        }
        // printf("%.2f %.2f %.2f\n", norm.x, norm.y, norm.z);
        p->setNormal(norm);
    }

    void computeForces(FluidParticle* p) {
        vec3 pos = p->getPosition();
        vec3 vel = p->getVelocity();
        float mass = p->getMass();
        float dens = p->getDensity();
        float pres = p->getPressure();
        vector<FluidParticle*>* neighbors = p->getNeighbors(true);
        // -----------------------
        // compute F_pressure
        // -----------------------
        vec3 dPres(0, 0, 0);
        if (SHOW_FPRES) {
            vec3 tpos = p->getPosition();
            printf("stats of %.0f-%.0f-%.0f:\tmass: %.2f\tdensity: %.2f\tpressure: %.2f\n",
                tpos.x, tpos.y, tpos.z, mass, dens, pres);
        }
#pragma omp parallel for
        for (int j = 0; j < neighbors->size(); j++) {
            FluidParticle* p_j = neighbors->at(j);
            float mass_j = p_j->getMass();
            float dens_j = p_j->getDensity();
            float pres_j = p_j->getPressure();
            vec3 dW_ij = dW(p, p_j);
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
        p->setFpres(Fpres);

        // -----------------------
        // compute F_viscosity
        // -----------------------
        //vec3 ddVisc(0, 0, 0);
        //if (SHOW_FVISC) {
        //    vec3 tpos = pos;
        //    printf("stats of %.0f-%.0f-%.0f:\tmass: %.2f\tvel: %.2f-%.2f-%.2f\n",
        //        tpos.x, tpos.y, tpos.z, mass, vel.x, vel.y, vel.z);
        //}
        //#pragma omp parallel for
        //for (int j = 0; j < neighbors->size(); j++) {
        //    FluidParticle* p_j = neighbors->at(j);
        //    float mass_j = p_j->getMass();
        //    float dens_j = p_j->getDensity();
        //    vec3 x_ij = pos - p_j->getPosition();
        //    vec3 v_ij = vel - p_j->getVelocity();
        //    vec3 dW_ij = dW(p, p_j);
        //    vec3 ddVisc_j = mass_j / dens_j * v_ij * (dot(x_ij, dW_ij) / (dot(x_ij, x_ij) + 0.01 * pow(smtRadius, 2)));
        //    ddVisc += ddVisc_j;

        //    if (SHOW_FVISC) {
        //        vec3 tpos = p_j->getPosition();
        //        printf("\tstats of %.0f-%.0f-%.0f:\tmass: %.2f  density: %.2f  dW: %.2f-%.2f-%.2f  ddVisc: %.2f-%.2f-%.2f\n",
        //            tpos.x, tpos.y, tpos.z,
        //            mass_j, dens_j,
        //            dW_ij.x, dW_ij.y, dW_ij.z,
        //            ddVisc_j.x, ddVisc_j.y, ddVisc_j.z
        //        );
        //    }
        //}
        //vec3 Fvisc = 2 * mass * viscosity * ddVisc;
        //if (SHOW_FVISC) {
        //    printf("ddVisc: %.4f %.4f %.4f\n", ddVisc.x, ddVisc.y, ddVisc.z);
        //    printf("Fvisc: %.4f %.4f %.4f\n", Fvisc.x, Fvisc.y, Fvisc.z);
        //}

        // -----------------------
        // compute F_cohesion F_curvature
        // -----------------------
        vec3 Fcohs(0, 0, 0), Fcurv(0, 0, 0);
        for (int j = 0; j < neighbors->size(); j++) {
            FluidParticle* p_j = neighbors->at(j);
            vec3 pos_j = p_j->getPosition();
            float mass_j = p_j->getMass();
            float dens_j = p_j->getDensity();

            float K_ij = 2 * density0 / (dens + dens_j);
            Fcohs += mass_j * C(length(pos - pos_j)) * normalize(pos - pos_j) * K_ij;
            Fcurv += (p->getNormal() - p_j->getNormal()) * K_ij;
        }
        Fcohs *= -surfaceTension * mass;
        Fcurv *= -surfaceTension * mass;
        p->setFvisc(Fcurv + Fcohs);
        p->setBuiltForcesFlag(true);
    }

    void computeDiffusePotential(FluidParticle* p) {
        vector<FluidParticle*>* neighbors = p->getNeighbors(true);
        vec3 x_i = p->getPosition();
        vec3 v_i = p->getVelocity();
        vec3 n_i = p->getNormal();

        float vi_diff = 0, ki = 0;
        // Energy
        float Eki = 0.5 * p->getMass() * dot(v_i, v_i);
        for (int j = 0 ; j < neighbors->size() ; j++) {
            FluidParticle* p_j = neighbors->at(j);
            vec3 x_j = p_j->getPosition();
            vec3 v_j = p_j->getVelocity();
            vec3 n_j = p_j->getNormal();

            vec3 v_ij = v_i - v_j, x_ij = x_i - x_j;
            float W_ij = Wdiffuse(length(x_ij), smtRadius);

            // Trapped Air
            vi_diff += length(v_ij) * (1.0f - dot(normalize(v_ij), normalize(x_ij))) * W_ij;

            // Wave Crest
            if (dot(-x_ij, n_i) < 0 && dot(normalize(v_i), normalize(n_i)) >= 0.6)
                ki += (1.0f - dot(normalize(n_i), normalize(n_j))) * W_ij;
        }
        p->setDiffusePotential(vi_diff, ki, Eki);
    }

    float clamp(float I, float Tmin, float Tmax) {
        return (std::min(I, Tmax) - std::min(I, Tmin)) / (Tmax - Tmin);
    }

    void computeDiffusePotentials() {
        #pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++)
            computeDiffusePotential(Ps[i]);

        //vec3 Imax(0, 0, 0);
        //for (int i = 0; i < Ps.size(); i++) {
        //    vec3 diffusePotential = Ps[i]->getDiffusePotential();
        //    for (int j = 0 ; j < 3 ; j++)
        //        Imax[j] = std::max(diffusePotential[j], Imax[j]);
        //}
        //printf("%.4f %.4f %.4f\n", Imax.x, Imax.y, Imax.z);


        for (int i = 0; i < Ps.size(); i++) {
            vec3 I = Ps[i]->getDiffusePotential();
            Ps[i]->setDiffusePotential(
                clamp(I.x, 0.1f * Imax.x, Imax.x),
                clamp(I.y, 0.1f * Imax.y, Imax.y),
                clamp(I.z, 0.1f * Imax.z, Imax.z)
            );
        }
    }



    void applyImplicitViscosity(FluidParticle* p) {
        vec3 vel = p->getVelocity();
        vector<FluidParticle*>* neighbors = p->getNeighbors(true);

        vec3 dVel(0, 0, 0);
#pragma omp parallel for
        for (int j = 0; j < neighbors->size(); j++) {
            FluidParticle* p_j = neighbors->at(j);
            vec3 vel_j = p_j->getVelocity();
            dVel += viscosity * W(p, p_j) * (vel_j - vel);
        }
        p->setVelocity(vel + dVel);
    }

    void applyImplicitViscosity() {
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++) applyImplicitViscosity(Ps[i]);
    }
    
    

    void maintainDiffuseParticles() {
        auto it = DPs.begin();
        while (it != DPs.end()) {
            if ((*it)->reduceLifetime()) it = DPs.erase(it); else it++;
        }
        for (int i = 0; i < Ps.size(); i++) {
            vec3 I = Ps[i]->getDiffusePotential();
            int n = dt* kdiffuse* I.z* (I.x + I.y);
            // printf("%.2f ", dt * kdiffuse * I.z * (I.x + I.y));
            for (int k = 0; k < n; k++) {
                vec3 pos = Ps[i]->getPosition() + 
                    vec3(0.02 * rand() / RAND_MAX, 0.02 * rand() / RAND_MAX, 0.02 * rand() / RAND_MAX);
                DiffuseParticle* p = new DiffuseParticle(pos, Ps[i]->getMass(), 0, 20);
                DPs.push_back(p);
            }
        }
        for (int i = 0; i < DPs.size(); i++) computeNeighborBlocks(DPs[i]);
#pragma omp parallel for
        for (int i = 0; i < DPs.size(); i++) computeNeighbors(DPs[i]);
#pragma omp parallel for
        for (int i = 0; i < DPs.size(); i++) computeDiffuseType(DPs[i]);
    }

    vec3 computeDiffuseNeighborAvgVel(DiffuseParticle* dp) {
        vec3 res(0, 0, 0);
        float totalW = 0.0f;
        vector<FluidParticle*>* neighbors = dp->getNeighbors(true);
        for (int i = 0; i < neighbors->size(); i++) {
            float Wi = W(dp, neighbors->at(i));
            res += neighbors->at(i)->getVelocity() * Wi;
            totalW += Wi;
        }
        return res / totalW;
    }

    vec3 computeDiffuseFexp(DiffuseParticle* dp) {
        return vec3(0, 0, 0);
    }

    void computeDiffuseType(DiffuseParticle* dp) {
        int numNeighbors = dp->getNeighbors(true)->size();
        if (numNeighbors < 6) 
            dp->setDiffuseType(DiffuseParticle::TYPE_SPRY);
        else if (numNeighbors > 20) 
            dp->setDiffuseType(DiffuseParticle::TYPE_BUBB);
        else 
            dp->setDiffuseType(DiffuseParticle::TYPE_FOAM);
    }

    void applyDiffuseForces(DiffuseParticle* dp) {
        vec3 acc(0, 0, 0);
        if (dp->getDiffuseType() == DiffuseParticle::TYPE_SPRY) {
            acc = (computeDiffuseFexp(dp) + vec3(0, -g, 0)) / dp->getMass();
            dp->applyTempForce(acc * dp->getMass());
        }
        else if (dp->getDiffuseType() == DiffuseParticle::TYPE_FOAM) {
            dp->setVelocity(computeDiffuseNeighborAvgVel(dp));
        }
        else if (dp->getDiffuseType() == DiffuseParticle::TYPE_BUBB) {
            acc = -kb * vec3(0, -g, 0) + kd * (computeDiffuseNeighborAvgVel(dp) - dp->getVelocity());
            dp->applyTempForce(acc * dp->getMass());
        }
    }



    void applyGForces() {
        vec3 g(0.0f, -g, 0.0f);
#pragma omp parallel for
        for (int i = 0 ; i < Ps.size() ; i++)
            Ps[i]->applyPermForce(Ps[i]->getMass() * g);
    }

    void applyTempForces(vec3& f) {
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++)
            Ps[i]->applyTempForce(f);
    }   

    void clearTempForces() {
#pragma omp parallel for
        for (int i = 0; i < Ps.size(); i++)
            Ps[i]->clearTempForce();
#pragma omp parallel for
        for (int i = 0; i < DPs.size(); i++)
            DPs[i]->clearTempForce();
    }

    void applySPHForces() {
#pragma omp parallel for
        for (int i = 0 ; i < Ps.size() ; i++) {
            Ps[i]->applyTempForce(Ps[i]->getFpres());
            Ps[i]->applyTempForce(Ps[i]->getFvisc());
        }
#pragma omp parallel for
        for (int i = 0; i < DPs.size(); i++)
            applyDiffuseForces(DPs[i]);
    }

    void setContainer(vec3 container_lb, vec3 container_ub) {
        this->container_lb = container_lb;
        this->container_ub = container_ub;
    }

    void applyPenaltyForces() {
        #pragma omp parallel for
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
        #pragma omp parallel for
        for (int i = 0; i < DPs.size(); i++) {
            vec3 pos = DPs[i]->getPosition();
            if (pos.y <= container_lb.y) {
                vec3 penalty_force(0, -penalty * (pos.y - container_lb.y), 0);
                DPs[i]->applyTempForce(penalty_force);
            }
            if (pos.y > container_ub.y) {
                vec3 penalty_force(0, -penalty * (pos.y - container_ub.y), 0);
                DPs[i]->applyTempForce(penalty_force);
            }
            if (pos.x <= container_lb.x) {
                vec3 penalty_force(-penalty * (pos.x - container_lb.x), 0, 0);
                DPs[i]->applyTempForce(penalty_force);
            }
            if (pos.x > container_ub.x) {
                vec3 penalty_force(-penalty * (pos.x - container_ub.x), 0, 0);
                DPs[i]->applyTempForce(penalty_force);
            }
            if (pos.z <= container_lb.z) {
                vec3 penalty_force(0, 0, -penalty * (pos.z - container_lb.z));
                DPs[i]->applyTempForce(penalty_force);
            }
            if (pos.z > container_ub.z) {
                vec3 penalty_force(0, 0, -penalty * (pos.z - container_ub.z));
                DPs[i]->applyTempForce(penalty_force);
            }
        }
    }

    vector<FluidParticle*>* getFluidParticles() { return &Ps; }

    vector<DiffuseParticle*>* getDiffuseParticles() { return &DPs;  }

    int getPositions(float* pos) {
#pragma omp parallel for
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i]->getPosition();
            for (int j = 0 ; j < 3 ; j++) pos[3 * i + j] = ppos[j];
        }
        return Ps.size();
    }

    void getVelocities(float* vel) {
#pragma omp parallel for
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 pvel = Ps[i]->getVelocity();
            for (int j = 0 ; j < 3 ; j++) vel[3 * i + j] = pvel[j];
        }
    }

    int getDiffusePositions(float* pos, int type) {
        int cnt = 0;
        for (int i = 0; i < DPs.size(); i++) {
            if (DPs[i]->getDiffuseType() != type) continue;
            vec3 ppos = DPs[i]->getPosition();
            for (int j = 0; j < 3; j++) pos[3 * cnt + j] = ppos[j];
            cnt++;
        }
        return cnt;
    }

    void getPosAcc(float* posacc) {
#pragma omp parallel for
        for (int i = 0 ; i < Ps.size() ; i++) {
            vec3 ppos = Ps[i]->getPosition();
            vec3 pacc = ppos + Ps[i]->computeTempAcceleration();
            for (int j = 0 ; j < 3 ; j++) posacc[6 * i + 0 + j] = ppos[j];
            for (int j = 0 ; j < 3 ; j++) posacc[6 * i + 3 + j] = pacc[j];
        }
    }

    int getTriangleFacets(float* pos, float isolevel) {
        return grid->computeTriangleFacets(pos, isolevel);
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