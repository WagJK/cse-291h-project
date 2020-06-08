#ifndef SPHGRID_H
#define SPHGRID_H

#include <omp.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath> 
#include <algorithm>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>

#include "sphparticle.h"
#include "sphmap.h"
#include "sphparticle.h"
#include "polygonoise.h"

#define PI 3.14159265358979323846264338327950288


class SpatialGrid {

private:
    int imax, jmax, kmax;
    vec3 lb, ub;
    float g_d, smtRadius;

    float* densityTable; vec3* normTable;
    Particle** gridPs;

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

    float W(Particle* x_i, Particle* x_j) {
        float q = distance(x_i->getPosition(), x_j->getPosition()) / smtRadius;
        return 1.0 / pow(smtRadius, 3) * f(q);
    }

    vec3 dW(Particle* x_i, Particle* x_j) {
        float q = distance(x_i->getPosition(), x_j->getPosition()) / smtRadius;
        vec3 norm = x_i->getPosition() - x_j->getPosition();
        norm = norm / sqrt(dot(norm, norm));
        return 1.0 / pow(smtRadius, 4) * df(q) * norm;
    }

    vec3 XYZ2Pos(int3 xyz) {
        return vec3(lb.x + xyz.x * g_d, lb.y + xyz.y * g_d, lb.z + xyz.z * g_d);
    }

    int3 index2XYZ(int index) {
        return int3(index / jmax / kmax, index / kmax % jmax, index % kmax);
    }

    vec3 index2Pos(int index) {
        return XYZ2Pos(index2XYZ(index));
    }

    int XYZ2Index(int3 xyz) {
        return xyz.x * jmax * kmax + xyz.y * kmax + xyz.z;
    }

public:
    SpatialGrid(vec3 lb, vec3 ub, float g_d, float smtRadius) {
        this->lb = lb;
        this->ub = ub;
        this->g_d = g_d;
        this->smtRadius = smtRadius;

        this->imax = (ub.x - lb.x) / g_d + 1;
        this->jmax = (ub.y - lb.y) / g_d + 1;
        this->kmax = (ub.z - lb.z) / g_d + 1;

        densityTable = new float[imax * jmax * kmax];
        normTable = new vec3[imax * jmax * kmax];
        gridPs = new Particle * [imax * jmax * kmax];
        for (int i = 0; i < imax * jmax * kmax; i++) {
            int3 xyz = index2XYZ(i);
            vec3 pos = XYZ2Pos(xyz);

            gridPs[i] = new Particle(pos, 0.0f);
            densityTable[i] = 0.0f;
            normTable[i] = vec3(0.0f, 0.0f, 0.0f);
        }
    }

    void computeDensityTable(SpatialHashTable* hashTable) {
        for (int i = 0; i < imax * jmax * kmax; i++) {
            hashTable->computeNeighborBlocks(gridPs[i]);
        }
        #pragma omp parallel for
        for (int i = 0; i < imax * jmax * kmax; i++) {
            hashTable->computeNeighbors(gridPs[i]);
            gridPs[i]->setBuiltNeighborsFlag(true);
        }
        #pragma omp parallel for
        for (int i = 0; i < imax * jmax * kmax; i++) {
            float dens = 0; vec3 norm(0, 0, 0);
            vector<Particle*>* neighbors = gridPs[i]->getNeighbors(true);
            #pragma omp parallel for reduction(+:dens)
            for (int j = 0; j < neighbors->size(); j++) {
                Particle* p_j = neighbors->at(j);
                dens += p_j->getMass() * W(gridPs[i], p_j);
                norm += p_j->getMass() * dW(gridPs[i], p_j);
            }
            densityTable[i] = dens;
            normTable[i] = -norm;
        }
    }

    int computeTriangleFacets(float* res, float isolevel) {

        vector<vec3> temp;
        int numTriangle = 0;
        int cube_edgeToVerts[12][2] = {
            {0,1}, {1,2}, {2,3}, {3,0},
            {4,5}, {5,6}, {6,7}, {7,4},
            {0,4}, {1,5}, {2,6}, {3,7},
        };

        for (int i = 0 ; i < imax * jmax * kmax ; i++) {
            int3 xyz = index2XYZ(i);
            vec3 pos = XYZ2Pos(xyz);
            if (xyz.x == imax - 1 || xyz.y == jmax - 1 || xyz.z == kmax - 1) continue;
            
            int3 indices[8] = {
                xyz + int3(0,0,0), xyz + int3(0,1,0), xyz + int3(1,1,0), xyz + int3(1,0,0),
                xyz + int3(0,0,1), xyz + int3(0,1,1), xyz + int3(1,1,1), xyz + int3(1,0,1)
            };
            float dens[8]; vec3 poss[8]; vec3 norm[8];
            for (int j = 0 ; j < 8 ; j++) {
                poss[j] = XYZ2Pos(indices[j]);
                dens[j] = densityTable[XYZ2Index(indices[j])];
                norm[j] = normTable[XYZ2Index(indices[j])];
            }
            numTriangle += Polygonise(poss, norm, dens, isolevel, &temp);
        }
        for (int i = 0; i < temp.size(); i++) {
            res[3 * i + 0] = temp[i].x;
            res[3 * i + 1] = temp[i].y;
            res[3 * i + 2] = temp[i].z;
        }
        return numTriangle;
    }
};
#endif