#ifndef SPHMAP_H
#define SPHMAP_H

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

#define LOG 0

using namespace glm;
using namespace std;


const int MAX_BID = 500;
const int LFT_BID = 10;

struct int3 {
    int x, y, z;

    int3(int x, int y, int z) {
        // if (abs(x) > MAX_BID || abs(y) > MAX_BID || abs(z) > MAX_BID) throw "int3 max size exceeded!";
        this->set(x, y, z);
    }

    void set(int x, int y, int z) {
        // if (abs(x) > MAX_BID || abs(y) > MAX_BID || abs(z) > MAX_BID) throw "int3 max size exceeded!";
        this->x = x;
        this->y = y;
        this->z = z;
    }

    int hash_func() const {
        return ((x + MAX_BID) << (2 * LFT_BID)) + ((y + MAX_BID) << LFT_BID) + (z + MAX_BID);
    }

    int3 operator+(int3 const& a) {
        return int3(x + a.x, y + a.y, z + a.z);
    }
};

class SpatialHashTable {

private:
    float m_d;
    unordered_map<int, vector<FluidParticle*>*> hash_map;

    int hash_func(vec3 pos) {
        int3 ret(pos.x / m_d, pos.y / m_d, pos.z / m_d);
        if (pos.x < 0.0) ret.x--;
        if (pos.y < 0.0) ret.y--;
        if (pos.z < 0.0) ret.z--;
        return ret.hash_func();
    }

    int3 int3_hash_func(vec3 pos) {
        int3 ret(pos.x / m_d, pos.y / m_d, pos.z / m_d);
        if (pos.x < 0.0) ret.x--;
        if (pos.y < 0.0) ret.y--;
        if (pos.z < 0.0) ret.z--;
        return ret;
    }

public:
    SpatialHashTable(float m_d) {
        this->m_d = m_d;
        this->hash_map.clear();
    }

    void clear() {
        hash_map.clear();
    }

    void build(vector<FluidParticle*> ps) {
        clear();
        for (FluidParticle* p : ps) {
            p->setBuiltBasicsFlag(false);
            p->setBuiltForcesFlag(false);
            p->setBuiltNeighborsFlag(false);
            int h = hash_func(p->getPosition());
            if (hash_map.find(h) == hash_map.end() || hash_map[h] == NULL)
                hash_map[h] = new vector<FluidParticle*>();
            hash_map[h]->push_back(p);
        }
    }

    void computeNeighborBlocks(FluidParticle* p) {
        auto res = p->getNeighborBlocks();
        res->clear();
        vec3 pos = p->getPosition();
        int3 h = int3_hash_func(pos);
        for (int ix = h.x - 1; ix < h.x + 2; ix++) {
            for (int iy = h.y - 1; iy < h.y + 2; iy++) {
                for (int iz = h.z - 1; iz < h.z + 2; iz++) {
                    int3 block_id(ix, iy, iz);
                    auto entry = hash_map[block_id.hash_func()];
                    if (entry != NULL) res->push_back(entry);
                }
            }
        }
    }

    void computeNeighbors(FluidParticle* p) {
        auto res = p->getNeighbors(false);
        res->clear();
        vec3 pos = p->getPosition();
        vector<vector<FluidParticle*>*>* blocks = p->getNeighborBlocks();
#pragma omp parallel for 
        for (int i = 0; i < blocks->size(); i++) {
            vector<FluidParticle*>* entry = blocks->at(i);
#pragma omp parallel for 
            for (int j = 0; j < entry->size(); j++) {
                if (entry->at(j) != p && distance(entry->at(j) ->getPosition(), pos) <= m_d)
#pragma omp critical
                    res->push_back(entry->at(j));
            }
        }
    }
};

#endif