#ifndef SPHMAP_H
#define SPHMAP_H

#include <vector>
#include <map>
#include <cmath> 
#include <algorithm>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/scalar_multiplication.hpp>
#include "sphparticle.h"

#define LOG 0
#define SHOW_FPRES 0
#define SHOW_FVISC 0

using namespace glm;
using namespace std;

const float PI = 3.14159265358979323846264338327950288;

class int3 {
private:
    const int MAX_BID = 500;
    const int LFT_BID = 10;
public:
    int x, y, z;

    int3(int x, int y, int z) {
        if (abs(x) > MAX_BID || abs(y) > MAX_BID || abs(z) > MAX_BID)
            throw "int3 max size exceeded!";
        this->set(x, y, z);
    }

    void set(int x, int y, int z) {
        if (abs(x) > MAX_BID || abs(y) > MAX_BID || abs(z) > MAX_BID)
            throw "int3 max size exceeded!";
        this->x = x;
        this->y = y;
        this->z = z;
    }

    int hash_func() const {
        return ((x + MAX_BID) << (2 * LFT_BID)) + ((y + MAX_BID) << LFT_BID) + (z + MAX_BID);
    }

    bool operator< (const int3& y) const {
        return this->hash_func() < y.hash_func();
    }

    bool operator== (const int3& y) const {
        return this->hash_func() == y.hash_func();
    }
};

class SpatialHashTable {

private:
    vec3 m_d;
    map<int3, vector<Particle*>> hash_map;

    int3 hash_func(Particle* p) {
        vec3 pos = p->getPosition();
        return hash_func(pos);
    }

    int3 hash_func(vec3 pos) {
        int3 ret(pos.x / m_d.x, pos.y / m_d.y, pos.z / m_d.z);
        if (pos.x < 0.0) ret.x--;
        if (pos.y < 0.0) ret.y--;
        if (pos.z < 0.0) ret.z--;
        return ret;
    }

public:
    SpatialHashTable(vec3 m_d) {
        this->m_d = m_d;
        this->hash_map.clear();
    }

    void clear() {
        hash_map.clear();
    }

    void build(vector<Particle*> ps) {
        clear();
        for (Particle* p : ps)
            hash_map[hash_func(p->getPosition())].push_back(p);
    }

    vector<Particle*> neighbors(Particle* p, float dist) {
        vec3 pos = p->getPosition();
        if (LOG) printf("\tneighbors of %.0f-%.0f-%.0f: ", pos.x, pos.y, pos.z);
        vector<Particle*> res;
        vec3 lb(pos.x - dist, pos.y - dist, pos.z - dist);
        vec3 ub(pos.x + dist, pos.y + dist, pos.z + dist);
        int3 lbh = hash_func(lb), ubh = hash_func(ub);
        for (int ix = lbh.x; ix < ubh.x + 1; ix++) {
            for (int iy = lbh.y; iy < ubh.y + 1; iy++) {
                for (int iz = lbh.z; iz < ubh.z + 1; iz++) {
                    int3 block_id(ix, iy, iz);
                    for (Particle* tp : hash_map[block_id]) {
                        if (tp != p && distance(tp->getPosition(), pos) <= dist) {
                            res.push_back(tp);
                            vec3 tp_pos = tp->getPosition();
                            if (LOG) printf("%.0f-%.0f-%.0f ", tp_pos.x, tp_pos.y, tp_pos.z);
                        }
                    }
                }
            }
        }
        if (LOG) printf("\n");
        return res;
    }

    vector<Particle*> neighbors(Particle* p) {
        vec3 pos = p->getPosition();
        if (LOG) printf("\tneighbors of %.0f-%.0f-%.0f: ", pos.x, pos.y, pos.z);
        vector<Particle*> res;
        int3 h = hash_func(pos);
        for (int ix = h.x - 1; ix < h.x + 2; ix++) {
            for (int iy = h.y - 1; iy < h.y + 2; iy++) {
                for (int iz = h.z - 1; iz < h.z + 2; iz++) {
                    int3 block_id(ix, iy, iz);
                    for (Particle* tp : hash_map[block_id]) {
                        if (tp != p && distance(tp->getPosition(), pos) <= m_d.x) {
                            res.push_back(tp);
                            if (LOG) {
                                vec3 tp_pos = tp->getPosition();
                                printf("%.0f-%.0f-%.0f ", tp_pos.x, tp_pos.y, tp_pos.z);
                            }
                        }
                    }
                }
            }
        }
        if (LOG) printf("\n");
        return res;
    }
};

#endif