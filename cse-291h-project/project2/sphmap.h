#ifndef SPHMAP_H
#define SPHMAP_H

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
};

class SpatialHashTable {

private:
    float m_d;
    unordered_map<int, vector<Particle*>*> hash_map;

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

    void build(vector<Particle*> ps) {
        clear();
        for (Particle* p : ps) {
            p->setBuiltBasicsFlag(false);
            p->setBuiltForcesFlag(false);
            p->setBuiltNeighborsFlag(false);
            int h = hash_func(p->getPosition());
            if (hash_map.find(h) == hash_map.end() || hash_map[h] == NULL)
                hash_map[h] = new vector<Particle*>();
            hash_map[h]->push_back(p);
        }
    }

    void neighbors(Particle* p, float dist, vector<Particle*>* res) {
        res->clear();
        vec3 pos = p->getPosition();
        if (LOG) printf("\tneighbors of %.0f-%.0f-%.0f: ", pos.x, pos.y, pos.z);
        vec3 lb(pos.x - dist, pos.y - dist, pos.z - dist);
        vec3 ub(pos.x + dist, pos.y + dist, pos.z + dist);
        int3 lbh = int3_hash_func(lb), ubh = int3_hash_func(ub);
        for (int ix = lbh.x; ix < ubh.x + 1; ix++) {
            for (int iy = lbh.y; iy < ubh.y + 1; iy++) {
                for (int iz = lbh.z; iz < ubh.z + 1; iz++) {
                    int3 block_id(ix, iy, iz);
                    vector<Particle*>* entry = hash_map[block_id.hash_func()];
                    if (entry == NULL) continue;
                    for (Particle* tp : *entry) {
                        if (tp != p && distance(tp->getPosition(), pos) <= dist) {
                            res->push_back(tp);
                            vec3 tp_pos = tp->getPosition();
                            if (LOG) printf("%.0f-%.0f-%.0f ", tp_pos.x, tp_pos.y, tp_pos.z);
                        }
                    }
                }
            }
        } if (LOG) printf("\n");
    }

    void neighbors(Particle* p, vector<Particle*>* res) {
        res->clear();
        vec3 pos = p->getPosition();
        if (LOG) printf("\tneighbors of %.0f-%.0f-%.0f: ", pos.x, pos.y, pos.z);
        int3 h = int3_hash_func(pos);
        for (int ix = h.x - 1; ix < h.x + 2; ix++) {
            for (int iy = h.y - 1; iy < h.y + 2; iy++) {
                for (int iz = h.z - 1; iz < h.z + 2; iz++) {
                    int3 block_id(ix, iy, iz);
                    vector<Particle*>* entry = hash_map[block_id.hash_func()];
                    if (entry == NULL) continue;
                    for (Particle* tp : *entry) {
                        if (tp != p && distance(tp->getPosition(), pos) <= m_d) {
                            res->push_back(tp);
                            if (LOG) {
                                vec3 tp_pos = tp->getPosition();
                                printf("%.0f-%.0f-%.0f ", tp_pos.x, tp_pos.y, tp_pos.z);
                            }
                        }
                    }
                }
            }
        } if (LOG) printf("\n");
    }
};

#endif