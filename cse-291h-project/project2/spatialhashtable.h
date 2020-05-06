#ifndef SPATIALHASHTABLE_H
#define SPATIALHASHTABLE_H

#include <vector>
#include <map>
#include <cmath> 
#include <algorithm>
#include "sphparticle.h"

using namespace std;
using namespace glm;

class int3{
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
};

class SpatialHashTable {
    
private:
    vec3 m_d;
    map<int3, vector<Particle>> hash_map;

    int3 hash_func(Particle p) {
        vec3 pos = p.getPosition();
        return hash_func(pos);
    }

    int3 hash_func(vec3 pos) {
        int3 ret(pos.x / m_d.x, pos.y / m_d.y, pos.z / m_d.z);
        if(pos.x < 0.0){ret.x--;}
        if(pos.y < 0.0){ret.y--;}
        if(pos.z < 0.0){ret.z--;}
        return ret;
    }

public:
    SpatialHashTable(vec3 m_d) {
        this->m_d = m_d;
        this->hash_map.clear();
    }

    void clear() {
        for (pair<int3, vector<Particle>> element : hash_map)
            element.second.clear();
        hash_map.clear();
    }

    void build(vector<Particle> ps) {
        for (Particle p : ps)
            hash_map[hash_func(p.getPosition())].push_back(p);
    }
    
    vector<Particle> neighbors(Particle p, float dist) {
        vec3 pos = p.getPosition();
        vector<Particle> res = neighbors(pos, dist);
        res.erase(find(res.begin(), res.end(), p));
        return res;
    }

    vector<Particle> neighbors(vec3 pos, float dist) {
        vector<Particle> res;
        vec3 lb(pos.x - dist, pos.y - dist, pos.z - dist);
        vec3 ub(pos.x + dist, pos.y + dist, pos.z + dist);
        int3 lbh = hash_func(lb), ubh = hash_func(ub);

        for (int ix = lbh.x ; ix < ubh.x + 1 ; ix++) {
            for (int iy = lbh.y ; iy < ubh.y + 1 ; iy++) {
                for (int iz = lbh.z ; iz < ubh.z + 1 ; iz++) {
                    int3 block_id(ix, iy, iz);
                    for (Particle p : hash_map[block_id])
                        if (distance(p.getPosition(), pos) <= dist)
                            res.push_back(p);
                }
            }
        }
    }
};

#endif