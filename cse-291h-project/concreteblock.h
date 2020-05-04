#ifndef CONCRETEBLOCK_H
#define CONCRETEBLOCK_H

#include <vector>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/scalar_multiplication.hpp>
#include "tetrahedron.h"

using namespace glm;
using namespace std;

class ConcreteBlockSystem {
public:
    ConcreteBlockSystem(vec3 pos, vec3 size, vec3 gap, float mass, float E, float v) :
        pos(pos), size(size), gap(gap), mass(mass), E(E), v(v)
    {

        for (int ix = 0; ix < size[0]; ix++) {
            Ps.push_back(vector<vector<Particle>>());
            for (int iy = 0; iy < size[1]; iy++) {
                Ps[ix].push_back(vector<Particle>());
                for (int iz = 0; iz < size[2]; iz++) {
                    
                    mat4 rotate_matrix(1.0f);
                    rotate_matrix = rotate(rotate_matrix, radians(0.0f), vec3(0, 0, 1));
                    rotate_matrix = rotate(rotate_matrix, radians(0.0f), vec3(1, 0, 0));
                    vec3 ppos(gap[0] * ix, gap[1] * iy, gap[2] * iz);
                    ppos = rotate_matrix * vec4(ppos, 1);
                    ppos = ppos + pos;

                    Particle P(ppos, 0); 
                    Ps[ix][iy].push_back(P);
                }
            }
        }
        for (int ix = 0; ix < size[0] - 1; ix++) {
            for (int iy = 0; iy < size[1] - 1; iy++) {
                for (int iz = 0; iz < size[2] - 1; iz++) {

                    // 0,0,0 - 1,0,0 - 0,1,0 - 0,0,1
                    Tetrahedron top1(&Ps[ix][iy][iz], &Ps[ix + 1][iy][iz], &Ps[ix][iy + 1][iz], &Ps[ix][iy][iz + 1], mass, E, v);

                    // 0,1,1 - 0,1,0 - 1,1,1 - 0,0,1
                    Tetrahedron top2(&Ps[ix][iy + 1][iz + 1], &Ps[ix][iy + 1][iz], &Ps[ix + 1][iy + 1][iz + 1], &Ps[ix][iy][iz + 1], mass, E, v);

                    // 1,1,0 - 1,1,1 - 0,1,0 - 1,0,0
                    Tetrahedron bot1(&Ps[ix + 1][iy + 1][iz], &Ps[ix + 1][iy + 1][iz + 1], &Ps[ix][iy + 1][iz], &Ps[ix + 1][iy][iz], mass, E, v);

                    // 1,0,1 - 0,0,1 - 1,1,1 - 1,0,0
                    Tetrahedron bot2(&Ps[ix + 1][iy][iz + 1], &Ps[ix][iy][iz + 1], &Ps[ix + 1][iy + 1][iz + 1], &Ps[ix + 1][iy][iz], mass, E, v);

                    // center: 1,0,0 - 0,1,0 - 0,0,1 - 1,1,1
                    Tetrahedron cntr(&Ps[ix + 1][iy][iz], &Ps[ix][iy + 1][iz], &Ps[ix][iy][iz + 1], &Ps[ix + 1][iy + 1][iz + 1], 2 * mass, E, v);

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

    void GetParticles(vector<Particle*>& ps) {
        ps.clear();
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    ps.push_back(&Ps[ix][iy][iz]);
    }

    int GetNumLines() {
        int cnt_lines = 0;
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    if (ix != size[0] - 1) cnt_lines++;
                    if (iy != size[1] - 1) cnt_lines++;
                    if (iz != size[2] - 1) cnt_lines++;

                    // if (ix != size[0] - 1 && iy != size[1] - 1) cnt_lines += 2;
                    // if (iy != size[1] - 1 && iz != size[2] - 1) cnt_lines += 2;
                    // if (ix != size[0] - 1 && iz != size[2] - 1) cnt_lines += 2;
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
                    /*if (ix != size[0] - 1 && iy != size[1] - 1) {
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
                    }*/
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

    void GetAccelations(float* accs) {
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    vec3 ppos = Ps[ix][iy][iz].GetPosition();
                    vec3 pacc = ppos + Ps[ix][iy][iz].ComputeTempAcceleration();
                    int idx = 6 * (ix * size[1] * size[2] + iy * size[2] + iz);
                    accs[idx + 0] = ppos[0];
                    accs[idx + 1] = ppos[1];
                    accs[idx + 2] = ppos[2];
                    accs[idx + 3] = pacc[0];
                    accs[idx + 4] = pacc[1];
                    accs[idx + 5] = pacc[2];
                }
            }
        }
    }

    void GetAccelations(vector<vec3>& accs) {
        accs.clear();
        for (int ix = 0; ix < size[0]; ix++) {
            for (int iy = 0; iy < size[1]; iy++) {
                for (int iz = 0; iz < size[2]; iz++) {
                    accs.push_back(Ps[ix][iy][iz].GetPosition());
                    accs.push_back(Ps[ix][iy][iz].GetPosition() + Ps[ix][iy][iz].ComputeAcceleration());
                }
            }
        }
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

    void ApplyGForces() {
        vec3 g(0.0f, -9.8f, 0.0f);
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].ApplyPermForce(Ps[ix][iy][iz].GetMass() * g);
    }

    void ApplyTempForces(vec3& f) {
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].ApplyTempForce(f);
    }

    void ClearTempForces() {
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    Ps[ix][iy][iz].ClearTempForce();
    }

    void ComputeForces() {
        for (auto tthd : Ts) tthd.ComputeForces();
    }

    void ComputeAccelerations(vector<vec3>& acc) {
        acc.clear();
        for (int ix = 0; ix < size[0]; ix++)
            for (int iy = 0; iy < size[1]; iy++)
                for (int iz = 0; iz < size[2]; iz++)
                    acc.push_back(Ps[ix][iy][iz].ComputeAcceleration());
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

#endif