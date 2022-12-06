#pragma once
#include "polyhedron.h"
#include "Poliline.h"
#include "icMatrix.H"
#include "drawUtil.h"
#include <cmath>
#include <tuple>
# define M_PI           3.14159265358979323846  /* pi */
# define MIN_K          0.05
# define STEP           0.03



extern Polyhedron* poly;
extern std::vector<POLYLINE> polylines;
extern bool singularity_IBVF;
double npn	=	256; //orig was 64 -> 256
double scale = 	4.0;
extern int win_width;
extern int win_height;
extern float tmax;
extern float dmax;


struct Singularity{
    int type = -1;                      // -1 non singular, 0 source, 1 sink, 2 saddle, 3 center, 4focus
    icVector3 p;                        // position
    icVector3 rgb = icVector3(0.);      // color
    icMatrix2x2 jacobi;
};
extern std::list<Singularity> singularities;

namespace ez{
    typedef icVector3 v3;
    typedef icVector2 v2;
    typedef double db;
    #define ENDL '\n'
}

void findMinMaxField(icVector3& min, icVector3& max){ //XXX: MOST LIKELY RIGHT
    min.x = poly->vlist[0]->x;
    min.y = poly->vlist[0]->y;
    min.z = poly->vlist[0]->z;
    max = min;

    for (int i=1; i < poly->nverts; ++i){
        if (min.x > poly->vlist[i]->x)
            min.x = poly->vlist[i]->x;
        if (min.y > poly->vlist[i]->y)
            min.y = poly->vlist[i]->y;
        if (min.z > poly->vlist[i]->z)
            min.z = poly->vlist[i]->z;

        if (max.x < poly->vlist[i]->x)
            max.x = poly->vlist[i]->x;
        if (max.y < poly->vlist[i]->y)
            max.y = poly->vlist[i]->y;
        if (max.z < poly->vlist[i]->z)
            max.z = poly->vlist[i]->z;
    }
}

//XXX: done
bool onQuad(const Quad* q, const icVector3& p){ 
    double v0x = q->verts[2]->x;
    double v0y = q->verts[2]->y;
    double v2x = q->verts[0]->x;
    double v2y = q->verts[0]->y;
    if (p.x >= v0x && p.x <= v2x && p.y >= v0y && p.y <= v2y)
        return true;
    else
        return false;
}

// XXX: done
Quad* findQuad(const icVector3& v){
    for (int i=0; i < poly->nquads; ++i){
        Quad* qtemp = poly->qlist[i];
        if (onQuad(qtemp, v))
            return qtemp;
    }
    return nullptr;
}


// Get vector from V field by bilinear interpolation
icVector3 getVector(Quad* q, const icVector3& p){ // XXX: possibly right
    using namespace ez;
    db x1 = q->verts[2]->x;
    db y1 = q->verts[2]->y;
    db x2 = q->verts[0]->x;
    db y2 = q->verts[0]->y;
    v3 v11(q->verts[2]->vx, q->verts[2]->vy, q->verts[2]->vz);
    v3 v21(q->verts[3]->vx, q->verts[3]->vy, q->verts[3]->vz);
    v3 v12(q->verts[1]->vx, q->verts[1]->vy, q->verts[1]->vz);
    v3 v22(q->verts[0]->vx, q->verts[0]->vy, q->verts[0]->vz);

    db x = (x2 - x1);
    db y = (y2 - y1);
    db x2x = (x2 - p.x);
    db y2y = (y2 - p.y);
    db xx1 = (p.x - x1);
    db yy1 = (p.y - y1);

    v3 nv = ((x2x/x)*(y2y/y))*v11 + ((xx1/x)*(y2y/y))*v21 + ((x2x/x)*(yy1/y))*v12 + ((xx1/x)*(yy1/y))*v22; // std::cout << "IV: " << interpol_vector.x << " " << interpol_vector.y << " " << interpol_vector.z << std::endl;
    // std::cout << "START vec: " << p.x << " " << p.y << " " << p.z << " next vec: " << nv.x << " " << nv.y << " " << nv.z << std::endl;
    return nv;
} 

bool sinp2Boundary(icVector3& np, const icVector3& min, const icVector3& max){ //XXX: done
    bool hitBound = false;
    if (np.x < min.x){
        hitBound = true;
        np.x = min.x;
    }else if (np.x > max.x){
        hitBound = true;
        np.x = max.x;
    }
    if (np.y < min.y){
        hitBound = true;
        np.y = min.y;
    }else if (np.y > max.y){
        hitBound = true;
        np.y = max.y;
    }
    return hitBound;
}


void streamlineTrace(Quad*& nextQuad, Quad* currQuad, icVector3 currPos, icVector3 currVec, double t,
    const icVector3& min, const icVector3& max){     // XXX: done
    using namespace ez;
    bool insideQuad = false;

    while(!insideQuad){
        //std::cout << "NOT INSIDE QUAD YET\n";
        // outside of field
        if (currPos.x < min.x || currPos.x > max.x ||
            currPos.y < min.y || currPos.y > max.y){
            //std::cout << "EXIT1\n";
            nextQuad = nullptr;
            return;
        }else if (onQuad(currQuad, currPos)){
            currPos = currPos + currVec * t;
            nextQuad = findQuad(currPos);
            return;
        }

        double t_ = INFINITY;
        Quad* nextQuad_ = nullptr;

        for (int ei = 0; ei < 4; ei++){
            Edge* edge = currQuad->edges[ei];
            Vertex* v0 = edge->verts[0];
            Vertex* v1 = edge->verts[1];
            double t_temp;

            if (std::abs(v0->x - v1->x) < EPSILON){
                t_temp = (v0->x - currPos.x) / currVec.x;
            }
            else{
                t_temp = (v0->y - currPos.y) / currVec.y;
            }

            if (t_temp > 0 && t_temp < t_){
                // get next quad
                t_ = t_temp;
                if (edge->quads[0] != currQuad && edge->quads[1] == currQuad)
                    nextQuad_ = edge->quads[0];
                else if (edge->quads[0] == currQuad && edge->quads[1] != currQuad)
                    nextQuad_ = edge->quads[1];
            }
        }

        if (nextQuad_ == nullptr){
            currPos = currPos + currVec * t;
            nextQuad = findQuad(currPos);
            return;
        }
        else{
            if (t_ >= t){
                insideQuad = true;
            }
            else{
                currQuad = nextQuad_;
                t = t - t_;
                currPos = currPos + currVec * t_;
            }
        }
    }
    nextQuad = currQuad;
}

void streamlineFB(POLYLINE& line, const icVector3& seed, const double& step, bool forward = true){ // XXX: done
    using namespace ez;
    line.vertices.push_back(seed);
    Quad* quad = findQuad(seed);
    v3 min, max;
    findMinMaxField(min, max);
    v3 currPos = seed;
    double coef = forward? 1.0 : -1.0;
    
    while(quad != nullptr){
        v3 currVec = getVector(quad, currPos);

        if (currVec.length() < EPSILON){
            break;}

        v3 nextPos = currPos + step * currVec * coef;
        
        // on boundary
        if (sinp2Boundary(nextPos, min, max)){
            line.vertices.push_back(nextPos);
            break;
        }

        // get next quad
        Quad* nextQuad = nullptr;
        streamlineTrace(nextQuad, quad, currPos, currVec * coef, step, min, max);

        // update quad, pos, vector
        quad = nextQuad;
        currPos = nextPos;
        line.vertices.push_back(currPos);
    }
}

// Create streanline forward and backward
void streamline(POLYLINE& line, const icVector3& seed, const double& step){     // XXX: done
    // forward
    streamlineFB(line, seed, step);
    
    // backward
    POLYLINE line_back;
    streamlineFB(line_back, seed, step, false);

    // merge lines
    line.merge(line_back);
}





bool quadricRoot(double& r0, double& r1, const double& a, const double&b, const double&c){
    using namespace ez;
    db m = b*b - 4 * a * c;
    if (m < 0) return false;

    r0 = (-b - std::sqrt(m)) / (2 * a);
    r1 = (-b + std::sqrt(m)) / (2 * a);
    return true;
}

bool singRoot(double &r0, double& r1, const double&a, const double&b, const double&c, const double& d){
    using namespace ez;
    db f0 = b - a - (c + d);
    db f1 = (c + d);
    db f2 = a;
    return quadricRoot(r0, r1, f0, f1, f2);
}

void createCircle(float radius, icVector3 point, const icVector3 color,  std::vector<std::tuple<icVector3, icVector3>>& seeds){
    icVector3 seed = point;
    // up
    seed.y += radius;
    seeds.push_back(std::make_tuple(seed, color));
    // down
    seed = point;
    seed.y -= radius;
    seeds.push_back(std::make_tuple(seed, color));
    // right
    seed = point;
    seed.x += radius;
    seeds.push_back(std::make_tuple(seed, color));
    // left
    seed = point;
    seed.x -= radius;
    seeds.push_back(std::make_tuple(seed, color));
    // up + right
    seed = point;
    seed.y += radius;
    seed.x += radius;
    seeds.push_back(std::make_tuple(seed, color));
    // up + left
    seed = point;
    seed.y += radius;
    seed.x -= radius;
    seeds.push_back(std::make_tuple(seed, color));
    // down + right
    seed = point;
    seed.y -= radius;
    seed.x += radius;
    seeds.push_back(std::make_tuple(seed, color));
    // down + left
    seed = point;
    seed.y -= radius;
    seed.x -= radius;
    seeds.push_back(std::make_tuple(seed, color));
}


void seedingTemplates(const Singularity& s, std::vector<std::tuple<icVector3, icVector3>>& seeds, const float offset){
    icVector3 color;
    if (s.type == 0 || s.type == 1){
        // source/sink: red/blue
        // one concentric circle
        color = icVector3(1, 0, 0);
        if (s.type == 1){
            color = icVector3(0, 0, 1);
        }
        std::cout << "SOURCE/SINK\n";
        createCircle(offset, s.p, color, seeds);

    }else if (s.type == 2){
        // saddle: green
        // one concentric circle with another inside and a seed on the center point
        color = icVector3(0, 1, 0);
        createCircle(offset, s.p, color, seeds);
        createCircle(offset*2, s.p, color, seeds);
        seeds.push_back(std::make_tuple(s.p, color));
        std::cout << "SADDLE\n";
    }else if (s.type == 3){
        // center: cyan
        // line with segments and seeds on the line
        color = icVector3(0, 1, 1);
        // go up
        for (int i=1; i <= 2; ++i){
            icVector3 seed = s.p;
            seed.y += (i*offset);
            seeds.push_back(std::make_tuple(seed, color));
        }

        // go down
        for (int i=1; i <= 2; ++i){
            icVector3 seed = s.p;
            seed.y -= (i*offset);
            seeds.push_back(std::make_tuple(seed, color));
        }
        std::cout << "CENTER\n";
    }else if (s.type == 4){
        // focus: yellow
        // start shaped lines and seeds on the line
        icVector3 seed = s.p;
        color = icVector3(1, 1, 0);
        std::cout << "FOCUS\n";

        seed.y += offset;
        seeds.push_back(std::make_tuple(seed, color));

        seed = s.p;
        seed.y -= offset;
        seeds.push_back(std::make_tuple(seed, color));

        seed = s.p;
        seed.x += offset;
        seeds.push_back(std::make_tuple(seed, color));
        
        seed = s.p;
        seed.x -= offset;
        seeds.push_back(std::make_tuple(seed, color));
    }

    
}

void extractSingularity(){
    using namespace ez;
    singularities.clear();
    // go through faces
    for (int i=0; i < poly->nquads; ++i){
        v3 vx1y1 = poly->qlist[i]->verts[2]->vec();
        v3 px1y1 = poly->qlist[i]->verts[2]->pos();
        v3 vx2y1 = poly->qlist[i]->verts[3]->vec();
        v3 px2y1 = poly->qlist[i]->verts[3]->pos();
        v3 vx2y2 = poly->qlist[i]->verts[0]->vec();
        v3 px2y2 = poly->qlist[i]->verts[0]->pos();
        v3 vx1y2 = poly->qlist[i]->verts[1]->vec();
        v3 px1y2 = poly->qlist[i]->verts[1]->pos();
        
        v3 pt(0.);
        db f11 = vx1y1.x;
        db f12 = vx1y2.x;
        db f21 = vx2y1.x;
        db f22 = vx2y2.x;

        db g11 = vx1y1.y;
        db g12 = vx1y2.y;
        db g21 = vx2y1.y;
        db g22 = vx2y2.y;

        db a00 = f11;
        db a10 = f21 - f11;
        db a01 = f12 - f11;
        db a11 = f11 - f21 - f12 + f22;
        db b00 = g11;
        db b10 = g21 - g11;
        db b01 = g12 - g11;
        db b11 = g11 - g21 - g12 + g22;
        db c00 = a11 * b00 - a00 * b11;
        db c10 = a11 * b10 - a10 * b11;
        db c01 = a11 * b01 - a01 * b11;

        if (c01 == 0.0){
            continue;
        }

        db div = c10/c01;
        db a = -a11 * div;
        db b = a10 - a01 * div;
        db c = a00 - (a01 + a11) * c00 / c01;
        db s[2];
        bool flag = quadricRoot(s[0], s[1], a, b, c);
        if (!flag)
            continue;
        db t[2];
        t[0] = -c00 / c01 - div * s[0];
        t[1] = -c00 / c01 - div * s[1];

        for (int j=0; j < 2; ++j){
            if (t[j] > -EPSILON && t[j] < 1 + EPSILON && s[j] > -EPSILON && s[j] < 1 + EPSILON){
                
                pt.x = s[j];
                pt.y = t[j];
                Singularity point;
                point.p = pt + px1y1;
                singularities.push_back(point);
            }
        }
    }
}

void extractSeparatrix(){
    using namespace ez;
    for (auto & s : singularities){
        if (s.type == 2){
            // saddle
            db a = s.jacobi.entry[0][0];
            db b = s.jacobi.entry[0][1];
            db c = s.jacobi.entry[1][0];
            db d = s.jacobi.entry[1][1];

            db Yd = (a + d)/2;
            db Yr = (c - b)/2;
            db Ys = std::sqrt((a-d) * (a-d) + (b + c) * (b + c)) / 2;
            db theta = std::atan2(b+c, a-d);
            db phi = std::atan(Yr/Ys);
            db a_cos = std::cos(theta/2);
            db b_sin = -std::sin(theta/2);
            db c_sin = std::sin(theta/2);
            db d_cos = std::cos(theta/2); //QUESTION: SAME AS A_COS??????
            db a_sin_phi = std::sqrt(std::sin(phi + M_PI/4));
            db a_cos_phi = std::sqrt(std::cos(phi + M_PI/4));
            v3 maj_v(0.0), min_v(0.0);
            maj_v.x = a_cos * (a_sin_phi + a_cos_phi) + b_sin * (a_sin_phi - a_cos_phi);
            maj_v.y = c_sin * (a_sin_phi + a_cos_phi) + d_cos * (a_sin_phi - a_cos_phi);
            min_v.x = a_cos * (a_sin_phi - a_cos_phi) + b_sin * (a_sin_phi + a_cos_phi);
            min_v.y = c_sin * (a_sin_phi - a_cos_phi) + d_cos * (a_sin_phi + a_cos_phi);

            db k = MIN_K / maj_v.length();

            // outgoing p + kv
            POLYLINE separatrix;
            streamlineFB(separatrix, s.p + k * maj_v, STEP);
            separatrix.rgb = v3(1, 0 ,0);
            polylines.push_back(separatrix);

            // outgoing p - kv
            separatrix.clear();
            streamlineFB(separatrix, s.p - k * maj_v, STEP);
            separatrix.rgb = v3(1, 0, 0);
            polylines.push_back(separatrix);

            // incoming p + kw
            k = MIN_K / min_v.length();
            separatrix.clear();
            streamlineFB(separatrix, s.p + k * min_v, STEP, false);
            separatrix.rgb = v3(0, 0, 1);
            polylines.push_back(separatrix);

            // incoming p - kw
            separatrix.clear();
            streamlineFB(separatrix, s.p - k * min_v, STEP, false);
            separatrix.rgb = v3(0, 0, 1);
            polylines.push_back(separatrix);
        }
    }
}

void classifySingularityByWinding(){
    using namespace ez;
    for (auto&s : singularities){
        v3 posn = s.p;
        Quad* quad = findQuad(posn);

        //winding number
        db winding_angle = 0;
        db angles[4];

        for (int i=0; i < 4; ++i){
            auto vi = quad->verts[i]->vec();
            db angle = atan2(vi.entry[1], vi.entry[0]);
            if (angle < 0)
                angle += 2 * M_PI;
            angles[i] = angle;
        }

        for (int i=0; i < 4; ++i){
            int nexti = i+1;
            if (nexti > 3)
                nexti = 0;
            db diff = angles[nexti] - angles[i];
            if (diff < -M_PI)
                diff += 2 * M_PI;
            if (diff > M_PI)
                diff -= 2 * M_PI;
            winding_angle += diff;
        }

        if (std::abs(winding_angle) < EPSILON){
            // 0 = non singularity
            s.type = -1;
        }
        else if (winding_angle < 0){
            // saddle
            s.type = 2;
            s.rgb = v3(0, 1, 0);
        }else if (winding_angle > 0){
            // soruce, sink, center, focus
            s.type = 0;
            s.rgb = v3(1, 0, 0);
        }
    }


}

void classifySingularity(){
    using namespace ez;
    for (auto& s : singularities){
        v3 pos = s.p;
        Quad* quad = findQuad(pos);

        int R[4] = {2, 3, 0, 1};
        const v2 min(quad->verts[R[0]]->x, quad->verts[R[0]]->y);
        const v2 max(quad->verts[R[2]]->x, quad->verts[R[2]]->y);
        db xlen = max.x - min.x;
        db ylen = max.y - min.y;

        db nix = -1/xlen;
        db ix = 1/xlen;

        db niy = -1/ylen;
        db iy = 1/ylen;

        db qv0x = quad->verts[R[0]]->vx;
        db qv1x = quad->verts[R[1]]->vx;
        db qv2x = quad->verts[R[2]]->vx;
        db qv3x = quad->verts[R[3]]->vx;

        db qv0y = quad->verts[R[0]]->vy;
        db qv1y = quad->verts[R[1]]->vy;
        db qv2y = quad->verts[R[2]]->vy;
        db qv3y = quad->verts[R[3]]->vy;

        db mpx = max.x - pos.x;
        db mpy = max.y - pos.y;
        db pmx = pos.x - min.x;
        db pmy = pos.y - min.y;

        // calculate jacobian
        db dfdx = (nix * mpy * qv0x) +  (ix * mpy * qv1x)    + (nix * pmy * qv3x)    + (ix * pmy * qv2x);
        db dfdy = (mpx * niy * qv0x) +  (pmx * niy * qv1x)   + (mpx * iy * qv3x)     + (pmx * iy * qv2x);
        db dgdx = (nix * mpy * qv0y) +  (ix * mpy * qv1y)    + (nix * pmy * qv3y)    + (ix * pmy * qv2y);
        db dgdy = (mpx * niy * qv0y) +  (pmx * niy * qv1y)   + (mpx * iy * qv3y)     + (pmx * iy * qv2y);

        s.jacobi.entry[0][0] = dfdx;
        s.jacobi.entry[0][1] = dfdy;
        s.jacobi.entry[1][0] = dgdx;
        s.jacobi.entry[1][1] = dgdy;

        db tr = dfdx + dgdy;
        db det = dfdx * dgdy - dfdy * dgdx;

        db delta = tr * tr - 4 * det;
        if (delta >= 0){
            db r1 = 0.5 * (tr + sqrt(delta));
            db r2 = 0.5 * (tr - sqrt(delta));
            if (r1 == 0 && r2 == 0){
                // unknown
                s.type = -1;
            }else if (r1 > 0 && r2 > 0){
                // source
                s.type = 0;
                s.rgb = v3(1, 0, 0);
            }else if (r1 < 0 && r2 < 0){
                // sink
                s.type = 1;
                s.rgb = v3(0, 0, 1);
            }else if (r1 > 0 && r2 < 0 || r1 < 0 && r2 > 0){
                // saddle
                s.type = 2;
                s.rgb = v3(0, 1, 0);
            }else{
                s.type = -1;
            }
        }else{
            if (tr == 0){
                // center
                s.type = 3;
                s.rgb = v3(0, 1, 1);
            }else{
                // focus
                s.type = 4;
                s.rgb = v3(1, 1, 0);
            }
        }
    }
}

void drawSingularities(){
    for (auto& s : singularities){
        drawDot(s.p.x, s.p.y, s.p.z, 0.15, s.rgb.x, s.rgb.y, s.rgb.z);
    }
}

