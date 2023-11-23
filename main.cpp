#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <strstream>
#include <algorithm>
#include <list>
#include "SDL2/SDL.h"
using namespace std;

struct vec3d {
    float x=0; float y=0; float z=0;
    float w=1;
};

struct triangle {
    vec3d p[3];
    int color[3];
};

struct mat4x4 {
    float m[4][4] = {0};
};

struct mesh {
    vector<triangle> tris;
};



bool loadMeshFromObjFile(string filename, mesh &mesh) {
    ifstream f(filename);
    if (!f.is_open()) return false;

    vector<vec3d> verts;
    while (!f.eof()) {
        char line[128];
        f.getline(line, 128);
        strstream s; s << line;

        char ident;
        if (line[0] == 'v') {
            vec3d v;
            s >> ident >> v.x >> v.y >> v.z;
            verts.push_back(v);
        }

        if (line[0] == 'f') {
            int f[3];
            s >> ident >> f[0] >> f[1] >> f[2];
            mesh.tris.push_back({verts[f[0]-1], verts[f[1]-1], verts[f[2]-1] });
        }
    }
    return true;
}
vec3d vectorAdd(vec3d &v1, vec3d &v2) {
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

vec3d vectorSub(vec3d &v1, vec3d &v2) {
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

vec3d vectorMul(vec3d &v, float s) {
    return { v.x * s, v.y * s, v.z * s};
}


vec3d vectorDiv(vec3d &v, float s) {
    return { v.x / s, v.y / s, v.z / s};
}

float vectorDotProd(vec3d &v1, vec3d &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
vec3d vectorCrossProd(vec3d &v1, vec3d &v2) {
    vec3d v;
    v.x = v1.y * v2.z - v1.z * v2.y;
    v.y = v1.z * v2.x - v1.x * v2.z;
    v.z = v1.x * v2.y - v1.y * v2.x;
    return v;
}

float vectorLength(vec3d& v) {
    return sqrtf(vectorDotProd(v,v));
}

vec3d vectorNormalise(vec3d &v) {
    float l = vectorLength(v);
    return {v.x/l, v.y/l, v.z/l};
}

vec3d vectorIntersect(vec3d &planeP, vec3d &planeN, vec3d &lineStart, vec3d &lineEnd) {
    planeN = vectorNormalise(planeN);
    float planeD = -vectorDotProd(planeN, planeP);
    float ad = vectorDotProd(lineStart, planeN);
    float bd = vectorDotProd(lineEnd, planeN);
    float t = (-planeD - ad) / (bd - ad);
    vec3d startToEnd = vectorSub(lineEnd, lineStart);
    vec3d lineIntersect = vectorMul(startToEnd, t);
    return vectorAdd(lineStart, lineIntersect);
}

vec3d matrixMultiplyVector(mat4x4 &m, vec3d &i) {
    vec3d v;
    v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
    return v;
}

mat4x4 matrixInitIdentity(void) {
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f; matrix.m[1][1] = 1.0f; matrix.m[2][2] = 1.0f; matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixRotX(float angleRad)
{
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = cosf(angleRad);
    matrix.m[1][2] = sinf(angleRad);
    matrix.m[2][1] = -sinf(angleRad);
    matrix.m[2][2] = cosf(angleRad);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixRotY(float angleRad)
{
    mat4x4 matrix;
    matrix.m[0][0] = cosf(angleRad);
    matrix.m[0][2] = sinf(angleRad);
    matrix.m[2][0] = -sinf(angleRad);
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = cosf(angleRad);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixRotZ(float angleRad)
{
    mat4x4 matrix;
    matrix.m[0][0] = cosf(angleRad);
    matrix.m[0][1] = sinf(angleRad);
    matrix.m[1][0] = -sinf(angleRad);
    matrix.m[1][1] = cosf(angleRad);
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixTranslate(float x, float y, float z) {
    mat4x4 matrix;
    matrix.m[0][0] = 1.0f; matrix.m[1][1] = 1.0f; matrix.m[2][2] = 1.0f; matrix.m[3][3] = 1.0f;
    matrix.m[3][0] = x;
    matrix.m[3][1] = y;
    matrix.m[3][2] = z;
    return matrix; 
}

mat4x4 matrixProject(float fovDeg, float aspectRatio, float near, float far) {
    float fovRad = 1.0f / tanf(fovDeg * 0.5f / 180.0f * 3.14159f);
    mat4x4 matrix;
    matrix.m[0][0] = aspectRatio * fovRad;
    matrix.m[1][1] = fovRad;
    matrix.m[2][2] = far / (far - near);
    matrix.m[3][2] = (-far * near) / (far - near);
    matrix.m[2][3] = 1.0f;
    matrix.m[3][3] = 0.0f;
    return matrix;
}

mat4x4 matrixMultiplyMatrix(mat4x4 &m1, mat4x4 &m2) {
    mat4x4 matrix;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
    return matrix;
}

mat4x4 matrixPoint(vec3d &pos, vec3d &target, vec3d &up) {
    vec3d newForward = vectorSub(target, pos);
    newForward = vectorNormalise(newForward);

    vec3d a = vectorMul(newForward, vectorDotProd(up, newForward));
    vec3d newUp = vectorSub(up, a);
    newUp = vectorNormalise(newUp);

    vec3d newRight = vectorCrossProd(newUp, newForward);

    mat4x4 matrix;
    matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 matrixQuickInverse(mat4x4 &m) {
    mat4x4 matrix;
    matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
    matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
    matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

int triangleClip(vec3d planeP, vec3d planeN, triangle &in, triangle &out1, triangle &out2) {
    planeN = vectorNormalise(planeN);

    auto dist = [&](vec3d &p) {
        vec3d n = vectorNormalise(p);
        return (planeN.x*p.x  + planeN.y*p.y + planeN.z*p.z - vectorDotProd(planeN, planeP));
    };

    vec3d* insidePoints[3]; int insidePointCount = 0;
    vec3d* outsidePoints[3]; int outsidePointCount = 0;

    float d0 = dist(in.p[0]); float d1 = dist(in.p[1]);float d2 = dist(in.p[2]);

    if (d0 >=0) { insidePoints[insidePointCount++] = &in.p[0];}
    else { outsidePoints[outsidePointCount++] = &in.p[0]; }
    if (d1 >=0) { insidePoints[insidePointCount++] = &in.p[1];}
    else { outsidePoints[outsidePointCount++] = &in.p[1]; }
    if (d2 >=0) { insidePoints[insidePointCount++] = &in.p[2];}
    else { outsidePoints[outsidePointCount++] = &in.p[2]; }

    switch (insidePointCount) {
        case 0:{
            return 0;
        }
        case 3: {
            out1 = in;
            return 1;
        } 
        case 1: {
            if (outsidePointCount == 2) {
                out1.p[0] = *insidePoints[0];
                out1.p[1] = vectorIntersect(planeP, planeN, *insidePoints[0], *outsidePoints[0]);
                out1.p[2] = vectorIntersect(planeP, planeN, *insidePoints[0], *outsidePoints[1]);

                return 1;
            }
        }
        case 2: {
            if (outsidePointCount == 1) {
                out1.p[0] = *insidePoints[0];
                out1.p[1] = *insidePoints[1];
                out1.p[2] = vectorIntersect(planeP, planeN, *insidePoints[0], *outsidePoints[0]);

                out2.p[0] = *insidePoints[1];
                out2.p[1] = out1.p[2];
                out2.p[2] = vectorIntersect(planeP, planeN, *insidePoints[1], *outsidePoints[0]);

                return 2;
            }
        }
    }
}


bool spin = false;
bool update(SDL_Window* window, SDL_Renderer *renderer, mesh &m, double delta, vec3d camera, float yaw, vec3d lookDir, float theta, mat4x4 matProj) {
    SDL_RenderClear(renderer);
    SDL_SetRenderDrawColor(renderer, 0,0,0,255);

    mat4x4 matTrans = matrixTranslate(0.0f, 0.0f, 5.0f);
    mat4x4 matWorld = matrixInitIdentity();
    if (spin) {
        mat4x4 matRotZ, matRotX;
        matRotZ = matrixRotZ(theta * 0.5f);
        matRotX = matrixRotX(theta);
        matWorld = matrixMultiplyMatrix(matRotZ, matRotX);
    }
    matWorld = matrixMultiplyMatrix(matWorld, matTrans);

    vec3d up = {0,1,0};
    vec3d target  = {0,0,1};
    mat4x4 cameraRot = matrixRotY(yaw);
    lookDir = matrixMultiplyVector(cameraRot, target);
    target = vectorAdd(camera, lookDir);
    mat4x4 matCamera = matrixPoint(camera, target, up);
    mat4x4 matView = matrixQuickInverse(matCamera);

    std::vector<triangle> trianglesToRaster;

    for (auto tri: m.tris) {
        triangle triProjected, triTransformed, triViewed;

        triTransformed.p[0] = matrixMultiplyVector(matWorld, tri.p[0]);
        triTransformed.p[1] = matrixMultiplyVector(matWorld, tri.p[1]);
        triTransformed.p[2] = matrixMultiplyVector(matWorld, tri.p[2]);

        vec3d normal, line1, line2;

        line1 = vectorSub(triTransformed.p[1], triTransformed.p[0]);
        line2 = vectorSub(triTransformed.p[2], triTransformed.p[0]);

        normal = vectorCrossProd(line1, line2);
        normal = vectorNormalise(normal);

        vec3d cameraRay = vectorSub(triTransformed.p[0], camera);

        if (vectorDotProd(normal, cameraRay) < 0.0f) {
            vec3d lightDirection = { 0.0f, 1.0f, -1.0f };
            lightDirection = vectorNormalise(lightDirection);
            float dp = max(0.1f, vectorDotProd(lightDirection, normal));
            triViewed.color[0] = triViewed.color[1] = triViewed.color[2] = 255*dp;
            triViewed.p[0] = matrixMultiplyVector(matView, triTransformed.p[0]);
            triViewed.p[1] = matrixMultiplyVector(matView, triTransformed.p[1]);
            triViewed.p[2] = matrixMultiplyVector(matView, triTransformed.p[2]);
        
            int clippedTrangles = 0;
            triangle clipped[2];
            clippedTrangles = triangleClip({0.0f, 0.0f, 0.1f}, {0.0f, 0.0f, 1.0f}, triViewed, clipped[0], clipped[1]);

            for (int n = 0; n < clippedTrangles; n++) {
                triProjected.color[0] = triViewed.color[0];
                triProjected.color[1] = triViewed.color[1];
                triProjected.color[2] = triViewed.color[2];
                triProjected.p[0] = matrixMultiplyVector(matProj, clipped[n].p[0]);
                triProjected.p[1] = matrixMultiplyVector(matProj, clipped[n].p[1]);
                triProjected.p[2] = matrixMultiplyVector(matProj, clipped[n].p[2]);

                triProjected.p[0] = vectorDiv(triProjected.p[0], triProjected.p[0].w);
                triProjected.p[1] = vectorDiv(triProjected.p[1], triProjected.p[1].w);
                triProjected.p[2] = vectorDiv(triProjected.p[2], triProjected.p[2].w);

                triProjected.p[0].x *= -1.0f;
                triProjected.p[1].x *= -1.0f;
                triProjected.p[2].x *= -1.0f;
                triProjected.p[0].y *= -1.0f;
                triProjected.p[1].y *= -1.0f;
                triProjected.p[2].y *= -1.0f;

                vec3d offsetView = {1,1,0};

                triProjected.p[0] = vectorAdd(triProjected.p[0], offsetView);
                triProjected.p[1] = vectorAdd(triProjected.p[1], offsetView);
                triProjected.p[2] = vectorAdd(triProjected.p[2], offsetView);
                triProjected.p[0].x *=0.5f * (float)640;
                triProjected.p[0].y *=0.5f * (float)480;
                triProjected.p[1].x *=0.5f * (float)640;
                triProjected.p[1].y *=0.5f * (float)480;
                triProjected.p[2].x *=0.5f * (float)640;
                triProjected.p[2].y *=0.5f * (float)480;

                trianglesToRaster.push_back(triProjected);

            }
        }   
    }

    sort(trianglesToRaster.begin(), trianglesToRaster.end(), [](triangle &t1, triangle &t2) {
        float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
        float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
        return z1 > z2;
    });

    for (auto &triToRaster: trianglesToRaster) {
        triangle clipped[2];
        list<triangle> triList;
        triList.push_back(triToRaster);
        int newTris = 1;
        for (int p = 0; p < 4; p++)  {
            int trisToAdd = 0;
            while (newTris > 0) {
                triangle test = triList.front();
                triList.pop_front();
                newTris--;

                switch(p) {
                    case 0: trisToAdd = triangleClip({ 0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, test, clipped[0], clipped[1]); break;
                    case 1: trisToAdd = triangleClip({ 0.0f, (float)480 - 1, 0.0f}, { 0.0f, -1.0f, 0.0f}, test, clipped[0], clipped[1]); break;
                    case 2: trisToAdd = triangleClip({ 0.0f, 0.0f, 0.0f}, { 1.0f, 0.0f, 0.0f}, test, clipped[0], clipped[1]); break;
                    case 3: trisToAdd = triangleClip({ (float)640 - 1, 0.0f, 0.0f}, {-1.0f, 0.0f, 0.0f}, test, clipped[0], clipped[1]); break;
                }

                for (int w = 0; w < trisToAdd; w++) {
                    triList.push_back(clipped[w]);
                }
            }
            newTris = triList.size();
        }
        for (auto &t : triList) {
            SDL_Vertex vert1 = {{t.p[0].x, t.p[0].y}, {t.color[0],t.color[0],t.color[0], 255}, {1,1}};
            SDL_Vertex vert2 = {{t.p[1].x, t.p[1].y}, {t.color[0],t.color[0],t.color[0], 255}, {1,1}};
            SDL_Vertex vert3 = {{t.p[2].x, t.p[2].y}, {t.color[0],t.color[0],t.color[0], 255}, {1,1}};
            SDL_Vertex verts[] = { vert1, vert2, vert3};
            SDL_RenderGeometry(renderer, nullptr, verts, 3, NULL, 0);
        }
    }
    SDL_RenderPresent(renderer);
    return true;
}

int main(void) {
    if(SDL_Init(SDL_INIT_VIDEO) < 0)
        std::cout << "SDL could not be initialized: " << SDL_GetError();
    else
        std::cout << "SDL video system is ready to go\n";
    
    SDL_Window* window = SDL_CreateWindow("Atel",SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 640, 480, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    vec3d camera = {0,0,0};
    vec3d lookDir;
    float yaw;
    float theta;
    mesh mesh1;
    mat4x4 matProj = matrixProject(90.0f, ((float)640/(float)480),0.1f, 1000.0f);;
    loadMeshFromObjFile("hand.obj", mesh1);

    bool isRunning = true;
    Uint64 now = SDL_GetPerformanceCounter();
    Uint64 prev = 0;
    double deltaTime = 0;

    while(isRunning) {
        theta += 5.0f * (deltaTime * .0001f);
        SDL_Event event;
        while(SDL_PollEvent(&event)) {
            switch(event.type) {
                case SDL_QUIT:
                isRunning = false;
                break;   
                case SDL_KEYDOWN:
                    switch(event.key.keysym.sym) {
                        case SDLK_RIGHT:
                            camera.x +=1.0f * (deltaTime * .001f);
                            break;
                        case SDLK_LEFT:
                            camera.x -=1.0f * (deltaTime * .001f);
                            break;
                        case SDLK_UP:
                            camera.y +=1.0f * (deltaTime * .001f);
                            break;
                        case SDLK_DOWN:
                            camera.y -=1.0f * (deltaTime * .001f);
                            break;
                        case SDLK_w:
                            camera.z +=1.0f * (deltaTime * .001f);
                            break;
                        case SDLK_s:
                            std::cout << spin << std::endl;
                            spin = true;
                            break;
                        case SDLK_o:
                            spin = false;
                            break;
                    }
            }
        }
        prev = now;
        now = SDL_GetPerformanceCounter();
        deltaTime = ((double)((now - prev)*1000 / (double)SDL_GetPerformanceFrequency()));
        update(window, renderer, mesh1, deltaTime, camera, yaw, lookDir, theta, matProj);
    } 

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
