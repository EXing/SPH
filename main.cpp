#include <iostream>
#include <GL/glut.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <memory.h>
#include "Vector2.h"

using namespace std;

#define ParticleCount 1000
#define Pi 3.1415926535f

#define ParticleRadius 0.05
#define Gamma 7
#define Eta 0.01
// alpha is between 0.08~0.5
#define Alpha 0.08
#define Cs 1500
#define Density0 1000
#define Epsilon 0.01
#define Gravity 9.81
#define TimeStep 4.52E-4
#define H (6 * ParticleRadius)
#define Viscocity 0.5
#define Mass (Density0 * 4 / 3 * Pi * pow(ParticleRadius, 3))
#define MaxNeighborCount 64

#define ScreenWidth 500
#define ScreenHeight 500
#define ViewWidth 10.0f
#define ViewHeight (ScreenHeight * ViewWidth / ScreenWidth)
#define CellSize H

#define SubSteps 100
#define ShowStep SubSteps * TimeStep

class Particle
{
public:
    Vector2 pos;
    Vector2 v;
    Vector2 a;
    double density;
    double nearDensity;
    double pressure;
    double nearPressure;
    Particle* next;
};

// c is a parameter for collision
struct Wall
{
    Wall(double x, double y, double c) : x(x), y(y), c(c) {}
    double x, y, c;
};

struct Neighbors
{
    const Particle* particles[MaxNeighborCount];
    double dis[MaxNeighborCount];
    int count;
};

Particle particles[ParticleCount];
Neighbors neighbors[ParticleCount];
Vector2 prePosition[ParticleCount];
Vector2 relaxedPosition[ParticleCount];

Wall walls[4] = {
        Wall(1, 0, 0),
        Wall(0, 1, 0),
        Wall(-1, 0, -ViewWidth),
        Wall(0, -1, -ViewHeight)
};

const int GridWidth = (int)(ViewWidth / CellSize);
const int GridHeight = (int)(ViewHeight / CellSize);
const int GridCellCount = GridHeight * GridWidth;
Particle* grid[GridCellCount];
int gridCoords[ParticleCount * 2];

void updateGrid()
{
    // Clear grid
    memset(grid, 0, GridCellCount * sizeof(Particle*));

    // Add particles to grid
    for(int i = 0; i < ParticleCount; i++)
    {
        Particle &p = particles[i];
        int x = (int) (p.pos.x / CellSize);
        int y = (int) (p.pos.y / CellSize);

        // TODO : maybe wrong
        if(x < 1)
            x = 1;
        else if (x > GridWidth - 2 )
            x = GridWidth - 2;

        if(y < 1)
            y = 1;
        else if (y > GridHeight - 2)
            y = GridHeight - 2;

        // maybe there are many particles in one grid, so I use a linked list to store them
        p.next = grid[x + y * GridWidth];
        grid[x + y * GridWidth] = &p;

        // the coordinates of the ith particle
        gridCoords[2 * i] = x;
        gridCoords[2 * i + 1] = y;
    }
}

void gravityForce()
{
    for(int i = 0; i < ParticleCount; i++)
        particles[i].v.y -= Gravity * TimeStep;
}

void advance()
{
    for(int i = 0; i < ParticleCount; i++)
    {
        // preserve current position
        prePosition[i] = particles[i].pos;

        particles[i].pos = particles[i].pos + (particles[i].v * TimeStep);
    }
}

double W(double x)
{
    double h = 1.0;
    double sigma3 = 2 / (3 * h);
    if(2 > x)
        return 0;
    if(1 < x && 2 >= x)
        return sigma3 / 4 * pow((2 - x), 3);
    if(0 <= x && 1 >= x)
        return sigma3 * (1 - 1.5 * x * x * (1 - x / 2));
}

void calculatePressure()
{
    double B = Density0 * Cs * Cs / Gamma;
    for(int i = 0; i < ParticleCount; i++)
    {
        Particle &pi = particles[i];
        int gi = gridCoords[2 * i];
        int gj = gridCoords[2 * i + 1] * GridWidth;

        neighbors[i].count = 0;

        double density = 0;
        for(int ni = gi - 1; ni <= gi + 1; ni++)
        {
            for(int nj = gj - GridWidth; nj <= gj + GridWidth; nj+= GridWidth)
            {
                for(Particle* pj = grid[ni + nj]; NULL != pj; pj = pj->next)
                {
                    double r = (pj->pos - pi.pos).Sqrt();
                    if(r < sqrt(Epsilon) || r > H)
                        continue;

                    density += Mass * W(r);
                    if(neighbors[i].count < MaxNeighborCount)
                    {
                        neighbors[i].particles[neighbors[i].count] = pj;
                        neighbors[i].dis[neighbors[i].count] = r;
                        neighbors[i].count ++;
                    }
                }
            }
        }
        pi.density = density;
        pi.pressure = B * ( pow(pi.density / Density0 , Gamma) - 1);
    }
}

void momentumEquation()
{
    for(int i = 0; i < ParticleCount; i++)
    {
        const Particle &p = particles[i];
        Vector2 deltaV(0,0);

        for(int j;j<neighbors[i].count;j++)
        {
            const Particle& pj = *neighbors[i].particles[j];
            double dis = neighbors[i].dis[j];
        }
    }
}

void Render()
{
    glClearColor(0.02f, 0.01f, 0.01f, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, ViewWidth, 0, ViewHeight, 0, 1);

    glEnable(GL_POINT_SMOOTH);

    glEnableClientState(GL_VERTEX_ARRAY);
    glPointSize(2.5f*ParticleRadius*ScreenWidth/ViewWidth);

    glVertexPointer(2, GL_DOUBLE, sizeof(Particle), particles);
    glDrawArrays(GL_POINTS, 0, ParticleCount);

    glDisableClientState(GL_VERTEX_ARRAY);

    glutSwapBuffers();
}

void Update()
{
   /* for (size_t step=0; step<SubSteps; ++step)
    {
        EmitParticles();

        gravityForce();
        advance();
        updateGrid();
        calculatePressure();
        calculateRelaxedPositions();
        moveToRelaxedPositions();
        updateGrid();
        resolveCollisions();
    }*/

    glutPostRedisplay();
}

int main(int argc, char** argv) {
    glutInitWindowSize(ScreenWidth, ScreenHeight);
    glutInit(&argc, argv);
    glutInitDisplayString("samples stencil>=3 rgb double depth");
    glutCreateWindow("SPH");
    glutDisplayFunc(Render);
    glutIdleFunc(Update);

    //INIT
    memset(particles, 0, ParticleCount*sizeof(Particle));
    UpdateGrid();

    glutMainLoop();

    return 0;
}

