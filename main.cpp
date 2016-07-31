#include <iostream>
#include <GL/glut.h>
#include <cmath>
#include <algorithm>
#include <cstring>

using namespace std;

#define ParticleCount 1000
#define Pi 3.1415926535f
#define ParticleRadius 0.05
#define Gamma 7
// #define Eta 0.01
// alpha is between 0.08~0.5
#define Alpha 0.5
#define Cs 88.5
#define Density0 1000.0
#define Epsilon 0.01
#define Gravity 9.81
#define TimeStep 0.001
#define H (6 * ParticleRadius)
#define Mass (Density0 * 2 * ParticleRadius * 2 * ParticleRadius)
#define MaxNeighborCount 64
#define ScreenWidth 500
#define ScreenHeight 500
#define ViewWidth 10.0f
#define ViewHeight (ScreenHeight * ViewWidth / ScreenWidth)
#define CellSize H
#define SubSteps 100
#define h ParticleRadius
#define sigma3 ( 2 / (3 * h))
#define MatrixRow 25
#define MatrixCol (ParticleCount/MatrixRow)

class Vector2 {
public:
    Vector2() {}

    Vector2(double x, double y) : x(x), y(y) {}

    double x;
    double y;

    Vector2 operator+(Vector2 x) {
        return Vector2(x.x + this->x, x.y + this->y);
    }

    Vector2 operator-(Vector2 x) {
        return Vector2(this->x - x.x, this->y - x.y);
    }

    double operator*(Vector2 x) {
        return x.x * this->x + x.y * this->y;
    }

    Vector2 operator*(double x) {
        return Vector2(x * this->x, x * this->y);
    }

    double Sqrt() {
        return sqrt(this->x * this->x + this->y * this->y);
    }
};

class Particle {
public:
    Vector2 pos;
    Vector2 v;
    double density;
    double pressure;
    Particle *next;
};

//neighbors[i] used to record the ith particle's neighbors in its centered grid.
struct Neighbors {
     Particle *particles[MaxNeighborCount];
    double dis[MaxNeighborCount];
    int count;
};

Particle particles[ParticleCount];
Neighbors neighbors[ParticleCount];
Vector2 prePosition[ParticleCount];

const int GridWidth = (int) (ViewWidth / CellSize);
const int GridHeight = (int) (ViewHeight / CellSize);
const int GridCellCount = GridHeight * GridWidth;

//grid[GridCellCount] point to the header of list that particles in the same grid
Particle *grid[GridCellCount];

//index of grid[]
int gridCoords[ParticleCount * 2];

void updateGrid() {
    // Clear grid
    memset(grid, 0, GridCellCount * sizeof(Particle *));

    // Add particles to grid
    for (int i = 0; i < ParticleCount; i++) {
        Particle &p = particles[i];
        int x = (int) (p.pos.x / CellSize);
        int y = (int) (p.pos.y / CellSize);

        // TODO : maybe wrong
        if (x < 1)
            x = 1;
        else if (x > GridWidth - 2)
            x = GridWidth - 2;

        if (y < 1)
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

void gravityForces() {
    for (int i = 0; i < ParticleCount; i++)
        particles[i].v.y -= Gravity * TimeStep;
}

void advance() {
    for (int i = 0; i < ParticleCount; i++) {
        // preserve current position
        prePosition[i] = particles[i].pos;

        particles[i].pos = particles[i].pos + (particles[i].v * TimeStep);
    }
}

double W(double x) {
    if (2 < x)
        return 0;
    else if (1 < x)
        return sigma3 / 4 * pow((2 - x), 3);
    else
        return sigma3 * (1 - 1.5 * x * x * (1 - x / 2.0));
}

Vector2 deltaW(Vector2 a) {
    double x = a.Sqrt();
    if (2 < x)
        return Vector2(0, 0);
    else if (1 < x) {
        double k = -3.0 * sigma3 / 4.0 * pow((2 - x), 2);
        return a * -k * TimeStep;
    }
    else {
        double k = sigma3 * (9.0 / 4.0 * x * x - 3.0 * x);
        return a * -k* TimeStep;
    }
}

void calculatePressure() {
    double B = Density0 * Cs * Cs / Gamma /100;

    for (int i = 0; i < ParticleCount; i++) {
        Particle &pi = particles[i];
        int gi = gridCoords[2 * i];
        int gj = gridCoords[2 * i + 1] * GridWidth;

        neighbors[i].count = 0;

        double density = 0;
        for (int ni = gi - 1; ni <= gi + 1; ni++) {
            for (int nj = gj - GridWidth; nj <= gj + GridWidth; nj += GridWidth) {
                for (Particle *pj = grid[ni + nj]; pj != NULL; pj = pj->next) {
                    double r = (pj->pos - pi.pos).Sqrt();
                    if (/*r < sqrt(Epsilon) || */r > H)
                        continue;

                    density += Mass * W(r);
                    if (neighbors[i].count < MaxNeighborCount) {
                        neighbors[i].particles[neighbors[i].count] = pj;
                        neighbors[i].dis[neighbors[i].count] = r;
                        neighbors[i].count++;
                    }
                }
            }
        }
        pi.density = density;
        pi.pressure = B * (pow(pi.density / Density0, Gamma) - 1);
        //printf("%f %f\n",pi.density,pi.pressure);
    }
}

inline void momentumEquation() {
    for (int i = 0; i < ParticleCount; i++) {
        Particle &p = particles[i];
        Vector2 delta(0, 0);

        for (int j = 0; j < neighbors[i].count; j++) {
            Particle &pj = *neighbors[i].particles[j];
            double dis = neighbors[i].dis[j];

            delta = delta -
                    deltaW(pj.pos - p.pos) * Mass * (p.pressure / p.density / p.density + pj.pressure / pj.density / pj.density);
        }
        p.v = p.v + delta * TimeStep;
    }
}

inline void viscosityEquation() {
    for (int i = 0; i < ParticleCount; i++) {
        Particle &p = particles[i];
        Vector2 delta(0, 0);

        for (int j = 0; j < neighbors[i].count; j++) {
             Particle &pj = *neighbors[i].particles[j];
            double dis = neighbors[i].dis[j];

            double dianji = (p.v - pj.v) * (p.pos - pj.pos);
            if (dianji < 0) {
                delta = delta - deltaW(pj.pos - p.pos) * Mass * ((-2 * Alpha * ParticleRadius * Cs / (p.density + pj.density)) *
                                                      (dianji /
                                                       (dis * dis + Epsilon * ParticleRadius * ParticleRadius)));
            }
        }
        p.v = p.v + delta * TimeStep;
    }
}

void collisions() {
    for (int i = 0; i < ParticleCount; i++) {
        Particle &p = particles[i];
        if (p.pos.x >= 10 || p.pos.x <= 0)
            p.v.x = -0.2 * p.v.x;
        if (p.pos.y >= 10 || p.pos.y <= 0)
            p.v.y = -0.2 * p.v.y;
    }
}

/*double random(double x, double y) {
    return x + (y - x) * ((double) rand() / (double) (RAND_MAX - 1));
}*/

void generateParticles() {
    for (int i = 0; i < MatrixRow; i++) {
        for (int j = 0; j < MatrixCol; j++) {
            particles[i*MatrixCol+j].pos = Vector2(h+j*2*h, h+i*2*h);
            particles[i*MatrixCol+j].v = Vector2(0, 0);
        }
    }
}

void Render() {
    glClearColor(0.02f, 0.01f, 0.01f, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, ViewWidth, 0, ViewHeight, 0, 1);

    glEnable(GL_POINT_SMOOTH);

    glEnableClientState(GL_VERTEX_ARRAY);
    glPointSize((GLfloat) (2.5f * ParticleRadius * ScreenWidth / ViewWidth));

    glVertexPointer(2, GL_DOUBLE, sizeof(Particle), particles);
    glDrawArrays(GL_POINTS, 0, ParticleCount);

    glDisableClientState(GL_VERTEX_ARRAY);

    glutSwapBuffers();
}

void update() {
    for (int step = 0; step < SubSteps; step++) {
        calculatePressure();
        gravityForces();
        momentumEquation();
        viscosityEquation();
        advance();
        collisions();
        updateGrid();
    }
    glutPostRedisplay();
}

int main(int argc, char **argv) {
    glutInitWindowSize(ScreenWidth, ScreenHeight);
    glutInit(&argc, argv);
    glutCreateWindow("SPH");
    glutDisplayFunc(Render);
    glutIdleFunc(update);

    //INIT
    memset(particles, 0, ParticleCount * sizeof(Particle));
    generateParticles();
    updateGrid();

    glutMainLoop();

    return 0;
}