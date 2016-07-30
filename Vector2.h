//
// Created by kevin on 7/30/16.
//

#ifndef SPH_VECTOR2_H
#define SPH_VECTOR2_H

class Vector2
{
    Vector2() {}

    Vector2(double x, double y) : x(x) ,y(y) {}

public:
    double x;
    double y;

    Vector2 operator + (Vector2 x) {
        return Vector2(x.x + this->x, x.y + this->y);
    }

    Vector2 operator - (Vector2 x) {
        return Vector2(this->x - x.x, this->y - x.y);
    }

    double operator * (Vector2 x) {
        return x.x * this->x + x.y * this->y;
    }

    Vector2 operator * (double x) {
        return Vector2(x * this->x, x * this->y);
    }

    double Sqrt() {
        return sqrt(this->x * this->x + this->y * this->y);
    }
};

#endif //SPH_VECTOR2_H
