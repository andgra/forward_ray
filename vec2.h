#ifndef __VEC2_H__
#define __VEC2_H__

#include <math.h>

template <class T>
class vec2 {
public:
    T X, Y;

    vec2() :X(0), Y(0) {}
    vec2(T X, T Y) : X(X), Y(Y) {}
    vec2(const vec2& v) : X(v.X), Y(v.Y) {}

    vec2& operator=(const vec2& v) {
        X = v.X;
        Y = v.Y;
        return *this;
    }

    vec2 operator+(vec2& v) {
        return vec2(X + v.X, Y + v.Y);
    }

    vec2 operator+(vec2 v) {
        return vec2(X + v.X, Y + v.Y);
    }
    vec2 operator-(vec2& v) {
        return vec2(X - v.X, Y - v.Y);
    }
    vec2 operator-(vec2 v) {
        return vec2(X - v.X, Y - v.Y);
    }

    vec2& operator+=(vec2& v) {
        X += v.X;
        Y += v.Y;
        return *this;
    }
    vec2& operator-=(vec2& v) {
        X -= v.X;
        Y -= v.Y;
        return *this;
    }

    vec2 operator+(double s) {
        return vec2(X + s, Y + s);
    }
    vec2 operator-(double s) {
        return vec2(X - s, Y - s);
    }
    vec2 operator*(double s) {
        return vec2(X * s, Y * s);
    }
    vec2 operator/(double s) {
        return vec2(X / s, Y / s);
    }


    vec2& operator+=(double s) {
        X += s;
        Y += s;
        return *this;
    }
    vec2& operator-=(double s) {
        X -= s;
        Y -= s;
        return *this;
    }
    vec2& operator*=(double s) {
        X *= s;
        Y *= s;
        return *this;
    }
    vec2& operator/=(double s) {
        X /= s;
        Y /= s;
        return *this;
    }

    void set(T X, T Y) {
        this->X = X;
        this->Y = Y;
    }

    void rotate(double deg) {
        double theta = deg / 180.0 * M_PI;
        double c = cos(theta);
        double s = sin(theta);
        double tX = X * c - Y * s;
        double tY = X * s + Y * c;
        X = tX;
        Y = tY;
    }

    vec2& normalize() {
        if (length() == 0) return *this;
        *this *= (1.0 / length());
        return *this;
    }

    vec2& reflect(vec2 n) {
        vec2 r = n * dot(n, *this) * 2.0;
        *this -= r;
        return *this;
    }

    double dist(vec2 v) const {
        vec2 d(v.X - X, v.Y - Y);
        return d.length();
    }
    double length() const {
        return std::sqrt(X * X + Y * Y);
    }
    void truncate(double length) {
        double angle = atan2f(Y, X);
        X = length * cos(angle);
        Y = length * sin(angle);
    }

    vec2 ortho() const {
        return vec2(Y, -X);
    }

    static double dot(vec2 v1, vec2 v2) {
        return v1.X * v2.X + v1.Y * v2.Y;
    }
    static vec2 normalize(vec2 v) {
        return v.normalize();
    }
    static vec2 reflect(vec2 v, vec2 n) {
        return v.reflect(n);
    }
    static double cross(vec2 v1, vec2 v2) {
        return (v1.X * v2.Y) - (v1.Y * v2.X);
    }

};

typedef vec2<float> vec2f;
typedef vec2<double> vec2d;
typedef vec2<int> vec2i;


#endif