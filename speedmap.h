//
// Created by ANDGRA on 06.05.2019.
//

#ifndef STRAIGHT_RAY_SPEEDMAP_H
#define STRAIGHT_RAY_SPEEDMAP_H

#include "lib/bitmap_image.hpp"
#include <vector>
#include "figure.h"


class speedmap {
public:
    speedmap(const vector<figure>& _fc) {
        fc = _fc;
    }

    speedmap draw() {
        image = bitmap_image(static_cast<const unsigned int>(fc[0].right - fc[0].left),
                static_cast<const unsigned int>(fc[0].up - fc[0].bottom));

        // set background to env
        rgb_t envBg = speedToRGB(fc[0].speed);
        image.set_all_channels(envBg.red, envBg.green, envBg.blue);

        for (unsigned int i = 0; i < image.width(); i++) {
            for (unsigned int j = 0; j < image.height(); j++) {
                for (figure f: fc) {
                    if (f.index == 0) continue; // фоновую пропускаем
                    if (f.ContainsPoint(pointF(i, j))) {
                        image.set_pixel(i, j, speedToRGB(f.speed));
                    }
                }
            }
        }

        drawn = true;

        return *this;
    }

    bool save(const string &path) {
        if (drawn) {
            image.save_image(path);
            return true;
        }
        return false;
    }

private:
    bool drawn;
    vector<figure> fc;
    bitmap_image image;
    int maxSpeed = 16777215;
    int minSpeed = 0;
    rgb_t speedToRGB(int speed) {
        rgb_t result;
        if (speed > maxSpeed) {
            throw string("speed is much higher");
        }
        if (speed < minSpeed) {
            throw string("speed is much lower");
        }
        result.blue = static_cast<unsigned char>(speed % 256);
        speed -= result.blue;
        speed /= 256;
        result.green = static_cast<unsigned char>(speed % 256);
        speed -= result.green;
        speed /= 256;
        result.red = static_cast<unsigned char>(speed);
        return result;
    }
};
#endif //STRAIGHT_RAY_SPEEDMAP_H
