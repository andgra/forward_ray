//
// Created by ANDGRA on 20.05.2018.
//

#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <fstream>
#include "arrayOfReceiversTransmitters.h"
#include <sstream>  //for std::istringstream
#include <iterator> //for std::istream_iterator
#include <clocale>
#include <omp.h>
#include "beam.h"
#include "data.h"
#include <queue>
#include <string>
#include <vector>
#include <math.h>



using std::string;
using std::vector;
using std::stringstream;
using std::getline;
using std::ifstream;
using std::istream_iterator;
using std::to_string;
using std::cout;

string basePath( string const pathname )
{
    return pathname.substr(pathname.find_last_of("/\\") + 1);
}

string fileNameWoExt( string const pathname )
{
    string fileName = basePath(pathname);
    string::size_type const p(fileName.find_last_of('.'));
    return fileName.substr(0, p);
}

vector<string> split(const string s, char delim)
{
    vector<string> elems;
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

class algo {
public:
    static void OpenFile(string path) {
        fileName = fileNameWoExt(path);
        vector<string> lines;

        ifstream file(path);
        copy(istream_iterator<string>(file),
             istream_iterator<string>(),
             back_inserter(lines));

        int cntLines = lines.size();


        auto size = split(lines[0], ';');//нулевая строка - размеры
        width = stoi(size[0]);
        height = stoi(size[1]);


        vector<string> bckgrnd = split(lines[1], ';');//первая строка - параметры фона(среды, не принадлежащей ни к одному из заданных объектов)

        background = figure(vector<edge>{edge(pointI(0, 0), pointI(width - 1, 0), 0, 0),
                                            edge(pointI(0, height - 1), pointI(width - 1, height - 1), 1, 0)},
                            stoi(bckgrnd[0]), stod(bckgrnd[1]), stod(bckgrnd[2]), 0);
        //background.edges.Add();//грань с приёмниками
        //background.edges.Add();//противоположная грань
        figureCollection = vector<figure>();
        figureCollection.push_back(background);
        for (int j = 2; j < cntLines; j++)//начиная с второй строки идут описания фигур:
        {
            string figureStr = lines[j];
            vector<string> figureLine = split(figureStr, ';');
            vector<string> coordinates = split(figureLine.back(), ':');
            int cntCoord = coordinates.size();
            auto figureEdges = vector<edge>();
            pointI lastPoint, tempPoint, startPoint;

            auto xy = split(coordinates[0], '.');
            int x = stoi(xy[0]);
            int y = stoi(xy[1]);
            startPoint = pointI(x, y);
            lastPoint = pointI(x, y);

            for (int k = 1; k < cntCoord; k++) {
                xy = split(coordinates[k], '.');
                x = stoi(xy[0]);
                y = stoi(xy[1]);
                tempPoint = pointI(x, y);
                figureEdges.emplace_back(edge(lastPoint, tempPoint, k - 1, j - 1));
                lastPoint = tempPoint;
            }
            figureEdges.emplace_back(edge(lastPoint, startPoint, cntCoord - 1, j - 1));

            //var id = figureLine[0];//по сути, нам не нужно
            int speed = stoi(figureLine[1]);
            double density = stod(figureLine[2]);
            double absorption_c = stod(figureLine[3]);

            figureCollection.emplace_back(figure(figureEdges, speed, density, absorption_c, j - 1));
        }
        //maxTime = (int)(dY * height / (figureCollection.Min(x => x.speed)) / dt * 2);
        maxTime = (int) ((sqrt(pow(dY * height, 2) + pow(dX * width / 2, 2)) /
                          minSpeedCollection() / dt * 2));
    }


    static vector<vec2d> GenerateArrayOfVectors(int treshold, double step, bool onlyUp = false) {
        vector<vec2d> temp;
        double normal = 90;
        int size;

        if (onlyUp)
            size = 0;
        else
            size = (int) ((normal - treshold) / step);

        vector<vec2d> directions = vector<vec2d>(size + 3);

        directions[0] = vec2d(0, 1);//всегда добавляется вертикальный вектор
        directions[1] = vec2d(1, 0);//и прямая волна
        directions[2] = vec2d(-1, 0);

        for (int i = 3; i < size + 3; i += 2) {
            directions[i] = vec2d((double) cos((normal - (i - 2) * step) * degreeToRadians),
                                        (double) sin((normal - (i - 2) * step) * degreeToRadians));
            directions[i + 1] = vec2d((double) cos((normal + (i - 2) * step) * degreeToRadians),
                                            (double) sin((normal + (i - 2) * step) * degreeToRadians));
        }
        return directions;
    }

    static void MainAlgorithm() {
        int step = 2;//расстояние между датчиками
        vector<int> coordinatesOfTransmitters = vector<int>(width / step);
        int cntCoord = coordinatesOfTransmitters.size();
        for (int i = 0; i < cntCoord; i++) {
            coordinatesOfTransmitters[i] = i * step;
        }

        vector<vec2d> directions = GenerateArrayOfVectors(30, 0.01);//0.0001);

#pragma omp parallel for
        for (int i = 0; i < cntCoord; i++) {
            arrayOfReceiversTransmitters receivers = arrayOfReceiversTransmitters(coordinatesOfTransmitters, directions,
                                                                                  1, maxTime);//запускаем новую фиксацию
            queue<beam> calculationQueue = queue<beam>();
            beam current;

            for (vec2d d: directions)//добавляем все направления расчёта луча из данной точки в очередь
            {
                calculationQueue.push(beam(pointD(coordinatesOfTransmitters[i], 0), d, 0, 1, 0));

                while (calculationQueue.size() != 0) {
                    current = calculationQueue.front();
                    calculationQueue.pop();
                    current.direction = vec2d::normalize(current.direction);

                    try {
                        receivers.SendBeam(current, calculationQueue);
                    }
                    catch (const std::exception& ex) {
                        cout << "!" << ex.what();
                    }
                }


            }
            receivers.ProcessDifs();//Отрисовываем зафиксированные диф.

            SaveData("L" + fileName + "(withoutAbs.data)"/* + "-" + to_string(i + 1)*/ + ".data" + to_string(i),
                     receivers.convolvedData(), coordinatesOfTransmitters.size());//записываем результаты данной фиксации

        }


        SaveInfo("L" + fileName + "(withoutAbs.data)", maxTime, coordinatesOfTransmitters.size(), step);
        //SaveInfoInUniversalFormat("L3-0.5_2D-z80_150.data", maxTime, coordinatesOfTransmitters.Length, step);
        cout << "Finished";
    }

private:

    static void SaveInfo(string path, int maxtime, int cnt, int step)//, int impulseLen)
    {
        auto dxmm = 1000 * dX;
        auto dymm = 1000 * dY;
        auto df = (int) (1 / dt);
        ofstream myFile (path + ".info", ios::out | ios::binary);
        myFile.write(reinterpret_cast<char*>(&cnt),sizeof(cnt));
        myFile.write(reinterpret_cast<char*>(&step),sizeof(step));
        myFile.write(reinterpret_cast<char*>(&maxtime),sizeof(maxtime));
        myFile.write(reinterpret_cast<char*>(&dxmm),sizeof(dxmm));
        myFile.write(reinterpret_cast<char*>(&dymm),sizeof(dymm));
        myFile.write(reinterpret_cast<char*>(&df),sizeof(df));
    }

    /*static void SaveInfoInUniversalFormat(string path, int maxtime, int cnt, int step) {
        using (BinaryWriter
        writer = new BinaryWriter(File.Open(path + ".info", FileMode.Create)))
        {
            writer.Write(width);
            writer.Write(height);
            writer.Write(maxtime);
            writer.Write(cnt);//количество трансмиттеров
            writer.Write((double) (dX));//мм в 1 пикселе по dx
            writer.Write((double) (dY));//мм в 1 пикселе по dy
            writer.Write((int) (1 / dt));//частота дискретизации
        }
    }*/

    static void SaveData(string path, vector<vector<double>> data, int cnt )//, int impulseLen)
    {
        ofstream myFile (path, ios::out | ios::binary);
        for(auto v: data) {
            for(auto d: v) {
                myFile.write(reinterpret_cast<char*>(&d),sizeof(d));
            }
        }

        cout << basePath(path) << " saved;" << endl;
    }
};

#endif //ALGO_H
