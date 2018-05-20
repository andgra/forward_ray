//
// Created by ANDGRA on 20.05.2018.
//

#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "intersection.h"
#include "beam.h"
#include "arrayOfReceiversTransmitters.h"
#include "data.h"
#include <sstream>  //for std::istringstream
#include <iterator> //for std::istream_iterator
#include <queue>  // подключили библиотеку queue
#include <assert.h>


using std::string;
using std::vector;
using std::count;
using std::stringstream;
using std::getline;
using std::ifstream;
using std::istream_iterator;
using std::to_string;
using std::cout;

string fileNameWoExt( string const& pathname )
{
    string fileName = pathname.substr(pathname.find_last_of("/\\") + 1);
    string::size_type const p(fileName.find_last_of('.'));
    return fileName.substr(0, p);
}

vector<string> split(const string &s, char delim)
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

        int cntLines = vcount(lines);


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
            int cntCoord = vcount(coordinates);
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

    static void SendBeam(beam curBeam, arrayOfReceiversTransmitters receivers,
                         queue <beam> calculationQueue)//луч - начальная точка, угол, амплитуда
    {
        double distance;
        double time;
        double abs;

        if (fabs(curBeam.value) < EPS)//луч слишком угас, дальше не рассматриваем
            return;

        intersection closestIntersection = receivers.FindClosestIntersection(curBeam, false);//находим точку ближайшего пересечения луча с границей

        if (!closestIntersection.existed) {
            return;//ничего не нашлось
        }

        if (closestIntersection.figureIndex == 0)//это не фигура, а один из концов пространства
        {
            distance = closestIntersection.distance;
            time = (double) (curBeam.time + ((distance * dX / figureCollection[curBeam.figureIndex].speed) / dt));
            if (time > maxTime)
                return;//луч превысил отведенное время, дальше не рассматриваем

            abs = GetAbsorption(figureCollection[0].absorption_c, distance);
            if (closestIntersection.intersectionPoint.Y <
                1)//попали в сторону с приёмниками. Отражение без преломления и фиксация при попадании в приемник
            {
                //считаем, что среда до приёмника - фоновая(логично)
                receivers.CheckAndRecord(closestIntersection.intersectionPoint.X, (int) round(time),
                                         curBeam.value * abs);
                if (curBeam.direction.Y != 0)
                    calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint,
                                                      vec2d::reflect(curBeam.direction, vec2d(0, 1)), time,
                                                      -curBeam.value * abs, 0));
                else
                    calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint,
                                                      vec2d::reflect(curBeam.direction, vec2d(0, 1)), time,
                                                      curBeam.value * abs,
                                                      0));//для горизонтальной волны не должен меняться знак - там нет отражений
            } else//попали в противоположную сторону. Отражение без преломления.
            {
                calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint,
                                                  vec2d::reflect(curBeam.direction, vec2d(0, -1)), time,
                                                  -curBeam.value * abs, 0));
            }
            return;
        }


        vec2d dir = curBeam.direction;

        distance = closestIntersection.distance;
        time = curBeam.time + (double) ((distance * dX / figureCollection[curBeam.figureIndex].speed) / dt);
        if (time > maxTime)
            return;//луч превысил отведенное время, дальше не рассматриваем
        abs = GetAbsorption(figureCollection[curBeam.figureIndex].absorption_c, distance);

        //дальше проверяем, не является ли эта находка касанием
        if (GetNextFigure(closestIntersection.intersectionPoint, dir) == curBeam.figureIndex ||
            figureCollection[closestIntersection.figureIndex].isItDifObject()) {
            if (figureCollection[closestIntersection.figureIndex].IsPointDif(closestIntersection.intersectionPoint))
                receivers.RecordDif(closestIntersection.intersectionPoint, curBeam.value * abs, (int) time,
                                    closestIntersection.figureIndex);
            //////////////просто обрабатываем срабатывание касания и посылаем продолжения этого луча с данной точки
            calculationQueue.push(
                    beam(curBeam, closestIntersection.intersectionPoint, curBeam.direction, time, curBeam.value * abs,
                             curBeam.figureIndex));
            return;
        }


        double rcoef;
        double pcoef;

        int targetFigureIndex;//номер той фигуры, в которую попадает преломленный луч
        int sourceFigureIndex;

        targetFigureIndex = GetNextFigure(closestIntersection.intersectionPoint,
                                          curBeam.direction);//спрашиваем, куда луч проходит после прохождение границы, направление при преломлении изменится, но останется в пределах данного прогноза

        int v1 = figureCollection[curBeam.figureIndex].speed;
        int v2 = figureCollection[targetFigureIndex].speed;


        if (figureCollection[closestIntersection.figureIndex].IsPointDif(
                closestIntersection.intersectionPoint))//isPointDif(closestIntersection.intersectionPoint, closestIntersection.figureIndex))
        {
            receivers.RecordDif(closestIntersection.intersectionPoint, curBeam.value * abs, (int) time,
                                closestIntersection.figureIndex);
        }

        bool difP = false;
        vec2d reflectionDir = vec2d::reflect(closestIntersection.direction, closestIntersection.normalVector);
        sourceFigureIndex = GetNextFigure(closestIntersection.intersectionPoint, reflectionDir);
        if (curBeam.figureIndex != sourceFigureIndex) {
            //Debug.WriteLine("Аномалия отражения в точке " + Math.Round(closestIntersection.intersectionPoint.X) + ", " + Math.Round(closestIntersection.intersectionPoint.Y));
            difP = true;//если возникла такая ситуация, когда при грубом подсчёте луч отражается не в ту среду, из которой пришёл, значит мы на внутреннем угле
        }
        sourceFigureIndex = curBeam.figureIndex;

        if (v1 * sin(closestIntersection.angle * degreeToRadians) >=
            v2)//при таком раскладе должно быть полное внутреннее отражение
        {
            calculationQueue.push(
                    beam(curBeam, closestIntersection.intersectionPoint, reflectionDir, time, curBeam.value * abs,
                             sourceFigureIndex));
            //Console.WriteLine("полное внутреннее отражение");
        } else if ((v1 <= v2) && (closestIntersection.angle >= asin((double) v1 / v2) /
                                                               degreeToRadians))//проходящая волна скользит по поверхности вместо прохождения насквозь - мы не считаем эту волну
        {
            double cosReflectionAngle = cos(closestIntersection.angle * degreeToRadians);
            rcoef = (figureCollection[targetFigureIndex].gamma / cosReflectionAngle -
                     figureCollection[sourceFigureIndex].gamma / cosReflectionAngle) /
                    (figureCollection[targetFigureIndex].gamma / cosReflectionAngle +
                     figureCollection[sourceFigureIndex].gamma / cosReflectionAngle);
            calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint, reflectionDir, time,
                                              (double) (curBeam.value * rcoef) * abs, sourceFigureIndex));
            //Console.WriteLine("закритический угол " + closestIntersection.angle + " при критическом " + Math.Asin((double)v1 / v2) / degreeToRadians);
        } else {
            double temp = vec2d::dot(curBeam.direction * v1, closestIntersection.normalVector);

            double k = (double) (sqrt((v2 * v2 - v1 * v1) / (temp * temp) + 1) - 1);

            vec2d refractionDir = curBeam.direction * v1 + closestIntersection.normalVector * temp * k;

            targetFigureIndex = GetNextFigure(closestIntersection.intersectionPoint,
                                              refractionDir);//поздновато, но если мы выяснили, что дальше всё не так, как мы считали, это самое ранее, где мы можем поправить итог

            double cosReflectionAngle = cos(closestIntersection.angle * degreeToRadians);
            rcoef = (figureCollection[targetFigureIndex].gamma / cosReflectionAngle -
                     figureCollection[sourceFigureIndex].gamma / cosReflectionAngle) /
                    (figureCollection[targetFigureIndex].gamma / cosReflectionAngle +
                     figureCollection[sourceFigureIndex].gamma / cosReflectionAngle);

            double cosRefractionAngle =
                    vec2d::dot(curBeam.direction, refractionDir) / (curBeam.direction.length() * refractionDir.length());
            pcoef = (2 * figureCollection[sourceFigureIndex].gamma / cosReflectionAngle) /
                    (figureCollection[targetFigureIndex].gamma / cosRefractionAngle +
                     figureCollection[sourceFigureIndex].gamma / cosReflectionAngle);


            if (isPointDif(closestIntersection.intersectionPoint, closestIntersection.figureIndex) || difP) {
                receivers.RecordDif(closestIntersection.intersectionPoint, curBeam.value * abs, (int) time,
                                    closestIntersection.figureIndex);
            }

            calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint, reflectionDir, time,
                                              (double) (curBeam.value * rcoef) * abs, sourceFigureIndex));
            calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint, refractionDir, time,
                                              (double) (curBeam.value * pcoef) * abs, targetFigureIndex));
        }
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
        int cntCoord = vcount(coordinatesOfTransmitters);
        for (int i = 0; i < width / step; i++) {
            coordinatesOfTransmitters[i] = i * step;
        }

        vector<vec2d> directions = GenerateArrayOfVectors(30, 0.01);//0.0001);

//        ParallelOptions options = new ParallelOptions();
//        options.MaxDegreeOfParallelism = 7;

#pragma omp parallel for
        for (int i = 0; i < cntCoord; i++) {


//        Parallel.For(0, coordinatesOfTransmitters.Length, options, i = >
//        {
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
                        SendBeam(current, receivers, calculationQueue);
                    }
                    catch (const std::exception& ex) {
                        cout << "!" << ex.what();
                    }
                }


            }
            receivers.ProcessDifs();//Отрисовываем зафиксированные диф.

//            SaveData("L" + fileName + "(withoutAbs.data)" + "-" + to_string(i + 1) + ".data" + to_string(i),
//                     receivers.convolvedData, vcount(coordinatesOfTransmitters));//записываем результаты данной фиксации

        }


//        SaveInfo("L" + fileName + "(withoutAbs.data)", maxTime, vcount(coordinatesOfTransmitters), step);
        //SaveInfoInUniversalFormat("L3-0.5_2D-z80_150.data", maxTime, coordinatesOfTransmitters.Length, step);
        cout << "Finished";
    }

private:

    /*static void SaveInfo(string path, int maxtime, int cnt, int step)//, int impulseLen)
    {
        using (BinaryWriter
        writer = new BinaryWriter(File.Open(path + ".info", FileMode.Create)))
        {
            writer.Write(cnt);//количество трансмиттеров
            writer.Write(step);//шаг, с которым они расположены, начиная с 0
            writer.Write(maxtime*//*impulseLen *//*);
            writer.Write((double) (1000 * dX));//мм в 1 пикселе по dx
            writer.Write((double) (1000 * dY));//мм в 1 пикселе по dy
            writer.Write((int) (1 / dt));//частота дискретизации
            //writer.Write(difRadius);
        }
    }

    static void SaveInfoInUniversalFormat(string path, int maxtime, int cnt, int step) {
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
    }

    static void SaveData(string path, vector<vector<double>> data,
    int cnt
    )//, int impulseLen)
    {
        using (BinaryWriter
        writer = new BinaryWriter(File.Create(path)))
        {
            for (int x = 0; x < cnt; x++) {
                for (int t = 0; t < maxTime; t++) {
                    writer.Write(data[x, t]);
                }
            }
        }

        Console.WriteLine(Path.GetFileName(path) + " сохранён");
    }*/
};

#endif //ALGO_H
