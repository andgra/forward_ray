//
// Created by ANDGRA on 20.05.2018.
//

#ifndef ARRAY_OF_RECEIVERS_TRANSMITTERS_H
#define ARRAY_OF_RECEIVERS_TRANSMITTERS_H

#include "point.h"
#include "vec2.h"
#include "beam.h"
#include "edge.h"
#include "figure.h"
#include "intersection.h"
#include "data.h"
#include <complex>
#include "dft.h"
#include <unordered_map>
#include <queue>
#include <math.h>

typedef complex<double> comp;

using namespace std;

typedef unordered_map<int, double> diffMap;
typedef unordered_map<pointD*, diffMap> diffsMap;

class arrayOfReceiversTransmitters {
public:
    vector<int> coordinates;//список координат приёмников
    int radius;//радиус луча, т.е. насколько близкое к приемнику попадание засчитывается
    vector<vector<double>> recordedData;

    vector<vector<double>> convolvedData() {
        return MakeConvolution();
    }

    int maxTime;
    vector<vec2d> directions;//список направлений лучей, посылаемых из каждого приемника


    arrayOfReceiversTransmitters(vector<int> coordinates_i, vector<vec2d> directons_i, int radius_i, int maxTime_i) {
        coordinates = coordinates_i;
        radius = radius_i;
        maxTime = maxTime_i;
        recordedData = vector<vector<double>>(coordinates.size());
        for(int i = 0; i<coordinates.size(); i++) {
            recordedData[i] = vector<double>(maxTime + 1);
        }
        directions = directons_i;
        step = coordinates[1] - coordinates[0];
        if (radius > step)
            throw string("Overlapping areas of reception");
        //diffactions = new List<DifEffectInstance>();
        difs = diffsMap();
    }

    void CheckAndRecord(double x, int time, double value)//, bool dif)
    {
        int pos = Check(x);
        if (pos != -1)
            if (fabs(float(recordedData[pos][time])) <
                fabs(float(value)))//мы не должны суммировать сигналы лучей, пришедшие в один приемник в одно и то же время
                recordedData[pos][time] = value;//просто берём самый сильный сигнал
        //Debug.WriteLine(i + " " + time);

        return;//в два приемника одновременно попасть не можем
    }

    int Next(double x)//следующий от этой точки вправо приемник
    {
        int sourceIndex = Check(x);
        if (sourceIndex + 1 < coordinates.size())
            return coordinates[sourceIndex + 1];
        else
            return -1;
    }

    int Previous(double x)//следующий от этой точки влево приемник
    {
        int sourceIndex = Check(x);
        if (sourceIndex - 1 >= 0)
            return coordinates[sourceIndex - 1];
        else
            return -1;
    }

    void RecordDif(pointD p, double val, int prevTime, int fIndex)//, int fIndex)
    {
        auto point = figureCollection[fIndex].GetNearestDifPoint(p);
        bool existed = false;
        bool existedDif = false;
        for(auto it1 : difs) {
            if(*(it1.first) == point) {
                existed = true;
                for(auto it2 : it1.second) {
                    if(it2.first == prevTime) {
                        if (it2.second < val)
                            difs[it1.first][it2.first] = val;
                        existedDif = true;

                        break;
                    }
                }
                if(!existedDif) {
                    difs[it1.first].insert(diffMap::value_type(prevTime, val));
                }
                break;
            }
        }
        if (!existed) {
            diffMap tempDif;
            tempDif.insert(diffMap::value_type(prevTime, val));
            difs.insert(diffsMap::value_type(&point, tempDif));
        }
        //diffactions.Add(new DifEffectInstance(figureCollection[fIndex].GetNearestDifPoint(p), prevTime, val));//записываем точную точку диф, а не то, где сработал луч
    }

    void ProcessDifs() {
        //нам нужно попробовать добиться того, чтобы на одну точку в один момент времени отрисовывалась только одна диф

        //var t = diffactions.GroupBy(e => e.difPoint);
        //foreach(var t0 in t)
        //{
        //    var t1 = t0.GroupBy(e => (int)e.time);
        //    foreach(var t2 in t1)
        //    {
        //        var val = t2.Max(e => e.val);
        //        DrawDif(t0.Key, val, t2.Key);
        //    }
        //}

        for(auto it1 : difs) {
            for(auto it2 : it1.second) {
                DrawDif(*(it1.first), it2.second, it2.first);
            }
        }
        /*for (auto iter = difs.begin(); iter != difs.end(); ++iter) {
            pointD p = iter->first;
            auto difsP = iter->second;
            for (auto iter2 = difsP.begin(); iter2 != difsP.end(); ++iter2) {
                int time = iter2->first;
                DrawDif(p, iter2->second, time);
            }
        }*/
    }

    intersection FindClosestIntersection(beam beam, bool excludingTouches) {
        intersection closestIntersection;

        if (beam.direction.Y == 0)//прямая волна
        {
            int nextReceiverCoord;
            if (beam.direction.X > 0) {
                nextReceiverCoord = this->Next(beam.startPoint.X);
            } else {
                nextReceiverCoord = this->Previous(beam.startPoint.X);
            }
            if (nextReceiverCoord != -1)//если дошли до края, то дальше луч ни с чем не пересекается
                closestIntersection = intersection(pointD(nextReceiverCoord, 0), beam.direction,
                                                   beam.figureIndex, 0, vec2d(0, 1),
                                                   GetDistance(beam.startPoint,
                                                               pointD(nextReceiverCoord, 0)));//проверить!
        } else {
            double closestIntersectionDistance = numeric_limits<double>::max();

            double A = beam.startPoint.Y - (beam.startPoint.Y + beam.direction.Y);
            double B = (beam.direction.X + beam.startPoint.X) - beam.startPoint.X;
            double C = beam.startPoint.X * (beam.startPoint.Y + beam.direction.Y) -
                       (beam.direction.X + beam.startPoint.X) * beam.startPoint.Y;

            double currentDistance;

            for(figure f: figureCollection)//опрашиваем все сегменты фигур на пересечение с данным лучом
            {
                for(edge l: f.edges)
                {
                    double zn = det(l.A, l.B, A, B);
                    if (fabs(zn) < EPS) continue;//прямые не пересекаются

                    double x = (-det(l.C, l.B, C, B) / zn);
                    double y = (-det(l.A, l.C, A, C) / zn);//нашли точку пересечения

                    if (IsPointOnEdge(beam.startPoint, l))//если мы рассматриваем ту самую границу, на которой находимся, она нам не интересна
                        continue;

                    if (x >= l.xStart && x <= l.xEnd && y >= l.yStart && y <= l.yEnd)//точка лежит на отрезке
                    {
                        //проверяем, с нужной ли стороны луча получилось пересечение
                        if ((beam.direction.X * (x - beam.startPoint.X) >= 0) &&
                            (beam.direction.Y * (y - beam.startPoint.Y) >= 0)) {
                            currentDistance = GetDistance(pointD(x, y), beam.startPoint);
                            if (currentDistance < closestIntersectionDistance) {
                                vec2d dir = beam.direction;
                                int fIndex = InWhichSmallestFigureIsPoint(pointD(x + dir.X, y + dir.Y));
                                if (excludingTouches && (((fIndex == beam.figureIndex) && (f.index != 0)) ||
                                                         f.isItDifObject())) { ;// Console.WriteLine(x + ", " + y + " - касание из " + beam.startPoint.X + ", " + beam.startPoint.Y);
                                } else {
                                    closestIntersection = intersection(pointD(x, y), beam.direction,
                                                                       l.figureIndex, l.edgeIndex, l.normalVector(),
                                                                       currentDistance);
                                    closestIntersectionDistance = currentDistance;
                                }
                            }
                        }
                    }
                }
            }
        }
        return closestIntersection;
    }

    double GetDistance(pointD a, pointD b) {
        return sqrt(pow(a.X - b.X, 2) * pow(dX, 2) + pow(a.Y - b.Y, 2) * pow(dY, 2)) / dX;
    }

    double det(double a, double b, double c, double d) {
        return a * d - b * c;
    }

    bool IsPointOnEdge(pointD p, edge e) {
        double k, c;

        if (e.p2.X == e.p1.X) {
            return (p.X == e.p1.X && p.Y >= min(e.p1.Y, e.p2.Y) && p.Y <= max(e.p1.Y, e.p2.Y));
        }

        k = (double(e.p2.Y) - e.p1.Y) / (e.p2.X - e.p1.X);

        if (k == 0) {
            return (p.Y == e.p1.Y && p.X >= min(e.p1.X, e.p2.X) && p.X <= max(e.p1.X, e.p2.X));
        }

        c = e.p1.Y - k * e.p1.X;

        return p.Y == p.X * k + c;
    }

    void SendBeam(beam curBeam, queue <beam>& calculationQueue)//луч - начальная точка, угол, амплитуда
    {
        double distance;
        double time;
        double abs;

        if (fabs(curBeam.value) < EPS)//луч слишком угас, дальше не рассматриваем
            return;

        intersection closestIntersection = FindClosestIntersection(curBeam, false);//находим точку ближайшего пересечения луча с границей

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
                CheckAndRecord(closestIntersection.intersectionPoint.X, (int) round(time),
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
                RecordDif(closestIntersection.intersectionPoint, curBeam.value * abs, (int) time,
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
            RecordDif(closestIntersection.intersectionPoint, curBeam.value * abs, (int) time,
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
                RecordDif(closestIntersection.intersectionPoint, curBeam.value * abs, (int) time,
                                    closestIntersection.figureIndex);
            }

            calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint, reflectionDir, time,
                                       (double) (curBeam.value * rcoef) * abs, sourceFigureIndex));
            calculationQueue.push(beam(curBeam, closestIntersection.intersectionPoint, refractionDir, time,
                                       (double) (curBeam.value * pcoef) * abs, targetFigureIndex));
        }
    }

private:
//private double[] signal = { 0, 0.688605F, 0.0746371F, -0.320453F, -0.0703439F, 0.145268F, 0.0493075F, -0.0639651F, -0.0304584F, 0.0272173F, 0.0174822F, 0, 0, 0 };
//private double[] signal = { 0.0f, 0.049901336421413582f, 0.012533323356430454f, -0.12769734259472953f, -0.0323296853314312f, 0.12363734711836992f, 0.058899928429548734f, -0.14477232839456303f, -0.086715661338308936f, 0.16042230584538295f, 0.13519060802726859f, -0.19262831069394754f, -0.20536413177860582f, 0.61609239533582116f, 0.7319875806369972f, -0.58778525229247813f, -0.84432792550201252f, 0.48175367410171815f, 0.90482705246601725f, -0.33131209741621637f, -0.808398038850879f, 0.14921393229891752f, 0.49114362536434342f, -0.025066646712862559f, -0.29940801852848131f, 0.0000000000000025482654181230304f, 0.29940801852848165f, 0.037599970069288363f, -0.19645745014573834f, -0.024868988716484131f, 0.19021130325903149f, 0.036812455268466777f };
//private double[] signal = { 0f, 0.01154333f, 0.02304818f, 0.03447618f, 0.04578925f, 0.05694965f, 0.06792018f, 0.07866427f, 0.08914609f, 0.09933069f, 0.1091841f, 0.1186735f, 0.1277673f, 0.1364351f, 0.144648f, 0.1523787f, 0.1596013f, 0.1662918f, 0.1724279f, 0.1779892f, 0.182957f, 0.1873148f, 0.1910482f, 0.1941445f, 0.1965936f, 0.1983873f, 0.1995195f, 0.1999866f, 0.1997869f, 0.1989211f, 0.1973921f, 0.195205f, 0.1923671f, 0.1888878f, 0.1847788f, 0.1800537f, 0.1747284f, 0.1688205f, 0.1623497f, 0.1553377f, 0.1478078f, 0.1397851f, 0.1312963f, 0.1223698f, 0.1130353f, 0.103324f, 0.09326819f, 0.08290142f, 0.07225826f, 0.06137419f, 0.0502855f, 0.03902915f, 0.02764269f, 0.01616407f, 0.00463155f, -0.007105037f, -0.0194472f, -0.03235147f, -0.04577045f, -0.05965291f, -0.07394405f, -0.08858567f, -0.1035164f, -0.118672f, -0.1339855f, -0.1493877f, -0.1648072f, -0.1801708f, -0.1954039f, -0.2104307f, -0.2251749f, -0.2395592f, -0.2535067f, -0.2669406f, -0.2797845f, -0.2919632f, -0.3034028f, -0.314031f, -0.3237776f, -0.3325747f, -0.3403572f, -0.3470632f, -0.3526339f, -0.3570144f, -0.3601539f, -0.3620057f, -0.3625278f, -0.361683f, -0.3594393f, -0.3557697f, -0.3506531f, -0.344074f, -0.3360224f, -0.3264948f, -0.3154936f, -0.3030272f, -0.2891104f, -0.2737646f, -0.257017f, -0.2389015f, -0.2194581f, -0.1987329f, -0.1767784f, -0.1536529f, -0.1294204f, -0.1041509f, -0.07791957f, -0.05080708f, -0.02289898f, 0.005714433f, 0.03520138f, 0.06563723f, 0.09691449f, 0.1289193f, 0.1615318f, 0.1946269f, 0.2280741f, 0.2617388f, 0.2954821f, 0.3291618f, 0.3626329f, 0.395748f, 0.428358f, 0.4603131f, 0.4914628f, 0.5216569f, 0.550746f, 0.5785827f, 0.6050213f, 0.6299194f, 0.6531378f, 0.6745418f, 0.6940014f, 0.711392f, 0.7265953f, 0.7394994f, 0.75f, 0.7580002f, 0.7634116f, 0.7661549f, 0.7661598f, 0.7633657f, 0.7577223f, 0.7491899f, 0.7377395f, 0.7233532f, 0.7060245f, 0.6857586f, 0.6625723f, 0.6364943f, 0.6075652f, 0.5758376f, 0.541376f, 0.5042565f, 0.4645673f, 0.4224078f, 0.3778889f, 0.3311324f, 0.282271f, 0.2314475f, 0.1788149f, 0.1245355f, 0.06878076f, 0.01173044f, -0.04575089f, -0.1024507f, -0.1581827f, -0.2127659f, -0.2660252f, -0.3177919f, -0.3679044f, -0.4162085f, -0.4625582f, -0.5068156f, -0.5488517f, -0.5885467f, -0.6257905f, -0.6604825f, -0.6925323f, -0.7218598f, -0.7483954f, -0.7720798f, -0.7928648f, -0.8107126f, -0.8255965f, -0.8375003f, -0.8464185f, -0.8523564f, -0.8553296f, -0.855364f, -0.8524956f, -0.8467703f, -0.8382437f, -0.8269802f, -0.8130537f, -0.7965463f, -0.7775484f, -0.7561581f, -0.7324808f, -0.7066289f, -0.678721f, -0.6488816f, -0.6172405f, -0.5839324f, -0.5490962f, -0.5128743f, -0.4754126f, -0.4368594f, -0.3973649f, -0.3570809f, -0.3161598f, -0.2747546f, -0.2330177f, -0.191101f, -0.1491545f, -0.1073269f, -0.06576395f, -0.02460869f, 0.01599937f, 0.05551623f, 0.09363686f, 0.1302418f, 0.1652208f, 0.1984727f, 0.2299064f, 0.2594404f, 0.2870035f, 0.3125347f, 0.3359834f, 0.3573093f, 0.3764825f, 0.3934838f, 0.4083039f, 0.4209439f, 0.4314148f, 0.4397374f, 0.4459422f, 0.4500687f, 0.4521653f, 0.4522891f, 0.4505053f, 0.4468866f, 0.4415132f, 0.4344718f, 0.4258555f, 0.415763f, 0.4042981f, 0.391569f, 0.377688f, 0.3627708f, 0.3469356f, 0.3303029f, 0.3129944f, 0.2951329f, 0.2768413f, 0.258242f };

    vector<double> signal = {0.0000000e+00f, 2.4705567e-05f, 2.0080991e-04f, 6.7034760e-04f, 1.5285072e-03f, 2.7887512e-03f,
                       4.3621813e-03f, 6.0582370e-03f, 7.6099872e-03f, 8.9451009e-03f, 1.0357220e-02f, 1.1854203e-02f,
                       1.3434635e-02f, 1.5096378e-02f, 1.6836528e-02f, 1.8651363e-02f, 2.0536290e-02f, 2.2485809e-02f,
                       2.4493465e-02f, 2.6551804e-02f, 2.8652346e-02f, 3.0785553e-02f, 3.2940805e-02f, 3.5106368e-02f,
                       3.7269410e-02f, 3.9415959e-02f, 4.1530937e-02f, 4.3598153e-02f, 4.5600336e-02f, 4.7519144e-02f,
                       4.9335226e-02f, 5.1028263e-02f, 5.2577026e-02f, 5.3959444e-02f, 5.5152707e-02f, 5.6133337e-02f,
                       5.6877315e-02f, 5.7360191e-02f, 5.7557207e-02f, 5.7443462e-02f, 5.6994051e-02f, 5.6184221e-02f,
                       5.4989576e-02f, 5.3386237e-02f, 5.1351044e-02f, 4.8861768e-02f, 4.5897320e-02f, 4.2437952e-02f,
                       3.8465515e-02f, 3.3963643e-02f, 2.8918024e-02f, 2.3316601e-02f, 1.7149819e-02f, 1.0410842e-02f,
                       3.0957828e-03f, -4.7960808e-03f, -1.3262098e-02f, -2.2296047e-02f, -3.1887945e-02f,
                       -4.2023882e-02f,
                       -5.2685857e-02f, -6.3851640e-02f, -7.5494669e-02f, -8.7583944e-02f, -1.0008395e-01f,
                       -1.1295464e-01f,
                       -1.2615137e-01f, -1.3962500e-01f, -1.5332183e-01f, -1.6718377e-01f, -1.8114841e-01f,
                       -1.9514917e-01f,
                       -2.0911551e-01f, -2.2297308e-01f, -2.3664409e-01f, -2.5004750e-01f, -2.6309937e-01f,
                       -2.7571326e-01f,
                       -2.8780061e-01f, -2.9927114e-01f, -3.1003338e-01f, -3.1999511e-01f, -3.2906386e-01f,
                       -3.3714759e-01f,
                       -3.4415504e-01f, -3.4999654e-01f, -3.5458449e-01f, -3.5783401e-01f, -3.5966358e-01f,
                       -3.5999569e-01f,
                       -3.5875738e-01f, -3.5588104e-01f, -3.5130486e-01f, -3.4497359e-01f, -3.3683902e-01f,
                       -3.2686061e-01f,
                       -3.1500605e-01f, -3.0125168e-01f, -2.8558314e-01f, -2.6799560e-01f, -2.4849430e-01f,
                       -2.2709483e-01f,
                       -2.0382345e-01f, -1.7871724e-01f, -1.5182438e-01f, -1.2320417e-01f, -9.2927061e-02f,
                       -6.1074678e-02f,
                       -2.7739663e-02f, 6.9744838e-03f, 4.2953670e-02f, 8.0073528e-02f, 1.1819984e-01f, 1.5718900e-01f,
                       1.9688855e-01f, 2.3713785e-01f, 2.7776876e-01f, 3.1860632e-01f, 3.5946965e-01f, 4.0017286e-01f,
                       4.4052586e-01f, 4.8033547e-01f, 5.1940638e-01f, 5.5754209e-01f, 5.9454638e-01f, 6.3022399e-01f,
                       6.6438192e-01f, 6.9683081e-01f, 7.2738564e-01f, 7.5586730e-01f, 7.8210342e-01f, 8.0592954e-01f,
                       8.2719040e-01f, 8.4574068e-01f, 8.6144620e-01f, 8.7418485e-01f, 8.8384730e-01f, 8.9033806e-01f,
                       8.9357620e-01f, 8.9349574e-01f, 8.9004660e-01f, 8.8319486e-01f, 8.7292337e-01f, 8.5923189e-01f,
                       8.4213722e-01f, 8.2167357e-01f, 7.9789233e-01f, 7.7086210e-01f, 7.4066848e-01f, 7.0741349e-01f,
                       6.7121553e-01f, 6.3220865e-01f, 5.9054178e-01f, 5.4637843e-01f, 4.9989539e-01f, 4.5128208e-01f,
                       4.0073961e-01f, 3.4847954e-01f, 2.9472288e-01f, 2.3969890e-01f, 1.8364376e-01f, 1.2679933e-01f,
                       6.9411822e-02f, 1.1730438e-02f, -4.5994017e-02f, -1.0351054e-01f, -1.6056928e-01f,
                       -2.1692297e-01f,
                       -2.7232826e-01f, -3.2654709e-01f, -3.7934792e-01f, -4.3050709e-01f, -4.7981000e-01f,
                       -5.2705216e-01f,
                       -5.7204050e-01f, -6.1459404e-01f, -6.5454525e-01f, -6.9174051e-01f, -7.2604096e-01f,
                       -7.5732338e-01f,
                       -7.8548044e-01f, -8.1042135e-01f, -8.3207208e-01f, -8.5037583e-01f, -8.6529285e-01f,
                       -8.7680072e-01f,
                       -8.8489413e-01f, -8.8958496e-01f, -8.9090151e-01f, -8.8888872e-01f, -8.8360721e-01f,
                       -8.7513298e-01f,
                       -8.6355674e-01f, -8.4898311e-01f, -8.3152980e-01f, -8.1132692e-01f, -7.8851581e-01f,
                       -7.6324826e-01f,
                       -7.3568535e-01f, -7.0599639e-01f, -6.7435795e-01f, -6.4095265e-01f, -6.0596794e-01f,
                       -5.6959516e-01f,
                       -5.3202814e-01f, -4.9346226e-01f, -4.5409340e-01f, -4.1411650e-01f, -3.7372491e-01f,
                       -3.3310896e-01f,
                       -2.9245535e-01f, -2.5194588e-01f, -2.1175675e-01f, -1.7205767e-01f, -1.3301112e-01f,
                       -9.4771549e-02f,
                       -5.7484843e-02f, -2.1287683e-02f, 1.3692919e-02f, 4.7340058e-02f, 7.9547286e-02f, 1.1021888e-01f,
                       1.3927005e-01f, 1.6662708e-01f, 1.9222736e-01f, 2.1601939e-01f, 2.3796269e-01f, 2.5802764e-01f,
                       2.7619535e-01f, 2.9245722e-01f, 3.0681485e-01f, 3.1927946e-01f, 3.2987154e-01f, 3.3862048e-01f,
                       3.4556398e-01f, 3.5074741e-01f, 3.5422349e-01f, 3.5605150e-01f, 3.5629687e-01f, 3.5503030e-01f,
                       3.5232738e-01f, 3.4826782e-01f, 3.4293488e-01f, 3.3641466e-01f, 3.2879567e-01f, 3.2016793e-01f,
                       3.1062266e-01f, 3.0025154e-01f, 2.8914624e-01f, 2.7739778e-01f, 2.6509622e-01f, 2.5232998e-01f,
                       2.3918559e-01f, 2.2574714e-01f, 2.1209598e-01f, 1.9831035e-01f, 1.8446516e-01f, 1.7063159e-01f,
                       1.5687700e-01f, 1.4326461e-01f, 1.2985344e-01f, 1.1669817e-01f, 1.0384899e-01f, 9.1351636e-02f,
                       7.9247288e-02f, 6.7572616e-02f, 5.6359809e-02f, 4.5636635e-02f, 3.5426527e-02f, 2.5748687e-02f,
                       1.6618228e-02f, 8.0463104e-03f, 4.0317649e-05f, -7.3959655e-03f, -1.4262161e-02f,
                       -2.0561092e-02f,
                       -2.6298566e-02f, -3.1483162e-02f, -3.6126003e-02f, -4.0240526e-02f, -4.3842256e-02f,
                       -4.6948578e-02f,
                       -4.9578514e-02f, -5.1752489e-02f, -5.3492129e-02f, -5.4820038e-02f, -5.5759601e-02f,
                       -5.6334782e-02f,
                       -5.6569938e-02f, -5.6489646e-02f, -5.6118533e-02f, -5.5481113e-02f, -5.4601658e-02f,
                       -5.3504050e-02f,
                       -5.2211665e-02f, -5.0747264e-02f, -4.9132895e-02f, -4.7389805e-02f, -4.5538362e-02f,
                       -4.3598000e-02f,
                       -4.1587163e-02f, -3.9523251e-02f, -3.7422612e-02f, -3.5300497e-02f, -3.3171058e-02f,
                       -3.1047333e-02f,
                       -2.8941268e-02f, -2.6863704e-02f, -2.4824409e-02f, -2.2832101e-02f, -2.0894464e-02f,
                       -1.9018197e-02f,
                       -1.7209040e-02f, -1.5471818e-02f, -1.3810488e-02f, -1.2228184e-02f, -1.0727265e-02f,
                       -9.3093645e-03f,
                       -7.9664271e-03f, -6.3887839e-03f, -4.6444070e-03f, -3.0082264e-03f, -1.6805470e-03f,
                       -7.6002488e-04f,
                       -2.4136060e-04f, -3.4723744e-05f, -0.0000000e+00f};
    int step;
// List<DifEffectInstance> diffactions;

    diffsMap difs;

    int Check(double x)//находим, в зону приёма какого приёмнка попадает луч, попавший в точку x
    {
        for (int i = 0; i < coordinates.size(); i++)//для каждого элемента
        {
            if ((x >= coordinates[i] - radius) && (x <= coordinates[i] +
                                                        radius))//длина площадки сенсора radius*2+1, при 0 засчитывается только прямое попадание
            {
                return i;
            }
        }
        return -1;
    }

    vector<vector<double>> MakeConvolution() {
        auto convolvedData = vector<vector<double>>(coordinates.size())/*new double[coordinates.Length, ]*/;
        for(int i = 0; i<coordinates.size(); i++) {
            convolvedData[i] = vector<double>(maxTime + 1);
        }

        comp* tempColumn;
        auto impulseSpec = new comp[4096];
        for (int i = 0; i < signal.size(); i++) {
            //if (i < impulse.Length)
            impulseSpec[i].real(signal[i]);
        }

        fourier_transform(impulseSpec, 4096);
        for (int x = 0; x < width / 2; x++) {
            tempColumn = new comp[4096];
            for (int i = 0; i < min(4096, maxTime); i++) {
                tempColumn[i].real(recordedData[x][i]);
            }
            Convolution(tempColumn, impulseSpec);

            for (int i = 0; i <= maxTime; i++)//обрезаем то, что вылезло за границы
            {
                if (i < 4096) {
                    convolvedData[x][i] = tempColumn[i].real() / 4096;
                    //Math.Sqrt(tempColumn[i].Re * tempColumn[i].Re + tempColumn[i].Im * tempColumn[i].Im) / 4096;//(4096 * 4096);
                }else {
                    convolvedData[x][i] = 0;
                    //if (convolvedData[x, i] == 0) convolvedData[x, i] = (double)Math.Pow(10, -20);//&
                }
            }
            delete [] tempColumn;
        }
        delete [] impulseSpec;
        return convolvedData;
    }

    void Convolution(comp* tempColumn, comp* impulseSpec)//спектр импульса уже умноженный на себя
    {
        fourier_transform(tempColumn, 4096);
        for (int i = 0;i < 4096; i++) {
            tempColumn[i] *= impulseSpec[i];
        }
        inverse_fourier_transform(tempColumn, 4096);
    }

    void DrawDif(pointD p, double val, double prevTime) {
        int fIndex = GetNextFigure(p);//получаем скорость среды, по которой сигнал пойдёт вниз!
        int speed = figureCollection[fIndex].speed;

        int position = Check(p.X);


        double curTime = 0;
        double curVal = 1;
        int curFigureIndex = 0;
        vec2d vUp = vec2d(0, 1);

        //делаем хитрый ход, считаем луч не сверху вниз, а снизу вверх - разницы нет, а это поможет избежать громоздкого механизма

        beam curBeam = beam(pointD(p.X, 0), vUp, 0, 1, 0);//исходный луч
        intersection intersection = FindClosestIntersection(curBeam, false);//отправили луч с приёмника

        while (intersection.intersectionPoint.Y <= p.Y) {
            curTime = (double) ((intersection.distance * dX / figureCollection[curBeam.figureIndex].speed) / dt);//считаем время отрезка
            curVal = GetAbsorption(figureCollection[curBeam.figureIndex].absorption_c,
                                        intersection.distance);//считаем поглощение отрезка
            curFigureIndex = GetNextFigure(intersection.intersectionPoint, vUp);
            curBeam = beam(curBeam, intersection.intersectionPoint, vUp, curBeam.time + curTime,
                               curBeam.value * curVal,
                               curFigureIndex);
            intersection = FindClosestIntersection(curBeam, false);
        }

        curTime = curBeam.time;
        curVal = curBeam.value;

        double delta = (double) (p.Y * dX / speed / dt -
                                 curTime);//разность между расчетным и реальным временем

        double time = (double) (p.Y * dX / speed / dt) - delta + prevTime;
        if (time > maxTime)
            return;
        recordedData[position][(int) round(time)] += val * curVal / 2;
        //Console.WriteLine(position + "," + (int)Math.Round(time));

        bool stop = false;

        double newTime;

        for (int i = 1; stop != true; i++) {
            stop = true;
            newTime = (double) ((sqrt(step * step * i * i + p.Y * p.Y) * dX / speed) / dt) - delta +
                      prevTime;
            if (newTime > maxTime)
                return;
            if (position + i < coordinates.size()) {
                recordedData[position + i][(int) round(newTime)] += val * curVal / i;
                stop = false;
            }
            if (position - i >= 0) {
                recordedData[position - i][(int) round(newTime)] += val * curVal / i;
                stop = false;
            }
        }
    }
};


#endif //ARRAY_OF_RECEIVERS_TRANSMITTERS_H
