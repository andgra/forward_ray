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
#include "mpi.h"

#if defined(_WIN32)

#include <direct.h>
#include <windows.h>
#include <conio.h>

#else
#include <sys/stat.h>
#include <stdlib.h>
#endif


using std::string;
using std::vector;
using std::stringstream;
using std::getline;
using std::ifstream;
using std::istream_iterator;
using std::to_string;
using std::cout;
using std::exception;
using std::endl;
using std::ios;
using std::ofstream;

string basePath(string const pathname) {
    return pathname.substr(pathname.find_last_of("/\\") + 1);
}

string dirnameOf(const string fname) {
    size_t pos = fname.find_last_of("\\/");
    return (string::npos == pos)
           ? ""
           : fname.substr(0, pos);
}

string fileNameWoExt(string const pathname) {
    string fileName = basePath(pathname);
    string::size_type const p(fileName.find_last_of('.'));
    return fileName.substr(0, p);
}


#if defined(_WIN32)

int DeleteDirectory(const std::string &refcstrRootDirectory,
                    bool bDeleteSubdirectories = true) {
    bool bSubdirectory = false;       // Flag, indicating whether
    // subdirectories have been found
    HANDLE hFile;                       // Handle to directory
    std::string strFilePath;                 // Filepath
    std::string strPattern;                  // Pattern
    WIN32_FIND_DATA FileInformation;             // File information


    strPattern = refcstrRootDirectory + "\\*.*";
    hFile = ::FindFirstFile(strPattern.c_str(), &FileInformation);
    if (hFile != INVALID_HANDLE_VALUE) {
        do {
            if (FileInformation.cFileName[0] != '.') {
                strFilePath.erase();
                strFilePath = refcstrRootDirectory + "\\" + FileInformation.cFileName;

                if (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                    if (bDeleteSubdirectories) {
                        // Delete subdirectory
                        int iRC = DeleteDirectory(strFilePath, bDeleteSubdirectories);
                        if (iRC)
                            return iRC;
                    } else
                        bSubdirectory = true;
                } else {
                    // Set file attributes
                    if (::SetFileAttributes(strFilePath.c_str(),
                                            FILE_ATTRIBUTE_NORMAL) == FALSE)
                        return ::GetLastError();

                    // Delete file
                    if (::DeleteFile(strFilePath.c_str()) == FALSE)
                        return ::GetLastError();
                }
            }
        } while (::FindNextFile(hFile, &FileInformation) == TRUE);

        // Close handle
        ::FindClose(hFile);

        DWORD dwError = ::GetLastError();
        if (dwError != ERROR_NO_MORE_FILES)
            return dwError;
        else {
            if (!bSubdirectory) {
                // Set directory attributes
                if (::SetFileAttributes(refcstrRootDirectory.c_str(),
                                        FILE_ATTRIBUTE_NORMAL) == FALSE)
                    return ::GetLastError();

                // Delete directory
                if (::RemoveDirectory(refcstrRootDirectory.c_str()) == FALSE)
                    return ::GetLastError();
            }
        }
    }

    return 0;
}

#else

int DeleteDirectory(const std::string &refcstrRootDirectory)
{

system(("rm -r "+refcstrRootDirectory).c_str());
}
#endif

void RefreshDir(const string path) {
    DeleteDirectory(path);
#if defined(_WIN32)
//        _rmdir (path.c_str());
    Sleep(100);
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0777);
#endif
}


class algo {
public:
    static void OpenFile(const string path) {
        filePath = path;
        fileDir = dirnameOf(path);
        fileName = fileNameWoExt(path);
        vector<string> lines;

        ifstream file(path);
        if (!file) {
            throw string("file was not found");
        }
        copy(istream_iterator<string>(file),
             istream_iterator<string>(),
             back_inserter(lines));

        int cntLines = lines.size();

        auto size = split(lines[0], ';');//нулевая строка - размеры
        width = stoi(size[0]);
        height = stoi(size[1]);


        vector<string> bckgrnd = split(lines[1],
                                       ';');//первая строка - параметры фона(среды, не принадлежащей ни к одному из заданных объектов)

        background = figure(vector<edge>{edge(pointI(0, 0), pointI(width - 1, 0), 0, 0),
                                         edge(pointI(0, height - 1), pointI(width - 1, height - 1), 1, 0)},
                            stoi(bckgrnd[0]), stod_c(bckgrnd[1]), stod_c(bckgrnd[2]), 0);
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
        maxTime = (int)(dY * height / minSpeedCollection() / dt * 2);
//        maxTime = (int) ((sqrt(pow(dY * height, 2) + pow(dX * width / 2, 2)) /
//                          minSpeedCollection() / dt * 2));
    }


    static vector<vec2f> GenerateArrayOfVectors(int treshold, double step, bool onlyUp = false) {
        vector<vec2f> temp;
        double normal = 90;
        int size;

        if (onlyUp)
            size = 0;
        else
            size = (int) ((normal - treshold) / step);

        vector<vec2f> directions = vector<vec2f>(size + 1);

        directions[0] = vec2f(0, 1);//всегда добавляется вертикальный вектор
        directions[1] = vec2f(1, 0);//и прямая волна
        directions[2] = vec2f(-1, 0);

        for (int i = 3; i < size + 1; i++) {
            if (i % 2 == 0) {
                directions[i] = vec2f((float) cos((normal + (i - 2) * step) * degreeToRadians),
                                      (float) sin((normal + (i - 2) * step) * degreeToRadians));
            } else {
                directions[i] = vec2f((float) cos((normal - (i - 2) * step) * degreeToRadians),
                                      (float) sin((normal - (i - 2) * step) * degreeToRadians));
            }
        }
        return directions;
    }

    static void MainAlgorithm() {
//        int hlf_thr = omp_get_max_threads() / 2;
//        cout << "threads: " << omp_get_max_threads() << endl;
//        omp_set_nested(1);
//        cout << "nested: " << omp_get_nested() << endl;
        int step = 2;//расстояние между датчиками
        vector<int> coordinatesOfTransmitters = vector<int>(width / step);
        int cntCoord = coordinatesOfTransmitters.size();
        for (int i = 0; i < cntCoord; i++) {
            coordinatesOfTransmitters[i] = i * step;
        }
        Timer mtmr;
        double mt1 = mtmr.elapsed();

        int sectorThreshold = 30;
        double beamGradStep = 0.1;
        vector<vec2f> directions = GenerateArrayOfVectors(sectorThreshold, beamGradStep);//0.0001);

        int done = 0;
        int busy_threads = 0;
        int busy_transmitters = 0;



        int rank, size;
        int buff[1];
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank == 0) {
            string outDir = getOutPath();
            cout << "dir: " << outDir << endl;
            RefreshDir(outDir);

            vector<string> headers = vector<string>();

            headers.push_back("MPI ranks: " + to_string(size));
            headers.push_back("OpenMP threads: " + to_string(omp_get_max_threads()));
            headers.push_back("Sector threshold: " + to_string(sectorThreshold));
            headers.push_back("Beam angle step: " + to_string(beamGradStep));

            unsigned int maxHeaderSize = 0;
            for (auto h: headers) {
                unsigned int headerSize = h.length();
                maxHeaderSize = max(maxHeaderSize, headerSize);
            }

            char borderSymbol = '$';
            unsigned int borderCnt = 3;

            string verticalBorder = string(borderCnt * 2 + 2 + maxHeaderSize, borderSymbol);
            string horizontalBorder = string(borderCnt, borderSymbol);

            string header;
            header += "\n";
            header += verticalBorder + "\n";
            for (auto h: headers) {
                header += horizontalBorder + " " + h + string(maxHeaderSize - h.length(), ' ') + " " + horizontalBorder + "\n";
            }
            header += verticalBorder + "\n";


            cout << header << endl;

            // после подготовки папки, можем начинать работу в других процессах
            for (int r = 1; r < size; r++) {
                MPI_Send(buff, 1, MPI_INT, r, 1, MPI_COMM_WORLD);
            }
        } else {
            // ждем готовность папки для вывода
            MPI_Status sC;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &sC);
        }

        bool is_parallel = size > 1;

//        std::cout << "rank: " << rank << std::endl;
//        return;
        //не создавать разные варианты генерации с одним и тем же именем!
//#pragma omp parallel for schedule(dynamic, 3) num_threads(hlf_thr)
        for (int i = rank; i < cntCoord; i += size) {
//#pragma omp critical
            {
                busy_transmitters++;
                busy_threads++;
            };
            //если файл с данным именем уже существует, мы считаем, что он создан раньше и уже посчитан
            if (!ifstream("L" + fileName + ".data" + to_string(i)))
            {
                arrayOfReceiversTransmitters receivers = arrayOfReceiversTransmitters(coordinatesOfTransmitters,
                                                                                      directions,
                                                                                      1,
                                                                                      maxTime);//запускаем новую фиксацию

                std::cout  << endl << "-----" << endl << "started " << i << endl << "-----" << endl;
//            int j = 0;
                int dSize = directions.size();
                Timer tmr;
                double t1 = tmr.elapsed();
//                int remains = cntCoord - 2 - done;
//                cout << "remains: " << remains << endl;
//                remains = remains < 1 ? 1 : remains;
//                int using_thr = max(hlf_thr / remains, 2);
//                cout << "using threads: " << using_thr << endl;
//#pragma omp parallel for schedule(dynamic)
                for (int j = 0; j < dSize; j++)
                { //добавляем все направления расчёта луча из данной точки в очередь

#pragma omp critical
                    {
                        busy_threads++;
                    };
//                std::cout << "thr " << omp_get_num_threads() << std::endl;
                    vec2f d = directions[j];

                    queue<beam> calculationQueue = queue<beam>();
                    beam current;
//                tmr.reset();
                    calculationQueue.push(beam(pointF(coordinatesOfTransmitters[i], 0), d, 0, 1, 0));

                    while (!calculationQueue.empty()) {
                        current = calculationQueue.front();
                        calculationQueue.pop();
                        current.direction = vec2f::normalize(current.direction);

                        try {
//                        Timer tmrb;
//                        double tb0 = tmrb.elapsed();
                            receivers.SendBeam(current, calculationQueue);
//                        double tb = tmrb.elapsed() - tb0;
//                        if(tb != 0)
//                        cout << "b " << tb << endl;
                        }
                        catch (const exception &ex) {
                            cout << "!" << ex.what();
                        }
                    }
//                auto c = 1;
//                j++;

#pragma omp critical
                    {
                        busy_threads--;
                    };
                }
                double t = tmr.elapsed() - t1;
                std::cout << "-----" << endl << "ended " << i << "; time " << t << endl << "-----";
                receivers.ProcessDifs();//Отрисовываем зафиксированные диф.

                auto convolved = receivers.convolvedData();

//#pragma omp critical
                {
                    //записываем результаты данной фиксации
                    SaveData(
                        getOutPath() + "L" + fileName + ".data" + to_string(i),
                        convolved, coordinatesOfTransmitters.size()
                    );
                }
            }
            done++;

//#pragma omp critical
            {
                busy_transmitters--;
                busy_threads--;
            };
        }
        std::cout << std::endl;


        if (rank == 0) {
            int waiting = size - 1;
            while(waiting > 0) {
                MPI_Status sR;
                MPI_Recv(buff, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &sR);
                waiting--;
                cout << "taken signal from " << sR.MPI_SOURCE << " source" << endl;
                cout << "waiting " << waiting << " more" << endl << endl;
            }
            SaveInfo(getOutPath() + "L" + fileName + ".data", maxTime, coordinatesOfTransmitters.size(), step);
            //SaveInfoInUniversalFormat("L3-0.5_2D-z80_150.data", maxTime, coordinatesOfTransmitters.Length, step);
            cout << "Finished" << endl;
            double mt = mtmr.elapsed() - mt1;
            std::cout << "total time " << mt << std::endl;
        } else {
            buff[0] = rank;
            MPI_Send(buff, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

private:

    static void SaveInfo(string path, int maxtime, int cnt, int step)//, int impulseLen)
    {
        float dxmm = 1000 * dX;
        float dymm = 1000 * dY;
        float df = (int) (1 / dt);
        ofstream myFile(path + ".info", ios::out | ios::binary);
        myFile.write(reinterpret_cast<char *>(&cnt), sizeof(cnt));
        myFile.write(reinterpret_cast<char *>(&step), sizeof(step));
        myFile.write(reinterpret_cast<char *>(&maxtime), sizeof(maxtime));
        myFile.write(reinterpret_cast<char *>(&dxmm), sizeof(dxmm));
        myFile.write(reinterpret_cast<char *>(&dymm), sizeof(dymm));
        myFile.write(reinterpret_cast<char *>(&df), sizeof(df));

        cout << basePath(path + ".info") << " saved;" << endl << endl;
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

    static void SaveData(string path, vector<vector<float>> data, int cnt)//, int impulseLen)
    {
        //SaveSpeedMap();
        ofstream myFile(path, ios::out | ios::binary);
        for (auto v: data) {
            for (auto f: v) {
                myFile.write(reinterpret_cast<char *>(&f), sizeof(f));
            }
        }

        cout << "-----" << endl << basePath(path) << " saved;" << endl << "-----";
    }
};

#endif //ALGO_H
