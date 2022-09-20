#ifndef _PGM_H_
#define _PGM_H_
#include <omp.h>
#include <stdio.h>
#include <time.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bitset> /*std::bitset<32>*/
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#define THREADS 16
using namespace std;
using namespace std::chrono;

enum PFM_endianness { BIG, LITTLE, ERROR };

class PFM {
   public:
    PFM() = default;
    inline bool is_little_big_endianness_swap() {
        if (this->endianess == 0.f) {
            std::cerr << "this-> endianness is not assigned yet!\n";
            exit(0);
        } else {
            uint32_t endianness = 0xdeadbeef;
            unsigned char* temp = (unsigned char*)&endianness;
            PFM_endianness endianType_ = ((*temp) ^ 0xef == 0     ? LITTLE
                                          : (*temp) ^ (0xde) == 0 ? BIG
                                                                  : ERROR);
            return ((BIG == endianType_) && (this->endianess < 0.f)) ||
                   ((LITTLE == endianType_) && (this->endianess > 0.f));
        }
    }

    template <typename T>
    T* read_pfm(const std::string& filename) {
        FILE* pFile;
        pFile = fopen(filename.c_str(), "rb");
        char c[100];
        if (pFile != NULL) {
            fscanf(pFile, "%s", c);
            // strcmp() returns 0 if they are equal.
            if (!strcmp(c, "Pf")) {
                fscanf(pFile, "%s", c);
                this->width = atoi(c);
                fscanf(pFile, "%s", c);
                this->height = atoi(c);
                int length_ = this->width * this->height;
                fscanf(pFile, "%s", c);
                this->endianess = atof(c);

                fseek(pFile, 0, SEEK_END);
                long lSize = ftell(pFile);
                long pos = lSize - this->width * this->height * sizeof(T);
                fseek(pFile, pos, SEEK_SET);

                T* img = new T[length_];
                // cout << "sizeof(T) = " << sizeof(T);
                fread(img, sizeof(T), length_, pFile);
                fclose(pFile);

                T* tbimg = (T*)malloc(length_ * sizeof(T));  // top-to-bottom.
                // PFM SPEC image stored bottom -> top reversing image
                for (int i = 0; i < this->height; i++) {
                    memcpy(&tbimg[(this->height - i - 1) * (this->width)],
                           &img[(i * (this->width))], (this->width) * sizeof(T));
                }

                if (this->is_little_big_endianness_swap()) {
                    std::cout << "little-big endianness transformation "
                                 "is needed.\n";
                    // little-big endianness transformation is needed.
                    union {
                        T f;
                        unsigned char u8[sizeof(T)];
                    } source, dest;

                    for (int i = 0; i < length_; ++i) {
                        source.f = tbimg[i];
                        for (unsigned int k = 0, s_T = sizeof(T); k < s_T; k++)
                            dest.u8[k] = source.u8[s_T - k - 1];
                        tbimg[i] = dest.f;
                        // cout << dest.f << ", ";
                    }
                }
                delete[] img;
                return tbimg;

            } else {
                std::cout << "Invalid magic number!"
                          << " No Pf (meaning grayscale pfm) is missing!!\n";
                fclose(pFile);
                exit(0);
            }

        } else {
            std::cout << "Cannot open file " << filename << ", or it does not exist!\n";
            fclose(pFile);
            exit(0);
        }
    }

    template <typename T>
    void write_pfm(const std::string& filename, const T* imgbuffer, const float& endianess_) {
        std::ofstream ofs(filename.c_str(), std::ifstream::binary);

        ofs << "Pf\n" << this->width << " " << this->height << "\n" << endianess_ << "\n";
        int length_ = this->width * this->height;
        this->endianess = endianess_;
        T* tbimg = (T*)malloc(length_ * sizeof(T));
        // PFM SPEC image stored bottom -> top reversing image
        for (int i = 0; i < this->height; i++) {
            memcpy(&tbimg[(this->height - i - 1) * this->width], &imgbuffer[(i * this->width)],
                   this->width * sizeof(T));
        }

        if (this->is_little_big_endianness_swap()) {
            std::cout << "little-big endianness transformation is needed.\n";
            // little-big endianness transformation is needed.
            union {
                T f;
                unsigned char u8[sizeof(T)];
            } source, dest;

            for (int i = 0; i < length_; ++i) {
                source.f = tbimg[i];
                for (size_t k = 0, s_T = sizeof(T); k < s_T; k++)
                    dest.u8[k] = source.u8[s_T - k - 1];
                tbimg[i] = dest.f;
                // cout << dest.f << ", ";
            }
        }

        ofs.write((char*)tbimg, this->width * this->height * sizeof(T));
        ofs.close();
        free(tbimg);
    }

    inline float getEndianess() { return endianess; }
    inline int getHeight(void) { return height; }
    inline int getWidth(void) { return width; }
    inline void setHeight(const int& h) { height = h; }
    inline void setWidth(const int& w) { width = w; }

   private:
    int height;
    int width;
    float endianess;
};

int Next8[8][2] = {{0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}};
int Next4[4][2] = {{0, 1}, {-1, 0}, {0, -1}, {1, 0}};
#endif /* PGM_H_ */

struct Coordinate {
    int X, Y;

    bool operator==(const Coordinate& that) const { return X == that.X && Y == that.Y; }
};

struct CoordinateHash {
    size_t operator()(const Coordinate& that) const { return that.X * size_t(9875321) + that.Y; }
};

/*
complex 64 64 8

*/

struct Solution {
    PFM pfm_rw;
    string path_in = "/home/tianxiangw/Desktop/CG/CPJFA";
    string path_out_pjfa = "/home/tianxiangw/Desktop/CG/CPJFA/large_out_pjfa.pfm";
    string path_out_jfa = "/home/tianxiangw/Desktop/CG/CPJFA/large_out_jfa.pfm";
    string path_out_BF = "/home/tianxiangw/Desktop/CG/CPJFA/large_out_BF.pfm";

    float* input = pfm_rw.read_pfm<float>(path_in);  //
    float* output = NULL;
    
    const int imgH = pfm_rw.getHeight();
    const int imgW = pfm_rw.getWidth();
    int size = imgH * imgW;
    float norm = sqrt(pow(imgH, 2) + pow(imgW, 2));

    void ResetOut() {
        delete output;
        output = new float[imgH * imgW * sizeof(float)];
        for (int i = 0; i < imgH * imgW; i++) {
            output[i] = 1.0;
        }
    }

    void Display(float* in) {
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                printf("%-3.2lf ", in[i * imgW + j]);
            }
            cout << endl;
        }
    }

    void InitialInput() {
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                if (input[i * imgW + j] < 0.8) {
                    input[i * imgW + j] = 0;
                } else {
                    input[i * imgW + j] = 1.0;
                }
            }
        }
    }

    void Recover(float*& In) {
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                if (In[i * imgW + j] != 0) {
                    In[i * imgW + j] = 1;
                }
            }
        }
    }

    void Solve() {
        freopen("/home/tianxiangw/Desktop/CG/CPJFA/report.txt", "w", stdout);
        printf("Input: %s\n", path_in.c_str());

        printf("Method: Brute Force, no parallel\n");
        ResetOut();
        auto start1 = high_resolution_clock::now();
        Cal1();
        pfm_rw.write_pfm<float>(path_out_BF, output, -1.0f);
        auto end1 = high_resolution_clock::now();
        auto ans1 = duration_cast<microseconds>(end1 - start1);
        printf("Run Time : %d microseconds\n", ans1);

        printf("Method: Parallel JFA through cpu\n");
        ResetOut();
        auto start2 = high_resolution_clock::now();
        PJFA();
        pfm_rw.write_pfm<float>(path_out_pjfa, output, -1.0f);
        auto end2 = high_resolution_clock::now();
        auto ans2 = duration_cast<microseconds>(end2 - start2);
        printf("Run Time : %d microseconds\n", ans2);

        printf("Method: None Parallel JFA through cpu\n");
        ResetOut();
        auto start3 = high_resolution_clock::now();
        JFA();
        pfm_rw.write_pfm<float>(path_out_jfa, output, -1.0f);
        auto end3 = high_resolution_clock::now();
        auto ans3 = duration_cast<microseconds>(end3 - start3);
        printf("Run Time : %d microseconds\n", ans3);
        printf("\n");

        fclose(stdout);
    }

    void Cal1() {
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                for (int x = 0; x < imgH; x++) {
                    for (int y = 0; y < imgW; y++) {
                        if (input[x * imgW + y] == 0) {
                            float dis = sqrt(pow(i - x, 2) + pow(j - y, 2)) / norm;
                            if (dis < output[i * imgW + j]) {
                                output[i * imgW + j] = dis;
                            }
                        }
                    }
                }
            }
        }
    }

    float Distance(const Coordinate& a, const Coordinate& b) {
        return sqrt(pow(a.X - b.X, 2) + pow(a.Y - b.Y, 2)) / norm;
    }

    void PJFA() {
        vector<pair<Coordinate, float>> InfoBest;
        InfoBest.resize(imgH * imgW);
        bool q[imgH * imgW];
        for (int i = 0; i < imgH * imgW; i++) {
            q[i] = false;
        }
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                if (input[i * imgW + j] == 0) {
                    q[i * imgW + j] = true;
                    InfoBest[i * imgW + j] = make_pair(Coordinate{i, j}, 0.0);
                } else {
                    InfoBest[i * imgW + j] =
                        make_pair(Coordinate{i, j}, numeric_limits<float>::max());
                }
            }
        }

        int stepLength = max(imgH, imgW) / 2;
        while (stepLength != 0) {
            bool newq[imgH * imgW];
            for (int i = 0; i < imgH * imgW; i++) {
                newq[i] = false;
            }
#pragma omp parallel for schedule(guided) num_threads(16)
            for (int x = 0; x < imgH; x++) {
#pragma omp parallel for schedule(guided) num_threads(16)
                for (int y = 0; y < imgW; y++) {
                    if (q[x * imgW + y]) {
                        newq[x * imgW + y] = true;
#pragma omp parallel for schedule(guided) num_threads(8)
                        for (int i = 0; i < 8; i++) {
                            Coordinate newC = {x + Next8[i][0] * stepLength,
                                               y + Next8[i][1] * stepLength};
                            if (newC.X < 0 || newC.X >= imgH) continue;
                            if (newC.Y < 0 || newC.Y >= imgW) continue;
#pragma omp atomic write
                            newq[newC.X * imgW + newC.Y] = true;
                            float dis = Distance(InfoBest[x * imgW + y].first, newC);
                            if (dis < InfoBest[newC.X * imgW + newC.Y].second) {
                                InfoBest[newC.X * imgW + newC.Y].first =
                                    InfoBest[x * imgW + y].first;
                                InfoBest[newC.X * imgW + newC.Y].second = dis;
                            }
                        }
                    }
                }
            }
            stepLength /= 2;
            for (int i = 0; i < imgH * imgW; i++) {
                q[i] = newq[i];
            }
        }

        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                output[i * imgW + j] = InfoBest[i * imgW + j].second;
            }
        }
    }

    void JFA() {
        vector<pair<Coordinate, float>> InfoBest;
        InfoBest.resize(imgH * imgW);
        bool q[imgH * imgW];
        for (int i = 0; i < imgH * imgW; i++) {
            q[i] = false;
        }
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                if (input[i * imgW + j] == 0) {
                    q[i * imgW + j] = true;
                    InfoBest[i * imgW + j] = make_pair(Coordinate{i, j}, 0.0);
                } else {
                    InfoBest[i * imgW + j] =
                        make_pair(Coordinate{i, j}, numeric_limits<float>::max());
                }
            }
        }

        int stepLength = max(imgH, imgW) / 2;
        while (stepLength != 0) {
            bool newq[imgH * imgW];
            for (int i = 0; i < imgH * imgW; i++) {
                newq[i] = false;
            }
            for (int x = 0; x < imgH; x++) {
                for (int y = 0; y < imgW; y++) {
                    if (q[x * imgW + y]) {
                        newq[x * imgW + y] = true;
                        for (int i = 0; i < 8; i++) {
                            Coordinate newC = {x + Next8[i][0] * stepLength,
                                               y + Next8[i][1] * stepLength};
                            if (newC.X < 0 || newC.X >= imgH) continue;
                            if (newC.Y < 0 || newC.Y >= imgW) continue;
                            newq[newC.X * imgW + newC.Y] = true;
                            float dis = Distance(InfoBest[x * imgW + y].first, newC);
                            if (dis < InfoBest[newC.X * imgW + newC.Y].second) {
                                InfoBest[newC.X * imgW + newC.Y].first =
                                    InfoBest[x * imgW + y].first;
                                InfoBest[newC.X * imgW + newC.Y].second = dis;
                            }
                        }
                    }
                }
            }
            stepLength /= 2;
            for (int i = 0; i < imgH * imgW; i++) {
                q[i] = newq[i];
            }
        }

        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                output[i * imgW + j] = InfoBest[i * imgW + j].second;
            }
        }
    }
};

int main() {
    Solution().Solve();
    return 1;
}
