#ifndef _PGM_H_
#define _PGM_H_
#include <stdio.h>

#include <algorithm>
#include <bitset> /*std::bitset<32>*/
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

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

struct Solution {
    PFM pfm_rw;
    string path_in = "/Users/wayne_tx/Desktop/CG/singlesharp_in.pfm";
    string path_out_C1 = "/Users/wayne_tx/Desktop/CG/onesquare_out_m1.pfm";
    string path_out_C2 = "/Users/wayne_tx/Desktop/CG/onesquare_out_m2.pfm";
    string path_out_C2_r2 = "/Users/wayne_tx/Desktop/CG/onesquare_out_m2_r2.pfm";
    string path_out_C2_r3 = "/Users/wayne_tx/Desktop/CG/onesquare_out_m2_r3.pfm";
    string path_out_jfa = "/Users/wayne_tx/Desktop/CG/Complex/singlesharp_out_jfa.pfm";
    string path_in_win = "C:/CG/singlesharp_in.pfm";
    string path_out_win = "C:/CG/singlesharp_in.pfm";
    float* input = pfm_rw.read_pfm<float>(path_in_win);
    float* output = NULL;
    int imgH = pfm_rw.getHeight();
    int imgW = pfm_rw.getWidth();
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

    void Solve() {
        // freopen("/Users/wayne_tx/Desktop/CG/ot.txt", "w", stdout);
        // cout << " path in " << path_in << endl;
        // cout << "Height " << imgH << " Width " << imgW << endl;
        InitialInput();
        pfm_rw.write_pfm<float>(path_in_win, input, -1.0f);
        // ResetOut();
        // Cal1();
        // pfm_rw.write_pfm<float>(path_out_C1, output, -1.0f);
        // ResetOut();
        // Cal2(1);
        // pfm_rw.write_pfm<float>(path_out_C2, output, -1.0f);
        // ResetOut();
        // Cal2(2);
        // pfm_rw.write_pfm<float>(path_out_C2_r2, output, -1.0f);
        // ResetOut();
        // Cal2(3);
        // pfm_rw.write_pfm<float>(path_out_C2_r3, output, -1.0f);
        ResetOut();
        JFA();
        pfm_rw.write_pfm<float>(path_out_win, output, -1.0f);
        // fclose(stdout);
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

    void FloodFill(int x, int y, int r) {
        vector<bool> visited;
        visited.resize(imgH * imgW);
        queue<Coordinate> q;
        q.push(Coordinate{x, y});
        visited[x * imgW + y] = true;
        int reachtime = 0;

        while (true) {
            queue<Coordinate> newq;
            bool reach = false;
            while (!q.empty()) {
                Coordinate front = q.front();
                q.pop();
                if (!input[front.X * imgW + front.Y]) {
                    reach = true;
                    float dis = sqrt(pow(front.X - x, 2) + pow(front.Y - y, 2)) / norm;
                    if (output[x * imgW + y] > dis) {
                        output[x * imgW + y] = dis;
                    }
                }
                for (int i = 0; i < 8; i++) {  // 加入新的点到newq
                    int newX = front.X + Next8[i][0];
                    if (newX < 0 || newX >= imgH) continue;
                    int newY = front.Y + Next8[i][1];
                    if (newY < 0 || newY >= imgW) continue;
                    if (visited[newX * imgW + newY]) continue;
                    newq.push(Coordinate{newX, newY});
                    visited[newX * imgW + newY] = true;
                }
            }
            if (reach) reachtime++;
            if (reachtime == r) {
                return;
            }
            swap(q, newq);
        }
        return;
    }

    void Cal2(int r) {
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                FloodFill(i, j, r);
            }
        }
    }

    float Distance(Coordinate a, Coordinate b) {
        return sqrt(pow(a.X - b.X, 2) + pow(a.Y - b.Y, 2)) / norm;
    }

    void OneSeedFlood(Coordinate start, vector<pair<Coordinate, float>>& infobest) {
        unordered_set<Coordinate, CoordinateHash> q;
        q.emplace(start);
        int maxSize = max(imgH, imgW);
        int stepLength = max(imgH, imgW) / 2;
        while (stepLength != 0) {
            unordered_set<Coordinate, CoordinateHash> newQ;
            newQ.emplace(start);
            for (const auto& c : q) {
                for (int i = 0; i < 8; i++) {
                    Coordinate newC = {c.X + Next8[i][0] * stepLength,
                                       c.Y + Next8[i][1] * stepLength};
                    newQ.emplace(newC);
                    if (newC.X < 0 || newC.X >= imgH) continue;
                    if (newC.Y < 0 || newC.Y >= imgW) continue;
                    float dis = Distance(start, newC);
                    if (dis < infobest[newC.X * imgW + newC.Y].second) {
                        infobest[newC.X * imgW + newC.Y].first = newC;
                        infobest[newC.X * imgW + newC.Y].second = dis;
                    }
                }
            }
            stepLength /= 2;
            swap(q, newQ);
        }
    }

    void JFA() {
        vector<pair<Coordinate, float>> InfoBest;
        queue<Coordinate> q;
        vector<Coordinate> Seeds;
        int maxSize = max(imgH, imgW);
        InfoBest.resize(maxSize * maxSize);
        for (int i = 0; i < imgH; i++) {
            for (int j = 0; j < imgW; j++) {
                if (input[i * imgW + j] == 0) {
                    InfoBest[i * imgW + j] = make_pair(Coordinate{i, j}, 0.0);
                    Seeds.emplace_back(Coordinate{i, j});
                } else {
                    InfoBest[i * imgW + j] = make_pair(Coordinate{i, j}, numeric_limits<float>::max());
                }
            }
        }
        for (int i = 0; i < Seeds.size(); i++) {
            OneSeedFlood(Seeds[i], InfoBest);
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
