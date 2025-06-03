#pragma once
#include <vector>

extern int Depth;
extern double Point;
extern int Backstep;
extern double Deviation;
extern int Digits;
extern int limit;

extern std::vector<double> sampleBuf;
extern std::vector<double> highestBuf;
extern std::vector<double> lowestBuf;

constexpr double INF = std::numeric_limits<double>::infinity();
const std::tm ZTIME = { 0,0,0,2,3,126,0,0,-1 };
const double PI = 3.14159265358979323846;
const double RADIAN2DEGREE = 180.0 / PI;
const double DEGREE2RADIAN = PI / 180.0;
double action = 0;


struct Statis {
    int N = 0;
    double Mean = 0.0;
    double Stdev = 0.0;
    double Min = INF;
    double Max = -INF;
    double RegAngle = 0.0;
    double RegAr = 0.0;
    double RegBr = 0.0;
};

extern "C" __declspec(dllexport)
Statis RunANN(const long* timeArr, const double* openArr, const double* highArr, const double* lowArr, const double* closeArr, int size,int digits, double point, int depth = 14);

int InitializeAll(int Bars, int depth, int digits, double point);

int FindFirstExtremum(double& curlow, double& curhigh, int& whatlookfor, int& limit);

int Lowest(const std::vector<double>& priceArray, int count, int start);

int Highest(const std::vector<double>& priceArray, int count, int start);

int LawnPeak(double& lasthigh, double& lastlow,
    const std::vector<double>& high,
    const std::vector<double>& low, int limit);

int FinalCutting(double& curlow, double& curhigh,
    double& lastlow, double& lasthigh,
    int& whatlookfor,
    const std::vector<double>& high,
    const std::vector<double>& low, int limit);

int CalculateLawnPeak(const std::vector<double>& high, const std::vector<double>& low);

double QuantStatsNoneZero(const std::vector<double>& fx, Statis& stc,int offset, int zlimit, double empty_value = 0.0, int direction = -1);

double DateTimeStruct(const long* times, int size);
