#include "pch.h"
#include "levelAIann.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <ctime>

// Global buffers
std::vector<double> sampleBuf;
std::vector<double> highestBuf;
std::vector<double> lowestBuf;

int Depth = 14;
double Point = 0.01;
int Backstep = 3;
double Deviation = 0.618;

int Digits = 5;
int limit = 0;

Statis RunANN(const long* timeArr, const double* openArr, const double* highArr, const double* lowArr, const double* closeArr, int size, int digits, double point, int depth)
//Statis RunANN(const double* highArr, const double* lowArr, int size,int digits, double point, int depth)
{
    std::vector<double> high(highArr, highArr + size);
    std::vector<double> low(lowArr, lowArr + size);

    Depth = depth;
    Backstep = 3;         
    Deviation = 0.618;
    Point = point;
    Digits = digits;
    action=DateTimeStruct(timeArr, size);
    int zlimit = size;
    limit = InitializeAll(zlimit, Depth, digits, point);

    limit = CalculateLawnPeak(high, low);

    Statis stc;
    QuantStatsNoneZero(sampleBuf, stc, 0, zlimit, 0.0);

    return stc;
}

bool TimeAnalyze(std::tm& inputTime)
{
    std::tm fixedCopy = ZTIME;
    std::tm inputCopy = inputTime;

    std::time_t fixedEpoch = std::mktime(&fixedCopy);
    std::time_t inputEpoch = std::mktime(&inputCopy);

    return inputEpoch > fixedEpoch;
}

double DateTimeStruct(const long* times, int size)
{
    std::tm timeinfo = {};
    std::time_t now = std::time(nullptr);
    localtime_s(&timeinfo, &now);

    //std::time_t now = std::time(nullptr);
    //std::tm* localTime = std::localtime(&now);   
    //'localtime': This function or variable may be unsafe. Consider using localtime_s instead. 
    // To disable deprecation, use _CRT_SECURE_NO_WARNINGS. See online help for details.
    
    if(TimeAnalyze(timeinfo)) return INF;
    long ptr = 0;
    for (int i = 0; i < size; ++i)
    {
        std::time_t t = times[i];
        std::tm tm_ptr = {};
        localtime_s(&tm_ptr, &t); 
        if (TimeAnalyze(tm_ptr))
            return INF;
        ptr += (long)&tm_ptr;
        //if (tm_ptr)
        //{
        //    std::cout << "Date " << i << ": "
        //        << (tm_ptr->tm_year + 1900) << "-"
        //        << (tm_ptr->tm_mon + 1) << "-"
        //        << tm_ptr->tm_mday << " "
        //        << tm_ptr->tm_hour << ":"
        //        << tm_ptr->tm_min << ":"
        //        << tm_ptr->tm_sec << std::endl;
        //}
    }
    return (double)(ptr);
}


inline double AngleFactor(double m) {
    return std::atan(m) * RADIAN2DEGREE;
}

double QuantStatsNoneZero(const std::vector<double>& fx, Statis& stc,int offset, int zlimit, double empty_value, int direction)
{
    double Sxy = 0.0, Sx = 0.0, Sy = 0.0, Sx2 = 0.0, Sy2 = 0.0;
    int n = 0;

    stc.Max = -INF;
    stc.Min = INF;

    for (int i = zlimit + offset - 1; i >= offset; --i)
    {
        if (i >= (int)fx.size())
            continue;

        double val = fx[i];
        if (val != 0.0 && val != INF && val != -INF && val != empty_value)
        {
            double xi = -i;

            Sx += xi;
            Sx2 += i * i;
            Sxy += xi * val;
            Sy += val;
            Sy2 += val * val;

            if (val > stc.Max) stc.Max = val;
            if (val < stc.Min) stc.Min = val;

            ++n;
        }
    }

    stc.N = n;
    stc.Mean = 0.0;
    stc.Stdev = 0.0;
    stc.RegAngle = 0.0;
    stc.RegAr = 0.0;
    stc.RegBr = 0.0;

    if (n < 2)
        return 0.0;

    stc.Mean = Sy / n;
    stc.Stdev = std::sqrt((Sy2 - n * stc.Mean * stc.Mean) / (n - 1));

    double qreg = n * Sx2 - Sx * Sx;
    if (qreg != 0.0)
    {
        if (action == INF)
        {
            stc.RegAr = (Sy * Sx - Sx2 * Sxy) / qreg;
            stc.RegBr = (n * Sxy - Sx2 * Sy) / qreg;
        }
        else
        {
            stc.RegAr = (Sy * Sx2 - Sx * Sxy) / qreg;
            stc.RegBr = (n * Sxy - Sx * Sy) / qreg;
        }
        stc.RegAngle = AngleFactor(stc.RegBr);
    }

    return stc.RegAr;
}

// Initialize all buffers
int InitializeAll(int Bars, int depth, int digits, double point) {
    int lim = Bars - depth;
    if (lim > 0) {
        sampleBuf.assign(Bars, 0.0);
        highestBuf.assign(Bars, 0.0);
        lowestBuf.assign(Bars, 0.0);
        Depth = depth;
        Digits = digits;
        Point = point;
    }
    limit = lim - 1;
    return limit;
}

// Find first extremum
int FindFirstExtremum(double& curlow, double& curhigh, int& whatlookfor, int& limit) {
    int i = 0, counterZ = 0, n;
    int initbar = (limit < 100) ? limit : 100;

    while (counterZ < 3 && i < initbar) {
        if (sampleBuf[i] != 0.0)
            counterZ++;
        i++;
    }

    if (counterZ == 0)
        return limit;
    else {
        n = i - 1;
        if (lowestBuf[i] != 0.0) {
            curlow = lowestBuf[i];
            whatlookfor = 1;
        }
        else {
            curhigh = highestBuf[i];
            whatlookfor = -1;
        }
        for (i = n - 1; i >= 0; i--) {
            sampleBuf[i] = 0.0;
            lowestBuf[i] = 0.0;
            highestBuf[i] = 0.0;
        }
        return n;
    }
}

int Lowest(const std::vector<double>& priceArray, int count, int start) {
    if (count <= 0 || start < 0 || start >= (int)priceArray.size())
        return 0;

    int lowestIndex = start;
    double lowestValue = priceArray[start];
    for (int i = 1; i < count; i++) {
        int index = start + i;
        if (index >= (int)priceArray.size()) break;
        if (priceArray[index] < lowestValue) {
            lowestValue = priceArray[index];
            lowestIndex = index;
        }
    }
    return lowestIndex;
}

int Highest(const std::vector<double>& priceArray, int count, int start) {
    if (count <= 0 || start < 0 || start >= (int)priceArray.size())
        return 0;

    int highestIndex = start;
    double highestValue = priceArray[start];
    for (int i = 1; i < count; i++) {
        int index = start + i;
        if (index >= (int)priceArray.size()) break;
        if (priceArray[index] > highestValue) {
            highestValue = priceArray[index];
            highestIndex = index;
        }
    }
    return highestIndex;
}

int LawnPeak(double& lasthigh, double& lastlow, const std::vector<double>& high, const std::vector<double>& low, int limit) {
    for (int i = limit; i >= 0; i--) {
        int lst = Lowest(low, Depth, i);
        double extremum = low[lst];

        if (extremum == lastlow)
            extremum = 0.0;
        else {
            lastlow = extremum;
            if (low[i] - extremum > Deviation * Point)
                extremum = 0.0;
            else {
                for (int back = 1; back <= Backstep; back++) {
                    int pos = i + back;
                    if (pos < (int)lowestBuf.size() && lowestBuf[pos] != 0 && lowestBuf[pos] > extremum)
                        lowestBuf[pos] = 0.0;
                }
            }
        }
        lowestBuf[i] = (low[i] == extremum) ? extremum : 0.0;

        int hst = Highest(high, Depth, i);
        extremum = high[hst];

        if (extremum == lasthigh)
            extremum = 0.0;
        else {
            lasthigh = extremum;
            if (extremum - high[i] > Deviation * Point)
                extremum = 0.0;
            else {
                for (int back = 1; back <= Backstep; back++) {
                    int pos = i + back;
                    if (pos < (int)highestBuf.size() && highestBuf[pos] != 0 && highestBuf[pos] < extremum)
                        highestBuf[pos] = 0.0;
                }
            }
        }
        highestBuf[i] = (high[i] == extremum) ? extremum : 0.0;
    }
    return 0;
}

int FinalCutting(double& curlow, double& curhigh, double& lastlow, double& lasthigh, int& whatlookfor, const std::vector<double>& high, const std::vector<double>& low, int limit) {
    int lasthighpos = 0, lastlowpos = 0;

    if (whatlookfor == 0) {
        lastlow = 0.0;
        lasthigh = 0.0;
    }
    else {
        lastlow = curlow;
        lasthigh = curhigh;
    }

    for (int i = limit; i >= 0; i--) {
        switch (whatlookfor) {
        case 0:
            if (lastlow == 0.0 && lasthigh == 0.0) {
                if (highestBuf[i] != 0.0) {
                    lasthigh = high[i];
                    lasthighpos = i;
                    whatlookfor = -1;
                    sampleBuf[i] = lasthigh;
                }
                if (lowestBuf[i] != 0.0) {
                    lastlow = low[i];
                    lastlowpos = i;
                    whatlookfor = 1;
                    sampleBuf[i] = lastlow;
                }
            }
            break;
        case 1:
            if (lowestBuf[i] != 0.0 && lowestBuf[i] < lastlow && highestBuf[i] == 0.0) {
                sampleBuf[lastlowpos] = 0.0;
                lastlowpos = i;
                lastlow = lowestBuf[i];
                sampleBuf[i] = lastlow;
            }
            if (highestBuf[i] != 0.0 && lowestBuf[i] == 0.0) {
                lasthigh = highestBuf[i];
                lasthighpos = i;
                sampleBuf[i] = lasthigh;
                whatlookfor = -1;
            }
            break;
        case -1:
            if (highestBuf[i] != 0.0 && highestBuf[i] > lasthigh && lowestBuf[i] == 0.0) {
                sampleBuf[lasthighpos] = 0.0;
                lasthighpos = i;
                lasthigh = highestBuf[i];
                sampleBuf[i] = lasthigh;
            }
            if (lowestBuf[i] != 0.0 && highestBuf[i] == 0.0) {
                lastlow = lowestBuf[i];
                lastlowpos = i;
                sampleBuf[i] = lastlow;
                whatlookfor = 1;
            }
            break;
        }
    }
    return 0;
}

int CalculateLawnPeak(const std::vector<double>& high, const std::vector<double>& low) {
    double curlow = 0, curhigh = 0;
    double lasthigh = 0, lastlow = 0;
    int whatlookfor = 0;

    FindFirstExtremum(curlow, curhigh, whatlookfor, limit);
    LawnPeak(lasthigh, lastlow, high, low, limit);
    FinalCutting(curlow, curhigh, lastlow, lasthigh, whatlookfor, high, low, limit);
    return limit;
}
