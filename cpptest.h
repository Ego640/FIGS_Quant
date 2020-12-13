//
// Created by 齐瑜鹏 on 2020/11/25.
//
using namespace std;
#ifndef CPPTEST_CPPTEST_H
#define CPPTEST_CPPTEST_H

#endif //CPPTEST_CPPTEST_H
struct item{
    double Coeff;
    int Scan;
    string Sequence;
    int Charge;
    double PrecursorMZ;
    double RetentionTime;
    double Correlation;
};
struct partOfItem{
    string Sequence;
    int Charge;
    bool operator==(const struct partOfItem& p1) const{
        return (Sequence==p1.Sequence&&Charge==p1.Charge);
    }
};
struct partOfItem1{
    int Scan;
    double RetentionTime;
    double Coeff;
    double PrecursorMZ;
    double Correlation;
    bool operator==(const struct partOfItem1& p1) const{
        return (RetentionTime==p1.RetentionTime);
    }
};
struct partOfItem2{
    int Scan;
    double RetentionTime;
    double Coeff;
    double PrecursorMZ;
    bool operator==(const struct partOfItem2& p2) const{
        return (RetentionTime==p2.RetentionTime);
    }
};
struct headItem{
    double PrecursorMZ;
    double retentionTime;
    int seqNum;
    int msLevel;
};
struct partOfHeader{
    int seqNum;
    double retentionTime;
    double Coeff;
};
struct Quants{
    string Sequence;
    int Charge;
    double Quantity;
    double LBPval;
    double SNR;
    double MaxCoeffTime;
    double MaxCoeff;
    double PeakWidth;
    double LBPvalOnUniformGrid;
    double PeakVariance;
    double PeakSkewness;
    double PeakKutosis;
    double PeakStart;
    double PeakEnd;
    double Correlation;
    bool isTypeDecoy=false;
};
struct resQuants{
    string Sequence;
    int Charge;
    bool isTypeDecoy;
    double score;
    double Quantity;
    double PeakStart;
    double PeakEnd;
};
struct peak{
    double maxValue;
    int left;
    int top;
    int right;
};