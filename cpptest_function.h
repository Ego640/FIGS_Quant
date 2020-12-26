//
// Created by 齐瑜鹏 on 2020/11/26.
//
#include"cpptest.h"
#include<numeric>
#include<string>
#include<cassert>
#include<vector>
#include<fstream>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<boost/math/distributions/chi_squared.hpp>
using namespace std;
#define MIN_NUM -1000000
#ifndef CPPTEST_CPPTEST_FUNCTION_H
#define CPPTEST_CPPTEST_FUNCTION_H

#endif //CPPTEST_CPPTEST_FUNCTION_H
bool cmp(const struct item &i1,const struct item &i2){//Scan,Sequence,Charge
    if(i1.Sequence<i2.Sequence){
        return true;
    }
    else if(i1.Sequence==i2.Sequence){
        if(i1.Charge<i2.Charge){
            return true;
        }
        else if(i1.Charge==i2.Charge){
            if(i1.Scan<i2.Scan){
                return true;
            }
            else{
                return false;
            }
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}
void SparkCoeffs(string rawCoeffsFilePath,vector<struct item>& DIA_RefSpectraCoeffs){
    ifstream _csvInput(rawCoeffsFilePath,ios::in);
    string _Oneline;//用于存每一行的string
    vector<string> eachRow;
    const char* split=",";//分隔符
    char* p;
    //将每一行读出来的数据按逗号分隔，存在eachRow中
    double MaxRetentionTime=MIN_NUM;
    while(getline(_csvInput,_Oneline)){
        eachRow.clear();
        char* str=(char*)_Oneline.c_str();
        p=strtok(str,split);
        while(p!=NULL){
            eachRow.push_back(p);
            p=strtok(NULL,split);
        }
        if(eachRow.size()==8){//处理不完整数据，不完善的数据就不存在vector中
            if(MaxRetentionTime<stod(eachRow[5])){//对于RententionTime是否最大值是否小于600
                // 判断是在删除所有Coeffs只为0元素操作之前完成的
                // 所以我们对Coeffs=0的元素的RetentionTime也要进行比较从而找到其最大值
                MaxRetentionTime=stod(eachRow[5]);
            }
            if(stod(eachRow[0])!=0){//将Coeffs为0的不写入vector中
                struct item ii;//我们将上面存在eachRow中的值分别赋给ii的每个成员
                ii.Coeff=stod(eachRow[0]);
                /////目前会将其小数点之后的部分忽略掉，后面我们来解决这个问题
                ii.Scan=stoi(eachRow[1])+1;
                ii.Sequence=eachRow[2];
                ii.Charge=stoi(eachRow[3]);
                ii.PrecursorMZ=stod(eachRow[4]);
                ii.RetentionTime=stod(eachRow[5]);
                ii.Correlation=stod(eachRow[7]);
                DIA_RefSpectraCoeffs.push_back(ii);
            }
        }
    }
    sort(DIA_RefSpectraCoeffs.begin(),DIA_RefSpectraCoeffs.end(),cmp);//对eachRow中的结构体按照cmp函数进行排序
    if(MaxRetentionTime<600){
        for(int i=0;i<DIA_RefSpectraCoeffs.size();i++){
            DIA_RefSpectraCoeffs[i].RetentionTime*=60;
        }
    }
    /////SparkCoeff函数调试部分
    //cout<<DIA_RefSpectraCoeffs[49].Scan<<" "<<DIA_RefSpectraCoeffs[49].Coeff<<" "<<DIA_RefSpectraCoeffs[49].Sequence<<" "<<DIA_RefSpectraCoeffs[49].Charge<<" "<<DIA_RefSpectraCoeffs[49].PrecursorMZ<<" "<<DIA_RefSpectraCoeffs[49].RetentionTime<<" "<<DIA_RefSpectraCoeffs[49].Correlation<<endl;
}
bool cmppeak(struct peak p1,struct peak p2){
    return p1.maxValue>p2.maxValue;
}
void allmomemt(vector<double> &PeakCentralMoments,vector<double> in){
    double exp=0;
    for(int i=0;i<in.size();i++){
        exp+=in[i];
    }
    exp/=(1.0*in.size());//期望
    //计算第2矩
    double var=0;
    for(int i=0;i<in.size();i++){
        var+=(in[i]-exp)*(in[i]-exp);
    }
    var/=(1.0*in.size());
    //计算第3矩
    double skew=0;
    for(int i=0;i<in.size();i++){
        skew+=pow((in[i]-exp),3);
    }
    skew/=(1.0*in.size());
    //计算第4矩
    double kurt=0;
    for(int i=0;i<in.size();i++){
        kurt+=pow((in[i]-exp),4);
    }
    kurt/=(1.0*in.size());
    PeakCentralMoments[0]=var;
    PeakCentralMoments[1]=skew;
    PeakCentralMoments[2]=kurt;
}
void kz(vector<double> &x,int m,int k){//x是待滤波的数组，m为核的大小，k为迭代次数
    vector<double> kztemp(x.size());
    for(int i=0;i<k;i++){
        for(int j=0;j<x.size();j++){
            double posi_val=0;
            if(j<(m-1)/2){
                for(int r=0;r<=j+(m-1)/2;r++){
                    posi_val+=x[r];
                }
                posi_val/=(j+(m-1)/2+1);
                kztemp[j]=posi_val;
            }
            else if(j+(m-1)/2>=x.size()){
                for(int r=j-(m-1)/2;r<x.size();r++){
                    posi_val+=x[r];
                }
                posi_val/=(x.size()-j+(m-1)/2);
                kztemp[j]=posi_val;
            }
            else{
                for(int r=j-(m-1)/2;r<=j+(m-1)/2;r++){
                    posi_val+=x[r];
                }
                posi_val/=m;
                kztemp[j]=posi_val;
            }
        }
        for(int j=0;j<x.size();j++){
            x[j]=kztemp[j];
        }
    }
}
/*double boxtest(vector<double> x,string dirPath){
    //写输入文件
    ofstream outFile;
    string inFilePath=dirPath+"/in.csv";
    string RFilePath=dirPath+"/BoxTest.R";
    outFile.open(inFilePath,ios::out);
    outFile<<"0"<<endl;
    for(int i=0;i<x.size();i++){
        outFile<<to_string(x[i])<<endl;
    }
    outFile.close();
    //system函数运行R脚本
    string cmdstr="Rscript "+RFilePath+" "+dirPath;
    system(cmdstr.c_str());
    //读取r脚本运行结果
    string outFilePath=dirPath+"/res.txt";
    ifstream infile;
    infile.open(outFilePath.data());
    assert(infile.is_open());
    string s;
    getline(infile,s);
    double p=stod(s);
    infile.close();
    //判断in.csv和res.txt文件是否存在，若存在，删除两个文件
    ifstream fin(inFilePath.c_str());
    if(fin.good()){
        remove(inFilePath.c_str());
    }
    ifstream fres(outFilePath.c_str());
    if(fres.good()){
        remove(outFilePath.c_str());
    }
    return p;
}*/
double boxtest(vector<double> x){//默认滞后数lag=1
    double xsum=accumulate(x.begin(),x.end(),0);
    double xmeans=xsum/x.size();
    double fenzi=0,fenmu=0;
    for(int i=0;i<x.size();i++){
        fenmu+=((x[i]-xmeans)*(x[i]-xmeans));
    }
    for(int i=0;i<x.size()-1;i++){
        fenzi+=((x[i]-xmeans)*(x[i+1]-xmeans));
    }
    double p=fenzi/fenmu;
    double Q=1.0*x.size()*(x.size()+2)*p*p/(x.size()-1);
    //cout<<Q<<endl;
    boost::math::chi_squared mydist(1);
    return 1-boost::math::cdf(mydist,Q);
}
double trapz(vector<double> x,vector<double> y){//求积分，梯形法
    double res=0;
    for(int i=1;i<x.size();i++){
        res+=(y[i-1]+y[i])*(x[i]-x[i-1])/2;
    }
    return res;
}
double sd(vector<double> a){//求标准差函数
    double sum=accumulate(a.begin(),a.end(),0);
    double mean=sum/a.size();
    double accum=0;
    for(int i=0;i<a.size();i++){
        accum+=(a[i]-mean)*(a[i]-mean);
    }
    double stdev=sqrt(accum/(1.0*(a.size()-1)));
    return stdev;
}
void scale(vector<double> &a){//标准化
    double stdev=0;
    for(int i=0;i<a.size();i++){
        stdev+=(a[i]*a[i]);
    }
    stdev/=(1.0*(a.size()-1));
    stdev=sqrt(stdev);
    for(int i=0;i<a.size();i++){
        a[i]/=stdev;
    }
}
void findPeaks(vector<double> d,vector<struct peak> &f,int leftstep,int rightstep,int peaknum){//对应R语言包中的findpeaks函数
    bool up=false;
    int l=0,t=0;
    for(int r=1;r<d.size();r++){
        if(!up){
            if(d[r]>d[r-1]){
                t=r;
                if(t-l>=leftstep)  up=true;
            }
            else{
                l=r;t=r;
            }
        }
        else{
            if(d[r]>d[r-1]){
                if(r-t==1)  t=r;
                else if(r-t<=rightstep){
                    l=r-1; t=r;
                    up=false;
                }
                else{
                    struct peak p1;
                    p1.maxValue=d[t];
                    p1.left=l+1;
                    p1.top=t+1;
                    p1.right=r;
                    f.push_back(p1);
                    l=r-1; t=r;
                    up=false;
                }
            }
            else if(d[r]==d[r-1]){
                if(r-t==1)  t=r;
                else if(r-t<=rightstep){
                    l=r; t=r;
                    up=false;
                }
                else{
                    struct peak p1;
                    p1.maxValue=d[t];
                    p1.left=l+1;
                    p1.top=t+1;
                    p1.right=r;
                    f.push_back(p1);
                    l=r;t=r;
                    up=false;
                }
            }
            else{
                if(r==d.size()-1&&r-t>=2){
                    struct peak p1;
                    p1.maxValue=d[t];
                    p1.left=l+1;
                    p1.top=t+1;
                    p1.right=r+1;
                    f.push_back(p1);
                }
            }
        }
    }
    sort(f.begin(),f.end(),cmppeak);
    for(int j=f.size();j>=peaknum;j--){
        f.erase(f.begin()+j-1);
    }
    /////
    /*cout<<f.size()<<endl;
    for(int j=0;j<f.size();j++){
        cout<<f[j].maxValue<<" "<<f[j].left<<" "<<f[j].top<<" "<<f[j].right<<endl;
    }*/
}
void FindPeaksInData(vector<struct peak> &PeptidePeaks,vector<struct item> Data,vector<struct partOfItem> Identifiers,int i,vector<struct headItem> header,bool Multiplexed=false,int IntensityCutoff=0,int QuantileCutoff=0,bool RTWindow=true,string smooth="rollmean",int FilterWindow=3,int KZiters=3){
    vector<struct partOfItem2> PeptideData;
    for(int j=0;j<Data.size();j++){
        if(Data[j].Sequence==Identifiers[i].Sequence&&Data[j].Charge==Identifiers[i].Charge){
            struct partOfItem2 pp2;
            pp2.RetentionTime=Data[j].RetentionTime;
            pp2.Scan=Data[j].Scan;
            pp2.PrecursorMZ=Data[j].PrecursorMZ;
            pp2.Coeff=Data[j].Coeff;
            PeptideData.push_back(pp2);
        }
    }
    if(RTWindow){
        /////待填入
    }
    PeptideData.erase(unique(PeptideData.begin(),PeptideData.end()), PeptideData.end());//去除重复元素
    //////
    /*cout<<PeptideData.size()<<endl;
    cout<<PeptideData[34].Scan<<" "<<PeptideData[34].RetentionTime<<" "<<PeptideData[34].Coeff<<" "<<PeptideData[34].PrecursorMZ<<endl;*/
    int count=0;
    int maxScan=PeptideData[0].Scan,minScan=PeptideData[0].Scan;
    for(int j=0;j<PeptideData.size();j++){
        if(PeptideData[j].Coeff<1){
            PeptideData[j].Coeff=0;
        }
        if(PeptideData[j].Coeff>0){
            count++;
        }
        if(PeptideData[j].Scan>maxScan){//找到PeptideData中Scan最大值
            maxScan=PeptideData[j].Scan;
        }
        if(PeptideData[j].Scan<minScan){//找到PeptideData中Scan最小值
            minScan=PeptideData[j].Scan;
        }
    }
    ////PeptidePeaks=NULL这个也是必要的，因为不一定会进入下面的if分支，所以没进去的时候返回null
    ////所以一开始传进来的PeptidePeaks就是null
    //vector<struct partOfHeader> DataOnuniformDrid;
    vector<double> DataOnUniformGridCoeff;
    vector<int> floorPeptideDataPrecursorMZ;//floor(unique(PeptideData$PrecursorMZ))但此处记录的不是unique的
    for(int j=0;j<PeptideData.size();j++){
        floorPeptideDataPrecursorMZ.push_back(int(PeptideData[j].PrecursorMZ));
    }
    sort(floorPeptideDataPrecursorMZ.begin(),floorPeptideDataPrecursorMZ.end());
    floorPeptideDataPrecursorMZ.erase(unique(floorPeptideDataPrecursorMZ.begin(),floorPeptideDataPrecursorMZ.end()), floorPeptideDataPrecursorMZ.end());
    /////
    /*for(int j=0;j<floorPeptideDataPrecursorMZ.size();j++){
        cout<<floorPeptideDataPrecursorMZ[j]<<" ";
    }
    cout<<endl;*/
    bool flag_break=false;
    if(count>3){
        //这其中的操作需要调用包中的函数，如findpeaks和kz
        for(int j=0;j<header.size();j++){
            if(header[j].seqNum>=minScan&&header[j].seqNum<=maxScan){
                if(header[j].msLevel==2){
                    vector<int>::iterator ret;
                    ret = find(floorPeptideDataPrecursorMZ.begin(), floorPeptideDataPrecursorMZ.end(), int(header[j].PrecursorMZ));
                    if(ret!=floorPeptideDataPrecursorMZ.end()){//即header[j].PrecursorMZ是在floorPeptidedataPrecursorMZ之中的
                        /*struct partOfHeader ph;
                        ph.seqNum=header[j].seqNum;
                        ph.retentionTime=header[j].retentionTime;
                        ph.Coeff=0;
                        for(int k=0;k<PeptideData.size();k++){
                            if(header[j].seqNum==PeptideData[k].Scan){
                                ph.Coeff=PeptideData[k].Coeff;
                                break;
                            }
                        }
                        DataOnuniformDrid.push_back(ph);*/
                        flag_break=false;
                        for(int k=0;k<PeptideData.size();k++){
                            if(header[j].seqNum==PeptideData[k].Scan){
                                DataOnUniformGridCoeff.push_back(PeptideData[k].Coeff);
                                flag_break=true;
                            }
                        }
                        if(!flag_break){
                            DataOnUniformGridCoeff.push_back(0);
                        }
                    }
                }
            }
        }
        //////
        /*for(int j=0;j<DataOnUniformGridCoeff.size();j++){
            cout<<DataOnUniformGridCoeff[j]<<" ";
        }
        cout<<endl;*/
        /*for(int j=0;j<DataOnuniformDrid.size();j++){
            DataOnUniformGridCoeff.push_back(DataOnuniformDrid[j].Coeff);
        }*/
        //下面进行寻找序列峰的操作
        findPeaks(DataOnUniformGridCoeff,PeptidePeaks,2,2,10);
        //cout<<PeptidePeaks.size()<<endl;
        if(smooth=="rollmean"){
            /////有关kz函数的操作
            vector<double> kzData(DataOnUniformGridCoeff);//因为后面还需要DataOnUniformGridCoeff的值进行判别，所以这里copy一份
            kz(kzData,FilterWindow,KZiters);
            ///////
            /*for(int j=0;j<DataOnUniformGridCoeff.size();j++){
                cout<<kzData[j]<<" ";
            }
            cout<<endl;*/
            for(int j=0;j<DataOnUniformGridCoeff.size();j++){
                if(DataOnUniformGridCoeff[j]<1){
                    kzData[j]=0;
                }
            }
            ////////
            /*for(int j=0;j<DataOnUniformGridCoeff.size();j++){
                cout<<kzData[j]<<" ";
            }
            cout<<endl;*/
            PeptidePeaks.clear();
            findPeaks(kzData,PeptidePeaks,2,2,10);
            //cout<<PeptidePeaks.size()<<endl;
        }
    }
}
void QuantifyPeptides(vector<double> &result,vector<struct item> Data,vector<struct partOfItem> Identifiers,vector<struct headItem> header,int i,int IntensityCutoff=0,int QuantileCutoff=0,bool RTWindow=true,string smooth="rollmean",int FilterWindow=3,int KZiters=3){
    vector<struct partOfItem1> PeptideData;
    vector<double> PeptideDataCoeff;
    for(int j=0;j<Data.size();j++){
        if(Data[j].Sequence==Identifiers[i].Sequence&&Data[j].Charge==Identifiers[i].Charge){
            struct partOfItem1 pp1;
            pp1.Coeff=Data[j].Coeff;
            pp1.PrecursorMZ=Data[j].PrecursorMZ;
            pp1.Scan=Data[j].Scan;
            pp1.RetentionTime=Data[j].RetentionTime;
            pp1.Correlation=Data[j].Correlation;
            PeptideData.push_back(pp1);
        }
    }
    if(RTWindow){
        //////待填入
    }
    ///////
    /*cout<<PeptideData.size()<<endl;
    cout<<PeptideData[13].Scan<<" "<<PeptideData[13].RetentionTime<<" "<<PeptideData[13].Coeff<<" "<<PeptideData[13].PrecursorMZ<<endl;*/
    PeptideData.erase(unique(PeptideData.begin(),PeptideData.end()), PeptideData.end());//去除重复元素
    //////
    /*cout<<PeptideData.size()<<endl;
    cout<<PeptideData[26].Scan<<" "<<PeptideData[26].RetentionTime<<" "<<PeptideData[26].Coeff<<" "<<PeptideData[26].PrecursorMZ<<endl;*/
    int count=0;
    int maxScan=PeptideData[0].Scan,minScan=PeptideData[0].Scan;
    for(int j=0;j<PeptideData.size();j++){
        if(PeptideData[j].Coeff<1){
            PeptideData[j].Coeff=0;
        }
        //在此处记录PeptideDataCoeff得vector
        PeptideDataCoeff.push_back(PeptideData[j].Coeff);
        if(PeptideData[j].Coeff>1){
            count++;
        }
        if(PeptideData[j].Scan>maxScan){//找到PeptideData中Scan最大值
            maxScan=PeptideData[j].Scan;
        }
        if(PeptideData[j].Scan<minScan){//找到PeptideData中Scan最小值
            minScan=PeptideData[j].Scan;
        }
    }
    vector<struct peak> RawPeaks,SmoothPeaks;
    FindPeaksInData(RawPeaks,Data,Identifiers,i,header,false,IntensityCutoff,QuantileCutoff,RTWindow,"none",FilterWindow);
    FindPeaksInData(SmoothPeaks,Data,Identifiers,i,header,false,IntensityCutoff,QuantileCutoff,RTWindow,smooth,FilterWindow,KZiters);

    double area=0,BoxTestPval=1,BoxTestPvalOnGrid=1,MaxCoeff=0,PeakWidth=0,TimeAtTopOfPeak=0,snr=0,Variance=0,Skewness=0,Kurtosis=0,PeakStart=0,PeakEnd=0,Correlation=0;
    result={area,BoxTestPval,snr,TimeAtTopOfPeak,MaxCoeff,PeakWidth,BoxTestPvalOnGrid,Variance,Skewness,Kurtosis,PeakStart,PeakEnd,Correlation};
    vector<struct partOfHeader> DataOnuniformGrid;
    vector<int> floorPeptideDataPrecursorMZ;
    for(int j=0;j<PeptideData.size();j++){
        floorPeptideDataPrecursorMZ.push_back(int(PeptideData[j].PrecursorMZ));
    }
    sort(floorPeptideDataPrecursorMZ.begin(),floorPeptideDataPrecursorMZ.end());
    floorPeptideDataPrecursorMZ.erase(unique(floorPeptideDataPrecursorMZ.begin(),floorPeptideDataPrecursorMZ.end()), floorPeptideDataPrecursorMZ.end());
    /////
    /*for(int j=0;j<floorPeptideDataPrecursorMZ.size();j++){
        cout<<floorPeptideDataPrecursorMZ[j]<<" ";
    }
    cout<<endl;*/
    if(count>3&&!SmoothPeaks.empty()){
        for(int j=0;j<header.size();j++){
            if(header[j].seqNum>=minScan&&header[j].seqNum<=maxScan){
                if(header[j].msLevel==2){
                    vector<int>::iterator ret;
                    ret = find(floorPeptideDataPrecursorMZ.begin(), floorPeptideDataPrecursorMZ.end(), int(header[j].PrecursorMZ));
                    if(ret!=floorPeptideDataPrecursorMZ.end()){//即header[j].PrecursorMZ是在floorPeptidedataPrecursorMZ之中的
                        struct partOfHeader ph;
                        ph.seqNum=header[j].seqNum;
                        ph.retentionTime=header[j].retentionTime;
                        ph.Coeff=0;
                        for(int k=0;k<PeptideData.size();k++){
                            if(header[j].seqNum==PeptideData[k].Scan){
                                ph.Coeff=PeptideData[k].Coeff;
                                break;
                            }
                        }
                        DataOnuniformGrid.push_back(ph);
                    }
                }
            }
        }
        vector<double> DataOnUniformGridCoeff;
        for(int k=0;k<DataOnuniformGrid.size();k++){
            DataOnUniformGridCoeff.push_back(DataOnuniformGrid[k].Coeff);
        }
        int start=SmoothPeaks[0].left;
        int end=SmoothPeaks[0].right;
        vector<struct partOfHeader> PeakDataOnUniformGrid;
        vector<double> PeakDataOnUniformGridRetentionTime,PeakDataOnUniformGridCoeff;
        for(int k=start-1;k<end;k++){
            PeakDataOnUniformGrid.push_back(DataOnuniformGrid[k]);
            PeakDataOnUniformGridRetentionTime.push_back(DataOnuniformGrid[k].retentionTime);
            PeakDataOnUniformGridCoeff.push_back(DataOnuniformGrid[k].Coeff);
        }
        area=trapz(PeakDataOnUniformGridRetentionTime,PeakDataOnUniformGridCoeff);
        TimeAtTopOfPeak=DataOnuniformGrid[SmoothPeaks[0].top-1].retentionTime;//////记得减去1
        MaxCoeff=DataOnuniformGrid[SmoothPeaks[0].top-1].Coeff;

        //PeakCentralMoments = all.moments(scale(PeakDataOnUniformGrid$Coeff,center=FALSE),order.max = 4,central = TRUE)[3:5]
        vector<double> PeakCentralMoments={0,0,0};
        //////
        /*for(int k=0;k<PeakDataOnUniformGridCoeff.size();k++){
            cout<<PeakDataOnUniformGridCoeff[k]<<" ";
        }
        cout<<endl;*/
        scale(PeakDataOnUniformGridCoeff);
        //////
        /*for(int k=0;k<PeakDataOnUniformGridCoeff.size();k++){
            cout<<PeakDataOnUniformGridCoeff[k]<<" ";
        }
        cout<<endl;*/
        allmomemt(PeakCentralMoments,PeakDataOnUniformGridCoeff);////这里需要进行归一化
        Variance=PeakCentralMoments[0];
        Skewness=PeakCentralMoments[1];
        Kurtosis=PeakCentralMoments[2];

        //string dirPath="/Users/patrickqi/Desktop/Figs_cpp";//这个地方后续还需要更改，需要将dirPath作为参数传入
        BoxTestPval=boxtest(PeptideDataCoeff);
        BoxTestPvalOnGrid=boxtest(DataOnUniformGridCoeff);
        cout<<i<<endl;
        PeakWidth=DataOnuniformGrid[end-1].retentionTime-DataOnuniformGrid[start-1].retentionTime;
        PeakStart=DataOnuniformGrid[start-1].retentionTime;
        PeakEnd=DataOnuniformGrid[end-1].retentionTime;

        DataOnUniformGridCoeff.erase(DataOnUniformGridCoeff.begin()+start-1,DataOnUniformGridCoeff.begin()+end);//为后面的计算snr做好准备，去除该向量中从start到end位置的元素
        //////
        /*for(int k=0;k<DataOnUniformGridCoeff.size();k++){
            cout<<DataOnUniformGridCoeff[k]<<" ";
        }
        cout<<endl;*/
        /*cout<<DataOnUniformGridCoeff.size()<<endl;
        cout<<sd(DataOnUniformGridCoeff)<<endl;*/
        if(sd(DataOnUniformGridCoeff)!=0)
            snr=area/sd(DataOnUniformGridCoeff);
        else
            snr=-1;
        Correlation=0;
        for(int k=0;k<PeptideData.size();k++){
            Correlation+=PeptideData[k].Correlation;
        }
        Correlation/=PeptideData.size();
        result.clear();
        result={area,BoxTestPval,snr,TimeAtTopOfPeak,MaxCoeff,PeakWidth,BoxTestPvalOnGrid,Variance,Skewness,Kurtosis,PeakStart,PeakEnd,Correlation};
        ////result={上面计算的一堆}，返回result
        ////实现的时候将result当作参数，引用传入
    }
}

bool EnoughData(vector<struct item> Data,const struct partOfItem &Ii,bool RTWindow=true){
    vector<struct item> PeptideData;
    for(int j=0;j<Data.size();j++){
        if(Data[j].Sequence==Ii.Sequence&&Data[j].Charge==Ii.Charge){
            PeptideData.push_back(Data[j]);
        }
    }
    if(RTWindow){
        //////
    }
    if(!PeptideData.empty())    return true;
    else    return false;
}

void f(vector<double> &x,int k,bool IsThereEnoughData,vector<struct item> FIGSData,vector<struct partOfItem> IDs,vector<struct headItem>header){
    if(IsThereEnoughData){
        QuantifyPeptides(x,FIGSData,IDs,header,k);
    }
    else{
        x={0,1,0,0,0,0,1,0,0,0,0,0,0};
    }
}
void QuantifyAllFromCoeffs(vector<struct Quants> &quants,vector<struct item> FIGSData,vector<struct headItem>header){//返回的vector还没加进去
    vector<struct partOfItem> IDs;
    bool IsThereEnoughData;
    ////////这里有所改动，如有问题需要注意，这里后续还需要进行进一步的检查，避免别的地方需要使用在此处简化的部分
    vector<vector<double> > Quant;
    vector<double> temp;
    for(int i=0;i<FIGSData.size();i++){
        struct partOfItem pp;
        pp.Sequence=FIGSData[i].Sequence;
        pp.Charge=FIGSData[i].Charge;
        IDs.push_back(pp);
    }
    IDs.erase(unique(IDs.begin(),IDs.end()), IDs.end());//根据已经重载的==号，来将IDs中的相同结构体清除
    /*cout<<IDs[171].Sequence<<" "<<IDs[171].Charge<<endl;
    cout<<IDs.size()<<endl;*/
    for(int i=0;i<IDs.size();i++){////////////i<IDs.size()
        IsThereEnoughData=EnoughData(FIGSData,IDs[i],false);
        f(temp,i,IsThereEnoughData,FIGSData,IDs,header);//IsThereEnoughData和R中有所不同，这里不存所有的bool值进一个数组了，而是运用当前计算所得的bool值
        //但是当前索引i还是要传进f函数的，因为f函数可能要调用Quantifypeptide函数，其中是需要索引的
        Quant.push_back(temp);
    }
    //////待填入关于Quants的元素扩展操作
    //quants为需要返回的vector
    struct Quants q1;
    for(int i=0;i<IDs.size();i++){
        q1.Sequence=IDs[i].Sequence;
        q1.Charge=IDs[i].Charge;
        q1.Quantity=Quant[i][0];
        q1.LBPval=Quant[i][1];
        q1.SNR=Quant[i][2];
        q1.MaxCoeffTime=Quant[i][3];
        q1.MaxCoeff=Quant[i][4];
        q1.PeakWidth=Quant[i][5];
        q1.LBPvalOnUniformGrid=Quant[i][6];
        q1.PeakVariance=Quant[i][7];
        q1.PeakSkewness=Quant[i][8];
        q1.PeakKutosis=Quant[i][9];
        q1.PeakStart=Quant[i][10];
        q1.PeakEnd=Quant[i][11];
        q1.Correlation=Quant[i][12];
        quants.push_back(q1);
    }
    /*for(int i=0;i<10;i++){
        q1.Sequence=IDs[i+10].Sequence;
        q1.Charge=IDs[i+10].Charge;
        q1.Quantity=Quant[i][0];
        q1.LBPval=Quant[i][1];
        q1.SNR=Quant[i][2];
        q1.MaxCoeffTime=Quant[i][3];
        q1.MaxCoeff=Quant[i][4];
        q1.PeakWidth=Quant[i][5];
        q1.LBPvalOnUniformGrid=Quant[i][6];
        q1.PeakVariance=Quant[i][7];
        q1.PeakSkewness=Quant[i][8];
        q1.PeakKutosis=Quant[i][9];
        q1.PeakStart=Quant[i][10];
        q1.PeakEnd=Quant[i][11];
        q1.Correlation=Quant[i][12];
        quants.push_back(q1);
    }*/
    // cout<<quants.size()<<endl;
    /*for(int ii=0;ii<quants.size();ii++){
        cout<<quants[ii].Sequence<<" "<<quants[ii].Charge<<" "<<quants[ii].Quantity<<" "<<quants[ii].LBPval<<" ";
        cout<<quants[ii].SNR<<" "<<quants[ii].MaxCoeffTime<<" "<<quants[ii].MaxCoeff<<" "<<quants[ii].PeakWidth<<" ";
        cout<<quants[ii].LBPvalOnUniformGrid<<" "<<quants[ii].PeakVariance<<" "<<quants[ii].PeakSkewness<<" ";
        cout<<quants[ii].PeakKutosis<<" "<<quants[ii].PeakStart<<" "<<quants[ii].PeakEnd<<" ";
        cout<<quants[ii].Correlation<<endl;
    }*/
}