#include<iostream>
#include<vector>
#include<string>
#include<cstring>
#include<algorithm>
#include<fstream>
#include<numeric>
#include<cassert>
#include<cmath>
#include<opencv2/opencv.hpp>
#include<algorithm>
using namespace std;
using namespace cv;
#define MIN_NUM -1000000
//vector<bool> IsThereEnoughData;//后期可能需要修改
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

vector<struct headItem> h;
vector<struct item> results;
vector<struct item> resultsWithDecoys;
vector<struct Quants> quants;
vector<struct Quants> quantsWithDecoys;
vector<struct resQuants> D;
//vector<bool> isTypeDecoy;
double irisdata[100000][4];
vector<int> Type;
vector<double> fdr;
vector<string> pasteSequenceCharge;
vector<double> quantsScore;

bool cmppeak(struct peak p1,struct peak p2){
    return p1.maxValue>p2.maxValue;
}
bool cmp(struct item i1,struct item i2){//Scan,Sequence,Charge
    //用于对结构体排序的函数
    if(i1.Scan<i2.Scan){
        return true;
    }
    else if(i1.Scan==i2.Scan){
        if(i1.Sequence.compare(i2.Sequence)==-1){
            return true;
        }
        else if(i1.Sequence.compare(i2.Sequence)==0){
            if(i1.Charge<i2.Charge){
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
double boxtest(vector<double> x,string dirPath){
    //写输入文件
    ofstream outFile;
    string inFilePath=dirPath+"/in.csv";
    string RFilePath=dirPath+"/test.R";
    outFile.open(inFilePath,ios::out);
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
    double stdev=sqrt(accum/a.size());
    return stdev;
}
void scale(vector<double> &a){//标准化
    double stdev=sd(a);
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
                }
            }
            else if(d[r]==d[r-1]){
                if(r-t==1)  t=r;
                else if(r-t<=rightstep){
                    l=r-1; t=r-1;
                    up=false;
                }
                else{
                    struct peak p1;
                    p1.maxValue=d[t];
                    p1.left=l+1;
                    p1.top=t+1;
                    p1.right=r;
                    f.push_back(p1);
                }
            }
        }
    }
    sort(f.begin(),f.end(),cmppeak);
    for(int j=f.size();j>=peaknum;j--){
        f.erase(f.begin()+j-1);
    }
}
void SparkCoeffs(string rawCoeffsFilePath,vector<struct item>& DIA_RefSpectraCoeffs){
    ifstream _csvInput(rawCoeffsFilePath);
    string _Oneline;//用于存每一行的string
    vector<string> eachRow;
    const char* split=",";//分隔符
    char* p;
    //将每一行读出来的数据按逗号分隔，存在eachRow中
    double MaxRetentionTime=MIN_NUM;
    while(getline(_csvInput,_Oneline)){
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
                ii.Coeff=stod(eachRow[1]);
                ii.Scan=stoi(eachRow[0])+1;
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
    vector<struct partOfHeader> DataOnuniformDrid;
    vector<int> floorPeptideDataPrecursorMZ;//floor(unique(PeptideData$PrecursorMZ))但此处记录的不是unique的
    for(int j=0;j<PeptideData.size();j++){
        floorPeptideDataPrecursorMZ.push_back(int(PeptideData[j].PrecursorMZ));
        //floorPeptideDataPrecursorMZ.push_back(PeptideData[j].Scan);
    }
    if(count>3){
        //这其中的操作需要调用包中的函数，如findpeaks和kz
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
                        DataOnuniformDrid.push_back(ph);
                    }
                }
            }
        }
        vector<double> DataOnUniformGridCoeff;
        for(int j=0;j<DataOnuniformDrid.size();j++){
            DataOnUniformGridCoeff.push_back(DataOnuniformDrid[j].Coeff);
        }
        //下面进行寻找序列峰的操作
        findPeaks(DataOnUniformGridCoeff,PeptidePeaks,2,2,10);
        if(smooth=="rollmean"){
            /////有关kz函数的操作
            vector<double> kzData(DataOnUniformGridCoeff);
            kz(kzData,FilterWindow,KZiters);
            for(int j=0;j<DataOnUniformGridCoeff.size();j++){
                if(DataOnUniformGridCoeff[j]<1){
                    kzData[j]=0;
                }
            }
            PeptidePeaks.clear();
            findPeaks(kzData,PeptidePeaks,2,2,10);
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
    PeptideData.erase(unique(PeptideData.begin(),PeptideData.end()), PeptideData.end());//去除重复元素
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
    vector<int> floorPeptideDataPrecursorMZ;//floor(unique(PeptideData$PrecursorMZ))但此处记录的不是unique的
    for(int j=0;j<PeptideData.size();j++){
        floorPeptideDataPrecursorMZ.push_back(int(PeptideData[j].PrecursorMZ));
    }
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
        DataOnUniformGridCoeff.erase(DataOnUniformGridCoeff.begin()+start-1,DataOnUniformGridCoeff.begin()+end-1);//为后面的计算做好准备，去除该向量中从start到end位置的元素
        vector<struct partOfHeader> PeakDataOnUniformGrid;
        vector<double> PeakDataOnUniformGridRetentionTime,PeakDataOnUniformGridCoeff;
        for(int k=start-1;k<end;k++){
            PeakDataOnUniformGrid.push_back(DataOnuniformGrid[k]);
            PeakDataOnUniformGridRetentionTime.push_back(DataOnuniformGrid[k].retentionTime);
            PeakDataOnUniformGridCoeff.push_back(DataOnuniformGrid[k].Coeff);
        }
        area=trapz(PeakDataOnUniformGridRetentionTime,PeakDataOnUniformGridCoeff);
        TimeAtTopOfPeak=DataOnuniformGrid[SmoothPeaks[0].top].retentionTime;
        MaxCoeff=DataOnuniformGrid[SmoothPeaks[0].top].Coeff;

        //PeakCentralMoments = all.moments(scale(PeakDataOnUniformGrid$Coeff,center=FALSE),order.max = 4,central = TRUE)[3:5]
        vector<double> PeakCentralMoments={0,0,0};
        allmomemt(PeakCentralMoments,PeakDataOnUniformGridCoeff);

        string dirPath="/Users/patrickqi/Desktop/Figs_C++";//这个地方后续还需要更改，需要将dirPath作为参数传入
        BoxTestPval=boxtest(PeptideDataCoeff,dirPath);
        BoxTestPvalOnGrid=boxtest(DataOnUniformGridCoeff,dirPath);
        snr=area/sd(DataOnUniformGridCoeff);
        PeakWidth=DataOnuniformGrid[end-1].retentionTime-DataOnuniformGrid[start-1].retentionTime;
        Variance=PeakCentralMoments[0];
        Skewness=PeakCentralMoments[1];
        Kurtosis=PeakCentralMoments[2];
        PeakStart=DataOnuniformGrid[start-1].retentionTime;
        PeakEnd=DataOnuniformGrid[end-1].retentionTime;
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

bool EnoughData(vector<struct item> Data,vector<struct partOfItem> Identifiers,int i,bool RTWindow=true){
    vector<struct item> PeptideData;
    for(int j=0;j<Data.size();j++){
        if(Data[j].Sequence==Identifiers[i].Sequence&&Data[j].Charge==Identifiers[i].Charge){
            PeptideData.push_back(Data[j]);
        }
    }
    if(RTWindow){

    }
    if(PeptideData.size()>0)    return true;
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

void QuantifyAllFromCoeffs(vector<struct Quants> &quants,vector<struct item> FIGSData,vector<struct headItem>header){
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
    for(int i=0;i<IDs.size();i++){
        IsThereEnoughData=EnoughData(FIGSData,IDs,i,false);
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
}

int main(int argc,char** argv){
    string dirPath;
    cout<<"helloworld"<<endl;
    if(argc==2){
        dirPath=argv[1];//只有一个命令行参数，传的是目录路径
    }
    char* str=(char*)dirPath.c_str();
    char* p=strtok(str,"/");
    vector<string> partOfPath;
    while(p!=NULL){
        partOfPath.push_back(p);
        p=strtok(NULL,"/");
    }
    string headerPath="";
    if(!partOfPath.empty()){
        headerPath+=partOfPath[0];
    }
    for(int i=1;i<partOfPath.size()-1;i++){
        headerPath=headerPath+"/"+partOfPath[i];
    }
    headerPath=headerPath+"/header.csv";
    string resultPath=dirPath+"/Coeffs.csv";
    string decoyResultPath=dirPath+"DecoyCoeffs.csv";
    string outputPath=dirPath+"/Quants.csv";

    ifstream _csvinput(headerPath,ios::in);
    if(!_csvInput){
        cout<<"打开文件夹失败"<<endl;
        exit(1);
    }
    string _Oneline;
    vector<string> eachRow;
    struct headItem h1;
    while(getline(_csvinput,_Oneline)){
        str=(char*)_Oneline.c_str();
        p=strtok(str,",");
        while(p!=NULL){
            eachRow.push_back(p);
            p=strtok(NULL,",");
        }
        h1.PrecursorMZ=stod(eachRow[0]);
        h1.retentionTime=stod(eachRow[1])*60;
        h1.seqNum=stoi(eachRow[2])+1;
        h1.msLevel=2;
        h.push_back(h1);
    }

    SparkCoeffs(resultPath,results);
    SparkCoeffs(decoyResultPath,resultsWithDecoys);

    QuantifyAllFromCoeffs(quants,results,h);
    QuantifyAllFromCoeffs(quantsWithDecoys,resultsWithDecoys,h);
    quantsWithDecoys.insert(quantsWithDecoys.end(),quants.begin(),quants.end());
    string quantsWithDecoysSequence;
    for(int i=0;i<quantsWithDecoys.size();i++){
        quantsWithDecoysSequence=quantsWithDecoys[i].Sequence;
        string::size_type position;
        position=quantsWithDecoysSequence.find("DECOY");
        if(position!=quantsWithDecoysSequence.npos){
            //isTypeDecoy.push_back(true);
            quantsWithDecoys[i].isTypeDecoy=true;
        }
        /*else{
            //isTypeDecoy.push_back(false);
        }*/
    }
    bool flag=false;
    for(vector<struct Quants>::iterator it=quantsWithDecoys.begin();it!=quantsWithDecoys.end();){
        struct Quants temp=*it;
        if(temp.Quantity<=0){
            it=quantsWithDecoys.erase(it);//删除该结构体，并将it指向已删除元素的下一个位置
        }
        else{
            if(it->isTypeDecoy) flag=true;
            ++it;
        }
    }
    if(flag){
        for(int i=0;i<quantsWithDecoys.size();i++){
            irisdata[i][0]=quantsWithDecoys[i].PeakVariance;
            irisdata[i][1]=quantsWithDecoys[i].PeakSkewness;
            irisdata[i][2]=quantsWithDecoys[i].PeakKutosis;
            irisdata[i][3]=quantsWithDecoys[i].Correlation;
            if(quantsWithDecoys[i].isTypeDecoy){
                Type.push_back(0);
            }
            else{
                Type.push_back(1);
            }
        }
        Mat mat=Mat(quantsWithDecoys.size(),4,CV_64FC1,irisdata);
        LDA lda=LDA(mat,Type,1);
        Mat eivector=lda.eigenvectors().clone();

        //针对两类分类问题，计算两个数据集的中心
        int classNum=2;
        vector<Mat> classmean(classNum);
        vector<int> setNum(classNum);

        for(int i=0;i<classNum;i++)
        {
            classmean[i]=Mat::zeros(1,mat.cols,mat.type());  //初始化类中均值为0
            setNum[i]=0;  //每一类中的条目数
        }

        Mat instance;
        for(int i=0;i<mat.rows;i++)
        {
            instance=mat.row(i);//获取第i行
            if(Type[i]==0)  //如果标签为0
            {
                add(classmean[0], instance, classmean[0]);  //矩阵相加
                setNum[0]++;  //数量相加
            }
            else if(Type[i]==1)  //对于第1类的处理
            {
                add(classmean[1], instance, classmean[1]);
                setNum[1]++;
            }
            else
            {}
        }
        for(int i=0;i<classNum;i++)   //计算每一类的均值
        {
            classmean[i].convertTo(classmean[i],CV_64FC1,1.0/static_cast<double>(setNum[i]));
        }
        double total_means=(lda.project(classmean[0].at<double>(0)).at<double>(0)*setNum[0]+lda.project(classmean[1].at<double>(0)).at<double>(0)*setNum[1])/(setNum[0]+setNum[1]);
        Mat matEivector=mat*eivector;
        for(int i=0;i<quantsWithDecoys.size();i++){
            struct resQuants rQ;
            rQ.Sequence=quantsWithDecoys[i].Sequence;
            rQ.Charge=quantsWithDecoys[i].Charge;
            rQ.Quantity=quantsWithDecoys[i].Quantity;
            rQ.PeakStart=quantsWithDecoys[i].PeakStart;
            rQ.PeakEnd=quantsWithDecoys[i].PeakEnd;
            rQ.score=matEivector.at<double>(i)-total_means;
            rQ.isTypeDecoy=quantsWithDecoys[i].isTypeDecoy;
            D.push_back(rQ);
        }
        double cutoff=1000;
        for(int i=0;i<D.size();i++){
            int fenzi,fenmu=0;
            for(int j=0;j<D.size();j++){
                if(D[j].score>D[i].score){
                    fenmu++;
                    if(D[j].isTypeDecoy){
                        fenzi++;
                    }
                }
            }
            fdr.push_back(fenzi*1.0/fenmu);
            if(fdr[i]<0.01){
                if(D[i].score<cutoff){
                    cutoff=D[i].score;
                }
            }
        }

        for(vector<struct resQuants>::iterator it=D.begin();it!=D.end();){
            struct resQuants temp=*it;
            if(temp.score<cutoff){
                it=D.erase(it);
            }
            else{
                ++it;
            }
        }
        for(int i=0;i<D.size();i++){
            pasteSequenceCharge.push_back(D[i].Sequence+to_string(D[i].Charge));
        }
        for(vector<struct Quants>::iterator it=quants.begin();it!=quants.end();){
            struct Quants temp1=*it;
            vector<string>::iterator ret;
            ret=find(pasteSequenceCharge.begin(),pasteSequenceCharge.end(),it->Sequence+to_string(it->Charge));
            if(ret==pasteSequenceCharge.end()){
                it=quants.erase(it);
            }
            else{
                quantsScore.push_back(D[ret-pasteSequenceCharge.begin()].score);
                //////上面一句还有待检查
                ++it;
            }
        }
    }
    else{
        for(int i=0;i<quants.size();i++){
            quantsScore.push_back(1);
        }
    }

    ofstream outFile;
    outFile.open(outputPath,ios::out);
    outFile<<"Sequence"<<','<<"Charge"<<','<<"Quantity"<<','<<"Score"<<','<<"PeakStart"<<','<<"PeakEnd"<<endl;
    for(int i=0;i<quants.size();i++){
        outFile<<quants[i].Sequence<<','<<quants[i].Charge<<','<<quants[i].Quantity<<','<<quantsScore[i]<<','<<quants[i].PeakStart<<','<<quants[i].PeakEnd<<endl;
    }
    outFile.close();
    //system("pause");
    cout<<"hello world"<<endl;
    return 0;
}
