#include<iostream>
#include<vector>
#include<string>
#include<string.h>
#include<algorithm>
#include<fstream>
#include<numeric>
#include<math.h>
#include<algorithm>
using namespace std;
#define MIN -1000000
vector<bool> IsThereEnoughData;//后期可能需要修改
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
    bool operator==(const struct partOfItem& p1){
        return (Sequence==p1.Sequence&&Charge==p1.Charge);
    }
};
struct partOfItem1{
    int Scan;
    double RetentionTime;
    double Coeff;
    double PrecursorMZ;
    double Correlation;
    bool operator==(const struct partOfItem1& p1){
        return (RetentionTime==p1.RetentionTime);
    }
};
struct partOfItem2{
    int Scan;
    double RetentionTime;
    double Coeff;
    double PrecursorMZ;
    bool operator==(const struct partOfItem2& p2){
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

struct peak{
    double maxValue;
    int left;
    int top;
    int right;
};
bool cmppeak(struct peak p1,struct peak p2){
    return p1.maxValue>p2.maxValue;
}
bool cmp(struct item i1,struct item i2){//Scan,Sequence,Charge
    //用于对结构体排序的函数
    if(i1.Scan<i2.Scan){
        return true;
    }
    else if(i1.Scan==i2.Scan){
        if(i1.Sequence<i2.Sequence){
            return true;
        }
        else if(i1.Sequence==i2.Sequence){
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
double scale(vector<double> &a){//标准化
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
void SparkCoeffs(string rawCoeffsFilePath,vector<struct item> DIA_RefSpectraCoeffs){
    ifstream _csvInput(rawCoeffsFilePath);
    string _Oneline;//用于存每一行的string
    vector<string> eachRow;
    const char* split=",";//分隔符
    char* p;
    //将每一行读出来的数据按逗号分隔，存在eachRow中
    double MaxRetentionTime=MIN;
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
    sort(eachRow.begin(),eachRow.end(),cmp);//对eachRow中的结构体按照cmp函数进行排序
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
    /////PeptidePeaks=NULL这个也是必要的，因为不一定会进入下面的if分支，所以没进去的时候返回null
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
        //vector<struct peak> PeptidePeaks;
        findPeaks(DataOnUniformGridCoeff,PeptidePeaks,2,2,10);
        if(smooth=="rollmean"){
            /////有关kz函数的操作
        }
    }
}

void QuantifyPeptides(vector<struct item> Data,vector<struct partOfItem> Identifiers,vector<struct headItem> header,int i,int IntensityCutoff=0,int QuantileCutoff=0,bool RTWindow=true,string smooth="rollmean",int FilterWindow=3,int KZiters=3){
    vector<struct partOfItem1> PeptideData;
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

    vector<double> result={0,1,0,0,0,0,1,0,0,0,0,0,0};
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
        double area=trapz(PeakDataOnUniformGridRetentionTime,PeakDataOnUniformGridCoeff);
        double TimeAtTopOfPeak=DataOnuniformGrid[SmoothPeaks[0].top].retentionTime;
        double MaxCoeff=DataOnuniformGrid[SmoothPeaks[0].top].Coeff;
        //PeakCentralMoments = all.moments(scale(PeakDataOnUniformGrid$Coeff,center=FALSE),order.max = 4,central = TRUE)[3:5]
        /*
         BoxTest = Box.test(PeptideData$Coeff,type="Ljung-Box")
         BoxTestPval = round(BoxTest$p.value,digits=5)

         BoxTestOnGrid = Box.test(DataOnUniformGrid$Coeff,type="Ljung-Box")
         BoxTestPvalOnGrid = round(BoxTestOnGrid$p.value,digits=5)
         */
        double snr=area/sd(DataOnUniformGridCoeff);
        double peakWidth=DataOnuniformGrid[end-1].retentionTime-DataOnuniformGrid[start-1].retentionTime;
        double Correlation=0;
        for(int k=0;k<PeptideData.size();k++){
            Correlation+=PeptideData[k].Correlation;
        }
        Correlation/=PeptideData.size();
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
    /*if (RTWindow)
        PeptideData = subset(PeptideData, abs(RetentionTime - RefRetentionTime) < 300)
        */
    if(PeptideData.size()>0)    return true;
    else    return false;
}

vector<int> f(int k,vector<bool> IsThereEnoughData2){
    if(IsThereEnoughData2[k]){
        ////////待填入
    }
    else{
        return {0,1,0,0,0,0,1,0,0,0,0,0,0};
    }
}

vector<vector<int> > QuantifyAllFromCoeffs(vector<struct item> FIGSData,vector<bool>& IsThereEnoughData1){
    vector<struct partOfItem> IDs;
    vector<vector<int> > Quant;
    for(int i=0;i<FIGSData.size();i++){
        struct partOfItem pp;
        pp.Sequence=FIGSData[i].Sequence;
        pp.Charge=FIGSData[i].Charge;
        IDs.push_back(pp);
    }
    IDs.erase(unique(IDs.begin(),IDs.end()), IDs.end());//根据已经重载的==号，来将IDs中的相同结构体清除
    for(int i=0;i<IDs.size();i++){
        IsThereEnoughData1.push_back(EnoughData(FIGSData,IDs,i,false));
        Quant.push_back(f(i,IsThereEnoughData1));
    }
    //////待填入关于Quants的元素扩展操作
    ////
}

int main(int argc,char** argv){
    vector<struct peak> p;
    vector<struct item> d;
    vector<struct partOfItem> I;
    int i=1;
    vector<headItem> header;
    FindPeaksInData(p,d,I,i,header);
    system("pause");
    return 0;
}
