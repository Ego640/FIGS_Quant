#include <iostream>
#include"cpptest_function.h"
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;

vector<headItem> h;
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

int main(int argc,char** argv) {
    //确定输入输出文件的路径
    string dirPath;
    if(argc==2){
        dirPath=argv[1];//只有一个命令行参数，传的是目录路径
    }
    string resultPath=dirPath+"/Coeffs.csv";
    string decoyResultPath=dirPath+"/DecoyCoeffs.csv";
    string outputPath=dirPath+"/Quants.csv";
    char* str=(char*)dirPath.c_str();
    char* p=strtok(str,"/");
    vector<string> partOfPath;
    while(p!=NULL){
        partOfPath.push_back(p);
        p=strtok(NULL,"/");
    }
    string headerPath="/";
    if(!partOfPath.empty()){
        headerPath+=partOfPath[0];
    }
    for(int i=1;i<partOfPath.size()-1;i++){
        headerPath=headerPath+"/"+partOfPath[i];
    }
    headerPath=headerPath+"/header.csv";

    //读取header.csv文件内容
    /////我们后面再回来考虑数据精度的问题
    ifstream _csvinput(headerPath,ios::in);
    string _Oneline;
    vector<string> eachRow;
    while(getline(_csvinput,_Oneline)){
        eachRow.clear();//因为我们在循环中重复利用eachRow这个vector，所以我们在循环的开始要将其清空
        struct headItem h1;
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
    quantsWithDecoys.insert(quantsWithDecoys.begin(),quants.begin(),quants.end());
    /*cout<<quantsWithDecoys.size()<<endl;
    for(int i=0;i<quantsWithDecoys.size();i++){
        cout<<quantsWithDecoys[i].Sequence<<" "<<quantsWithDecoys[i].Quantity<<endl;
    }*/
    string quantsWithDecoysSequence;
    for(int i=0;i<quantsWithDecoys.size();i++){
        quantsWithDecoysSequence=quantsWithDecoys[i].Sequence;
        string::size_type position;
        position=quantsWithDecoysSequence.find("DECOY");
        if(position!=quantsWithDecoysSequence.npos){
            //isTypeDecoy.push_back(true);
            quantsWithDecoys[i].isTypeDecoy=true;
        }
        else{
            //isTypeDecoy.push_back(false);
            quantsWithDecoys[i].isTypeDecoy=false;
        }
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
    /*ofstream outTestFile;
    string qdFile="/Users/patrickqi/Desktop/qd.csv";
    outTestFile.open(qdFile,ios::out);
    outTestFile<<"Sequence"<<','<<"Charge"<<','<<"PeakVariance"<<','<<"PeakSkewness"<<','<<"PeakKutosis"<<','<<"Correlation"<<','<<"Type"<<endl;
    for(int i=0;i<quantsWithDecoys.size();i++){
        outTestFile<<quantsWithDecoys[i].Sequence<<','<<quantsWithDecoys[i].Charge<<','<<quantsWithDecoys[i].PeakVariance<<','<<quantsWithDecoys[i].PeakSkewness<<','<<quantsWithDecoys[i].PeakKutosis<<','<<quantsWithDecoys[i].Correlation<<','<<quantsWithDecoys[i].isTypeDecoy<<endl;
    }
    outTestFile.close();*/
    /*cout<<quantsWithDecoys.size()<<endl;
    for(int i=0;i<quantsWithDecoys.size();i++){
        cout<<quantsWithDecoys[i].Sequence<<" "<<quantsWithDecoys[i].Quantity<<endl;
    }*/
    if(flag) cout<<"flag is true"<<endl;
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

        Mat mat=Mat(Type.size(),4,CV_64FC1,irisdata);
        LDA lda=LDA(mat,Type,1);
        Mat eivector=lda.eigenvectors();
        eivector=eivector*9.7237093;
        cout<<eivector<<endl;

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
        double total_means=(lda.project(classmean[0]).at<double>(0)*setNum[0]+lda.project(classmean[1]).at<double>(0)*setNum[1])/(1.0*(setNum[0]+setNum[1]));

        Mat matEivector=mat*eivector;

        //针对两类分类问题，计算两个数据集的中心
        /*int classNum=2;
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
        Mat matEivector=mat*eivector;*/
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
        /////
        /*cout<<D.size()<<endl;
        for(int i=0;i<D.size();i++){
            cout<<D[i].Sequence<<" "<<D[i].Charge<<" "<<D[i].Quantity<<" "<<D[i].PeakStart<<" "<<D[i].PeakEnd<<" "<<D[i].score<<" "<<D[i].isTypeDecoy<<endl;
        }*/
        double cutoff=1000;
        for(int i=0;i<D.size();i++){
            int fenzi=0,fenmu=0;
            for(int j=0;j<D.size();j++){
                if(D[j].score>D[i].score){
                    fenmu++;
                    if(D[j].isTypeDecoy){
                        fenzi++;
                    }
                }
            }
            fdr.push_back((1.0*fenzi)/(1.0*fenmu));
            if(fdr[i]<0.01){
                if(D[i].score<cutoff){
                    cutoff=D[i].score;
                }
            }
        }

        /*for(int i=0;i<fdr.size();i++){
            cout<<fdr[i]<<" ";
        }
        cout<<endl;
        cout<<cutoff<<endl;*/
        for(vector<struct resQuants>::iterator it=D.begin();it!=D.end();){
            struct resQuants temp=*it;
            if(temp.score<cutoff){
                it=D.erase(it);
            }
            else{
                ++it;
            }
        }
        /////
        /*cout<<D.size()<<endl;
        for(int i=0;i<D.size();i++){
            cout<<D[i].Sequence<<" "<<D[i].Charge<<" "<<D[i].Quantity<<" "<<D[i].PeakStart<<" "<<D[i].PeakEnd<<" "<<D[i].score<<" "<<D[i].isTypeDecoy<<endl;
        }*/
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
    /*cout<<quants.size()<<" "<<quantsScore.size()<<endl;
    for(int i=0;i<quants.size();i++){
        cout<<quants[i].Sequence<<" "<<quants[i].Charge<<endl;
    }
    for(int i=0;i<quantsScore.size();i++){
        cout<<quantsScore[i]<<" ";
    }
    cout<<endl;*/
    ofstream outFile;
    outFile.open(outputPath,ios::out);
    outFile<<"Sequence"<<','<<"Charge"<<','<<"Quantity"<<','<<"Score"<<','<<"PeakStart"<<','<<"PeakEnd"<<endl;
    for(int i=0;i<quants.size();i++){
        outFile<<quants[i].Sequence<<','<<quants[i].Charge<<','<<quants[i].Quantity<<','<<quantsScore[i]<<','<<quants[i].PeakStart<<','<<quants[i].PeakEnd<<endl;
    }
    outFile.close();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
