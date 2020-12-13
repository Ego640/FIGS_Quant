library("readr")
args=commandArgs(TRUE)
dirPath=args[1]
inputFilePath = file.path(dirPath, "in.csv")
x=read.csv(inputFilePath)
names(x)=c("input")
x$input=as.numeric(x$input)
x=x$input
BoxTest=Box.test(x,type="Ljung-Box")
BoxTestPval=round(BoxTest$p.value,digits=5)
r<-as.character(BoxTestPval)
outputFilePath=file.path(dirPath,"res.txt")
write_file(r,outputFilePath)
