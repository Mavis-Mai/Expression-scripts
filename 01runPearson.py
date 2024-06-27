import pandas as pd
import sys
import os
from scipy.stats import pearsonr
from scipy.stats import spearmanr


class countPearson():
    def fpkmtblPreTreatment(self,fpkmfile):
        #表达量处理
        fpkmdataset=pd.read_csv(fpkmfile,sep="\t",header=0,index_col=0)
        #fpkmdataset=fpkmdataset.loc[:,fpkmdataset.sum()>10]
        #fpkmdataset = fpkmdataset.loc[:, (fpkmdataset.sum()>10).sort_index()]
        fpkmdataset=fpkmdataset.T
        return fpkmdataset

    #读取TF
    def getSampleGroupt(self,TFgroups):
        Samplegroup=pd.read_csv(TFgroups,sep="\t",names=["sample","group"],index_col=False)
        return Samplegroup

    #获得基因对：
    def getGenepairs(self,list1,list2):
        genepairsList=[]
        for i in list1:
            for j in list2:
                if i !=j:
                    genepairsList.append((i,j))
        genepairsList=set(genepairsList)
        #print(genepairsList[:3])
        return genepairsList

    def checkDupRow(self,id,dataset):
        tmpdata=pd.DataFrame(dataset[id])
        try:
            if len(tmpdata.columns)>1:
                idflat=dataset[id].iloc[:,0].values.flatten()
            else:
                idflat=dataset[id].values.flatten()
        except AttributeError:
            print(id)
        return idflat


    #相关系数计算
    def corCounter(self,fpkmdataset,TFfpkmdataset,genepairsList,Samplegroup,output):
        x,y,yS,value,valueS,ave1,ave2,typpe=[],[],[],[],[],[],[],[]
        for genepairs in genepairsList:
            (idi,idj)=genepairs
            try:
                namei=Samplegroup.loc[Samplegroup["sample"]==idi]["group"][0]
            except:
                namei="Null"
            idiflat=self.checkDupRow(idi,TFfpkmdataset)
            idjflat=self.checkDupRow(idj,fpkmdataset)
            try:
                corr,p_value=pearsonr(idiflat,idjflat)
                corrS,p_valueS=spearmanr(idiflat,idjflat)
            except ValueError:
                print(idj)
                print(idjflat)

            ave1.append(round(pd.DataFrame(TFfpkmdataset[idi]).iloc[:,0].mean(),2))
            ave2.append(round(pd.DataFrame(fpkmdataset[idj]).iloc[:,0].mean(),2))

        # ave1.append(round(fpkmdataset[idi].iloc[:,0].mean(),2))
        # ave2.append(round(fpkmdataset[idj].iloc[:,0].mean(),2))
            try:
                x.append(idi+"-"+idj)
                y.append(corr)
                yS.append(corrS)
                value.append(str(p_value))
                valueS.append(str(p_valueS))
                typpe.append(namei)
            except TypeError:
                pass
    #AT1G75250ai测试
        #newdata=pd.DataFrame({"genepairs":x,
         #   "pearson_corr":y,"pearson_pvalue":value,
         #   "Spearman_corr":yS,"Spearman_pvalue":valueS,
         #   "ave1":ave1,"aver2":ave2,"pairtypes":typpe})
        #多样品的情况下，感觉表达量均值用途不大。干脆去掉
        newdata=pd.DataFrame({"genepairs":x,
            "pearson_corr":y,"pearson_pvalue":value,
            "Spearman_corr":yS,"Spearman_pvalue":valueS,
            "pairtypes":typpe})

        newdata["pearson_corr"]=newdata["pearson_corr"].astype("float")
        newdata["Spearman_corr"]=newdata["Spearman_corr"].astype("float")
        newdata.to_csv(output,sep="\t",quoting=False,index=False,header=True)

if __name__=="__main__":
    path=os.path.split(TFgroups)[0]
    fpkmfileName=os.path.split(fpkmfile)[1]
    if path=="":
        path="."
    output=path+"/"+fpkmfileName+"_Cor.txt"


    mycountPearson=countPearson()
    Samplegroup=mycountPearson.getSampleGroupt(TFgroups)
    fpkmdataset=mycountPearson.fpkmtblPreTreatment(fpkmfile)
    TFfpkmdatset=mycountPearson.fpkmtblPreTreatment(tffpkmfile)
    print(TFfpkmdatset.head())

    genepairsList=mycountPearson.getGenepairs(Samplegroup["sample"],fpkmdataset.columns)
    mycountPearson.corCounter(fpkmdataset,TFfpkmdatset,genepairsList,Samplegroup,output)

