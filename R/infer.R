#' Title
#'
#' @return
#' @export
#'
#' @examples
infer<-function(rnaseq,insituexpr,spatialposition,clsize=15,k1=10,l1=5,k2=10,l2=5,k3=10,l3=5,count=100,lambda=0.01){

  #使用Seurat包里的函数对RNAseq数据进行处理*****************************************************************************
  seqdata <- Seurat::CreateSeuratObject(counts = rnaseq,project = "atlasInf", min.cells = 3, min.features = 200)
  seqdata <- Seurat::NormalizeData(object = seqdata, normalization.method = "LogNormalize", scale.factor = 10000)
  #seqdata <- Seurat::FindVariableFeatures(object = seqdata,selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seqdata)
  seqdata <- Seurat::ScaleData(seqdata, features = all.genes)

  #single_cell是1297*2002单细胞数据矩阵，用Seurat包处理后的VariableFeature有2000个，84元基因集中有2个不在这2000个之中，所以出现2+2000=2002
  single_cell<-t(seqdata[["RNA"]]@scale.data)
  print("step 1/6 done")

  #构建metacell，并建3个knn网络*****************************************************************************************
  sc_knn<-get_affinity(k1,l1,single_cell)#单细胞的knn图

  lab<-nnit(single_cell,clsize,method="maxd")#执行平衡聚类，聚成87个类，得到类标签向量lab

  metacell<-create_metacell(single_cell,sc_knn,lab)#构建87个metacell

  mc_knn<-get_affinity(k2,l2,metacell)##metacell的knn图

  bulk_knn<-get_affinity(k3,l3,insituexpr)#bulk细胞的knn图
  print("step 2/6 done")

  #将三个knn网络中的bulkcell网络可视化**********************************************************************************
  plot_affinity(bulk_knn)
  print("step 3/6 done")

  #整合三个knn网络******************************************************************************************************
  mut_inf_mat1<-matrix(0,nrow(sc_knn),nrow(mc_knn))#互信息熵矩阵1初始化

  for(i in 1:nrow(mc_knn))
  {
    index<-which(lab==i)
    mut_inf_mat1[index,i]<-compute_mut_inf(single_cell[index,],metacell[i,])
  }

  mut_inf_mat2<-compute_mut_inf(metacell[,colnames(insituexpr)],insituexpr) #互信息熵矩阵2

  for(i in 1:nrow(mut_inf_mat2))#让第二个互信息熵矩阵稀疏一些
  {
    nearest<-order(mut_inf_mat2[i,],decreasing=T)[1:15]
    mut_inf_mat2[i,-nearest]<-0
  }

  graph<-network_integrate(sc_knn,mc_knn,bulk_knn,mut_inf_mat1,mut_inf_mat2)
  print("step 4/6 done")

  #计算扩散核矩阵******************************************************************************************************
  graph_kernel<-compute_graph_kernel(eps=1e-5,graph,count,lambda)#计算扩散核矩阵
  print("step 5/6 done")

  #用加权求和方式推断得到3039*8924维表达矩阵***************************************************************************
  data<-t(metacell)

  atlasinf<-matrix(0,nrow(insituexpr),nrow(rnaseq))#表示通过推断过程得到的3039*8924矩阵
  colnames(atlasinf)<-rownames(data)#列名为基因名
  gene<-vector('integer',nrow(single_cell))#计算过程中的辅助向量
  temp<-vector('integer',nrow(insituexpr))#计算过程中的辅助向量

  for(i in rownames(data))
  {
    gene <- data[i, ]
    temp <-
      as.matrix(graph_kernel)[(nrow(metacell) + 1):nrow(graph), 1:nrow(metacell)] %*% gene #加权求和
    atlasinf[, i] <- temp
  }
  print("step 6/6 done")
  atlasinf
}
