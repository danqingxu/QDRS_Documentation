load("/Volumes/Seagate Bac/Multiple Phenotypes/Package/QDRS_Analysis/data/ckd_phedat_test_set1.rda")
subcol = c("Race","Group","G_Staging","Min_Age","Max_Age","code_585.3","code_433.1","code_172.1")
example_scaled = head(ckd_phedat_test[!is.na(ckd_phedat_test$NKF_Stage),subcol])
example_unscaled = head(ckd_phedat_test_unscaled[!is.na(ckd_phedat_test$NKF_Stage),subcol])
rownames(example_scaled) = NULL
rownames(example_unscaled) = NULL
save(example_scaled,example_unscaled,file = "data/example_feature_matrix.rda")


library(readxl)
Case_defining_phecodes <- read_excel("data/Case_defining_phecodes.xlsx",
col_types = c("text", "text", "text",
"text"))
View(Case_defining_phecodes)
Case_defining_phecodes$Phecode[1:4] <- c("585.3","585.31","585.32","585.33")
Dementia_label_phecodes
Case_defining_phecodes$Phecode[c(21,25,27)] <- c("290.1","290.16","317.1")
save(Case_defining_phecodes,file="data/Case_defining_phecodes.rda")

load("/Volumes/Seagate Bac/Multiple Phenotypes/common/phecode_prevalence_eMERGE.rda")
save(prevalence.ckd,file="data/prevalence.ckd.rda")

load("/Volumes/Seagate Bac/Multiple Phenotypes/Package/QDRS_Analysis/data/ckd_phedat_training.rda")
ckd_phedat_training = ckd_phedat_training[,-1]
ckd_phedat_training_unscaled = ckd_phedat_training_unscaled[,-1]
save(ckd_phedat_training, ckd_phedat_training_unscaled, file = "/Volumes/Seagate Bac/Multiple Phenotypes/Package/QDRS_Analysis/data/ckd_phedat_training.rda")

load("/Volumes/Seagate Bac/Multiple Phenotypes/Package/QDRS_Analysis/data/ckd_phedat_test.rda")
ckd_phedat_test = ckd_phedat_test[,-1]
ckd_phedat_test_unscaled = ckd_phedat_test_unscaled[,-1]
save(ckd_phedat_test, ckd_phedat_test_unscaled, file = "/Volumes/Seagate Bac/Multiple Phenotypes/Package/QDRS_Analysis/data/ckd_phedat_test.rda")

fname.phenos = c("ckd","cad","T2D","hf","Dementia","gerd")
name.phenos = c("CKD","CAD","T2D","HF","Dementia","GERD")
num.phecodes = c(110,95,139,95,154,165)
add = c(200,1200)
Wilcoxon = NULL
method = c("PheRS","Eigen","PC1","PC2","LPC")
for (k in 1:6){
  wilcox.weight = matrix(rep(NA,length(method)*2),nc=length(method))
  print(fname.phenos[k])
  for(i in 1:2){
    print(i)
    n.add = add[i]
    load(paste0("/Volumes/Seagate Bac/Multiple Phenotypes/",fname.phenos[k],"/added_phecodes_",num.phecodes[k],"/weights_",num.phecodes[k],"phecodes_",n.add,".rda"))
    for (j in 1:length(method)){
      m = method[j]
      if (k==3&i==2&j==2){
        wilcox.weight[i,j] = NA
      } else {
        wilcox.weight[i,j] = wilcox.test(abs(Wt[Wt.adj[,paste0(fname.phenos[k],"_related")],m]),abs(Wt[!Wt.adj[,paste0(fname.phenos[k],"_related")],m]),alt="greater")$p.value
      }
    }
  }
  colnames(wilcox.weight) = method
  Wilcoxon = rbind(Wilcoxon, wilcox.weight)
}
Wilcoxon1 = formatC(Wilcoxon, digits = 2, format = "E")
wilcox.added = data.frame(Phenotype = rep(name.phenos,rep(2,6)), added = rep(c(200,1200),6),Wilcoxon1)
save(wilcox.added,file="data/wilcox_added_phecodes.rda")


load("/Volumes/Seagate Bac/Multiple Phenotypes/ckd/104_phecodes_only/ckd_phedat_training_gllvm.rda")
ckd_phedat_training = ckd_phedat_training[,-1]
ckd_phedat_training_unscaled = ckd_phedat_training_unscaled[,-1]
save(ckd_phedat_training, ckd_phedat_training_unscaled, file = "/Volumes/Seagate Bac/Multiple Phenotypes/Package/QDRS_Analysis/data/ckd_phedat_training_lvs.rda")


#=============== Debug Purpose ================

load("data/ckd_phedat_training.rda")
dim(ckd_phedat_training)
#[1] 25231  1874
dim(ckd_phedat_training_unscaled)
#[1] 25231  1874
load("data/ckd_phedat_test.rda")
dim(ckd_phedat_test)
#[1] 25232  1874
dim(ckd_phedat_test_unscaled)
#[1] 25232  1874

# load the list of relevant phecodes including the case defining phecodes
load("data/ckd_phecodes_list.rda")
length(ckd_phecodes)
#[1] 104
length(ckd_label_phecodes)
#[1] 6

load("data/prevalence.ckd.rda")
PheRS.res = PheRS(X = ckd_phedat_training_unscaled[,ckd_phecodes],
                  feature.prevalence = as.numeric(prevalence.ckd[,ckd_phecodes]))
PheRS.test = predictQDRS(Y = ckd_phedat_test_unscaled[,ckd_phecodes], 
                         weights = PheRS.res$weights)
ckd.QDRS.test = data.frame(Group = ckd_phedat_test$Group, 
                           G_Staging = ckd_phedat_test$G_Staging, 
                           PheRS = as.vector(PheRS.test))

phecode_m = ckd_phedat_training[,ckd_phecodes]
# this is the scaled training set (pre-processed the union of training and test set)
# so set training to all row numbers
# then the function returns weights and scores 
Eigen.res = eigen.score(X = phecode_m, 
                        training = 1:nrow(phecode_m), 
                        scale = FALSE)
Eigen.test = predictQDRS(Y = ckd_phedat_test[,ckd_phecodes], 
                         weights = Eigen.res$weights)
ckd.QDRS.test = data.frame(ckd.QDRS.test, 
                           Eigen = as.vector(Eigen.test))

PC.res = PC(X = phecode_m, 
            group = ckd_phedat_training$Group, 
            training = 1:nrow(phecode_m), 
            scale = FALSE, 
            pc.num = 1:2)
PC12.test = as.data.frame(predictQDRS(Y = ckd_phedat_test[,ckd_phecodes], 
                                      weights = PC.res$weights))
ckd.QDRS.test = data.frame(ckd.QDRS.test, 
                           PC12.test)

LPC.res = LPC(X = phecode_m, 
              group = ckd_phedat_training$Group,
              training = 1:nrow(phecode_m), 
              scale = FALSE)
LPC.test = as.data.frame(predictQDRS(Y = ckd_phedat_test[,ckd_phecodes], 
                                     weights = LPC.res$weights))
ckd.QDRS.test = data.frame(ckd.QDRS.test, 
                           LPC.test)
colnames(ckd.QDRS.test)[3:7] = c("PheRS", "Eigen", 
                                 "PC1","PC2","LPC")

#============ Debug Purpose LVS ============
# Issue with gllvm package
# Open issue with the unknown error
# Error in .Call("FreeADFunObject", ptr, PACKAGE = DLL) : "FreeADFunObject" not available for .Call() for package "gllvm"
# LVS
load("data/ckd_phedat_training_lvs.rda")
dim(ckd_phedat_training_unscaled)
#[1] 8410 1874
# this computation may take a long time 
family = "binomial";p.seed=3080;starting.choice="random"
X = ckd_phedat_training_unscaled[1:500,ckd_phecodes]
rowsum_u = rowSums(X)
colsum_u = colSums(X)
Xc = X[rowsum_u!=0,colsum_u!=0]
fit.VA <- gllvm(Xc, family = "binomial"("probit"), method = "VA", num.lv = 1, starting.val = starting.choice, seed = p.seed)
