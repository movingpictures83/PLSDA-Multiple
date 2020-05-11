require(graphics)
library("caret")
library("spls")
library(mlbench)

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  print("READING INPUT FILES...");
  t1 <<- read.table(toString(parameters["training",2]), sep = "\t", header =FALSE, stringsAsFactors=FALSE)#, nrow=20000)
  t2 <<- read.table(toString(parameters["clinical",2]), sep="\t", header = TRUE,  stringsAsFactors=FALSE)
  print("DONE");
  #prefix <<- toString(parameters["prefix", 2]);
#t1 <- read.table("ViralChallenge_training_EXPRESSION_RMA.tsv", sep = "\t", header =FALSE, stringsAsFactors=FALSE)
#t2 <- read.table("ViralChallenge_training_CLINICAL.tsv", sep="\t", header = TRUE,  stringsAsFactors=FALSE)


#t1 <- as.data.frame(t(t1), stringsAsFactors=FALSE)
#colnames(t1)[1] = "CEL"
# rownames(t1) <- substring(rownames(t1), 2, length(rownames(t1)))
# write.table(t1, file = "shortRMA.csv", sep = ",", row.names = FALSE, col.names = TRUE)
#x <- as.data.frame(merge(t1, t2, by ="CEL", stringsAsFactors=FALSE))

#studyID = unique(x$STUDYID)
}

run <- function() {

#rm(list=ls(all=TRUE))
#setwd("appData/Res_Challenge/Data")
#t1 = read.table("ViralChallenge_training_EXPRESSION_RMA.tsv", sep = "\t", header = FALSE, stringsAsFactors=FALSE)
#t2 <- read.table("ViralChallenge_training_CLINICAL.tsv", sep="\t", header = TRUE,  stringsAsFactors=FALSE)

t2[t2$STUDYID == "DEE4X H1N1",]$STUDYID="H1N1"
t2[t2$STUDYID == "DEE3 H1N1",]$STUDYID="H1N1"
t2[t2$STUDYID == "DEE2 H3N2",]$STUDYID="H3N2"
t2[t2$STUDYID == "DEE5 H3N2",]$STUDYID="H3N2"
t2[t2$STUDYID == "Rhinovirus Duke",]$STUDYID="Rhinovirus"
t2[t2$STUDYID == "Rhinovirus UVA",]$STUDYID="Rhinovirus"

t1 <<- as.data.frame(t(t1), stringsAsFactors=FALSE)
colnames(t1)[1] <<- "CEL"
# rownames(t1) <- substring(rownames(t1), 2, length(rownames(t1)))
# write.table(t1, file = "shortRMA.csv", sep = ",", row.names = FALSE, col.names = TRUE)
x <<- merge(t1, t2, by ="CEL", stringsAsFactors=FALSE)

studyID <<- unique(x$STUDYID)
}

output <- function(outputfile) {
for(virus in studyID){

  v1 = x[x[,"STUDYID"]==virus,]
  times = unique(v1$TIMEHOURS)
  for(t in times){
    t.x = data.matrix(v1[v1[,"TIMEHOURS"]==t,2:22278])
    t.y = as.factor(v1[v1[,"TIMEHOURS"]==t,22286])

    # t.x = data.matrix(v1[v1[,"TIMEHOURS"]==t,2:100])
    # t.y = as.factor(v1[v1[,"TIMEHOURS"]==t,108])
    my_pls1 = plsda(t.x, t.y)

    res = order(varImp(my_pls1), decreasing = TRUE)[1:1000]

    # resCon = sapply(res, function(x) varImp(my_pls1)[1][[1]][x] > 0.01)
    # fin = res[resCon]
    imp = varImp(my_pls1)[1][[1]][c(res)]


    df = data.frame(res, imp)
    write.csv(df, file = paste(outputfile,virus,"_",t,".csv", sep=""), row.names = FALSE)
    #hist(varImp(my_pls1)[1][[1]][c(res)])
    # break;
  }
  # break;
}
}
