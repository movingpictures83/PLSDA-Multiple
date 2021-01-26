
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

require(graphics)
library("caret")
library("spls")
library(mlbench)

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }
  rownames(parameters) <<- parameters[,1];
  print("READING INPUT FILES...");
  training_set <<- read.table(paste(pfix, toString(parameters["training",2]), sep=""), sep = "\t", header =FALSE, stringsAsFactors=FALSE)#, nrow=20000)
  test_set <<- read.table(paste(pfix, toString(parameters["clinical",2]), sep=""), sep="\t", header = TRUE,  stringsAsFactors=FALSE)
  joincol <<- toString(parameters["joinby", 2])
  id <<- toString(parameters["id", 2])
  xcoor <<- toString(parameters["x", 2])
  classcol <<- toString(parameters["classcol", 2])
  print("DONE");
}

run <- function() {
# transpose the matrix
training_set <<- as.data.frame(t(training_set), stringsAsFactors=FALSE)
train_set_size <<- ncol(training_set)
class_index <<- grep(classcol, colnames(test_set))
# provide a label for the first column
colnames(training_set)[1] <<- joincol
# merge test and training sets
merged_set <<- merge(training_set, test_set, by =joincol, stringsAsFactors=FALSE)
# obtain all unique STUDYIDs
studyID <<- unique(as.character(unlist(merged_set[id])))
#studyID <<- unique(merged_set$STUDYID)
}

output <- function(outputfile) {
# All viruses
for(virus in studyID){
  # Pick only the samples for this virus
  virus_set = merged_set[merged_set[,id]==virus,]
  # Get the different times
  #times = unique(virus_set$TIMEHOURS)
  times = unique(as.numeric(unlist(virus_set[xcoor])))
  # For every time
  for(t in times){
    # Expression data.  Column 1 is CEL, remove.  22278 is hardcoded, only works
    # for this case.  Where the expression data stops.
    #vv = virus_set[virus_set[,"TIMEHOURS"]==t, 2:22278]
    vv = virus_set[virus_set[,xcoor]==t, 2:train_set_size]
    t.x = data.matrix(matrix(as.numeric(unlist(vv)),nr=nrow(vv)))
    # Column SYMPTOMATIC_SC2, hardcoded.
    #t.y = as.factor(virus_set[virus_set[,"TIMEHOURS"]==t,22286])
    t.y = as.factor(virus_set[virus_set[,xcoor]==t,train_set_size+class_index])
    my_pls1 = plsda(t.x, t.y)

    res = order(varImp(my_pls1), decreasing = TRUE)[1:1000]
    imp = varImp(my_pls1)[1][[1]][c(res)]


    df = data.frame(res, imp)
    write.csv(df, file = paste(outputfile,virus,"_",t,".csv", sep=""), row.names = FALSE)
  }
}
}
