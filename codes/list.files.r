


FilesDir=file.path("C:/IMARES/FAMA/EU/SC05 NAFO Multispecies Model/Task 3/Subtask 3.2/MSY_ref/gadgetR-FLR_MSY_ref/MSY_truthplusnoise/gadgetR-MSY_ref_tpn_2/")
CombiDir=file.path(FilesDir, "paramfiles")

files=list.files(FilesDir)

ncombi=vector(length=length(files))
for(i in 5:length(files))
{
  filenchar=nchar(files[i])
  combinchar=nchar('results-combination')
  rdsnchar=nchar('.rds')
  
  ncombi[i]=as.integer(as.character(substring(files[i], combinchar+1, (filenchar-rdsnchar))))
  
}

ncombi=sort(ncombi)
ncombi=ncombi[-c(1:4)]

originalcomb=c(1:8000)
diffcombi=setdiff(originalcomb, ncombi)
length(diffcombi)

originalcombi=read.csv(paste(CombiDir, "/effort_combination.csv", sep=""))

newcombi=originalcombi[originalcombi$comb.N %in% diffcombi,]









