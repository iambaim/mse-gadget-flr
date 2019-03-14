
library(dplyr)

mydir ="/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/red/simulation/"
redcatchkgest_I=read.table(paste0(mydir,"redfish_I.catch.kg", comment.char=";"))
colnames(redcatchkgest_I)=c("year","step","area","fleet","biomass")

redcatchkgest_II=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/red/simulation/redfish_II.catch.kg",
                            comment.char=";")
colnames(redcatchkgest_II)=c("year","step","area","fleet","biomass")

redcatchkgest_III=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/red/simulation/shrimp.catch.kg",
                             comment.char=";")
colnames(redcatchkgest_III)=c("year","step","area","fleet","biomass")

redcatchkgest=rbind(redcatchkgest_I,redcatchkgest_II,redcatchkgest_III)
redcatchkgest=ddply(redcatchkgest, .(year),summarize, biomass=sum(biomass))
catchmse=as.vector(gadcapOut$red$stk@catch)
year=c(1988:2016)

## plotting 
plot(year, catchmse, type="l")
lines(year, redcatchkgest$biomass, col="red")



codcatchkgest_I=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/cod/simulation/trawl_I.catch.kg",comment.char=";")
colnames(codcatchkgest_I)=c("year","step","area","fleet","biomass")

codcatchkgest_II=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/cod/simulation/trawl_II.catch.kg",comment.char=";")
colnames(codcatchkgest_II)=c("year","step","area","fleet","biomass")

codcatchkgest_III=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/cod/simulation/gil.catch.kg",comment.char=";")
colnames(codcatchkgest_III)=c("year","step","area","fleet","biomass")

codcatchkgest_IV=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/cod/simulation/long.catch.kg",comment.char=";")
colnames(codcatchkgest_IV)=c("year","step","area","fleet","biomass")

codcatchkgest=rbind(codcatchkgest_I,codcatchkgest_II,codcatchkgest_III,codcatchkgest_IV)

codcatchkgest=ddply(codcatchkgest, .(year), summarize, biomass=sum(biomass))

catchmse=as.vector(gadcapOut$cod$stk@catch)
year=c(1988:2016)

plot(year,catchmse, type="l")
lines(year,codcatchkgest$biomass, col="red")





shrimpcatchkgest=read.table("/home/user/repos/MSE@JRC/2018-flemishcapGadget-mse/models/SC05_107/output/shrimp/simulation/shrimp.catch.kg",comment.char=";")
colnames(shrimpcatchkgest)=c("year","step","area","fleet","biomass")
catchmse=as.vector(gadcapOut$shrimp$stk@catch)
year=c(1988:2016)

plot(year,catchmse, type="l")
lines(year,shrimpcatchkgest$biomass, col="red")


