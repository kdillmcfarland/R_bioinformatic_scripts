#Save list of currently installed packages
tmp <- installed.packages()
installedpkgs <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpkgs, file="~/Desktop/installed_old.rda")

#Install R as usual
#Restart RStudio

#Re-install packages
##Load prev list
load("~/Desktop/installed_old.rda")

##Bioconductor
tmp <- installed.packages()
installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
missing <- setdiff(installedpkgs, installedpkgs.new)
install.packages(c("testthat","BiocManager")) #causes BiocManager to fail if not installed first

for (i in 1:length(missing)) {
  tryCatch({ BiocManager::install(missing[i], update=TRUE, ask=FALSE)
  }, error=function(e){})
}

##CRAN
tmp <- installed.packages()
installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
missing <- setdiff(installedpkgs, installedpkgs.new)

for (i in 1:length(missing)) {
  tryCatch({ install.packages(missing[i], dependencies=TRUE, quiet=TRUE, Ncpus = 4)
  }, error=function(e){})
}

update.packages()

#GitHub
devtools::install_github("BIGslu/BIGverse", upgrade = "always")

#Check
tmp <- installed.packages()
installedpkgs.new <- as.vector(tmp[is.na(tmp[,"Priority"]), 1])
setdiff(installedpkgs, installedpkgs.new)