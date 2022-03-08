# No Remotes ----
# Attachments ----
to_install <- c("MASS", "numDeriv")
no.packages = length(to_install)

for (i in 1:no.packages) {
    message(paste("looking for ", to_install[i]))
    message(paste("     installing", to_install[i] ))
    install.packages(to_install[i])
}


