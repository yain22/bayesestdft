# No Remotes ----
# Attachments ----
to_install <- c("mvtnorm", "truncdist", "nleqslv", "ggplot2", "compiler", "scales", "ggthemes", "svMisc")
for (i in to_install) {
    message(paste("looking for ", i))
    if (!requireNamespace(i)) {
        message(paste("     installing", i))
        install.packages(i)
    }
}
