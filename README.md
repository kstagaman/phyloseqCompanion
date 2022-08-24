# Installation

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
for (pkg in c("remotes", "dada2", "phyloseq", "ALDEx2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
remotes::install_github("ggloor/CoDaSeq/CoDaSeq")
remotes::install_github("kstagaman/phyloseqCompanion")
```
