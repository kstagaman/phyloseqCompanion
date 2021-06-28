# Installation

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
for (pkg in c("devtools", "dada2", "phyloseq", "ALDEx2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
devtools::install_github("ggloor/CoDaSeq/CoDaSeq")
devtools::install_github("kstagaman/phyloseqCompanion")
```
