# The Sharpton Lab dada2 pipeline

This package is designed to make it easy to replicably conduct the dada2 pipeline as run by the Sharpton lab at Oregon State University.

To begin, if you want to build phylogenetic trees, you need the [mothur](https://mothur.org/) and [FastTree](http://www.microbesonline.org/fasttree/) programs installed on the machine on which you will be conducting the pipeline.

Additionally, it relies on another R package that I've written and is only available on GitHub. You can install that package by running the following:

```
devtools::install_github("kstagaman/phyloseqCompanion")
```

This package can be installed by running the following

```
devtools::install_github("/kstagaman/sharpton-lab-dada2-pipeline")
```

To get started, first generate the template script by running in your working directory.

```
library(dada2.pipeline)
generate.template()
```

