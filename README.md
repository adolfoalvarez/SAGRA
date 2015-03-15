# SAGRA
Splittind And Group Recombining Algorithm

This is the implementation of the SAGRA algorithm, proposed by 
[Álvarez and Peña (2014)](http://e-archivo.uc3m.es/handle/10016/19233)

SAGRA is a new clustering methodology which, based in a strategy of splitting, 
cleaning and recombining, is able to detect groups inside a data sample. 
As splitting rule we use the discriminator function, where points 
sharing the same discriminator are classified into the same group, defining in this
way a partition of the sample. For cleaning we detect and purge the outliers of
each group, and finally for recombining, we propose the use of a Bayes factor
to weight the likelihood of the sample given two models: one where all data is
generated from a single distribution, and other when the distribution is a mixture
estimated from the obtained partition.

The main function of the package is:

* `sagra()` takes a data frame and return a clustering solution.

## Installation

**This is a development version, mainly intended for internal improvements** since many bugs are expected.
Nevertheles, this version can be installed using:

```R
# install.packages("devtools")
devtools::install_github("adolfoalvarez/SAGRA")
```
