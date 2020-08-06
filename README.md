# Phenotype

Remove Outliers from Phenotypic Data and Performs the Best Linear Unbiased Prediction

## Installation

``` r
devtools::install_github("biozhp/Phenotype")
```

## Usage

### Remove outliers from phenotypic data

``` r
data("wheatph")
inlier <- outlier(x = wheatph, col = 5, fold = 1.5)
```

### Calculate statistical indicators of phenotypic data

``` r
data("wheatph")
inlier <- outlier(x = wheatph, col = 5, fold = 1.5)
stat_out <- stat(x = inlier, col = 2)
```

### Histogram drawing

``` r
data("wheatph")
inlier <- outlier(x = wheatph, col = 5, fold = 1.5)
stat_out <- stat(x = inlier, col = 2)
histplot(x = stat_out$mean)
```

### Performs the Best Linear Unbiased Prediction (BLUP)

``` r
data("wheatph")
blup_out <- blup(x = wheatph, fold = 1.5)
```

## Contact
For any bugs/issues/suggestions, please send emails to: Peng Zhao pengzhao@nwafu.edu.cn.