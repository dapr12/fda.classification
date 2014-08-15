  `fda.classification`


## About
fda.classification This set of functions visualise Functional Data (data comming from functions) at a discrete time using R.
The Detailed proposal is available on [this blog](http://dapr12.wordpress.com/) 

## Install

### Install the development version using `install_github` within Hadley's [devtools](https://github.com/hadley/devtools) package.

```R
install.packages("devtools")
require(devtools)

devtools::install_github("dapr12/fda.classification")
library('fda.classification')
```

Note: 

Windows users have to first install [Rtools](http://cran.r-project.org/bin/windows/Rtools/).

### Packages `fda.classification` depends on
+ [fda] (http://cran.r-project.org/web/packages/fda/index.html)
+ [splines] (http://cran.r-project.org/web/packages/splines/index.html)
+ [MASS] (http://cran.r-project.org/web/packages/MASS/index.html)


For examples
+ [fdaExamples] (https://github.com/dapr12/fdaexamples)


### Functions currently available


#### FdaClass

```r
fdaclass(mdata, argval, rangeval) 
```


#### Outliergram

```r
outliergram(fdaobj)
```


#### smoothfda

```r
smoothfda(fdaobj , bandwidth, degree )
```

#### GCV for Choosing the Smoothing Parameter 

```r
gcsvc( fdaobj, norder, lambda, Lf, Intv  )
```


#### Smooth using B-Splines

```r
smoothbsplines( fdaobj, norder, lambda, Lf)
```


#### Functional Derivative

```r
fderiv(fdaobj,nderiv) 
```


#### Medianfd

```r
Medianfd(Smoothfda)
```

#### D-Plot
```r
dpout(Eigenvalues, plotting)
```


#### pcafd

```r
pcafd( fdaobj, nharm )
```


#### Density Scores

```r
densityScores(pcaobj2, 2)
```


#### Functional Data Density Estimation base on the Harmonics 

```r
fdensity( fdaobj, pcaobj, bandwith, plotting)
```


#### Functional Variance 

```r
varfd(Smoothfda)
```


#### Classification 

```r
classfd(Classlearn, train, test)
```


#### Near Perfect Classification for Functional Data  (Based on the Centroid Method)

```r
pdfclasf(data, test, indClass0, indClass1, indtest )
```


#### Simulate Functional Data  

```r
simulatefda( nsamples, ndrawn, rangeval, mean, sigma )
```

