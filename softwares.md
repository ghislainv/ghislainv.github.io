---
layout: page
title: Softwares
permalink: /softwares/
---

## `hSDM` R package

<div id="img_softwares">
  <img src="/images/softwares/Agrandidieri.jpg" alt="Adansonia grandidieri"/>
</div>

`hSDM` is an R package for hierarchical Bayesian species distribution models. It includes functions to fit mixture models (site-occupancy, N-mixture, ZIB and ZIP models) accounting for imperfect detection, excess of zeros in the observations and spatial autocorrelation (through an intrinsic CAR process). The functions uses an adaptive Metropolis within Gibbs algorithm written in C code. This makes parameter inference faster than with softwares commonly used to fit such models (such as JAGS, WinBUGS or OpenBUGS) and allows analyzing very large data-sets (typically with more than tens of thousands grid cells).

### Related publications

**Latimer A. M., Wu S. S., Gelfand A. E. and Silander J. A.** 2006\. Building statistical models to analyze species distributions. _Ecological Applications_. **16**(1): 33-50\. [![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](/publications/biblio/Latimer2006-EcologicalApplications.pdf)

**MacKenzie D. I., Nichols J. D., Lachman G. B., Droege S., Royle J. A. and Langtimm C. A.** 2002\. Estimating site occupancy rates when detection probabilities are less than one. _Ecology_. **83**: 2248-2255\. [![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](/publications/biblio/MacKenzie2002-Ecology.pdf)

**Royle, J. A.** 2004\. N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. _Biometrics_. **60**: 108-115\. [![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](/publications/biblio/Royle2004-Biometrics.pdf)

### Links

`hSDM` website: [http://hSDM.sf.net](http://hSDM.sf.net "hSDM web site")

## Contribution to `MCMCpack`

![MCMCpack](/images/softwares/MCMCpack.jpg)

[`MCMCpack`](http://mcmcpack.wustl.edu/ "MCMCpack at wustl") (Markov chain Monte Carlo Package) is an R package which contains functions to perform Bayesian inference using posterior simulation for a number of statistical models. Most simulation is done in compiled C++ written with the [`Scythe`](http://scythe.wustl.edu/) Statistical Library Version 1.0.2. Authors: Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park.

Since version 1.1-1 (July 2011), `MCMCpack` includes additional functions for generalized linear mixed models (glmm): `MCMChregress()` for Gaussian models, `MCMChlogit()` for Bernoulli models (logit link function) and `MCMChpoisson()` for Poisson models (log link function). Author: Ghislain Vieilledent. 

### Links

Last version of `MCMCpack` is available from the Comprehensive R Archive Network at [http://cran.r-project.org/package=MCMCpack](http://CRAN.R-project.org/package=MCMCpack "MCMCpack at CRAN")

### Related publications

**Martin A. D., Quinn K. M. and Jong Hee Park** 2011. `MCMCpack`: Markov Chain Monte Carlo in R. _Journal of Statistical Software_. **42**(9): 1-21. [![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](/publications/biblio/Martin2011-JSS.pdf) 

**Pemstein D., Quinn K. M. and Martin A. D.** 2011\. The `Scythe` Statistical Library: An Open Source C++ Library for Statistical Computation. _Journal of Statistical Software_. **42**(12):1-26. [![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](/publications/biblio/Pemstein2011-JSS.pdf)

## `phcfM` R package

<div id="img_softwares">
  <img src="/images/softwares/deforestation.jpg" alt="deforestation"/>
</div>

`phcfM` is an R package for modelling anthropogenic deforestation. It was initially developed to obtain REDD+ baseline scenarios of deforestation for the _programme holistique de conservation des forêts à Madagascar_ (from which the package was named after). It includes two main functions:

1. `demography()`, to model the population growth with time in a hierarchical Bayesian framework using population census data and Gaussian linear mixed models.
2. `deforestation()`, to model the deforestation process in a hierarchical Bayesian framework using land-cover change data and Binomial logistic regression models with variable time-intervals between land-cover observations.

### Code and manual

The last stable version of the `phcfM` R package is officially available for several operating systems (Unix, Windows and Mac OSX) on the Comprehensive R Archive Network ([CRAN](http://cran.r-project.org/web/packages/phcfM/index.html)).

### Example

A GRASS location (`phcfM_SM`) and two mapsets with geographical data layers (`PERMANENT` and `study_area_4`) are available to illustrate the use of the `phcfM` R package. Associated to the GRASS location, a directory (`./scripts`) includes the data and the R/GRASS scripts used for the demographic and deforestation models.

- GRASS location: [phcfM_SM.zip](http://sourceforge.net/projects/phcfm/files/phcfM_SM.zip/download)
- `./scripts`: [scripts.zip](http://sourceforge.net/projects/phcfm/files/scripts.zip/download)

### Related publication

**<span style="text-decoration: underline;">Vieilledent G.</span>, C. Grinand and R. Vaudry.** 2013. Forecasting deforestation and carbon emissions in tropical developing countries facing demographic expansion: a case study in Madagascar. _Ecology and Evolution_. **3**:1702-1716.
\[doi:[10.1002/ece3.550](http://dx.doi.org/10.1002/ece3.550)\].
[![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](http://onlinelibrary.wiley.com/doi/10.1002/ece3.550/pdf) /
Supplementary materials [![supplements](/images/logos/logo-pdf.png "supplements")](http://onlinelibrary.wiley.com/store/10.1002/ece3.550/asset/supinfo/ece3550-sup-0001-AppendixS1-S8.pdf?v=1&s=60392770526cfbe1f23b0aa1cc92f3b0e922e267)

## `twoe` (2e) software:

![twoe](/images/softwares/twoe.png)

`twoe` (2e) is a software which aims first, at estimating the demographic parameters of tropical tree species from permanent forest plot data (through an R package) and second, at simulating forest dynamics (through a [Capsis](http://capsis.cirad.fr/) module). Authors: Ghislain Vieilledent, François de Coligny.

### Related publication

**Laurans M., B. Hérault, <span style="text-decoration: underline;">G. Vieilledent</span> and G. Vincent.** 2014. Vertical stratification reduces competition for light in dense tropical forests. _Forest Ecology and Management_. **329**: 79-88.
\[doi:[10.1016/j.foreco.2014.05.059](http://dx.doi.org/10.1016/j.foreco.2014.05.059)\].
[![manuscript in pdf](/images/logos/logo-pdf.png "manuscript in pdf")](/publications/Laurans2014-FEM.pdf) /
Supplementary materials [![supplements](/images/logos/logo-zip.png "supplements")](/publications/Laurans2014-FEM-SM.zip)

### Link

Last versions of the R package and of the Capsis module are available on the `twoe` website: [http://twoe.sf.net](http://twoe.sf.net "twoe").