# Méthodes de Monte Carlo

## Plan du cours

* Methods for sampling random variables: inverse of the cdf, rejection method, sampling of the normal distribution.
* Sampling the Brownian motion: direct method or Brownian bridges.
* Variance reduction methods: antithetic variates, control variates, conditioning approach, stratified sampling and importance sampling.
* Computing sensitivities.
* Discretization of Stochastic Differential Equation: the Euler scheme.
* Multi-level Monte Carlo.

## TP

The practical sessions heavily rely on the PNL library [https://pnlnum.github.io/pnl/](https://pnlnum.github.io/pnl/) (Download tab). Make sure to install the latest stable version from [GitHub](https://github.com/pnlnum/pnl/releases).

* Séance 1 ([sujet](document/tp1.pdf), solution sous la forme d'un notebook python [hitting-time.ipynb](document/hitting-time.ipynb))
* Séance 2 ([sujet](document/tp2.pdf), [squelette](document/skel-is-toy.zip))
* Séance 3 ([sujet](document/tp3.pdf), [squelette](document/skel-tp3.zip))
* Séance 4 ([sujet](document/tp4.pdf))
* Séance 5 ([sujet](document/tp5.pdf))
* Séance 6 ([sujet](document/tp6.pdf), [squelette](document/skel-tp6.zip))
* Séance 7 ([sujet](document/tp7.pdf), [squelette](document/skel-tp7.zip))
* Séance 8 ([sujet](document/tp8.pdf), [squelette](document/skel-tp8.zip))

## Exercices

* Simulation de variables aléatoires : [exo_simul.pdf](document/exo_simul.pdf)

## References on C++/.NET

### C++/CLI

* [Presentation slides](document/PresentationCppCli.pdf)
* [MSDN presentation of C++/CLI](http://msdn.microsoft.com/en-us/library/ms379617%28VS.80%29.aspx)
* [Using C++/CLI with Visual Studio (MSDN)](http://msdn.microsoft.com/fr-fr/library/bb384847.aspx)
* [Nishant Sivakumar's blog](http://blog.voidnish.com/?cat=2)
* [A C++/CLI course on Developpez.com (french)](http://dotnet.developpez.com/faq/cppcli/)
* [Another C++/CLI course on FunctionX.com](http://www.functionx.com/cppcli/index.htm)

### Graphical libraries

* [FusionCharts Free official page](http://www.fusioncharts.com/free/)
* [Zedgraph official page](http://zedgraph.org/wiki/index.php?title=Main_Page)
* Open Flash Chart (open-source, not yet tested)
  * [official page](http://teethgrinder.co.uk/open-flash-chart-2/)
  * [An example](http://monstersgotmy.net/post/OpenFlash-Chart-Fixing-Open-Source.aspx)
* Microsoft Chart Control (.NET 3.5)
  * [Components](http://www.microsoft.com/downloads/details.aspx?FamilyId=130F7986-BF49-4FE5-9CA8-910AE6EA442C&displaylang=en)
  * [Add-ons for VS2008](http://www.microsoft.com/downloads/details.aspx?familyid=1D69CE13-E1E5-4315-825C-F14D33A303E9&displaylang=en)
  * [Documentation](http://www.microsoft.com/downloads/details.aspx?familyid=EE8F6F35-B087-4324-9DBA-6DD5E844FD9F&displaylang=en)