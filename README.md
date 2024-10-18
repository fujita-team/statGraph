# StatGraph

Contains statistical methods to analyze graphs, such as graph parameter estimation, model selection based on the GIC (Graph Information Criterion), statistical tests to discriminate two or more populations of graphs (ANOGVA - Analysis of Graph Variability), correlation between graphs, and clustering of graphs. 

##### *References*:

* Takahashi et al. (2012) <doi:10.1371/journal.pone.0049949>, 
* Futija et al. (2017) <doi:10.3389/fnins.2017.00066>, 
* Fujita et al. (2017) <doi:10.1016/j.csda.2016.11.016>, 
* Tang et al. (2017) <doi:10.3150/15-BEJ789>, 
* Tang et al. (2017) <doi:10.1080/10618600.2016.1193505>, 
* Ghoshdastidar et al. (2017) <arXiv:1705.06168>,  
* Ghoshdastidar et al. (2017) <arXiv:1707.00833>, 
* Cerqueira et al. (2017) <doi:10.1109/TNSE.2017.2674026>, 
* Fraiman and Fraiman (2018) <doi:10.1038/s41598-018-23152-5>, 
* Fujita et al. (2019) <doi:10.1093/comnet/cnz028>.

## For Developers:

### Building:

Before attempting to build this project, make sure that you have all of the requirements needed by running:

> `apt install -y r-base libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev r-cran-tinytex texlive texinfo texlive-fonts-extra curl libcurl4 libharfbuzz-dev libfribidi-dev libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev libgit2-dev`

Once the requirements are installed, you can build the project by running:

    ./build.r


or

    Rscript build.r

This will generate a statGraph_{*version*}.pdf file containing the documentation and a statGraph_{*version*}.tar.gz file containing the finished library.

If you are running the [Nix Package Manager](https://nixos.org/guides/how-nix-works/), you can just run `nix develop` inside this directory and get a shell with all the requirements for building the library.

### Testing

Test files are located in directory statGraph/tests/testthat/. 

In order to perform tests, run the following:

    ./test.r 

Alternative, open the file `statGraph.Rproj` in RStudio and simply click on `Test`. This will run all tests. To run individual test files, open them in RStudio and click on `Run Tests`.

### Structure of this project

This project uses the standard structure of a R project described in https://r-pkgs.org/whole-game.html

In short: Methods should be added to .R files located in statGraph/R/. To make them public, they have to be exported by adding the corresponding @export tag before their constructors.

More specific guidelines for the structure of this codebase can be found on statGraph/R/README.md and statGraph/tests/testthat/README.md


