<head>
  <meta name="google-site-verification" content="8RxjfD3zAKng0Ht73dgzJk9f9ddsO9Do3k_d8CXpxpM" />
</head>

# Package "ULT" [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18914.svg)](http://dx.doi.org/10.5281/zenodo.18914)
# ULT stands for "Useful Little Things"
## This is a package designed to create some useful functions; no peculiar topic, only things that I found difficult to do and I tried to make easier

### Here are the main topics of my function :
#### - a particularily useful function, extract.from.factor, allowing one to "extract" values from a factor vector (or data associated to a factor vector, in a data frame for instance) according to given values
#### - some graphical useful functions, as contrasting.palette (a palette with 22 contrasting colors) and transparent.colors (I felt like, when I was plotting things with transparent colors, it was quite difficult to have the RGB values of these colors to re-use them, for instance in the legend or other things)
#### - a useful function to create a tree stepwise (work in progress), by specifying first the terminal units, then to link them one by one, specify the node ages etc ; I felt that, up to ~10 taxa, it was OK to write a tree in NEWICK format manually, but with a lot of taxa, I was struggling and doing mistakes
#### - some functions related to a multivariate extension of the Kruskal-Wallis test written by Fanyin He in her PhD; she published a work using this but did not published the R code, I improved and fitted it a bit, but all the credits go to her

### I plan to build other functions, especially a function to draw a tree and be completely parametrable

### Feel free to use, please just cite this package if used (and the references associated to the used functions, if relevant) and let me know it (jacob.maugoust@umontpellier.fr)
