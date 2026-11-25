# Reference manual

The supplied CRAN-style PDF is compiled from the standalone LaTeX file:

```sh
pdflatex PORTO-reference.tex
makeindex PORTO-reference.idx
pdflatex PORTO-reference.tex
pdflatex PORTO-reference.tex
```

For the native R-generated manual, run `R CMD Rd2pdf --output=PORTO-reference-native.pdf PORTO` from the parent directory. The `man/` files are the corresponding native Rd documentation. The supplied PDF is not an official CRAN publication or a native Rd2pdf rendering.
