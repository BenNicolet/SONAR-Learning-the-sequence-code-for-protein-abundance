Remove hybrids and pre-processing parallel reporter assay
================
Kaspar Bresser

- [Check SAM files](#check-sam-files)

## Check SAM files

Import SAM files

``` r
list.files("./Output/aligned/") %>% 
  map(~read_tsv(paste0("./Output/aligned/",.), , skip = 469, col_names = F)) %>% 
  set_names(list.files("./Output/aligned/")) -> tables
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

``` r
tables
```

    ## $`1A_R2.sam`
    ## # A tibble: 234,364 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 GCGG… GGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGGG… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    42 41M2S *         0     0 TCCT… 11FB… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TGTT… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TCTC… GGGG… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 GTAT… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGGG… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TCCT… GEHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GTCT… GGFH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TGGG… GGHH… AS:i… XN:i…
    ## # ℹ 234,354 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`1B_R2.sam`
    ## # A tibble: 197,209 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 ACCT… EGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GEHB… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 GCGG… F13@… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TCCT… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CACA… EGGF… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHG… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 ATAA… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 GCTT… A2EF… AS:i… XN:i…
    ## # ℹ 197,199 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`2A_R2.sam`
    ## # A tibble: 222,828 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GHHG… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M1S *         0     0 CATC… 5HHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TTGC… HFGF… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 ACCT… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GAGG… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 AAAT… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GCCC… GHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GGTG… GGHH… AS:i… XN:i…
    ## 10 M008…     4 *         0     0 *     *         0     0 CACA… 5D45… YT:Z… <NA> 
    ## # ℹ 222,818 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`2B_R2.sam`
    ## # A tibble: 156,033 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    23    41 3S36… *         0     0 CCTT… B255… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TCCT… FFHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    21    44 42M   *         0     0 CCCT… GHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 ATCC… GGGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CCGC… GGHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CTGC… GGGG… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TCGC… GGHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GACA… HHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 GTCT… HHHH… AS:i… XN:i…
    ## # ℹ 156,023 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`3A_R2.sam`
    ## # A tibble: 194,255 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 TCTC… GGEF… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 ATGC… GGGG… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 GTCT… HGHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 GTCT… HHGH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGCT… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TGCA… GHHG… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CTTC… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GTAA… 12FF… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GGGC… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    11 42M1S *         0     0 CCTT… A333… AS:i… XS:i…
    ## # ℹ 194,245 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`3B_R2.sam`
    ## # A tibble: 194,868 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CCCT… ?13G… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CCCT… ?E5G… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TTGC… FGGG… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGAA… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CACA… ADGB… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 ACCA… EE3A… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TGCT… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CCTT… EHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TGCA… HHHH… AS:i… XN:i…
    ## # ℹ 194,858 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`4A_R2.sam`
    ## # A tibble: 200,544 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 TCCC… HGHG… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 ATCC… GGGH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TACG… HGGH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CGAA… HHHH… AS:i… XN:i…
    ##  5 M008…     4 *         0     0 *     *         0     0 GCAG… 54AE… YT:Z… <NA> 
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 GGAT… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCCT… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TCGG… GGGH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GTGT… E1FG… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TCCC… GGGG… AS:i… XN:i…
    ## # ℹ 200,534 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`4B_R2.sam`
    ## # A tibble: 174,056 × 20
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    14 43M   *         0     0 GGCT… B23D… AS:i… XS:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCG… GGGG… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TTTT… HHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    25    42 5S38M *         0     0 CCTT… ?HHD… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TGGT… HHHG… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TCCA… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TGTT… HHHH… AS:i… XN:i…
    ##  8 M008…     4 *         0     0 *     *         0     0 GATT… HHHH… YT:Z… <NA> 
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CCTT… HHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 CAGT… FGDF… AS:i… XN:i…
    ## # ℹ 174,046 more rows
    ## # ℹ 7 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>, X20 <chr>
    ## 
    ## $`5A_R2.sam`
    ## # A tibble: 49,765 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 ACTT… HHHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    42 38M5S *         0     0 TTGC… GGHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    25    42 5S38M *         0     0 GCCC… AHHF… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 AGGT… HHGH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TACT… GGHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CATC… HHHH… AS:i… XN:i…
    ##  7 M008…     4 *         0     0 *     *         0     0 CAAA… HHHH… YT:Z… <NA> 
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 CTGC… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TCCC… GGGG… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ATCA… HHHH… AS:i… XN:i…
    ## # ℹ 49,755 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`5B_R2.sam`
    ## # A tibble: 197,927 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 GATC… HHHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TTTT… HFGH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 ATCC… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 GAAA… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCCA… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 ACGC… GGHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TTGC… EGHF… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 GAAA… HHHH… AS:i… XN:i…
    ## # ℹ 197,917 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`6A_R2.sam`
    ## # A tibble: 132,098 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 GACA… HHHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCA… EBHF… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 GCTG… HHHF… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TAAT… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 ACCT… 1EBE… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 CTCT… GFAG… AS:i… XN:i…
    ##  9 M008…    16 olig…    26    41 6S37M *         0     0 TAGG… EHGD… AS:i… XN:i…
    ## 10 M008…    16 olig…    26    41 6S37M *         0     0 TAGG… HHHG… AS:i… XN:i…
    ## # ℹ 132,088 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`6B_R2.sam`
    ## # A tibble: 158,772 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 41M2S *         0     0 TCCC… 31HF… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TCCC… HHHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 GCCT… HHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 CCGC… GGHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    14 43M   *         0     0 CCTT… 1255… AS:i… XS:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 GTGA… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    42 43M   *         0     0 GATC… //>/… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CTGC… >//>… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ACGA… GGHH… AS:i… XN:i…
    ## # ℹ 158,762 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`7A_R2.sam`
    ## # A tibble: 171,923 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CTCT… BGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCG… GGGH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CGCG… FFEB… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CGCG… EGEH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TCCT… GHGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TTTG… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCCT… E?FG… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TCCC… HGHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GCCC… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TTTC… HHHH… AS:i… XN:i…
    ## # ℹ 171,913 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $`7B_R2.sam`
    ## # A tibble: 166,673 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CCGC… GGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 AAAT… HHGF… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 ATTT… HFHG… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TCGC… GGGB… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGCC… GGGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 GGCT… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 42M2S *         0     0 TTTT… FHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    21    44 1S42M *         0     0 AGCC… 111B… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGGG… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ACTT… GHHG… AS:i… XN:i…
    ## # ℹ 166,663 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>

The reads contain the first 3 motifs of the synthetic UTR.

Define function that splits the motifs into columns and checks if all
motifs are identical.

``` r
get_motifs <- function(dat){
  dat %>% 
    mutate(read = X10) %>% 
 #   mutate(read = case_when(str_starts(X10, "T") ~ substring(X10, 2),
#                            TRUE ~ X10))  %>% 
    separate_wider_delim(read, "TATATCCGATCAT", names = c("motif1", "motif2"), too_few = "align_start") %>% 
    separate_wider_delim(motif2, "AACTGTACGCCT", names = c("motif2", "motif3"), too_few = "align_start") %>% 
    select(X1, X10, motif1, motif2, motif3) %>% 
    na.omit %>% 
    mutate(same = (motif1 == motif2 & motif2 == motif3 & motif3 == motif1))
}
```

Check the frequency of hybrids

``` r
tables %>% 
  map(get_motifs) %>% 
  map(pull, same) %>% 
  map(freq_table)
```

    ## $`1A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  31275  13.9
    ## 2 TRUE  193815  86.1
    ## 
    ## $`1B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  26940  14.2
    ## 2 TRUE  162552  85.8
    ## 
    ## $`2A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  29080  13.6
    ## 2 TRUE  184824  86.4
    ## 
    ## $`2B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  20214  13.5
    ## 2 TRUE  129285  86.5
    ## 
    ## $`3A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  23155  12.4
    ## 2 TRUE  163573  87.6
    ## 
    ## $`3B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  17890   9.6
    ## 2 TRUE  167830  90.4
    ## 
    ## $`4A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  21544  11.3
    ## 2 TRUE  169526  88.7
    ## 
    ## $`4B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  19318  11.7
    ## 2 TRUE  146046  88.3
    ## 
    ## $`5A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group     n  prop
    ##   <lgl> <int> <dbl>
    ## 1 FALSE  5342  11.2
    ## 2 TRUE  42304  88.8
    ## 
    ## $`5B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  22985  12.2
    ## 2 TRUE  165206  87.8
    ## 
    ## $`6A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  15895  12.6
    ## 2 TRUE  110630  87.4
    ## 
    ## $`6B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  19417  12.8
    ## 2 TRUE  132569  87.2
    ## 
    ## $`7A_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  20287  12.4
    ## 2 TRUE  143612  87.6
    ## 
    ## $`7B_R2.sam`
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  17169  10.8
    ## 2 TRUE  141110  89.2

Define function to remove hybrids. In this function, instead of
requiring exact matches, I allow 1 mismatch in the motifs to allow
sequencing errors. This is done using the `adist()` function.

``` r
remove_hybrids <- function(sam.table){
  sam.table %>% 
#       mutate(read = case_when(str_starts(X10, "T") ~ substring(X10, 2),
#                              TRUE ~ X10))  %>% 
    mutate(read = X10) %>% 
      separate_wider_delim(read, "TATATCCGATCAT", names = c("motif1", "motif2"), too_few = "align_start") %>% 
      separate_wider_delim(motif2, "AACTGTACGCCT", names = c("motif2", "motif3"), too_few = "align_start") %>% 
    na.omit() %>% 
    mutate(motif = map_chr(str_split(X3, "_"), 3)) %>% 
    mutate(across(c(motif1, motif2, motif3), ~map2_int(.x, motif, adist) )) %>% 
    filter(motif1 > 1 | motif2 > 1 | motif3 > 1) %>% 
    pull(X1) -> hybrids
  
  sam.table %>% 
    filter(!(X1 %in% hybrids)) %>% 
    na.omit()
}
```

Create a directory to put the filtered SAM files, and first transfer the
header sections to new files.

``` r
dir.create("./Output/aligned/filtered")

list.files("./Output/aligned/", pattern = ".sam" ) %>% 
  map(~read_lines(paste0("./Output/aligned/",.), n_max = 469)) %>% 
  map2(list.files("./Output/aligned/", pattern = ".sam" ), 
       ~write_lines(.x, paste0("./Output/aligned/filtered/filtered_",.y)))
```

Remove hybrids from the SAM tables and write the lines to the filtered
files.

``` r
tables %>% 
  map(remove_hybrids) %>% 
  map2(list.files("./Output/aligned/filtered/", pattern = ".sam" ), 
       ~write_tsv(.x, paste0("./Output/aligned/filtered/",.y), append = T, col_names = F))
```

Re-check the frequency of hybrids

``` r
list.files("./Output/aligned/filtered/") %>% 
  map(~read_tsv(paste0("./Output/aligned/filtered/",.), , skip = 469, col_names = F)) %>% 
  set_names(list.files("./Output/aligned/filtered/")) -> tables
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)
    ## One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

``` r
tables
```

    ## $filtered_1A_R2.sam
    ## # A tibble: 209,961 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 GCGG… GGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGGG… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TGTT… HHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TCTC… GGGG… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GTAT… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGGG… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCCT… GEHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GTCT… GGFH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TGGG… GGHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 GACA… HHHH… AS:i… XN:i…
    ## # ℹ 209,951 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_1B_R2.sam
    ## # A tibble: 176,141 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 ACCT… EGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GEHB… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 GCGG… F13@… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TCCT… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CACA… EGGF… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHG… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 ATAA… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 GCTT… A2EF… AS:i… XN:i…
    ## # ℹ 176,131 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_2A_R2.sam
    ## # A tibble: 200,483 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GHHG… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M1S *         0     0 CATC… 5HHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TTGC… HFGF… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 ACCT… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GAGG… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 AAAT… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GCCC… GHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GGTG… GGHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GHHH… AS:i… XN:i…
    ## # ℹ 200,473 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_2B_R2.sam
    ## # A tibble: 140,307 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 TCCT… FFHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    21    44 42M   *         0     0 CCCT… GHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 ATCC… GGGH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 CCGC… GGHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CTGC… GGGG… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCGC… GGHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GACA… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GTCT… HHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TCCT… GGHH… AS:i… XN:i…
    ## # ℹ 140,297 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_3A_R2.sam
    ## # A tibble: 176,592 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 TCTC… GGEF… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 ATGC… GGGG… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 GTCT… HGHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 GTCT… HHGH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGCT… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TGCA… GHHG… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CTTC… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GTAA… 12FF… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GGGC… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    11 42M1S *         0     0 CCTT… A333… AS:i… XS:i…
    ## # ℹ 176,582 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_3B_R2.sam
    ## # A tibble: 181,254 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CCCT… ?13G… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CCCT… ?E5G… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TTGC… FGGG… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGAA… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 CACA… ADGB… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 ACCA… EE3A… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TGCT… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CCTT… EHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TGCA… HHHH… AS:i… XN:i…
    ## # ℹ 181,244 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_4A_R2.sam
    ## # A tibble: 184,973 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 TCCC… HGHG… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 ATCC… GGGH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TACG… HGGH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CGAA… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGAT… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TCCT… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCGG… GGGH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 GTGT… E1FG… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TCCC… GGGG… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ATGC… HHHH… AS:i… XN:i…
    ## # ℹ 184,963 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_4B_R2.sam
    ## # A tibble: 2,245 × 20
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    14 43M   *         0     0 GGCT… B23D… AS:i… XS:i…
    ##  2 M008…    16 olig…    20    11 43M   *         0     0 ACGC… 103D… AS:i… XS:i…
    ##  3 M008…    16 olig…    20    11 43M   *         0     0 TAGT… 3B3A… AS:i… XS:i…
    ##  4 M008…    16 olig…    20    17 40M3S *         0     0 TACG… GGGH… AS:i… XS:i…
    ##  5 M008…    16 olig…    20    11 43M   *         0     0 TTCT… HECH… AS:i… XS:i…
    ##  6 M008…    16 olig…    22    11 2S41M *         0     0 CCCG… 0115… AS:i… XS:i…
    ##  7 M008…    16 olig…    20    11 43M   *         0     0 TCTT… F3HH… AS:i… XS:i…
    ##  8 M008…    16 olig…    20    14 43M   *         0     0 ACGC… ?FD5… AS:i… XS:i…
    ##  9 M008…    16 olig…    20    11 43M   *         0     0 CTGC… HHGH… AS:i… XS:i…
    ## 10 M008…    16 olig…    20    11 43M   *         0     0 CTGC… HHFH… AS:i… XS:i…
    ## # ℹ 2,235 more rows
    ## # ℹ 7 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>, X20 <chr>
    ## 
    ## $filtered_5A_R2.sam
    ## # A tibble: 45,801 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 ACTT… HHHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 AGGT… HHGH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TACT… GGHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CATC… HHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 CTGC… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TCCC… GGGG… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 ATCA… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 ATGT… HHHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 AGAT… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ACCT… HGHH… AS:i… XN:i…
    ## # ℹ 45,791 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_5B_R2.sam
    ## # A tibble: 180,191 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 GATC… HHHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 TTTT… HFGH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 ATCC… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 GAAA… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCCA… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 ACGC… GGHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TTGC… EGHF… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 GAAA… HHHH… AS:i… XN:i…
    ## # ℹ 180,181 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_6A_R2.sam
    ## # A tibble: 120,343 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 GACA… HHHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCA… EBHF… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CGCA… GHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 GCTG… HHHF… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TAAT… HHHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 ACCT… 1EBE… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 CCCT… GHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 CTCT… GFAG… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 ATAA… FA55… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TCCC… GHHH… AS:i… XN:i…
    ## # ℹ 120,333 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_6B_R2.sam
    ## # A tibble: 143,939 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 41M2S *         0     0 TCCC… 31HF… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 TCCC… HHHH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 GCCT… HHHH… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CCTT… GHHH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 CCGC… GGHH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    14 43M   *         0     0 CCTT… 1255… AS:i… XS:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 GTGA… HHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    42 43M   *         0     0 GATC… //>/… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 CTGC… >//>… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ACGA… GGHH… AS:i… XN:i…
    ## # ℹ 143,929 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_7A_R2.sam
    ## # A tibble: 156,960 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CTCT… BGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 CGCG… GGGH… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 CGCG… FFEB… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 CGCG… EGEH… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 TCCT… GHGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 TTTG… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 43M   *         0     0 TCCT… E?FG… AS:i… XN:i…
    ##  8 M008…    16 olig…    20    44 43M   *         0     0 TCCC… HGHH… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 GCCC… GHHH… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 TTTC… HHHH… AS:i… XN:i…
    ## # ℹ 156,950 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>
    ## 
    ## $filtered_7B_R2.sam
    ## # A tibble: 154,526 × 19
    ##    X1       X2 X3       X4    X5 X6    X7       X8    X9 X10   X11   X12   X13  
    ##    <chr> <dbl> <chr> <dbl> <dbl> <chr> <chr> <dbl> <dbl> <chr> <chr> <chr> <chr>
    ##  1 M008…    16 olig…    20    44 43M   *         0     0 CCGC… GGHH… AS:i… XN:i…
    ##  2 M008…    16 olig…    20    44 43M   *         0     0 AAAT… HHGF… AS:i… XN:i…
    ##  3 M008…    16 olig…    20    44 43M   *         0     0 ATTT… HFHG… AS:i… XN:i…
    ##  4 M008…    16 olig…    20    44 43M   *         0     0 TCGC… GGGB… AS:i… XN:i…
    ##  5 M008…    16 olig…    20    44 43M   *         0     0 GGCC… GGGH… AS:i… XN:i…
    ##  6 M008…    16 olig…    20    44 43M   *         0     0 GGCT… HHHH… AS:i… XN:i…
    ##  7 M008…    16 olig…    20    44 42M2S *         0     0 TTTT… FHHH… AS:i… XN:i…
    ##  8 M008…    16 olig…    21    44 1S42M *         0     0 AGCC… 111B… AS:i… XN:i…
    ##  9 M008…    16 olig…    20    44 43M   *         0     0 TTGC… GGGG… AS:i… XN:i…
    ## 10 M008…    16 olig…    20    44 43M   *         0     0 ACTT… GHHG… AS:i… XN:i…
    ## # ℹ 154,516 more rows
    ## # ℹ 6 more variables: X14 <chr>, X15 <chr>, X16 <chr>, X17 <chr>, X18 <chr>,
    ## #   X19 <chr>

``` r
tables %>% 
  map(get_motifs) %>% 
  map(pull, same) %>% 
  map(freq_table)
```

    ## $filtered_1A_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  11355   5.5
    ## 2 TRUE  193815  94.5
    ## 
    ## $filtered_1B_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   9340   5.4
    ## 2 TRUE  162552  94.6
    ## 
    ## $filtered_2A_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  10706   5.5
    ## 2 TRUE  184824  94.5
    ## 
    ## $filtered_2B_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   7390   5.4
    ## 2 TRUE  129285  94.6
    ## 
    ## $filtered_3A_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   8381   4.9
    ## 2 TRUE  163573  95.1
    ## 
    ## $filtered_3B_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   7528   4.3
    ## 2 TRUE  167830  95.7
    ## 
    ## $filtered_4A_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   9669   5.4
    ## 2 TRUE  169526  94.6
    ## 
    ## $filtered_4B_R2.sam
    ## # A tibble: 1 × 3
    ##   group     n  prop
    ##   <lgl> <int> <dbl>
    ## 1 FALSE  2209   100
    ## 
    ## $filtered_5A_R2.sam
    ## # A tibble: 2 × 3
    ##   group     n  prop
    ##   <lgl> <int> <dbl>
    ## 1 FALSE  2380   5.3
    ## 2 TRUE  42304  94.7
    ## 
    ## $filtered_5B_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE  10202   5.8
    ## 2 TRUE  165206  94.2
    ## 
    ## $filtered_6A_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   6699   5.7
    ## 2 TRUE  110630  94.3
    ## 
    ## $filtered_6B_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   7705   5.5
    ## 2 TRUE  132569  94.5
    ## 
    ## $filtered_7A_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   9026   5.9
    ## 2 TRUE  143612  94.1
    ## 
    ## $filtered_7B_R2.sam
    ## # A tibble: 2 × 3
    ##   group      n  prop
    ##   <lgl>  <int> <dbl>
    ## 1 FALSE   9071     6
    ## 2 TRUE  141110    94
