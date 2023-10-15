# README

## Refactor Notes:
- Check gene selection logic: sampling with replacement: https://www.statisticshowto.com/sampling-with-replacement-without
- 




### Assumptions:
- *where the gene representing each locus was chosen uniformly at random* -pg 8
- Assuming "uniformly at random" indicates the chosen gene from each of the 12 locus is chosen using the same randomization strategy.
        
        random.randrange(0,len(locus))
- 

### Permutation Test:
- [scipy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html)

- [concepts](https://towardsdatascience.com/how-to-use-permutation-tests-bacc79f45749)
