
a reproducible analysis of [haxby 2001 data](http://www.ncbi.nlm.nih.gov/pubmed/11577229)

download the subjects via 

`wget http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz`

where you can install wget on osx via homebrew (or use yum on some linux versions)

then run ... 

```
for each subject
  for each run 
    mask selects voxels ( either  haxby mask or sccan/other automated procedure ) from training data
    set current run as test data                            # matrix   y = 121 by size of mask 
    train on other runs                                          # matrix   x = 121*11 by size of mask 
    features <- svd( matrix x , n )                          # n selects number of components, e.g.  n = 100 
    trainingmodel <- svm( features | training class labels (rest, scissors, shoes ...  ) )
    Evaluate(  predict( y ) |  testing class labels )  # majority vote 
  done 
done
```
