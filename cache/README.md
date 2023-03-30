# cache

This directory should hold any temporary files that are created to cache results that are particularily slow to compute. The most versatile storage format for all kinds of data is the R Data Storage (`.rds`) file format. The following code example illustrates how to cache an object as an `.rds` file and read it back in from cache:

```
# x can be any variable name holding any type of data
x <- 1:5

# save to cache (change 'x.rds' to a filename reflecting the name of your variable)
readr::write_rds(x, path = file.path("cache", "x.rds"))

# read back from cache
x <- readr::read_rds(path = file.path("cache", "x.rds"))
```
