

#' @include registry-seed.R
NULL

## Register base seeding methods
# None: do nothing and return object unchanged
setKINOMOSeed('none', function(object, x, ...){object}, overwrite=TRUE)
# Random: use function rKINOMO
setKINOMOSeed('random', rKINOMO, overwrite=TRUE)
