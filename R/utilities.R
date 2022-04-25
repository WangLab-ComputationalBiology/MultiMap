# Modify parameters in calling environment
#
# Used exclusively for helper parameter validation functions
#
# @param param name of parameter to change
# @param value new value for parameter
#
ModifyParam <- function(param, value) {
  # modify in original function environment
  env1 <- sys.frame(which = length(x = sys.frames()) - 2)
  env1[[param]] <- value
  # also modify in validation function environment
  env2 <- sys.frame(which = length(x = sys.frames()) - 1)
  env2[[param]] <- value
}
