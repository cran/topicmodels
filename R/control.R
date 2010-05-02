##*******************************************************
## Classes TMcontrol, LDAcontrol, CTMcontrol
##
## control parameters for the lda and ctm functions
## + superclass (TMcontol)


##**********************************************************
## coercion
setAs("NULL", "CTM_VEMcontrol", function(from, to) new(to))
setAs("NULL", "LDA_VEMcontrol", function(from, to) new(to))
setAs("NULL", "LDA_Gibbscontrol", function(from, to) new(to))
setAs("NULL", "OPTcontrol", function(from, to) new(to))

setAs("list", "LDA_VEMcontrol", function(from, to) .list2VEMcontrol(from, to))
setAs("list", "LDA_Gibbscontrol", function(from, to) .list2object(from, to))
setAs("list", "CTM_VEMcontrol", function(from, to) .list2VEMcontrol(from, to))

.list2VEMcontrol <- function(from, to) {
  n = names(from)
  s = slotNames(to)
  p = pmatch(n, s)
  if(any(is.na(p)))
    stop(paste("\nInvalid slot name(s) for class",
               to, ":", paste(n[is.na(p)], collapse=" ")))
  slotTypes <- getClass(to)@slots[p]
  from <- lapply(seq_along(from), function(i) as(from[[i]], slotTypes[[i]]))
  names(from) = s[p]
  do.call("new", c(from, Class=to))
}

setAs("list", "OPTcontrol", function(from, to) .list2object(from, to))

##**********************************************************
## .list2object is copied from package flexmix 

.list2object = function(from, to){
    n = names(from)
    s = slotNames(to)
    p = pmatch(n, s)
    if(any(is.na(p)))
        stop(paste("\nInvalid slot name(s) for class",
                   to, ":", paste(n[is.na(p)], collapse=" ")))
    slotTypes <- getClass(to)@slots[p]
    from <- lapply(seq_along(from), function(i) as(from[[i]], slotTypes[[i]]))
    names(from) = s[p]
    do.call("new", c(from, Class=to))
}


