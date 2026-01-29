CSdataclean <- function(itime, delta) {
itorder <- order(itime)
Itime <- itime[itorder]
Delta <- delta[itorder]

begin <- which(Delta == 1)[1]  #### what if there is no yy == 1?
if( is.na(begin) ) stop("no delta = 1?")
end00 <- which(Delta == 0)   #### what if there is no yy == 0?
if( length(end00) == 0 ) stop("no delta = 0?")
Iend <- end00[length(end00)]
if( begin > Iend ) stop( "data = (zeros followed by ones)?" )
list(itime=Itime, delta=Delta, Istart=begin, Iend=Iend)
}
