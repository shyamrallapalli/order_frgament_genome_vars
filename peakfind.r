# peak finding in a linear data
# source of function http://stats.stackexchange.com/a/36326
argmax <- function(x, y, w=1, ...) {
	require(zoo)
	n <- length(y)
	y.smooth <- loess(y ~ x, ...)$fitted
	y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
	delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
	i.max <- which(delta <= 0) + w
	list(x=x[i.max], i=i.max, y.hat=y.smooth)
}
