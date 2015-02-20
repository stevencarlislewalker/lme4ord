data(fish)
data(limn)
Y <- as.matrix(fish)
n <- nrow(Y)
m <- ncol(Y)
x <- as.vector(scale(limn$pH))

dl <- data.list(Y = Y, x = x,
                dimids = c("sites", "species"))
summary(dl)

df <- as.data.frame(dims_to_vars(dl))

U <- try(matrix(0, nrow = dd[loadingsDim], ncol = latentDims))
if(inherits(U, "try-error")) stop("loadingsDim does not index a dimension of data")
nFreeLoadings <- (dd[loadingsDim] * latentDims) - choose(latentDims, 2)
U[lower.tri(U, TRUE)] <- 1:nFreeLoadings
latentVarNames <- paste(latentName, 1:latentDims, sep = "")
U <- setNames(as.data.frame(U), latentVarNames)
