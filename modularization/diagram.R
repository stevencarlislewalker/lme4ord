library(diagram)
library(Matrix)
nms <- c("structured\nGLMM",
         "formula\nparsing",
         "deviance\nfunction\nconstruction",
         "deviance\nfunction\noptimization",
         "output and\ninference",
         "random-effects\nstructures", "repeated\nsparse\nmatrices")
A <- sparseMatrix(i = c(2, 3, 4, 5, 6, 6, 7, 7, 7),
                  j = c(1, 1, 1, 1, 2, 5, 2, 5, 6),
                  x = 1,
                  dims = c(7, 7))
A <- as.matrix(A)
firstRow <- 0.5
secondRow <- seq(0.12, 0.88, length = 4)
thirdRow <- c(0.25, 0.75)
Apos <- matrix(c( firstRow[1], 0.8,
                 secondRow[1], 0.5,
                 secondRow[2], 0.5,
                 secondRow[3], 0.5,
                 secondRow[4], 0.5,
                  thirdRow[1], 0.2,
                  thirdRow[2], 0.2), 7, 2, byrow = TRUE)
dimnames(A) <- rep(list(nms), 2)

png("~/Documents/lme4git/lme4ord/modularization/diagram.png")
plotmat(A, cex.txt = 0, box.cex = 0.9,
        pos = Apos, curve = 0)
dev.off()


