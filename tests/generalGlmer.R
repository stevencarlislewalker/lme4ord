library(lme4ord)
library(lme4)

## test formula parsing
form <- mite ~ ddd + flexvar(WatrCont, a, phy = huh) + x + (z + w | g) + jjj + derek(gjg) + identity(1 + x + fhg | ghjk)

stopifnot(identical(nobars(form),
                    mite ~ ddd + flexvar(WatrCont, a, phy = huh) + x + jjj + derek(gjg)))
stopifnot(identical(noSpecials(form),
                    mite ~ ddd + x + (z + w | g) + jjj + derek(gjg)))
stopifnot(identical(nobars(noSpecials(form)),
                    noSpecials(nobars(form))))

