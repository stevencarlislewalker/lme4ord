## FIXME: both fail but shouldn't!
try(splitForm(incidence/size ~ period + lme4(1 | herd)))
try(splitForm(I(incidence/size) ~ period + lme4(1 | herd)))
