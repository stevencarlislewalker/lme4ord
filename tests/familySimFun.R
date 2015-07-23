library(lme4ord)
stopifnot(all.equal(familySimFun(gaussian),
                    familySimFun(gaussian())))          
stopifnot(all.equal(familySimFun("gaussian"),
                    familySimFun(gaussian())))                    


