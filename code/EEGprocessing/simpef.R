# "H:\mt_behavalluse\simpef.R"
simpef <- function(model, int, s1, s2) { # e.g. s1=action, s2=valence
  flab = names(fixef(model))
  cv = vector("numeric", length(fixef(model)))
  cv[flab == int] = 1
  cv[flab == s1]  = 1
  p1 = doBy::esticon(model, cv)
  cat("simple ", s1, " effect for ", s2, " == 1: est = ", p1$Estimate,
      ", X2(", p1$DF, ") = ", p1$X2.value, ", p = ", p1$Pr, "\n", sep="")
  cv[flab == s1]  = -1
  p2 = doBy::esticon(model, cv)
  cat("simple ", s1, " effect for ", s2, " == -1: est = ", p2$Estimate,
      ", X2(", p2$DF, ") = ", p2$X2.value, ", p = ", p2$Pr, "\n", sep="")
  c(p1$Std.Error, p2$Std.Error, p1$Estimate, p2$Estimate)
}
