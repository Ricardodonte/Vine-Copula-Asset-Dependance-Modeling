library(VineCopula)

FAMILIES <- 0:10

family_label <- function(f) {
  base <- c(`0` = "I", `1` = "N", `2` = "t", `3` = "C", `4` = "G", `5` = "F", `6` = "J",
            `7` = "BB1", `8` = "BB6", `9` = "BB7", `10` = "BB8")
  # Les dizaines indiquent la rotation dans les codes VineCopula.
  out <- character(length(f))
  for (i in seq_along(f)) {
    fi <- f[i]
    if (fi <= 10) { out[i] <- base[as.character(fi)]; next }
    base_code <- if (fi %in% 13:20) fi - 10 else if (fi %in% 23:30) fi - 20 else fi - 30
    deg <- if (fi %in% 13:20) "180" else if (fi %in% 23:30) "90" else "270"
    out[i] <- paste0(base[as.character(base_code)], deg)
  }
  out
}

# La version d’origine exclut les rotations et le test d’indépendance.
fit_vine <- function(U, original = FALSE, type = "RVine") {
  if (original) {
    RVineStructureSelect(U, familyset = FAMILIES, type = type, rotations = FALSE)
  } else {
    RVineStructureSelect(U, familyset = FAMILIES, type = type, selectioncrit = "AIC",
                         indeptest = TRUE, level = 0.05, rotations = TRUE)
  }
}

vine_table <- function(rvm) {
  M <- rvm$Matrix; d <- nrow(M); nm <- rvm$names
  td <- RVinePar2Tau(rvm)
  lt <- matrix(0, d, d); ut <- matrix(0, d, d)
  for (i in 2:d) for (j in 1:(i - 1)) {
    fam <- rvm$family[i, j]
    if (fam == 0) next
    tdep <- BiCopPar2TailDep(fam, rvm$par[i, j], rvm$par2[i, j])
    lt[i, j] <- tdep$lower; ut[i, j] <- tdep$upper
  }
  rows <- list()
  for (tree in 1:(d - 1)) {
    i <- d - tree + 1
    for (j in 1:(d - tree)) {
      a <- nm[M[j, j]]; b <- nm[M[i, j]]
      cond <- if (i < d) paste(nm[M[(i + 1):d, j]], collapse = ",") else ""
      rows[[length(rows) + 1]] <- data.frame(
        Arbre = tree,
        Paire = if (cond == "") paste0(a, "-", b) else paste0(a, "-", b, " | ", cond),
        Famille = family_label(rvm$family[i, j]),
        par1 = round(rvm$par[i, j], 3), par2 = round(rvm$par2[i, j], 3),
        tau = round(td[i, j], 3), LTD = round(lt[i, j], 3), UTD = round(ut[i, j], 3))
    }
  }
  out <- do.call(rbind, rows)
  attr(out, "fit") <- c(logLik = rvm$logLik, AIC = rvm$AIC, BIC = rvm$BIC)
  out
}

pair_dependence <- function(u, v, rotations = TRUE) {
  fams <- if (rotations) c(0:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40) else 0:10
  fit <- BiCopSelect(u, v, familyset = fams, selectioncrit = "AIC", indeptest = TRUE, level = 0.05)
  if (fit$family == 0) return(c(tau = 0, LTD = 0, UTD = 0, family = 0))
  td <- BiCopPar2TailDep(fit$family, fit$par, fit$par2)
  c(tau = fit$tau, LTD = td$lower, UTD = td$upper, family = fit$family)
}
