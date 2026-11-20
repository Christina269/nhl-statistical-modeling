#3
# NHL 2021-22 : ΜΟΝΤΕΛΟ ZERO-DEFLATED TRUNCATED SKELLAM
#
# Η διαφορά δύο ανεξάρτητων Poisson ακολουθεί κατανομή Skellam:
#
#   P(D = k) = exp(-(l1+l2)) * (l1/l2)^(k/2) *
#              I_|k|(2*sqrt(l1*l2))
#
# Στο ice hockey η ισοπαλία είναι αδύνατη, οπότε η μάζα στο D = 0
# αφαιρείται και αναδιανέμεται:
#
#   P(D = k | D != 0) = P(D = k) / (1 - P(D = 0))
#
# Το script προσαρμόζει δύο εκδοχές:
#   (α) post-hoc: τα lambda του Double Poisson
#   (β) ML: απευθείας μεγιστοποίηση της περικομμένης πιθανοφάνειας
#
# ΣΗΜΕΙΩΣΗ για την Αρνητική Διωνυμική:
# Δεν προέκυψε ένδειξη υπερδιασποράς που να δικαιολογεί την εφαρμογή
# Negative Binomial, επομένως δεν εξετάζεται περαιτέρω.
#
#
pal    <- c("#4C78A8", "#F58518", "#54A24B", "#B279A2")
accent <- "#E45756"
grey1  <- "#BAB0AC"

pngw <- 2100; pngh <- 1500; pngr <- 300

safe_write <- function(x, f) {
  ok <- try(write.csv(x, f, row.names = FALSE), silent = TRUE)
  
  if (inherits(ok, "try-error")) {
    stop("Failed to write ", f,
        "- πιθανόν είναι ανοιχτό σε άλλο πρόγραμμα\n")
  } else {
    cat("Αποθηκεύτηκε:", f, "\n")
  }
}

brier <- function(p, y) mean((p - y)^2)
acc   <- function(p, y) mean((p > 0.5) == (y == 1))

nhl <- read.csv("nhl_2021_22.csv", stringsAsFactors = FALSE)

nhl$date  <- as.Date(nhl$date)
nhl$phase <- factor(nhl$phase, levels = c("Regular", "Playoffs"))

# ---------------------------------------------------
# 1. Η συνάρτηση πιθανότητας της Skellam
# ---------------------------------------------------

dskel <- function(k, l1, l2) {
  x <- 2*sqrt(l1*l2)
  
  exp(-(l1 + l2) + x) *
    (l1/l2)^(k/2) *
    besselI(x, abs(k), expon.scaled = TRUE)
}

pskel_win <- function(l1, l2, kmax = 25) {
  s <- 0
  
  for (k in 1:kmax) {
    s <- s + dskel(k, l1, l2)
  }
  
  s
}

# ---------------------------------------------------
# 2. Δόμηση δεδομένων και Double Poisson
# ---------------------------------------------------

data0 <- nhl[nhl$phase == "Regular",
             c("team1", "team2", "goals1", "goals2", "neutral")]

n <- nrow(data0)

data       <- data.frame(goals = c(data0$goals1, data0$goals2))
data$games <- rep(1:n, 2)

data$home <- c(as.numeric(!data0$neutral), rep(0, n))

data$att <- factor(c(data0$team1, data0$team2))
data$def <- factor(c(data0$team2, data0$team1))

# Χρειάζεται μόνο για να ξεχωρίζουμε τις δύο γραμμές όταν υπάρχει neutral game.
data$side <- c(rep("home", n), rep("away", n))

data <- data[order(data$games), ]

nteams <- length(levels(data$att))

contrasts(data$att) <- contr.sum(nteams)
contrasts(data$def) <- contr.sum(nteams)

colnames(contrasts(data$att)) <- levels(data$att)[-nteams]
colnames(contrasts(data$def)) <- levels(data$def)[-nteams]

model <- glm(goals ~ home + att + def, family = poisson, data = data)

cat("Double Poisson: AIC =", round(AIC(model), 1), "\n")

Xall <- model.matrix(model)

ih <- data$side == "home"
ia <- data$side == "away"

X1 <- Xall[ih, ]
X2 <- Xall[ia, ]

dif <- data0$goals1 - data0$goals2

stopifnot(all(data$games[ih] == 1:n))
stopifnot(all(data$games[ia] == 1:n))

# ---------------------------------------------------
# 3. Εκδοχή (α): post-hoc περικοπή
# ---------------------------------------------------

b0 <- coef(model)

l1 <- as.numeric(exp(X1 %*% b0))
l2 <- as.numeric(exp(X2 %*% b0))

pw <- pskel_win(l1, l2)
p0 <- dskel(0, l1, l2)

cat("\nΜέση P(D = 0) στην κανονική περίοδο:",
    round(mean(p0), 4), "\n")

cat("Παρατηρούμενο ποσοστό ισοπαλιών:",
    round(mean(data0$goals1 - data0$goals2 == 0), 4),
    "(εξ ορισμού 0)\n")

# ---------------------------------------------------
# 4. Εκδοχή (β): ML
# ---------------------------------------------------

# Αριθμητικά ασφαλέστερη μορφή της ίδιας likelihood.
negll <- function(b) {
  eta1 <- as.numeric(X1 %*% b)
  eta2 <- as.numeric(X2 %*% b)
  
  if (any(!is.finite(c(eta1, eta2)))) return(1e100)
  
  if (any(eta1 > log(100)) || any(eta2 > log(100)) ||
      any(eta1 < log(1e-8)) || any(eta2 < log(1e-8))) return(1e100)
  
  la1 <- exp(eta1)
  la2 <- exp(eta2)
  
  p  <- dskel(dif, la1, la2)
  pz <- dskel(0, la1, la2)
  
  p  <- pmax(p, 1e-300)
  pz <- pmin(pmax(pz, 0), 1 - 1e-12)
  
  value <- -sum(log(p) - log1p(-pz))
  
  if (!is.finite(value)) return(1e100)
  
  value
}

cat("\nΈναρξη βελτιστοποίησης (BFGS) ...\n")

fit <- optim(
  b0, negll,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 1000, reltol = 1e-10)
)

cat("Σύγκλιση (0 = επιτυχία):", fit$convergence, "\n")
cat("Log-likelihood =", round(-fit$value, 2), "\n")
cat("Αριθμός παραμέτρων =", length(fit$par), "\n")
cat("AIC =", round(2*fit$value + 2*length(fit$par), 1), "\n")

bz <- fit$par

# Έλεγχος Hessian πριν από τα SE.
eigH <- eigen(fit$hessian, symmetric = TRUE, only.values = TRUE)$values

cat("Ελάχιστη ιδιοτιμή Hessian:", min(eigH), "\n")

if (any(!is.finite(eigH)) || min(eigH) <= 0) {
  warning("Το Hessian δεν είναι θετικά ορισμένο. Τα standard errors χρειάζονται προσοχή.")
  V  <- NULL
  se <- rep(NA_real_, length(bz))
} else {
  V  <- solve(fit$hessian)
  se <- sqrt(pmax(diag(V), 0))
}

cat("\nΠαράγοντας έδρας:\n")

cat("  Double Poisson :", round(b0["home"], 4), "\n")

cat("  ZDT-Skellam    :", round(bz[2], 4),
    " (SE", round(se[2], 4), ")\n")

# ---------------------------------------------------
# 5. Πίνακας ικανοτήτων ZDT
# ---------------------------------------------------

index <- 1:(nteams - 1)

ab <- matrix(NA_real_, nrow = nteams, ncol = 4)

ab[index, 1] <- bz[2 + index]
ab[index, 2] <- se[2 + index]

ab[index, 3] <- bz[nteams + 1 + index]
ab[index, 4] <- se[nteams + 1 + index]

ab[nteams, 1] <- -sum(bz[2 + index])
ab[nteams, 3] <- -sum(bz[nteams + 1 + index])

# Standard errors της τελευταίας STZ ομάδας.
if (!is.null(V)) {
  Latt <- numeric(length(bz))
  Latt[2 + index] <- -1
  
  ab[nteams, 2] <- sqrt(as.numeric(t(Latt) %*% V %*% Latt))
  
  Ldef <- numeric(length(bz))
  Ldef[nteams + 1 + index] <- -1
  
  ab[nteams, 4] <- sqrt(as.numeric(t(Ldef) %*% V %*% Ldef))
}

rownames(ab) <- levels(data$att)
colnames(ab) <- c("Att", "SD-Att", "Def", "SD-Def")

cat("\n---------- ΙΚΑΝΟΤΗΤΕΣ (ZDT-Skellam) ----------\n")

print(round(ab[order(-ab[, 1]), ], 3))

abdf <- as.data.frame(ab)
abdf$team <- rownames(ab)

safe_write(abdf, "abilities_zdt.csv")

# ---------------------------------------------------
# 6. Σύγκριση ικανοτήτων DP vs ZDT
# ---------------------------------------------------

dpa <- matrix(nrow = nteams, ncol = 2)

dpa[index, 1] <- b0[2 + index]
dpa[index, 2] <- b0[nteams + 1 + index]

dpa[nteams, 1] <- -sum(b0[2 + index])
dpa[nteams, 2] <- -sum(b0[nteams + 1 + index])

cat("\nΣυσχέτιση ικανοτήτων DP vs ZDT:\n")

cat("  επίθεση :", round(cor(dpa[, 1], ab[, 1]), 4), "\n")
cat("  άμυνα   :", round(cor(dpa[, 2], ab[, 3]), 4), "\n")

png("fig14_dp_vs_zdt_abilities.png",
    width = pngw, height = pngh, res = pngr)

par(mfrow = c(1, 2))

plot(dpa[, 1], ab[, 1],
     pch = 16, col = pal[1],
     xlab = "att (Double Poisson)",
     ylab = "att (ZDT-Skellam)",
     main = "Επιθετική ικανότητα")

abline(0, 1, lty = 2, col = grey1)

plot(dpa[, 2], ab[, 3],
     pch = 16, col = pal[2],
     xlab = "def (Double Poisson)",
     ylab = "def (ZDT-Skellam)",
     main = "Αμυντική ικανότητα")

abline(0, 1, lty = 2, col = grey1)

par(mfrow = c(1, 1))

dev.off()

# ---------------------------------------------------
# 7. Προσαρμογή στη διαφορά γκολ
# ---------------------------------------------------

zl1 <- as.numeric(exp(X1 %*% bz))
zl2 <- as.numeric(exp(X2 %*% bz))

zp0 <- dskel(0, zl1, zl2)

# Χρησιμοποιούμε όλο το πραγματικό εύρος της διαφοράς γκολ.
M  <- max(abs(dif))
ks <- -M:M

fit_k <- numeric(length(ks))

for (i in seq_along(ks)) {
  if (ks[i] == 0) {
    fit_k[i] <- 0
    next
  }
  
  fit_k[i] <- mean(dskel(ks[i], zl1, zl2)/(1 - zp0))
}

obs_k <- as.numeric(table(factor(dif, levels = ks)))/n

png("fig15_skellam_fit.png",
    width = pngw, height = pngh, res = pngr)

barplot(rbind(obs_k, fit_k),
        beside = TRUE,
        names.arg = ks,
        col = c(pal[1], pal[2]),
        border = "white",
        xlab = "Διαφορά γκολ (Γηπεδούχος - Φιλοξενούμενος)",
        ylab = "Σχετική συχνότητα",
        main = "Παρατηρούμενη κατανομή vs περικομμένη Skellam")

legend("topright",
       c("Παρατηρούμενα", "ZDT-Skellam"),
       fill = c(pal[1], pal[2]),
       bty = "n")

dev.off()

cat("\nΜέγιστη απόκλιση παρατηρούμενης - προβλεπόμενης συχνότητας:",
    round(max(abs(obs_k - fit_k)), 4), "\n")
obs_n <- as.numeric(table(factor(dif, levels = ks)))
exp_n <- fit_k * n

cat("\nΈλεγχος καλής προσαρμογής (κελιά με exp >= 5):\n")
keep <- exp_n >= 5
X2   <- sum((obs_n[keep] - exp_n[keep])^2 / exp_n[keep])
cat("  Pearson X^2 =", round(X2, 1), " κελιά =", sum(keep), "\n")

cat("\nΑθροιστικές πιθανότητες (εκεί που το μοντέλο στέκει):\n")
cat("  P(D > 0)  μοντέλο =", round(mean(pskel_win(zl1,zl2)/(1-zp0)), 4),
    " παρατηρούμενο =", round(mean(dif > 0), 4), "\n")
p2 <- (pskel_win(zl1,zl2) - dskel(1, zl1, zl2))/(1 - zp0)
cat("  P(D >= 2) μοντέλο =", round(mean(p2), 4),
    " παρατηρούμενο =", round(mean(dif >= 2), 4), "\n")

cat("\nΜη μονοτονία της παρατηρούμενης κατανομής:\n")
cat("  |D| = 2 :", sum(abs(dif) == 2), "αγώνες\n")
cat("  |D| = 3 :", sum(abs(dif) == 3), "αγώνες  <- μεγαλύτερο\n")

# ---------------------------------------------------
# 8. Πρόβλεψη στα playoffs
# ---------------------------------------------------

test <- nhl[nhl$phase == "Playoffs", ]

m <- nrow(test)

newdat <- data.frame(
  home = c(as.numeric(!test$neutral), rep(0, m)),
  att  = factor(c(test$team1, test$team2), levels = levels(data$att)),
  def  = factor(c(test$team2, test$team1), levels = levels(data$def))
)

Xnew <- model.matrix(
  delete.response(terms(model)),
  newdat,
  contrasts.arg = model$contrasts
)

X1n <- Xnew[1:m, ]
X2n <- Xnew[(m + 1):(2*m), ]

# (α) DP με αναλυτική περικοπή

a1 <- as.numeric(exp(X1n %*% b0))
a2 <- as.numeric(exp(X2n %*% b0))

pa  <- pskel_win(a1, a2)
pa0 <- dskel(0, a1, a2)

p_dp <- pa/(1 - pa0)

# (β) ZDT-Skellam

z1 <- as.numeric(exp(X1n %*% bz))
z2 <- as.numeric(exp(X2n %*% bz))

pz  <- pskel_win(z1, z2)
pz0 <- dskel(0, z1, z2)

p_zdt <- pz/(1 - pz0)

mkt <- (1/test$odds1)/((1/test$odds1) + (1/test$odds2))

market_idx <- test$market_ok & is.finite(mkt)

y <- test$win1

pred <- data.frame(
  date = test$date,
  team1 = test$team1,
  team2 = test$team2,
  goals1 = test$goals1,
  goals2 = test$goals2,
  lam1 = round(z1, 3),
  lam2 = round(z2, 3),
  pDP = round(p_dp, 3),
  pZDT = round(p_zdt, 3),
  market = round(mkt, 3),
  win1 = y
)

print(head(pred, 10))

safe_write(pred, "pred_zdt_playoffs.csv")

cat("\n---------- OUT-OF-SAMPLE (playoffs, n =", m, ") ----------\n")

cat("DP + περικοπή : Brier =", round(brier(p_dp, y), 4),
    " Accuracy =", round(100*acc(p_dp, y), 1), "%\n")

cat("ZDT-Skellam   : Brier =", round(brier(p_zdt, y), 4),
    " Accuracy =", round(100*acc(p_zdt, y), 1), "%\n")

cat("\nΚοινό δείγμα μοντέλων-αγοράς, n =", sum(market_idx), "\n")

cat("DP + περικοπή : Brier =",
    round(brier(p_dp[market_idx], y[market_idx]), 4),
    " Accuracy =", round(100*acc(p_dp[market_idx], y[market_idx]), 1), "%\n")

cat("ZDT-Skellam   : Brier =",
    round(brier(p_zdt[market_idx], y[market_idx]), 4),
    " Accuracy =", round(100*acc(p_zdt[market_idx], y[market_idx]), 1), "%\n")

cat("Αγορά         : Brier =",
    round(brier(mkt[market_idx], y[market_idx]), 4),
    " Accuracy =", round(100*acc(mkt[market_idx], y[market_idx]), 1), "%\n")

png("fig16_zdt_vs_market.png",
    width = pngw, height = pngh, res = pngr)

plot(mkt[market_idx], p_zdt[market_idx],
     pch = 16, col = pal[1],
     xlim = c(0,1), ylim = c(0,1),
     xlab = "Πιθανότητα αγοράς (fair)",
     ylab = "Πιθανότητα ZDT-Skellam",
     main = "ZDT-Skellam vs αγορά - αγώνες playoffs")

abline(0, 1, lty = 2, col = grey1)

dev.off()

# ---------------------------------------------------
# 9. Έλεγχος ευαισθησίας: χρονολογικό 70/30
# ---------------------------------------------------

reg <- nhl[nhl$phase == "Regular", ]
reg <- reg[order(reg$date), ]

k70 <- floor(0.7*nrow(reg))

tr <- reg[1:k70, ]
te <- reg[(k70 + 1):nrow(reg), ]

cat("\nΈλεγχος ευαισθησίας - εκπαίδευση:",
    nrow(tr), " έλεγχος:", nrow(te), "\n")

n2 <- nrow(tr)

d2       <- data.frame(goals = c(tr$goals1, tr$goals2))
d2$games <- rep(1:n2, 2)
d2$home  <- c(as.numeric(!tr$neutral), rep(0, n2))
d2$att   <- factor(c(tr$team1, tr$team2))
d2$def   <- factor(c(tr$team2, tr$team1))
d2$side  <- c(rep("home", n2), rep("away", n2))

d2 <- d2[order(d2$games), ]

nt2 <- length(levels(d2$att))

contrasts(d2$att) <- contr.sum(nt2)
contrasts(d2$def) <- contr.sum(nt2)

colnames(contrasts(d2$att)) <- levels(d2$att)[-nt2]
colnames(contrasts(d2$def)) <- levels(d2$def)[-nt2]

model2 <- glm(goals ~ home + att + def,
              family = poisson, data = d2)

Xa <- model.matrix(model2)

ih2 <- d2$side == "home"
ia2 <- d2$side == "away"

X1b <- Xa[ih2, ]
X2b <- Xa[ia2, ]

dif2 <- tr$goals1 - tr$goals2

negll2 <- function(b) {
  eta1 <- as.numeric(X1b %*% b)
  eta2 <- as.numeric(X2b %*% b)
  
  if (any(!is.finite(c(eta1, eta2)))) return(1e100)
  
  if (any(eta1 > log(100)) || any(eta2 > log(100)) ||
      any(eta1 < log(1e-8)) || any(eta2 < log(1e-8))) return(1e100)
  
  la1 <- exp(eta1)
  la2 <- exp(eta2)
  
  p  <- pmax(dskel(dif2, la1, la2), 1e-300)
  pz <- dskel(0, la1, la2)
  
  pz <- pmin(pmax(pz, 0), 1 - 1e-12)
  
  value <- -sum(log(p) - log1p(-pz))
  
  if (!is.finite(value)) return(1e100)
  
  value
}

fit2 <- optim(
  coef(model2),
  negll2,
  method = "BFGS",
  control = list(maxit = 1000, reltol = 1e-10)
)

cat("Σύγκλιση:", fit2$convergence, "\n")

m2 <- nrow(te)

nd2 <- data.frame(
  home = c(as.numeric(!te$neutral), rep(0, m2)),
  att  = factor(c(te$team1, te$team2), levels = levels(d2$att)),
  def  = factor(c(te$team2, te$team1), levels = levels(d2$def))
)

Xn2 <- model.matrix(
  delete.response(terms(model2)),
  nd2,
  contrasts.arg = model2$contrasts
)

q1 <- as.numeric(exp(Xn2[1:m2, ] %*% fit2$par))
q2 <- as.numeric(exp(Xn2[(m2 + 1):(2*m2), ] %*% fit2$par))

p_zdt2 <- pskel_win(q1, q2)/(1 - dskel(0, q1, q2))

r1 <- as.numeric(exp(Xn2[1:m2, ] %*% coef(model2)))
r2 <- as.numeric(exp(Xn2[(m2 + 1):(2*m2), ] %*% coef(model2)))

p_dp2 <- pskel_win(r1, r2)/(1 - dskel(0, r1, r2))

mkt2 <- (1/te$odds1)/((1/te$odds1) + (1/te$odds2))

market_idx2 <- te$market_ok & is.finite(mkt2)

cat("DP + περικοπή : Brier =",
    round(brier(p_dp2[market_idx2], te$win1[market_idx2]), 4),
    " Accuracy =",
    round(100*acc(p_dp2[market_idx2], te$win1[market_idx2]), 1), "%\n")

cat("ZDT-Skellam   : Brier =",
    round(brier(p_zdt2[market_idx2], te$win1[market_idx2]), 4),
    " Accuracy =",
    round(100*acc(p_zdt2[market_idx2], te$win1[market_idx2]), 1), "%\n")

cat("Αγορά         : Brier =",
    round(brier(mkt2[market_idx2], te$win1[market_idx2]), 4),
    " Accuracy =",
    round(100*acc(mkt2[market_idx2], te$win1[market_idx2]), 1), "%\n")

cat("\n===================================================\n")
cat(" Ολοκληρώθηκε το μοντέλο ZDT-Skellam\n")
cat("===================================================\n")