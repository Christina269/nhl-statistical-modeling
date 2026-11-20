
#4
#Nhl 04 value betting · R
#
#
# NHL 2021-22 : ΣΤΡΑΤΗΓΙΚΗ ΣΤΟΙΧΗΜΑΤΩΝ ΑΞΙΑΣ (VALUE BETTING)
#
# Λογική: τοποθετείται στοίχημα μίας μονάδας όταν η αναμενόμενη αξία
# είναι θετική, δηλαδή όταν
#
#     p_model * odds - 1 > threshold
#
# Προσοχή: η σύγκριση γίνεται με τις ΠΡΑΓΜΑΤΙΚΕΣ αποδόσεις (odds) και όχι
# με τις fair, γιατί αυτό που πληρώνεται στην πράξη είναι η απόδοση με το
# περιθώριο μέσα. Επειδή το overround είναι θετικό, δεν μπορούν να έχουν
# ταυτόχρονα θετική αξία και οι δύο πλευρές του ίδιου αγώνα.
#
# Η αξιολόγηση γίνεται ΜΟΝΟ out-of-sample, σε δύο χωρισμούς:
#   (Α) εκπαίδευση = κανονική περίοδος, έλεγχος = playoffs (89 αγώνες)
#   (Β) χρονολογικό 70/30 μέσα στην κανονική περίοδο (394 αγώνες)
#
# Είσοδος : nhl_2021_22.csv
# Έξοδος  : fig18 - fig20 png, value_bets.csv, value_summary.csv
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

nhl <- read.csv("nhl_2021_22.csv", stringsAsFactors = FALSE)

nhl$date  <- as.Date(nhl$date)
nhl$phase <- factor(nhl$phase, levels = c("Regular", "Playoffs"))

# ---------------------------------------------------
# 1. Συναρτήσεις Skellam (ίδιες με το script 03)
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
# 2. Προσαρμογή μοντέλων σε ένα σύνολο εκπαίδευσης
# ---------------------------------------------------
# Επιστρέφει ό,τι χρειάζεται για πρόβλεψη: το glm, τους συντελεστές του
# Double Poisson και τους συντελεστές του ZDT-Skellam.

fit_models <- function(tr) {
  n <- nrow(tr)
  
  d       <- data.frame(goals = c(tr$goals1, tr$goals2))
  d$games <- rep(1:n, 2)
  d$home  <- c(as.numeric(!tr$neutral), rep(0, n))
  d$att   <- factor(c(tr$team1, tr$team2))
  d$def   <- factor(c(tr$team2, tr$team1))
  d$side  <- c(rep("home", n), rep("away", n))
  
  d <- d[order(d$games), ]
  
  nt <- length(levels(d$att))
  
  contrasts(d$att) <- contr.sum(nt)
  contrasts(d$def) <- contr.sum(nt)
  
  colnames(contrasts(d$att)) <- levels(d$att)[-nt]
  colnames(contrasts(d$def)) <- levels(d$def)[-nt]
  
  mod <- glm(goals ~ home + att + def, family = poisson, data = d)
  
  X   <- model.matrix(mod)
  X1  <- X[d$side == "home", ]
  X2  <- X[d$side == "away", ]
  dif <- tr$goals1 - tr$goals2
  
  negll <- function(b) {
    eta1 <- as.numeric(X1 %*% b)
    eta2 <- as.numeric(X2 %*% b)
    
    if (any(!is.finite(c(eta1, eta2)))) return(1e100)
    
    if (any(eta1 > log(100)) || any(eta2 > log(100)) ||
        any(eta1 < log(1e-8)) || any(eta2 < log(1e-8))) return(1e100)
    
    la1 <- exp(eta1)
    la2 <- exp(eta2)
    
    p  <- pmax(dskel(dif, la1, la2), 1e-300)
    pz <- pmin(pmax(dskel(0, la1, la2), 0), 1 - 1e-12)
    
    v <- -sum(log(p) - log1p(-pz))
    
    if (!is.finite(v)) return(1e100)
    
    v
  }
  
  ft <- optim(coef(mod), negll, method = "BFGS",
              control = list(maxit = 1000, reltol = 1e-10))
  
  list(mod = mod, bdp = coef(mod), bzdt = ft$par,
       conv = ft$convergence, levs = levels(d$att))
}

# Πιθανότητα νίκης γηπεδούχου για ένα σύνολο ελέγχου
predict_p <- function(f, te, which = "zdt") {
  m <- nrow(te)
  
  nd <- data.frame(
    home = c(as.numeric(!te$neutral), rep(0, m)),
    att  = factor(c(te$team1, te$team2), levels = f$levs),
    def  = factor(c(te$team2, te$team1), levels = f$levs)
  )
  
  X <- model.matrix(delete.response(terms(f$mod)), nd,
                    contrasts.arg = f$mod$contrasts)
  
  b <- if (which == "zdt") f$bzdt else f$bdp
  
  l1 <- as.numeric(exp(X[1:m, ] %*% b))
  l2 <- as.numeric(exp(X[(m + 1):(2*m), ] %*% b))
  
  pskel_win(l1, l2)/(1 - dskel(0, l1, l2))
}

# ---------------------------------------------------
# 3. Η στρατηγική στοιχημάτων αξίας
# ---------------------------------------------------
# p    : πιθανότητα νίκης γηπεδούχου από το μοντέλο
# te   : το σύνολο ελέγχου (χρειάζεται odds1, odds2, win1, market_ok)
# thr  : ελάχιστη αναμενόμενη αξία για να μπει το στοίχημα

value_bets <- function(p, te, thr = 0) {
  ok <- te$market_ok & is.finite(te$odds1) & is.finite(te$odds2)
  
  ev1 <- p*te$odds1 - 1              # αξία στον γηπεδούχο
  ev2 <- (1 - p)*te$odds2 - 1        # αξία στον φιλοξενούμενο
  
  b1 <- ok & ev1 > thr
  b2 <- ok & ev2 > thr
  
  side   <- rep(NA_character_, nrow(te))
  side[b1] <- "Home"
  side[b2] <- "Away"
  
  bet <- !is.na(side)
  
  odds <- ifelse(side == "Home", te$odds1, te$odds2)
  won  <- ifelse(side == "Home", te$win1 == 1, te$win1 == 0)
  
  profit <- ifelse(won, odds - 1, -1)
  
  data.frame(
    date   = te$date[bet],
    team1  = te$team1[bet],
    team2  = te$team2[bet],
    side   = side[bet],
    p      = round(ifelse(side[bet] == "Home", p[bet], 1 - p[bet]), 4),
    odds   = round(odds[bet], 3),
    ev     = round(ifelse(side[bet] == "Home", ev1[bet], ev2[bet]), 4),
    won    = as.numeric(won[bet]),
    profit = round(profit[bet], 3),
    stringsAsFactors = FALSE
  )
}

# Σύνοψη μιας στρατηγικής
summarise_bets <- function(b, label) {
  if (nrow(b) == 0) {
    cat(label, ": κανένα στοίχημα\n")
    return(data.frame(strategy = label, bets = 0, hit = NA, roi = NA,
                      profit = NA, p.value = NA))
  }
  
  tt <- if (nrow(b) > 2) t.test(b$profit, mu = 0) else NULL
  
  cat(sprintf("%-26s n=%3d  hit=%5.1f%%  μ.απόδοση=%5.2f  κέρδος=%7.2f  ROI=%6.2f%%  p=%s\n",
              label, nrow(b), 100*mean(b$won), mean(b$odds),
              sum(b$profit), 100*mean(b$profit),
              if (is.null(tt)) "-" else format(round(tt$p.value, 3), nsmall = 3)))
  
  data.frame(strategy = label,
             bets     = nrow(b),
             hit      = 100*mean(b$won),
             roi      = 100*mean(b$profit),
             profit   = sum(b$profit),
             p.value  = if (is.null(tt)) NA else tt$p.value)
}

# ---------------------------------------------------
# 4. Χωρισμός (Α): κανονική περίοδος -> playoffs
# ---------------------------------------------------

trA <- nhl[nhl$phase == "Regular", ]
teA <- nhl[nhl$phase == "Playoffs", ]

cat("=== ΧΩΡΙΣΜΟΣ Α: εκπαίδευση", nrow(trA), " έλεγχος", nrow(teA), "===\n")

fA <- fit_models(trA)
cat("Σύγκλιση ZDT:", fA$conv, "\n\n")

pA_dp  <- predict_p(fA, teA, "dp")
pA_zdt <- predict_p(fA, teA, "zdt")

resA <- rbind(
  summarise_bets(value_bets(pA_dp,  teA, 0.00), "A / DP    / thr=0%"),
  summarise_bets(value_bets(pA_zdt, teA, 0.00), "A / ZDT   / thr=0%"),
  summarise_bets(value_bets(pA_zdt, teA, 0.02), "A / ZDT   / thr=2%"),
  summarise_bets(value_bets(pA_zdt, teA, 0.05), "A / ZDT   / thr=5%")
)

# ---------------------------------------------------
# 5. Χωρισμός (Β): χρονολογικό 70/30 στην κανονική περίοδο
# ---------------------------------------------------

reg <- nhl[nhl$phase == "Regular", ]
reg <- reg[order(reg$date), ]

k70 <- floor(0.7*nrow(reg))
trB <- reg[1:k70, ]
teB <- reg[(k70 + 1):nrow(reg), ]

cat("\n=== ΧΩΡΙΣΜΟΣ Β: εκπαίδευση", nrow(trB), " έλεγχος", nrow(teB), "===\n")

fB <- fit_models(trB)
cat("Σύγκλιση ZDT:", fB$conv, "\n\n")

pB_dp  <- predict_p(fB, teB, "dp")
pB_zdt <- predict_p(fB, teB, "zdt")

resB <- rbind(
  summarise_bets(value_bets(pB_dp,  teB, 0.00), "B / DP    / thr=0%"),
  summarise_bets(value_bets(pB_zdt, teB, 0.00), "B / ZDT   / thr=0%"),
  summarise_bets(value_bets(pB_zdt, teB, 0.02), "B / ZDT   / thr=2%"),
  summarise_bets(value_bets(pB_zdt, teB, 0.05), "B / ZDT   / thr=5%")
)

res <- rbind(resA, resB)

cat("\n---------- ΣΥΓΚΕΝΤΡΩΤΙΚΟΣ ΠΙΝΑΚΑΣ ----------\n")
print(res, digits = 3, row.names = FALSE)

safe_write(res, "value_summary.csv")

# ---------------------------------------------------
# 6. Ανάλυση της βασικής στρατηγικής (Β / ZDT / thr=0)
# ---------------------------------------------------
# Ο χωρισμός Β έχει το μεγαλύτερο σύνολο ελέγχου, οπότε είναι ο πιο
# αξιόπιστος για συμπεράσματα.

bets <- value_bets(pB_zdt, teB, 0.00)

cat("\n---------- ΑΝΑΛΥΣΗ ΒΑΣΙΚΗΣ ΣΤΡΑΤΗΓΙΚΗΣ ----------\n")
cat("Στοιχήματα:", nrow(bets), "από", nrow(teB), "αγώνες  (",
    round(100*nrow(bets)/nrow(teB), 1), "%)\n")

print(table(bets$side))

cat("Μέση απόδοση:", round(mean(bets$odds), 3), "\n")
cat("Μέση αναμενόμενη αξία:", round(100*mean(bets$ev), 2), "%\n")
cat("Ποσοστό επιτυχίας:", round(100*mean(bets$won), 1), "%\n")
cat("Συνολικό κέρδος:", round(sum(bets$profit), 2), "μονάδες\n")
cat("ROI:", round(100*mean(bets$profit), 2), "%\n")

tt <- t.test(bets$profit, mu = 0)

cat("\nΈλεγχος στατιστικής σημαντικότητας (one-sample t-test):\n")
cat("  t =", round(tt$statistic, 3),
    " df =", tt$parameter,
    " p-value =", round(tt$p.value, 4), "\n")
cat("  95% ΔΕ μέσου κέρδους ανά στοίχημα: [",
    round(tt$conf.int[1], 3), ",", round(tt$conf.int[2], 3), "]\n")

if (tt$p.value >= 0.05) {
  cat("  -> Το κέρδος ΔΕΝ είναι στατιστικά σημαντικό.\n")
} else {
  cat("  -> Το κέρδος είναι στατιστικά σημαντικό.\n")
}

safe_write(bets, "value_bets.csv")
cat("\n---------- ΑΝΑΛΥΣΗ ΑΝΑ ΕΥΡΟΣ ΑΠΟΔΟΣΕΩΝ ----------\n")

grp <- cut(bets$odds, breaks = c(1, 1.8, 2.2, 3.0, 100),
           labels = c("<1.80", "1.80-2.20", "2.20-3.00", ">3.00"))

brk <- data.frame(
  bucket   = levels(grp),
  n        = as.numeric(table(grp)),
  meanOdds = round(tapply(bets$odds,   grp, mean), 2),
  breakeven= round(100/tapply(bets$odds, grp, mean), 1),
  hit      = round(100*tapply(bets$won, grp, mean), 1),
  roi      = round(100*tapply(bets$profit, grp, mean), 1)
)
print(brk, row.names = FALSE)

cat("\nΜέση απόδοση κερδισμένων :", round(mean(bets$odds[bets$won == 1]), 3), "\n")
cat("Μέση απόδοση χαμένων     :", round(mean(bets$odds[bets$won == 0]), 3), "\n")

cat("\n---------- ΑΝΑΛΥΣΗ ΑΝΑ ΠΛΕΥΡΑ ----------\n")
sd2 <- data.frame(
  side = names(table(bets$side)),
  n    = as.numeric(table(bets$side)),
  hit  = round(100*tapply(bets$won, bets$side, mean), 1),
  roi  = round(100*tapply(bets$profit, bets$side, mean), 1)
)
print(sd2, row.names = FALSE)

# ---------------------------------------------------
# 7. Καμπύλη συσσωρευμένου κέρδους
# ---------------------------------------------------

cum <- cumsum(bets$profit)

png("fig18_pl_curve.png", width = pngw, height = pngh, res = pngr)

plot(1:length(cum), cum, type = "l", lwd = 2, col = pal[3],
     xlab = "Στοιχήματα (χρονολογική σειρά)",
     ylab = "Συσσωρευμένο καθαρό κέρδος (μονάδες)",
     main = "Καμπύλη P&L - στοιχήματα αξίας (ZDT-Skellam, out-of-sample)")

abline(h = 0, col = accent, lty = 2, lwd = 2)

dev.off()

# ---------------------------------------------------
# 8. Ευαισθησία ως προς το κατώφλι
# ---------------------------------------------------

thrs <- seq(0, 0.15, by = 0.01)

sens <- data.frame(thr = thrs, bets = NA, roi = NA, p.value = NA)

for (i in seq_along(thrs)) {
  b <- value_bets(pB_zdt, teB, thrs[i])
  
  sens$bets[i] <- nrow(b)
  
  if (nrow(b) > 2) {
    sens$roi[i]     <- 100*mean(b$profit)
    sens$p.value[i] <- t.test(b$profit, mu = 0)$p.value
  }
}

print(sens, digits = 3, row.names = FALSE)

png("fig19_threshold_sensitivity.png", width = pngw, height = pngh, res = pngr)

par(mfrow = c(1, 2))

plot(100*sens$thr, sens$bets, type = "b", pch = 16, col = pal[1],
     xlab = "Κατώφλι αναμενόμενης αξίας (%)",
     ylab = "Πλήθος στοιχημάτων",
     main = "Πόσα στοιχήματα επιλέγονται")

plot(100*sens$thr, sens$roi, type = "b", pch = 16, col = pal[2],
     xlab = "Κατώφλι αναμενόμενης αξίας (%)",
     ylab = "ROI (%)",
     main = "Απόδοση ανά κατώφλι")

abline(h = 0, col = accent, lty = 2, lwd = 2)

par(mfrow = c(1, 1))

dev.off()

# ---------------------------------------------------
# 9. Έλεγχος αναφοράς: τυχαία επιλογή στοιχημάτων
# ---------------------------------------------------
# Πόσο κέρδος θα έδινε μια στρατηγική που ποντάρει στον ίδιο αριθμό
# αγώνων αλλά χωρίς κριτήριο αξίας; Δείχνει αν η επιλογή του μοντέλου
# προσθέτει κάτι πέρα από την τύχη.

set.seed(1)

B    <- 2000
nbet <- nrow(bets)
ok   <- which(teB$market_ok)

sim <- numeric(B)

for (i in 1:B) {
  idx  <- sample(ok, nbet, replace = FALSE)
  side <- sample(c("Home", "Away"), nbet, replace = TRUE)
  
  o <- ifelse(side == "Home", teB$odds1[idx], teB$odds2[idx])
  w <- ifelse(side == "Home", teB$win1[idx] == 1, teB$win1[idx] == 0)
  
  sim[i] <- mean(ifelse(w, o - 1, -1))
}

cat("\nΤυχαία στρατηγική ίδιου μεγέθους (", B, "επαναλήψεις ):\n")
cat("  μέσο ROI =", round(100*mean(sim), 2),
    "%  τυπική απόκλιση =", round(100*sd(sim), 2), "%\n")
cat("  Ποσοστό τυχαίων στρατηγικών με ROI >= το δικό μας:",
    round(100*mean(sim >= mean(bets$profit)), 1), "%\n")

png("fig20_random_benchmark.png", width = pngw, height = pngh, res = pngr)

hist(100*sim, breaks = 40, col = pal[1], border = "white",
     xlab = "ROI (%)", ylab = "Συχνότητα",
     main = "Στρατηγική αξίας έναντι τυχαίας επιλογής")

abline(v = 100*mean(bets$profit), col = accent, lwd = 3)

legend("topright",
       c("Τυχαίες στρατηγικές", "Στρατηγική αξίας"),
       fill = c(pal[1], NA), border = c("white", NA),
       lty = c(NA, 1), lwd = c(NA, 3), col = c(NA, accent), bty = "n")

dev.off()

cat("\n===================================================\n")
cat(" Ολοκληρώθηκε η ανάλυση στοιχημάτων αξίας\n")
cat("===================================================\n")
cat("\n---------- ΑΠΛΟΪΚΕΣ ΣΤΡΑΤΗΓΙΚΕΣ ΑΝΑΦΟΡΑΣ ----------\n")

okB <- teB$market_ok
o1  <- teB$odds1[okB]; o2 <- teB$odds2[okB]; w <- teB$win1[okB]

naive <- function(side_home, label) {
  o <- ifelse(side_home, o1, o2)
  won <- ifelse(side_home, w == 1, w == 0)
  pr <- ifelse(won, o - 1, -1)
  cat(sprintf("%-22s n=%3d  hit=%5.1f%%  ROI=%6.2f%%\n",
              label, length(pr), 100*mean(won), 100*mean(pr)))
}

naive(rep(TRUE,  sum(okB)), "Πάντα γηπεδούχος")
naive(rep(FALSE, sum(okB)), "Πάντα φιλοξενούμενος")
naive(o1 <  o2, "Πάντα φαβορί")
naive(o1 >= o2, "Πάντα αουτσάιντερ")

# Bootstrap διάστημα για το ROI (πιο αξιόπιστο από το t-test
# επειδή η κατανομή του κέρδους είναι έντονα ασύμμετρη)
set.seed(2)
bs <- replicate(5000, mean(sample(bets$profit, length(bets$profit), replace = TRUE)))
cat("\nBootstrap 95% ΔΕ για το ROI: [",
    round(100*quantile(bs, 0.025), 2), ",",
    round(100*quantile(bs, 0.975), 2), "] %\n")
