#5
#
#
# NHL 2021-22 : ΣΥΓΚΡΙΤΙΚΗ ΑΞΙΟΛΟΓΗΣΗ ΚΑΙ ΣΥΜΠΕΡΑΣΜΑΤΑ
#
# Συγκεντρώνει σε έναν πίνακα όλα τα μοντέλα και την αγορά, με τρεις
# μετρικές και στους δύο χωρισμούς, in-sample και out-of-sample:
#
#   Brier Score : μέσο τετραγωνικό σφάλμα πιθανότητας (όσο μικρότερο τόσο
#                 καλύτερη η βαθμονόμηση)
#   Accuracy    : ποσοστό σωστής πρόβλεψης νικητή
#   Log Loss    : λογαριθμική ποινή, τιμωρεί έντονα τις σίγουρες λάθος
#                 προβλέψεις
#
# Μοντέλα που συγκρίνονται:
#   Αγορά                   : fair πιθανότητες από τις closing odds
#   ZDT-Skellam             : ML στην περικομμένη πιθανοφάνεια
#   Double Poisson (trunc)  : lambda του DP με αναλυτική περικοπή
#   Double Poisson (raw)    : χωρίς αναδιανομή της μάζας στο D = 0
#   Σταθερή 0.50            : μηδενικής πληροφορίας
#   Σταθερή = ποσοστό έδρας : μόνο πλεονέκτημα έδρας
#
# Είσοδος : nhl_2021_22.csv
# Έξοδος  : fig21 - fig22 png, comparison_table.csv
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
# 1. Συναρτήσεις Skellam και μετρικές
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

brier   <- function(p, y) mean((p - y)^2)
acc     <- function(p, y) mean((p > 0.5) == (y == 1))

logloss <- function(p, y) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  -mean(y*log(p) + (1 - y)*log(1 - p))
}

# ---------------------------------------------------
# 2. Προσαρμογή και πρόβλεψη (ίδιες με το script 04)
# ---------------------------------------------------

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

# Επιστρέφει και τις τρεις εκδοχές πιθανότητας νίκης γηπεδούχου
predict_all <- function(f, te) {
  m <- nrow(te)
  
  nd <- data.frame(
    home = c(as.numeric(!te$neutral), rep(0, m)),
    att  = factor(c(te$team1, te$team2), levels = f$levs),
    def  = factor(c(te$team2, te$team1), levels = f$levs)
  )
  
  X <- model.matrix(delete.response(terms(f$mod)), nd,
                    contrasts.arg = f$mod$contrasts)
  
  X1 <- X[1:m, ]
  X2 <- X[(m + 1):(2*m), ]
  
  a1 <- as.numeric(exp(X1 %*% f$bdp))
  a2 <- as.numeric(exp(X2 %*% f$bdp))
  
  z1 <- as.numeric(exp(X1 %*% f$bzdt))
  z2 <- as.numeric(exp(X2 %*% f$bzdt))
  
  list(
    dp_raw = pskel_win(a1, a2),
    dp_trc = pskel_win(a1, a2)/(1 - dskel(0, a1, a2)),
    zdt    = pskel_win(z1, z2)/(1 - dskel(0, z1, z2))
  )
}

# ---------------------------------------------------
# 3. Αξιολόγηση ενός συνόλου
# ---------------------------------------------------

evaluate <- function(f, dset, label, homerate) {
  ok <- dset$market_ok
  y  <- dset$win1
  p  <- predict_all(f, dset)
  
  mkt <- (1/dset$odds1)/((1/dset$odds1) + (1/dset$odds2))
  
  # Η σύγκριση με την αγορά γίνεται στο κοινό δείγμα.
  rows <- list(
    c("Αγορά (closing odds)",   brier(mkt[ok], y[ok]),
      100*acc(mkt[ok], y[ok]),  logloss(mkt[ok], y[ok])),
    c("ZDT-Skellam",            brier(p$zdt[ok], y[ok]),
      100*acc(p$zdt[ok], y[ok]), logloss(p$zdt[ok], y[ok])),
    c("Double Poisson (trunc)", brier(p$dp_trc[ok], y[ok]),
      100*acc(p$dp_trc[ok], y[ok]), logloss(p$dp_trc[ok], y[ok])),
    c("Double Poisson (raw)",   brier(p$dp_raw[ok], y[ok]),
      100*acc(p$dp_raw[ok], y[ok]), logloss(p$dp_raw[ok], y[ok])),
    c("Σταθερή 0.50",           brier(rep(0.5, sum(ok)), y[ok]),
      100*acc(rep(0.5, sum(ok)), y[ok]), logloss(rep(0.5, sum(ok)), y[ok])),
    c("Σταθερή = % έδρας",      brier(rep(homerate, sum(ok)), y[ok]),
      100*acc(rep(homerate, sum(ok)), y[ok]),
      logloss(rep(homerate, sum(ok)), y[ok]))
  )
  
  out <- data.frame(
    set     = label,
    n       = sum(ok),
    model   = sapply(rows, function(r) r[1]),
    Brier   = as.numeric(sapply(rows, function(r) r[2])),
    Accuracy= as.numeric(sapply(rows, function(r) r[3])),
    LogLoss = as.numeric(sapply(rows, function(r) r[4])),
    stringsAsFactors = FALSE
  )
  
  out
}

# ---------------------------------------------------
# 4. Χωρισμός Α: κανονική περίοδος -> playoffs
# ---------------------------------------------------

trA <- nhl[nhl$phase == "Regular", ]
teA <- nhl[nhl$phase == "Playoffs", ]

cat("=== ΧΩΡΙΣΜΟΣ Α ===\n")
cat("Εκπαίδευση:", nrow(trA), " Έλεγχος:", nrow(teA), "\n")

fA <- fit_models(trA)
cat("Σύγκλιση ZDT:", fA$conv, "\n")

hrA <- mean(trA$win1)   # ποσοστό έδρας από την εκπαίδευση

rA_in  <- evaluate(fA, trA, "A in-sample",      hrA)
rA_out <- evaluate(fA, teA, "A out-of-sample",  hrA)

# ---------------------------------------------------
# 5. Χωρισμός Β: χρονολογικό 70/30
# ---------------------------------------------------

reg <- nhl[nhl$phase == "Regular", ]
reg <- reg[order(reg$date), ]

k70 <- floor(0.7*nrow(reg))
trB <- reg[1:k70, ]
teB <- reg[(k70 + 1):nrow(reg), ]

cat("\n=== ΧΩΡΙΣΜΟΣ Β ===\n")
cat("Εκπαίδευση:", nrow(trB), " Έλεγχος:", nrow(teB), "\n")

fB <- fit_models(trB)
cat("Σύγκλιση ZDT:", fB$conv, "\n")

hrB <- mean(trB$win1)

rB_in  <- evaluate(fB, trB, "B in-sample",     hrB)
rB_out <- evaluate(fB, teB, "B out-of-sample", hrB)

# ---------------------------------------------------
# 6. Ο τελικός πίνακας
# ---------------------------------------------------

tab <- rbind(rA_in, rA_out, rB_in, rB_out)

tab$Brier    <- round(tab$Brier, 4)
tab$Accuracy <- round(tab$Accuracy, 1)
tab$LogLoss  <- round(tab$LogLoss, 4)

cat("\n================ ΤΕΛΙΚΟΣ ΠΙΝΑΚΑΣ ================\n")

for (s in unique(tab$set)) {
  cat("\n---", s, " (n =", tab$n[tab$set == s][1], ")---\n")
  print(tab[tab$set == s, c("model", "Brier", "Accuracy", "LogLoss")],
        row.names = FALSE)
}

safe_write(tab, "comparison_table.csv")

# ---------------------------------------------------
# 7. Κατάταξη των μοντέλων out-of-sample
# ---------------------------------------------------
# Κρατάμε μόνο τα ουσιαστικά μοντέλα (χωρίς τις σταθερές αναφοράς)
# και τα κατατάσσουμε με βάση το Brier στο μεγαλύτερο σύνολο ελέγχου.

main <- c("Αγορά (closing odds)", "ZDT-Skellam",
          "Double Poisson (trunc)", "Double Poisson (raw)")

rk <- tab[tab$set == "B out-of-sample" & tab$model %in% main, ]
rk <- rk[order(rk$Brier), ]

cat("\nΚατάταξη out-of-sample (χωρισμός Β, n =", rk$n[1], "):\n")
print(rk[, c("model", "Brier", "Accuracy", "LogLoss")], row.names = FALSE)

# ---------------------------------------------------
# 8. Γράφημα σύγκρισης μετρικών
# ---------------------------------------------------

sub <- tab[tab$set %in% c("A out-of-sample", "B out-of-sample") &
             tab$model %in% main, ]

mB <- matrix(sub$Brier, nrow = length(main),
             dimnames = list(main, c("Playoffs", "Κανονική 30%")))

png("fig21_model_comparison.png", width = pngw, height = pngh, res = pngr)

par(mar = c(5, 4, 4, 8))

b <- barplot(t(mB), beside = TRUE, col = c(pal[1], pal[2]), border = "white",
             ylim = c(0, 1.25*max(mB)),
             ylab = "Brier Score (χαμηλότερο = καλύτερο)",
             main = "Σύγκριση μοντέλων out-of-sample",
             names.arg = c("Αγορά", "ZDT", "DP trunc", "DP raw"))

legend("topleft", colnames(mB), fill = c(pal[1], pal[2]), bty = "n")

par(mar = c(5, 4, 4, 2))

dev.off()

# ---------------------------------------------------
# 9. Βαθμονόμηση: μοντέλο έναντι αγοράς
# ---------------------------------------------------

pB  <- predict_all(fB, teB)
okB <- teB$market_ok
yB  <- teB$win1
mkB <- (1/teB$odds1)/((1/teB$odds1) + (1/teB$odds2))

calib <- function(p, y, k = 6) {
  br <- quantile(p, probs = seq(0, 1, length.out = k + 1))
  br[1] <- br[1] - 1e-9
  
  g <- cut(p, breaks = br)
  
  data.frame(pred = as.numeric(tapply(p, g, mean)),
             obs  = as.numeric(tapply(y, g, mean)),
             n    = as.numeric(table(g)))
}

cM <- calib(mkB[okB],      yB[okB])
cZ <- calib(pB$zdt[okB],   yB[okB])

cat("\nΒαθμονόμηση αγοράς (ίσου πλήθους ομάδες):\n")
print(round(cM, 3), row.names = FALSE)

cat("\nΒαθμονόμηση ZDT-Skellam:\n")
print(round(cZ, 3), row.names = FALSE)

png("fig22_calibration_compare.png", width = 1800, height = 1600, res = pngr)

plot(cM$pred, cM$obs, type = "b", pch = 16, col = pal[1], lwd = 2,
     xlim = c(0.2, 0.85), ylim = c(0.2, 0.85),
     xlab = "Προβλεπόμενη πιθανότητα νίκης γηπεδούχου",
     ylab = "Παρατηρούμενο ποσοστό νικών",
     main = "Βαθμονόμηση: αγορά έναντι ZDT-Skellam")

lines(cZ$pred, cZ$obs, type = "b", pch = 17, col = pal[2], lwd = 2)

abline(0, 1, lty = 2, col = grey1)

legend("topleft", c("Αγορά", "ZDT-Skellam"),
       col = c(pal[1], pal[2]), pch = c(16, 17), lwd = 2, bty = "n")

dev.off()

# ---------------------------------------------------
# 10. Σύνοψη ευρημάτων
# ---------------------------------------------------

bm <- tab$Brier[tab$set == "B out-of-sample" & tab$model == "Αγορά (closing odds)"]
bz <- tab$Brier[tab$set == "B out-of-sample" & tab$model == "ZDT-Skellam"]
b0 <- tab$Brier[tab$set == "B out-of-sample" & tab$model == "Σταθερή 0.50"]

cat("\n================ ΣΥΝΟΨΗ ================\n")

cat("Βελτίωση έναντι της σταθερής 0.50 (skill score):\n")
cat("  Αγορά       :", round(100*(1 - bm/b0), 2), "%\n")
cat("  ZDT-Skellam :", round(100*(1 - bz/b0), 2), "%\n")

cat("\nΔιαφορά μοντέλου - αγοράς σε Brier:", round(bz - bm, 4),
    " (θετικό = η αγορά καλύτερη)\n")

cat("\nΤο μερίδιο της απόδοσης του μοντέλου που εξηγείται από την αγορά:\n")
cat("  ", round(100*(1 - bz/b0)/(1 - bm/b0), 1), "% της βελτίωσης της αγοράς\n")

cat("\n===================================================\n")
cat(" Ολοκληρώθηκε η συγκριτική αξιολόγηση\n")
cat("===================================================\n")