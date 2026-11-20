#2
# NHL 2021-22 : ΜΟΝΤΕΛΟ DOUBLE POISSON
#
# Ακολουθεί τη δομή των σημειώσεων του μαθήματος:
#   1_DP_EngData.R        -> δόμηση δεδομένων, STZ contrasts, glm, ικανότητες
#   2_ConfInt_Errorbars.r -> διαστήματα εμπιστοσύνης με glht, error bars
#   3_prediction.r        -> πρόβλεψη με προσομοίωση rpois και πιθανότητες
#
# Είσοδος : nhl_2021_22.csv
# Έξοδος  : fig10 - fig13 png, abilities_dp.csv, pred_dp_playoffs.csv
#
# Χωρισμός: εκπαίδευση = κανονική περίοδος
#           έλεγχος    = playoffs
#
#
library(multcomp)
library(ggplot2)

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
# 1. Δόμηση των δεδομένων
# ---------------------------------------------------

data0 <- nhl[nhl$phase == "Regular",
             c("team1", "team2", "goals1", "goals2", "neutral")]

n <- nrow(data0)

cat("Αγώνες εκπαίδευσης (κανονική περίοδος):", n, "\n")

data       <- data.frame(goals = c(data0$goals1, data0$goals2))
data$games <- rep(1:n, 2)

# Neutral αγώνας δεν παίρνει home effect.
data$home <- c(as.numeric(!data0$neutral), rep(0, n))

data$att <- c(data0$team1, data0$team2)
data$def <- c(data0$team2, data0$team1)

data <- data[, c(2, 4, 5, 3, 1)]
data <- data[order(data$games), ]

data$att <- factor(data$att)
data$def <- factor(data$def)

head(data)

nteams <- length(levels(data$att))

cat("Πλήθος ομάδων:", nteams, "\n")

# ---------------------------------------------------
# 2. Περιορισμοί STZ
# ---------------------------------------------------

contrasts(data$att) <- contr.sum(nteams)
contrasts(data$def) <- contr.sum(nteams)

colnames(contrasts(data$att)) <- levels(data$att)[-nteams]
colnames(contrasts(data$def)) <- levels(data$def)[-nteams]

# ---------------------------------------------------
# 3. Προσαρμογή Double Poisson
# ---------------------------------------------------

model <- glm(goals ~ home + att + def, family = poisson, data = data)

print(summary(model))
round(summary(model)$coef, 3)

cat("\nΠαράγοντας έδρας (home):",
    round(coef(model)["home"], 4),
    " -> πολλαπλασιαστικά",
    round(exp(coef(model)["home"]), 4), "\n")

# ---------------------------------------------------
# 4. Έλεγχος υπερδιασποράς
# ---------------------------------------------------

pearson <- sum(residuals(model, type = "pearson")^2)

cat("\nPearson chi^2 =", round(pearson, 1),
    " df =", model$df.residual,
    " λόγος =", round(pearson/model$df.residual, 3), "\n")

cat("Deviance =", round(model$deviance, 1),
    " λόγος =", round(model$deviance/model$df.residual, 3), "\n")

# ---------------------------------------------------
# 5. Πίνακας ικανοτήτων
# ---------------------------------------------------

abilities <- matrix(nrow = nteams, ncol = 4)

index <- 1:(nteams - 1)

abilities[index, 1:2] <- summary(model)$coef[2 + index, 1:2]
abilities[index, 3:4] <- summary(model)$coef[nteams + 1 + index, 1:2]

abilities[nteams, 1] <- -sum(summary(model)$coef[2 + index, 1])
abilities[nteams, 3] <- -sum(summary(model)$coef[nteams + 1 + index, 1])

rownames(abilities) <- levels(data$att)
colnames(abilities) <- c("Att", "SD-Att", "Def", "SD-Def")

# ---------------------------------------------------
# 6. Τυπικά σφάλματα της τελευταίας ομάδας
# ---------------------------------------------------

L <- numeric(length(model$coef))
L[2 + index] <- -1

attLast <- glht(model, linfct = rbind("missing_effect" = L))

abilities[nteams, 2] <- summary(attLast)$test$sigma

L <- numeric(length(model$coef))
L[nteams + 1 + index] <- -1

defLast <- glht(model, linfct = rbind("missing_effect" = L))

abilities[nteams, 4] <- summary(defLast)$test$sigma

cat("\n---------- ΙΚΑΝΟΤΗΤΕΣ ΟΜΑΔΩΝ ----------\n")

print(round(abilities[order(-abilities[, 1]), ], 3))

ab <- as.data.frame(abilities)
ab$team <- rownames(abilities)

safe_write(ab, "abilities_dp.csv")

# ---------------------------------------------------
# 7. Διάγραμμα επιθετικής vs αμυντικής ικανότητας
# ---------------------------------------------------

png("fig10_abilities_scatter.png", width = pngw, height = pngh, res = pngr)

plot(abilities[, 1], abilities[, 3],
     pch = 16, col = pal[1],
     xlab = "Επιθετική ικανότητα (att)",
     ylab = "Αμυντική ικανότητα (def)",
     main = "Ικανότητες ομάδων - Double Poisson (NHL 2021-22)")

a <- 0.015

text(abilities[, 1] - a, abilities[, 3] + a,
     labels = rownames(abilities), cex = 0.6)

abline(h = 0, col = grey1, lwd = 2)
abline(v = 0, col = grey1, lwd = 2)

dev.off()

# ---------------------------------------------------
# 8. Διαστήματα εμπιστοσύνης και error bars
# ---------------------------------------------------

att.ci <- matrix(NA, nteams, 3)

index0 <- 1:(nteams - 1)
ind <- 2 + index0

att.ci[index0, 1]   <- model$coef[ind]
att.ci[index0, 2:3] <- confint(model)[ind, ]

att.ci[nteams, 1]   <- -sum(att.ci[index0, 1])
att.ci[nteams, 2:3] <- confint(attLast)$confint[2:3]

rownames(att.ci) <- levels(data$att)
colnames(att.ci) <- c("Att", "LCI", "UCI")

def.ci <- matrix(NA, nteams, 3)

ind <- nteams + 1 + index0

def.ci[index0, 1]   <- model$coef[ind]
def.ci[index0, 2:3] <- confint(model)[ind, ]

def.ci[nteams, 1]   <- -sum(def.ci[index0, 1])
def.ci[nteams, 2:3] <- confint(defLast)$confint[2:3]

rownames(def.ci) <- levels(data$def)
colnames(def.ci) <- c("Def", "LCI", "UCI")

att.ci <- as.data.frame(att.ci)

att.ci$team <- levels(data$att)

att.ci <- att.ci[order(att.ci$Att), ]

att.ci$team <- factor(att.ci$team, levels = att.ci$team)

p <- ggplot(att.ci, aes(x = Att, y = team, xmin = LCI, xmax = UCI)) +
  geom_point(colour = pal[1]) +
  geom_errorbarh(height = .2, colour = pal[1]) +
  geom_vline(xintercept = 0, colour = accent, linewidth = 1) +
  labs(title = "Επιθετική ικανότητα με 95% διαστήματα εμπιστοσύνης",
       x = "att", y = NULL) +
  theme_minimal(base_size = 11)

ggsave("fig11_att_errorbars.png", p, width = 7, height = 8, dpi = 300)

def.ci <- as.data.frame(def.ci)

def.ci$team <- levels(data$def)

def.ci <- def.ci[order(def.ci$Def), ]

def.ci$team <- factor(def.ci$team, levels = def.ci$team)

p <- ggplot(def.ci, aes(x = Def, y = team, xmin = LCI, xmax = UCI)) +
  geom_point(colour = pal[2]) +
  geom_errorbarh(height = .2, colour = pal[2]) +
  geom_vline(xintercept = 0, colour = accent, linewidth = 1) +
  labs(title = "Αμυντική ικανότητα με 95% διαστήματα εμπιστοσύνης",
       x = "def", y = NULL) +
  theme_minimal(base_size = 11)

ggsave("fig12_def_errorbars.png", p, width = 7, height = 8, dpi = 300)

# ---------------------------------------------------
# 9. Πρόβλεψη στα playoffs
# ---------------------------------------------------

test <- nhl[nhl$phase == "Playoffs", ]

m <- nrow(test)

cat("\nΑγώνες ελέγχου (playoffs):", m, "\n")

newdat <- data.frame(
  home = c(as.numeric(!test$neutral), rep(0, m)),
  att  = factor(c(test$team1, test$team2), levels = levels(data$att)),
  def  = factor(c(test$team2, test$team1), levels = levels(data$def))
)

contrasts(newdat$att) <- contr.sum(nteams)
contrasts(newdat$def) <- contr.sum(nteams)

colnames(contrasts(newdat$att)) <- levels(data$att)[-nteams]
colnames(contrasts(newdat$def)) <- levels(data$def)[-nteams]

lambda <- exp(predict(model, newdat))

lam1 <- lambda[1:m]
lam2 <- lambda[(m + 1):(2 * m)]

B  <- 10000

p1 <- numeric(m)
p0 <- numeric(m)
p2 <- numeric(m)

set.seed(1)

for (j in 1:m) {
  g1 <- rpois(B, lam1[j])
  g2 <- rpois(B, lam2[j])
  
  d <- g1 - g2
  
  p1[j] <- mean(d > 0)
  p0[j] <- mean(d == 0)
  p2[j] <- mean(d < 0)
}

cat("Μέση πιθανότητα ισοπαλίας στην κανονική διάρκεια:",
    round(mean(p0), 3), "\n")

cat("(η μάζα αυτή πρέπει να αναδιανεμηθεί - βλ. επόμενο μοντέλο)\n")

p1c <- p1/(1 - p0)

# Πλήρης ακρίβεια στην πιθανότητα αγοράς.
mkt <- (1/test$odds1)/((1/test$odds1) + (1/test$odds2))

market_idx <- test$market_ok & is.finite(mkt)

pred <- data.frame(
  date = test$date,
  team1 = test$team1,
  team2 = test$team2,
  goals1 = test$goals1,
  goals2 = test$goals2,
  lam1 = round(lam1, 3),
  lam2 = round(lam2, 3),
  pHome = round(p1, 3),
  pDraw = round(p0, 3),
  pAway = round(p2, 3),
  pHomeCond = round(p1c, 3),
  market = round(mkt, 3),
  win1 = test$win1
)

print(head(pred, 10))

safe_write(pred, "pred_dp_playoffs.csv")

# ---------------------------------------------------
# 10. Αξιολόγηση έναντι της αγοράς
# ---------------------------------------------------

brier <- function(p, y) mean((p - y)^2)
acc   <- function(p, y) mean((p > 0.5) == (y == 1))

y <- test$win1

cat("\n---------- OUT-OF-SAMPLE (playoffs, n =", m, ") ----------\n")

cat("Double Poisson (raw)  : Brier =", round(brier(p1, y), 4),
    " Accuracy =", round(100*acc(p1, y), 1), "%\n")

cat("Double Poisson (cond) : Brier =", round(brier(p1c, y), 4),
    " Accuracy =", round(100*acc(p1c, y), 1), "%\n")

cat("\nΚοινό δείγμα μοντέλου-αγοράς, n =", sum(market_idx), "\n")

cat("Double Poisson (cond) : Brier =",
    round(brier(p1c[market_idx], y[market_idx]), 4),
    " Accuracy =", round(100*acc(p1c[market_idx], y[market_idx]), 1), "%\n")

cat("Αγορά                 : Brier =",
    round(brier(mkt[market_idx], y[market_idx]), 4),
    " Accuracy =", round(100*acc(mkt[market_idx], y[market_idx]), 1), "%\n")

png("fig13_model_vs_market.png", width = pngw, height = pngh, res = pngr)

plot(mkt[market_idx], p1c[market_idx],
     pch = 16, col = pal[1], xlim = c(0,1), ylim = c(0,1),
     xlab = "Πιθανότητα αγοράς (fair)",
     ylab = "Πιθανότητα μοντέλου (Double Poisson)",
     main = "Μοντέλο vs αγορά - αγώνες playoffs")

abline(0, 1, lty = 2, col = grey1)

dev.off()

# ---------------------------------------------------
# 11. Έλεγχος ευαισθησίας: χρονολογικό 70/30
# ---------------------------------------------------

reg <- nhl[nhl$phase == "Regular", ]
reg <- reg[order(reg$date), ]

k <- floor(0.7*nrow(reg))

tr <- reg[1:k, ]
te <- reg[(k + 1):nrow(reg), ]

cat("\nΈλεγχος ευαισθησίας - εκπαίδευση:",
    nrow(tr), " έλεγχος:", nrow(te), "\n")

n2 <- nrow(tr)

d2       <- data.frame(goals = c(tr$goals1, tr$goals2))
d2$games <- rep(1:n2, 2)
d2$home  <- c(as.numeric(!tr$neutral), rep(0, n2))
d2$att   <- factor(c(tr$team1, tr$team2))
d2$def   <- factor(c(tr$team2, tr$team1))

d2 <- d2[order(d2$games), ]

nt2 <- length(levels(d2$att))

contrasts(d2$att) <- contr.sum(nt2)
contrasts(d2$def) <- contr.sum(nt2)

colnames(contrasts(d2$att)) <- levels(d2$att)[-nt2]
colnames(contrasts(d2$def)) <- levels(d2$def)[-nt2]

model2 <- glm(goals ~ home + att + def, family = poisson, data = d2)

m2 <- nrow(te)

nd2 <- data.frame(
  home = c(as.numeric(!te$neutral), rep(0, m2)),
  att  = factor(c(te$team1, te$team2), levels = levels(d2$att)),
  def  = factor(c(te$team2, te$team1), levels = levels(d2$def))
)

contrasts(nd2$att) <- contr.sum(nt2)
contrasts(nd2$def) <- contr.sum(nt2)

colnames(contrasts(nd2$att)) <- levels(d2$att)[-nt2]
colnames(contrasts(nd2$def)) <- levels(d2$def)[-nt2]

lam <- exp(predict(model2, nd2))

la1 <- lam[1:m2]
la2 <- lam[(m2 + 1):(2*m2)]

q1 <- numeric(m2)
q0 <- numeric(m2)

set.seed(2)

for (j in 1:m2) {
  g1 <- rpois(B, la1[j])
  g2 <- rpois(B, la2[j])
  
  d <- g1 - g2
  
  q1[j] <- mean(d > 0)
  q0[j] <- mean(d == 0)
}

q1c <- q1/(1 - q0)

mkt2 <- (1/te$odds1)/((1/te$odds1) + (1/te$odds2))

market_idx2 <- te$market_ok & is.finite(mkt2)

cat("Double Poisson (cond) : Brier =",
    round(brier(q1c[market_idx2], te$win1[market_idx2]), 4),
    " Accuracy =",
    round(100*acc(q1c[market_idx2], te$win1[market_idx2]), 1), "%\n")

cat("Αγορά                 : Brier =",
    round(brier(mkt2[market_idx2], te$win1[market_idx2]), 4),
    " Accuracy =",
    round(100*acc(mkt2[market_idx2], te$win1[market_idx2]), 1), "%\n")

cat("\n===================================================\n")
cat(" Ολοκληρώθηκε το μοντέλο Double Poisson\n")
cat("===================================================\n")