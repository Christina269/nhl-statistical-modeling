#1
#ΣΤΑΤΙΣΤΙΚΗ ΑΝΑΛΥΣΗ ΚΑΙ ΠΡΟΓΝΩΣΤΙΚΗ ΜΟΝΤΕΛΟΠΟΙΗΣΗ ΣΤΟ ICE HOCKEY
# NHL 2021-22
#
# ΜΕΡΟΣ Α : Συλλογή, καθαρισμός και δόμηση των δεδομένων
# ΜΕΡΟΣ Β : Περιγραφική στατιστική και οπτικοποίηση
#
# Πηγή: Sportsbook Reviews Online - NHL scores & odds archive
# Δομή πηγής: δύο γραμμές ανά αγώνα
#   Date | Rot | VH | Team | 1st | 2nd | 3rd | Final | Open | Close |
#   PuckLine | PL odds | OpenOU | OU odds | CloseOU | OU odds
#   VH: V = φιλοξενούμενος, H = γηπεδούχος, N = ουδέτερη έδρα
#
# Επιλογή σεζόν: το αρχείο του SBR έχει αναρτημένη ολόκληρη τη σεζόν 2021-22
# (1312 αγώνες κανονικής περιόδου + 89 playoffs). Η σεζόν 2022-23 είναι
# αναρτημένη μόνο μέχρι τις 27/11/2022, γι' αυτό δεν χρησιμοποιείται.
#
# Παράγει: nhl_2021_22.csv, market_baseline.csv, team_summary_2021_22.csv
#          και τα γραφήματα fig1 - fig9 σε png
#
#
library(rvest)
library(moments)

# ---- Παράμετροι σεζόν -------------------------------------------------
season <- "2021-22"
yearA  <- 2021
yearB  <- 2022
cutoff <- as.Date("2022-05-02")

# ---- Ενιαία παλέτα γραφημάτων -----------------------------------------
pal    <- c("#4C78A8", "#F58518", "#54A24B", "#B279A2")
pal4   <- c("#C6DBEF", "#9ECAE1", "#4292C6", "#08519C")
accent <- "#E45756"
grey1  <- "#BAB0AC"
pngw <- 2100; pngh <- 1500; pngr <- 300

# ---- Ασφαλής εγγραφή αρχείων ------------------------------------------
safe_write <- function(x, f) {
  ok <- try(write.csv(x, f, row.names = FALSE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    stop("Failed to write ", f,
        "- πιθανόν είναι ανοιχτό σε άλλο πρόγραμμα\n")
  } else {
    cat("Αποθηκεύτηκε:", f, "\n")
  }
}


###########################################################################
#                                                                         #
#   ΜΕΡΟΣ Α :  ΔΕΔΟΜΕΝΑ                                                   #
#                                                                         #
###########################################################################

url  <- paste("https://www.sportsbookreviewsonline.com/scoresoddsarchives/nhl-odds-",
              season, "/", sep = "")
file <- paste("sbr_nhl_", season, ".html", sep = "")

# ---------------------------------------------------
# Α1. Κατέβασμα της σελίδας
# ---------------------------------------------------

if (!file.exists(file)) {
  download.file(url, destfile = file, quiet = FALSE, mode = "wb",
                headers = c("User-Agent" =
                              "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"))
}

cat("Μέγεθος αρχείου (bytes):", file.info(file)$size, "\n")

page <- read_html(file)

# ---------------------------------------------------
# Α2. Εντοπισμός των γραμμών δεδομένων
# ---------------------------------------------------

tabs <- html_table(page, fill = TRUE)
cat("Πίνακες HTML που βρέθηκαν:", length(tabs), "\n")

if (length(tabs) > 0) {
  for (i in 1:length(tabs)) {
    cat("   πίνακας", i, ":", nrow(tabs[[i]]), "γραμμές,",
        ncol(tabs[[i]]), "στήλες\n")
  }
}

big <- 0

if (length(tabs) > 0) {
  sizes <- sapply(tabs, nrow)
  if (max(sizes) >= 100) big <- which.max(sizes)
}

if (big > 0) {
  cat("Χρησιμοποιείται ο πίνακας", big, "\n")
  htab  <- as.data.frame(tabs[[big]])
  lines <- apply(htab, 1, function(r) paste(r, collapse = " "))
} else {
  cat("Δεν βρέθηκε πίνακας - διαβάζουμε το κείμενο της σελίδας\n")
  lines <- unlist(strsplit(html_text2(page), "\n"))
}

cat("Γραμμές προς έλεγχο:", length(lines), "\n")

# ---------------------------------------------------
# Α3. Ανάλυση μίας γραμμής σε πεδία
# ---------------------------------------------------

parse_row <- function(line) {
  line <- gsub("\u00a0", " ", line)
  line <- gsub("[[:space:]]+", " ", line)
  tk   <- strsplit(trimws(line), " ", fixed = TRUE)[[1]]
  tk   <- tk[tk != ""]
  
  if (length(tk) < 12) return(NULL)
  if (!grepl("^[0-9]{3,4}$", tk[1])) return(NULL)
  if (!(tk[3] %in% c("V", "H", "N"))) return(NULL)
  
  rest  <- tk[-(1:3)]
  isnum <- grepl("^-?[0-9]*\\.?[0-9]+$", rest)
  k     <- which(isnum)[1]
  
  if (is.na(k) || k == 1) return(NULL)
  
  team <- paste(rest[1:(k - 1)], collapse = "")
  v    <- rest[k:length(rest)]
  length(v) <- 12
  
  c(tk[1], tk[3], team, v)
}

out <- lapply(lines, parse_row)
out <- out[!sapply(out, is.null)]

cat("Γραμμές δεδομένων που αναγνωρίστηκαν:", length(out),
    "(αναμένονται 2 ανά αγώνα)\n")

if (length(out) == 0) {
  stop("Δεν αναγνωρίστηκε καμία γραμμή. Έλεγξε το ", file, ": ",
       "αν είναι μικρό (< 100 KB) ή σελίδα σφάλματος, σώσε το χειροκίνητα ",
       "από browser.")
}

M <- do.call(rbind, out)

data0 <- data.frame(M, stringsAsFactors = FALSE)

names(data0) <- c("dt", "vh", "team", "p1", "p2", "p3", "final",
                  "open", "close", "pline", "plodds",
                  "ouopen", "ouopenodds", "ouclose", "oucloseodds")

data0$team[data0$team == "SeattleKraken"] <- "Seattle"

head(data0)

# ---------------------------------------------------
# Α4. Μετατροπή τύπων και ημερομηνία
# ---------------------------------------------------

num <- c("p1", "p2", "p3", "final", "open", "close", "pline",
         "plodds", "ouopen", "ouopenodds", "ouclose", "oucloseodds")

for (v in num) data0[[v]] <- suppressWarnings(as.numeric(data0[[v]]))

data0$goalsreg <- data0$p1 + data0$p2 + data0$p3

dd    <- sprintf("%04d", as.integer(data0$dt))
month <- as.integer(substr(dd, 1, 2))
day   <- as.integer(substr(dd, 3, 4))
year  <- ifelse(month >= 9, yearA, yearB)

data0$date <- as.Date(sprintf("%d-%02d-%02d", year, month, day))

# ---------------------------------------------------
# Α5. Ζευγάρωμα των γραμμών σε αγώνες
# ---------------------------------------------------

nrows <- nrow(data0)

if (nrows %% 2 != 0) {
  stop("Περιττό πλήθος γραμμών - λείπει γραμμή κάποιου αγώνα")
}

i1 <- seq(1, nrows, by = 2)
i2 <- i1 + 1

A <- data0[i1, ]
B <- data0[i2, ]

bad_pair <- is.na(A$date) | is.na(B$date) | A$date != B$date
if (any(bad_pair)) stop("Invalid dates or unmatched rows: check the source file.")
valid_sides <- (A$vh == "V" & B$vh == "H") | (A$vh == "H" & B$vh == "V") | (A$vh == "N" & B$vh == "N")
if (any(!valid_sides)) stop("Invalid home/away pairing: check the source file.")
cat("Ζεύγη με ασυμφωνία ημερομηνίας:", sum(bad_pair), "\n")

homeB <- (B$vh == "H") | (A$vh == "V") | (A$vh == "N" & B$vh == "N")

hm <- function(v) ifelse(homeB, B[[v]], A[[v]])
aw <- function(v) ifelse(homeB, A[[v]], B[[v]])

am2dec <- function(x) {
  ifelse(!is.finite(x) | x == 0, NA_real_, ifelse(x > 0, 1 + x/100, 1 + 100/abs(x)))
}

dat <- data.frame(
  date        = A$date,
  neutral     = (A$vh == "N" & B$vh == "N"),
  team1       = hm("team"),
  team2       = aw("team"),
  goals1      = hm("final"),
  goals2      = aw("final"),
  goals1reg   = hm("goalsreg"),
  goals2reg   = aw("goalsreg"),
  odds1       = am2dec(hm("close")),
  odds2       = am2dec(aw("close")),
  odds1open   = am2dec(hm("open")),
  odds2open   = am2dec(aw("open")),
  totalline   = hm("ouclose"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------
# Α6. Μεταβλητές ανάλυσης
# ---------------------------------------------------

dat$phase <- factor(ifelse(dat$date >= cutoff, "Playoffs", "Regular"),
                    levels = c("Regular", "Playoffs"))

dat$total <- dat$goals1 + dat$goals2
dat$dif   <- dat$goals1 - dat$goals2
dat$win1  <- as.numeric(dat$dif > 0)

dat$extra   <- (dat$goals1 + dat$goals2) - (dat$goals1reg + dat$goals2reg)

dat$endtype <- factor(ifelse(dat$extra == 0, "REG", "OT/SO"),
                      levels = c("REG", "OT/SO"))

lev <- c("Away REG", "Away OT/SO", "Home OT/SO", "Home REG")

res <- rep(NA, nrow(dat))
res[dat$win1 == 0 & dat$endtype == "REG"]   <- lev[1]
res[dat$win1 == 0 & dat$endtype == "OT/SO"] <- lev[2]
res[dat$win1 == 1 & dat$endtype == "OT/SO"] <- lev[3]
res[dat$win1 == 1 & dat$endtype == "REG"]   <- lev[4]

dat$result  <- factor(res, levels = lev, ordered = TRUE)
dat$result4 <- as.numeric(dat$result)

dat$overround <- 1/dat$odds1 + 1/dat$odds2 - 1

# ---------------------------------------------------
# Α7. Έλεγχοι ποιότητας και καθαρισμός
# ---------------------------------------------------

cat("\n---------- ΕΛΕΓΧΟΙ ΠΟΙΟΤΗΤΑΣ ----------\n")
cat("Αγώνες πριν τον καθαρισμό:", nrow(dat), "\n")

e1 <- is.na(dat$extra) | !(dat$extra %in% c(0, 1))

cat("(α) Ασυνεπές άθροισμα περιόδων vs τελικό σκορ:", sum(e1), "\n")

if (sum(e1) > 0) {
  print(dat[e1, c("date","team2","team1","goals1","goals2","extra")])
}

e2 <- is.na(dat$dif) | dat$dif == 0

cat("(β) Ισοπαλίες (πρέπει 0):", sum(e2), "\n")

e3 <- is.na(dat$odds1) | is.na(dat$odds2) |
  dat$odds1 < 1.01 | dat$odds2 < 1.01 |
  dat$odds1 > 30   | dat$odds2 > 30

cat("(γ) Προβληματικές/ελλείπουσες αποδόσεις:", sum(e3), "\n")

if (sum(e3) > 0) {
  print(dat[e3, c("date","team2","team1","odds1","odds2")])
}

e4 <- is.na(dat$overround) | dat$overround < 0 | dat$overround > 0.20

cat("(δ) Ύποπτο overround:", sum(e4), "\n")


# -------------------- ΔΙΟΡΘΩΣΗ --------------------
# Οι e3/e4 αφορούν μόνο την αγορά και δεν πρέπει να αφαιρούν
# σωστό αγωνιστικό αποτέλεσμα από τα μοντέλα.

dat$market_ok <- !(e3 | e4)

nhl <- dat[!(e1 | e2), ]
nhl <- nhl[order(nhl$date), ]

nhl_market <- nhl[nhl$market_ok, ]

# ---------------------------------------------------

cat("\nΑγώνες μετά τον καθαρισμό:", nrow(nhl), "\n")

print(table(nhl$phase))
print(table(nhl$endtype))

cat("Αγώνες με έγκυρες αποδόσεις:", nrow(nhl_market), "\n")
cat("Neutral-site games:", sum(nhl$neutral), "\n")

cat("Ομάδες:", length(unique(c(nhl$team1, nhl$team2))), " (πρέπει 32)\n")

print(sort(unique(c(nhl$team1, nhl$team2))))

cat("Εύρος ημερομηνιών:", format(min(nhl$date)), "-", format(max(nhl$date)), "\n")

safe_write(nhl, "nhl_2021_22.csv")


###########################################################################
#                                                                         #
#   ΜΕΡΟΣ Β :  ΠΕΡΙΓΡΑΦΙΚΗ ΑΝΑΛΥΣΗ                                        #
#                                                                         #
###########################################################################

n <- nrow(nhl)

cat("\n===================================================\n")
cat(" NHL 2021-22 - ΠΕΡΙΓΡΑΦΙΚΗ ΑΝΑΛΥΣΗ\n")
cat("===================================================\n")

cat("Μέγεθος δείγματος n =", n, "αγώνες\n")
cat("Περίοδος:", format(min(nhl$date)), "-", format(max(nhl$date)), "\n")

print(table(nhl$phase))

# ---------------------------------------------------
# Β1. Η δυναμική της έδρας
# ---------------------------------------------------

cat("\n--- Β1. ΠΛΕΟΝΕΚΤΗΜΑ ΕΔΡΑΣ ---\n")

w <- table(factor(nhl$win1, levels = c(0,1),
                  labels = c("Away win","Home win")))

print(w)
print(round(100*prop.table(w), 1))

png("fig1_home_advantage.png", width = pngw, height = pngh, res = pngr)

pie(w,
    labels = paste(names(w), "\n", round(100*prop.table(w), 1), "%"),
    col = c(pal[2], pal[1]),
    main = "Ποσοστά νικών γηπεδούχων και φιλοξενούμενων (NHL 2021-22)")

dev.off()

# ---------------------------------------------------
# Β2. Η κατανομή της τελικής έκβασης
# ---------------------------------------------------

cat("\n--- Β2. ΤΕΛΙΚΗ ΕΚΒΑΣΗ ---\n")

tabres <- table(nhl$result)

print(tabres)
print(round(100*prop.table(tabres), 1))

cat("Αγώνες πέρα από την κανονική διάρκεια:",
    round(100*mean(nhl$endtype == "OT/SO"), 1), "%\n")

png("fig2_result_distribution.png", width = pngw, height = pngh, res = pngr)

b <- barplot(tabres, col = pal4, border = "white",
             ylim = c(0, 1.15*max(tabres)),
             xlab = "Έκβαση (από τη σκοπιά του γηπεδούχου)",
             ylab = "Αριθμός αγώνων",
             main = "Κατανομή της τελικής έκβασης (NHL 2021-22)")

text(b, tabres, labels = tabres, pos = 3, font = 2)

dev.off()

# ---------------------------------------------------
# Β3. Παραγωγικότητα σε γκολ
# ---------------------------------------------------

cat("\n--- Β3. ΣΥΝΟΛΙΚΑ ΓΚΟΛ ΑΝΑ ΑΓΩΝΑ ---\n")

tg <- nhl$total

print(round(c(
  mean = mean(tg),
  sd = sd(tg),
  min = min(tg),
  max = max(tg),
  median = median(tg),
  skewness = skewness(tg),
  kurtosis = kurtosis(tg) - 3
), 3))

png("fig3_total_goals.png", width = pngw, height = pngh, res = pngr)

hist(tg, breaks = seq(min(tg)-0.5, max(tg)+0.5, 1),
     col = pal[1], border = "white",
     xlab = "Συνολικά γκολ", ylab = "Αριθμός αγώνων",
     main = "Κατανομή των συνολικών γκολ ανά αγώνα")

abline(v = mean(tg), col = accent, lty = 2, lwd = 2)

legend("topright",
       legend = paste("Μέση τιμή:", round(mean(tg), 2)),
       col = accent, lty = 2, lwd = 2, bty = "n")

dev.off()

# ---------------------------------------------------
# Β4. Άρτια και περιττά σύνολα γκολ
# ---------------------------------------------------

cat("\n--- Β4. ΑΡΤΙΑ vs ΠΕΡΙΤΤΑ ΣΥΝΟΛΑ ---\n")

cat("Περιττό σύνολο γκολ - όλοι οι αγώνες:",
    round(100*mean(nhl$total %% 2 == 1), 1), "%\n")

cat("  μόνο REG  :",
    round(100*mean(nhl$total[nhl$endtype == "REG"] %% 2 == 1), 1), "%\n")

cat("  μόνο OT/SO:",
    round(100*mean(nhl$total[nhl$endtype == "OT/SO"] %% 2 == 1), 1),
    "%  (έλεγχος ορθότητας: πρέπει 100)\n")

png("fig4_total_goals_by_endtype.png", width = pngw, height = pngh, res = pngr)

barplot(t(table(nhl$total, nhl$endtype)),
        col = c(pal[1], pal[2]), border = "white",
        xlab = "Συνολικά γκολ", ylab = "Αριθμός αγώνων",
        main = "Συνολικά γκολ: κανονική διάρκεια vs OT/SO")

legend("topright", c("REG", "OT/SO"),
       fill = c(pal[1], pal[2]), bty = "n")

dev.off()

# ---------------------------------------------------
# Β5. Γκολ ανά ομάδα και υπόθεση Poisson
# ---------------------------------------------------

cat("\n--- Β5. ΓΚΟΛ ΑΝΑ ΟΜΑΔΑ ---\n")

cat("Γηπεδούχοι : μέση =", round(mean(nhl$goals1), 3),
    " διασπορά =", round(var(nhl$goals1), 3),
    " λόγος =", round(var(nhl$goals1)/mean(nhl$goals1), 3), "\n")

cat("Φιλοξ.     : μέση =", round(mean(nhl$goals2), 3),
    " διασπορά =", round(var(nhl$goals2), 3),
    " λόγος =", round(var(nhl$goals2)/mean(nhl$goals2), 3), "\n")

goals <- c(nhl$goals1, nhl$goals2)
lam   <- mean(goals)

obs  <- as.numeric(table(factor(goals, levels = 0:max(goals))))/length(goals)
expp <- dpois(0:max(goals), lam)

png("fig5_goals_vs_poisson.png", width = pngw, height = pngh, res = pngr)

barplot(rbind(obs, expp), beside = TRUE, names.arg = 0:max(goals),
        col = c(pal[1], pal[2]), border = "white",
        xlab = "Γκολ", ylab = "Σχετική συχνότητα",
        main = paste("Γκολ ανά ομάδα vs Poisson(", round(lam, 2), ")", sep = ""))

legend("topright", c("Παρατηρούμενα", "Poisson"),
       fill = c(pal[1], pal[2]), bty = "n")

dev.off()

# ---------------------------------------------------
# Β6. Η διαφορά γκολ
# ---------------------------------------------------

cat("\n--- Β6. ΔΙΑΦΟΡΑ ΓΚΟΛ ---\n")

cat("μέση =", round(mean(nhl$dif), 3),
    " τ.α. =", round(sd(nhl$dif), 3),
    " P(|dif| = 1) =", round(100*mean(abs(nhl$dif) == 1), 1), "%\n")

png("fig6_goal_difference.png", width = pngw, height = pngh, res = pngr)

barplot(table(nhl$dif), col = pal[1], border = "white",
        xlab = "Διαφορά γκολ (Γηπεδούχος - Φιλοξενούμενος)",
        ylab = "Αριθμός αγώνων",
        main = "Κατανομή της διαφοράς γκολ")

mtext("Η τιμή 0 είναι εξ ορισμού αδύνατη (zero deflation)",
      side = 3, cex = 0.8)

dev.off()

# ---------------------------------------------------
# Β7. Κανονική περίοδος vs playoffs
# ---------------------------------------------------

cat("\n--- Β7. ΚΑΝΟΝΙΚΗ ΠΕΡΙΟΔΟΣ vs PLAYOFFS ---\n")

f <- function(d) c(
  Agones = nrow(d),
  HomeWinPct = 100*mean(d$win1),
  BeyondRegPct = 100*mean(d$endtype == "OT/SO"),
  MeanTotal = mean(d$total),
  SdTotal = sd(d$total),
  MeanHome = mean(d$goals1),
  MeanAway = mean(d$goals2)
)

cmp <- rbind(
  Regular = f(nhl[nhl$phase == "Regular", ]),
  Playoffs = f(nhl[nhl$phase == "Playoffs", ])
)

print(round(cmp, 2))

png("fig7_phase_comparison.png", width = pngw, height = pngh, res = pngr)

par(mfrow = c(1, 2))

hw <- cmp[, "HomeWinPct"]

b <- barplot(hw, col = c(pal[1], pal[2]), border = "white",
             ylim = c(0, 1.2*max(hw)),
             ylab = "%", main = "Ποσοστό νικών έδρας ανά φάση")

text(b, hw, labels = paste(round(hw, 1), "%"), pos = 3, font = 2)

pt <- prop.table(table(nhl$phase, nhl$result), 1)

barplot(100*t(pt), beside = TRUE, col = pal4, border = "white",
        ylab = "% αγώνων", main = "Κατανομή έκβασης ανά φάση")

legend("topleft", lev, fill = pal4, bty = "n", cex = 0.7)

par(mfrow = c(1, 1))

dev.off()

# ---------------------------------------------------
# Β8. Οι αποδόσεις της αγοράς
# ---------------------------------------------------

# Χρησιμοποιούνται μόνο οι αγώνες με έγκυρες αποδόσεις.

nhl_all <- nhl
nhl <- nhl_market
n <- nrow(nhl)

cat("\n--- Β8. ΑΠΟΔΟΣΕΙΣ ΤΗΣ ΑΓΟΡΑΣ ---\n")

nhl$imp1  <- 1/nhl$odds1
nhl$imp2  <- 1/nhl$odds2
nhl$book  <- nhl$imp1 + nhl$imp2
nhl$fair1 <- nhl$imp1/nhl$book
nhl$fair2 <- nhl$imp2/nhl$book

cat("Μέσο άθροισμα υπονοούμενων πιθανοτήτων:",
    round(mean(nhl$book), 4), "\n")

cat("Μέσο περιθώριο (overround):",
    round(100*mean(nhl$book - 1), 2), "%\n")

cat("Διάμεσο:", round(100*median(nhl$book - 1), 2),
    "%  εύρος:", round(100*min(nhl$book - 1), 2), "-",
    round(100*max(nhl$book - 1), 2), "%\n")

print(tapply(100*(nhl$book - 1), nhl$phase, mean))

i <- which.min(nhl$odds1)

cat("\nΕνδεικτικό παράδειγμα (ο ισχυρότερος φαβορί του δείγματος):\n")

cat(format(nhl$date[i]), ":", nhl$team1[i],
    "(γηπ.) -", nhl$team2[i], "\n")

cat("  Αποδόσεις:", round(nhl$odds1[i], 2),
    "/", round(nhl$odds2[i], 2), "\n")

cat("  Υπονοούμενες:", round(100*nhl$imp1[i], 1),
    "% /", round(100*nhl$imp2[i], 1),
    "%  άθροισμα", round(nhl$book[i], 3),
    "-> overround", round(100*(nhl$book[i]-1), 2), "%\n")

cat("  Fair:", round(100*nhl$fair1[i], 1),
    "% /", round(100*nhl$fair2[i], 1), "%\n")

cat("  Πραγματικό σκορ:", nhl$goals1[i], "-", nhl$goals2[i], "\n")

cat("\nΗ αγορά έκανε φαβορί τον γηπεδούχο στο",
    round(100*mean(nhl$fair1 > 0.5), 1), "% των αγώνων\n")

cat("Μέση fair πιθανότητα νίκης έδρας:",
    round(100*mean(nhl$fair1), 1), "%\n")

cat("Πραγματικό ποσοστό νικών έδρας:",
    round(100*mean(nhl$win1), 1), "%\n")

png("fig8_overround.png", width = pngw, height = pngh, res = pngr)

hist(100*(nhl$book - 1), breaks = 30,
     col = pal[1], border = "white",
     xlab = "Overround (%)", ylab = "Αριθμός αγώνων",
     main = "Κατανομή του περιθωρίου του bookmaker")

abline(v = 100*mean(nhl$book - 1),
       col = accent, lty = 2, lwd = 2)

dev.off()

# ---------------------------------------------------
# Β9. Βαθμονόμηση και baseline της αγοράς
# ---------------------------------------------------

cat("\n--- Β9. ΒΑΘΜΟΝΟΜΗΣΗ ΚΑΙ BASELINE ΤΗΣ ΑΓΟΡΑΣ ---\n")

bin  <- cut(nhl$fair1, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
pred <- tapply(nhl$fair1, bin, mean)
obs2 <- tapply(nhl$win1, bin, mean)
cnt  <- as.numeric(table(bin))

print(round(cbind(n = cnt, pred = pred, obs = obs2), 3))

png("fig9_market_calibration.png", width = 1800, height = 1600, res = pngr)

plot(pred, obs2,
     xlim = c(0,1), ylim = c(0,1),
     pch = 19, col = pal[1],
     cex = 1 + 2*cnt/max(cnt, na.rm = TRUE),
     xlab = "Fair πιθανότητα νίκης γηπεδούχου (αγορά)",
     ylab = "Παρατηρούμενο ποσοστό νικών",
     main = "Βαθμονόμηση των πιθανοτήτων της αγοράς")

abline(0, 1, lty = 2, col = grey1)
lines(pred, obs2, col = pal[1])

dev.off()

brier <- function(p, y) mean((p - y)^2)
acc   <- function(p, y) mean((p > 0.5) == (y == 1))

ir <- nhl$phase == "Regular"
ip <- nhl$phase == "Playoffs"

cat("\nBrier Score (σύνολο)  :",
    round(brier(nhl$fair1, nhl$win1), 4), "\n")

cat("Accuracy    (σύνολο)  :",
    round(100*acc(nhl$fair1, nhl$win1), 1), "%\n")

cat("Brier Score (κανονική):",
    round(brier(nhl$fair1[ir], nhl$win1[ir]), 4), "\n")

cat("Brier Score (playoffs):",
    round(brier(nhl$fair1[ip], nhl$win1[ip]), 4), "\n")

cat("Σταθερή πρόβλεψη 0.50            -> Brier =",
    round(brier(rep(0.5, n), nhl$win1), 4), "\n")

cat("Σταθερή πρόβλεψη = ποσοστό έδρας -> Brier =",
    round(brier(rep(mean(nhl$win1), n), nhl$win1), 4), "\n")

safe_write(nhl[, c("date","phase","team1","team2","odds1","odds2",
                   "fair1","fair2","win1")],
           "market_baseline.csv")

# Επιστροφή στο πλήρες dataset.
nhl <- nhl_all
n <- nrow(nhl)
rm(nhl_all)

# ---------------------------------------------------
# Β10. Συγκεντρωτικός πίνακας ομάδων
# ---------------------------------------------------

cat("\n--- Β10. ΠΙΝΑΚΑΣ ΟΜΑΔΩΝ ---\n")

teams <- sort(unique(c(nhl$team1, nhl$team2)))

tt <- data.frame(team = teams, GP = 0, W = 0, GF = 0, GA = 0)

for (i in 1:length(teams)) {
  tm <- teams[i]
  
  ih <- nhl$team1 == tm
  ia <- nhl$team2 == tm
  
  tt$GP[i] <- sum(ih) + sum(ia)
  tt$W[i]  <- sum(nhl$win1[ih]) + sum(1 - nhl$win1[ia])
  
  tt$GF[i] <- (sum(nhl$goals1[ih]) + sum(nhl$goals2[ia])) / tt$GP[i]
  tt$GA[i] <- (sum(nhl$goals2[ih]) + sum(nhl$goals1[ia])) / tt$GP[i]
}

tt$WinPct <- 100*tt$W/tt$GP
tt$Dif    <- tt$GF - tt$GA

tt <- tt[order(-tt$WinPct), ]

print(tt, digits = 3, row.names = FALSE)

safe_write(tt, "team_summary_2021_22.csv")

cat("\n===================================================\n")
cat(" Ολοκληρώθηκε. Παρήχθησαν fig1 - fig9 png,\n")
cat(" nhl_2021_22.csv, market_baseline.csv, team_summary_2021_22.csv\n")
cat("===================================================\n")