require(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
devtools::install_github("rasigadelab/rlabutils")
devtools::install_github("rasigadelab/thd")
library(rlabutils)
library(readxl)
library(thd)
library(lmerTest)


# Distance matrix
d <- data.table(read_excel("distance_matrix.xlsx", sheet = 1, na = "NA"))
setnames(d, "ID", "id")
# Basic consistency check
stopifnot(all(d$id == names(d)[-1]))

# Metadata
m <- data.table(read_excel("meta.xlsx", sheet = 1, na = "NA"))

# Rename some variables
renames <- matrix(c(
  "strain ID", "id",
  "Country of Isolation", "country"
), ncol = 2, byrow = TRUE)

setnames(m, renames[,1], renames[,2])

# Names are OK, order is different
stopifnot(length(setdiff(m$id, d$id)) == 0)
stopifnot(length(setdiff(d$id, m$id)) == 0)

names(m)

# Align metadata
m <- merge(data.table(id = d$id), m, by = "id", sort = FALSE)
stopifnot(all(m$id == d$id))

# Hamming matrix
H <- as.matrix(d[, -c("id")])
stopifnot(all(is.finite(H)))
stopifnot(all(dim(H) == nrow(m)))

rownames(H) <- colnames(H)

# Distribution of SNP distances
hist(H[lower.tri(H)])

# THD parameters
effsize <- 3e6
mu <- 1e-6
timescale <- 10

# Compute THD; rescale x1000 for display
t <- thd(H, timescale, effsize, mu, scale = "density") * 1000
hist(t)
m[, thd := t]

# Define mercury as either COMER_merA or COMER_merB
m[, mer := COMER_merA | COMER_merB]
table(m$mer, useNA = "a")

m[, pvl := `lukF-PV` | `lukS-PV`]
table(m$pvl, useNA = "a")

m[, acme := `ACME_arcA` | `ACME_arcB`]
table(m$acme, useNA = "a")

##################################################################
# SANITIZE COUNTRIES

table(m$country)
m[ grepl("related", country),
   country := {
     x <- country
     x <- gsub("Switzerland \\(related ", "", x)
     x <- gsub("\\)", "", x)
     x
   }
   ]
table(m$country, useNA = "a")

m[, country_SA := {
  country %in% c(
    "Bolivia", "Brazil", "Colombia", "Ecuador", "Uruguay"
  )
}]

m[, clade := factor(Phylogenie_group)]
levels(m$clade) <- c(levels(m$clade), "Other")
m[ is.na(clade), clade := "Other"]

m[, display_mer_SA := {
  x <- paste(as.integer(mer), as.integer(country_SA), sep = "")
  y <- factor(x, levels = c("00", "10", "01", "11"))
  levels(y) <- c("\nOther country\nMer-", "\nOther country\nMer+", "\nContinental SA\nMer-", "\nContinental SA\nMer+")
  y
}]

svgf("thd_comer", 8, 6)
{
  par(mar = c(5,5,2,2))
  boxplot(thd ~ display_mer_SA, m, ylab = "Success index x1000", xlab = "",
          col = rep(stdcols$light[1:2], 2), outline = F)
  beeswarm(thd ~ display_mer_SA, m, add = T, pwcol = m$clade, method = "hex", pch = 19, cex = 0.75)
  legend("bottomright", bty = "n", legend = levels(m$clade), fill = 1:nlevels(m$clade))
}
dev.off()

########################
# Additional figure

# layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE))

svgf("thd_features", 7, 5)
{
  layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = T), widths = c(2,2,2.5), heights = c(3,1))
  
  blank_col <- rgb(0.5,0.5,0.5,0.2)
  
  plv_cols <- c(blank_col, stdcols$transp[1])
  mec_cols <- c(blank_col, stdcols$transp[2])
  mer_cols <- c(blank_col, stdcols$transp[3])
  acme_cols <- c(blank_col, stdcols$transp[4])
  
  par(mar = c(1,4,4,2))
  
  dat <- list(
    PVLneg = m[pvl == F]$thd,
    PVLpos = m[pvl == T]$thd
  )
  boxplot(dat, outline = F, xlab = "", ylab = "Success index x1000", ylim = c(0, 3), xaxt = "n", las = 1)
  beeswarm(dat, method = "hex", add = T, col = plv_cols, pch = 19)

  dat <- list(
    mecA_neg = m[pvl == T & mecA == F]$thd,
    mecA_pos = m[pvl == T & mecA == T]$thd
  )
  boxplot(dat, outline = F, xlab = "", ylab = "Success index x1000", ylim = c(0, 3), xaxt = "n", las = 1)
  beeswarm(dat, method = "hex", add = T, col = mec_cols, pch = 19)
  
  dat <- list(
    COMER_neg_ACME_neg = m[pvl & mecA & !mer & !acme]$thd,
    COMER_pos_ACME_neg = m[pvl & mecA &  mer & !acme]$thd,
    COMER_neg_ACME_pos = m[pvl & mecA & !mer &  acme]$thd
  )
  boxplot(dat, outline = F, xlab = "", ylab = "Success index x1000", ylim = c(0, 3), xaxt = "n", las = 1)
  beeswarm(dat, method = "hex", add = T, col = c(mec_cols[1], mer_cols[2], acme_cols[2]), pch = 19)  
  
  # Marks
  par(mar = c(1,4,0,2))
  par(xpd = NA)
  
  pplot <- function(range) {
    plot(0, bty = "n", type = "n", xlim = range, ylim = c(4, -1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(rep(0,4), 0:3, labels = c("PVL", "mecA", "COMER", "ACME"), adj = 0)
  }
  
  marksize <- 3
  
  pplot(c(.5, 2.5))
  text(c(1,2), c(0,0), labels = c("-", "+"), cex = 1.5, adj = 0.5)
  points(c(1,2), c(0,0), pch = 22, cex = marksize, bg = plv_cols)
  
  pplot(c(.5, 2.5))
  text(c(1,2), c(0,0), labels = c("+", "+"), cex = 1.5, adj = 0.5)
  points(c(1,2), c(0,0), pch = 22, cex = marksize, bg = plv_cols[c(2,2)])
  text(c(1,2), c(1,1), labels = c("-", "+"), cex = 1.5, adj = 0.5)
  points(c(1,2), c(1,1), pch = 22, cex = marksize, bg = mec_cols)
  
  pplot(c(.5, 3.5))
  text(1:3, rep(0,3), labels = rep("+", 3), cex = 1.5, adj = 0.5)
  points(1:3, rep(0,3), pch = 22, cex = marksize, bg = plv_cols[c(2, 2, 2)])
  text(1:3, rep(1,3), labels = rep("+", 3), cex = 1.5, adj = 0.5)
  points(1:3, rep(1,3), pch = 22, cex = marksize, bg = mec_cols[c(2, 2, 2)])
  text(1:3, rep(2,3), labels = c("-", "+", "-"), cex = 1.5, adj = 0.5)
  points(1:3, rep(2,3), pch = 22, cex = marksize, bg = mer_cols[c(1, 2, 1)])
  text(1:3, rep(3,3), labels = c("-", "-", "+"), cex = 1.5, adj = 0.5)
  points(1:3, rep(3,3), pch = 22, cex = marksize, bg = acme_cols[c(1, 1, 2)])

  

}
dev.off()

################################################################
# INTERACTION MODEL
ft <- with(m, pvl == TRUE & mecA == TRUE)
summary(lm(thd ~ mer*country_SA, m[ft]))
confint(lm(thd ~ mer*country_SA, m[ft]))

# Compare COMER effect in SA/not SA
mer_SA_levels <- matrix(c(
  "11", "COMER+\nin South Am.",
  "01", "COMER-\nin South Am.",
  "10", "COMER+\noutside South Am.",
  "00", "COMER-\noutside South Am."
), ncol = 2, byrow = TRUE)

m[, mer_SA := paste(as.integer(mer), as.integer(country_SA), sep = "")]
m[, mer_SA := factor(mer_SA, levels = c(mer_SA_levels[,1]))]
levels(m$mer_SA) <- mer_SA_levels[,2]

par(mar = c(5, 10, 2, 2))

dat <- list(
  `COMER+\nin South Am.` = m[pvl & mecA & mer & country_SA]$thd,
  `COMER-\nin South Am.` = m[pvl & mecA & !mer & country_SA]$thd,
  `COMER+\noutside South Am.` = m[pvl & mecA & mer & !country_SA]$thd,
  `COMER-\noutside South Am.` = m[pvl & mecA & !mer & !country_SA]$thd
)

svgf("comer_SA", 6, 3)
{
  par(mar = c(5, 10, 2, 2))
  boxplot(dat, horizontal = TRUE, las = 1, xlab = "Success index x1000", ylab = "", outline = FALSE)
  beeswarm(dat, vertical = FALSE, method = "hex", add = T, pch = 19, col = stdcols$transp)
}
dev.off()

wilcox.test(thd ~ mer, m[country_SA == FALSE])
