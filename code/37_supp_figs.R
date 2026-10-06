#!/usr/bin/env Rscript
# Supplementary figures assembled from existing results (no model refitting):
#   FigS_tas_clinical_immune : Spearman among TAS, age, stage, and the six immune cells of Figure 5B
#   FigS_cellchat            : the four CellChat panels combined (A-D)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(showtext); library(magick)
})
ROOT <- "/Users/xiaolin/WorkBuddy/2026-09-04-20-00-08/paper2_triaptosis"
OUT  <- file.path(ROOT, "figures_pub")
RES  <- file.path(ROOT, "results")
source(file.path(ROOT, "code", "lib", "viz_base.R"))

# ---- TAS x clinical x immune --------------------------------------------
sc <- read.csv(file.path(RES, "02_model/KIRC/score_table.csv"), row.names = 1)
sc$patientId <- substr(rownames(sc), 1, 12)
clin <- read.csv(file.path(ROOT, "data/external_validation/kirc_tcga_pan_can_atlas_2018_clin_patient.csv"))
st <- gsub("^STAGE\\s*", "", toupper(trimws(clin$AJCC_PATHOLOGIC_TUMOR_STAGE)))
clin$Stage <- match(st, c("I", "II", "III", "IV"))
clin$Age <- as.numeric(clin$AGE)
ssg <- read.csv(file.path(RES, "04_immune/ssgsea_per_sample.csv"), check.names = FALSE)
ssg$patientId <- substr(ssg$patientId, 1, 12)
cells <- c("M2 macrophage", "Mast cell", "Type 17 T helper cell",
           "CD56bright NK cell", "Memory B cell", "B cell")
d <- sc %>% select(patientId, TAS) %>%
  left_join(clin %>% select(patientId, Age, Stage), by = "patientId") %>%
  left_join(ssg %>% select(patientId, all_of(cells)), by = "patientId")
vars <- c("TAS", "Age", "Stage", cells)
cr <- expand.grid(a = vars, b = vars, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(ct = list(suppressWarnings(cor.test(d[[a]], d[[b]], method = "spearman", exact = FALSE))),
         rho = unname(ct$estimate), p = ct$p.value) %>%
  ungroup() %>% select(-ct)
off <- cr$a != cr$b
cr$fdr <- NA_real_
cr$fdr[off] <- p.adjust(cr$p[off], method = "BH")
cr <- cr %>% mutate(a = factor(a, vars), b = factor(b, rev(vars)),
                    lab = ifelse(!off, "", sprintf("%.2f%s", rho, ifelse(fdr < 0.05, "*", ""))))
write.csv(cr, file.path(RES, "04_immune/tas_clinical_immune_spearman.csv"), row.names = FALSE)
cat(sprintf("n=%d with stage; TAS-stage rho %.2f; TAS-age rho %.2f\n", sum(!is.na(d$Stage)),
            cr$rho[cr$a == "TAS" & cr$b == "Stage"], cr$rho[cr$a == "TAS" & cr$b == "Age"]))

p <- ggplot(cr, aes(a, b, fill = ifelse(off, rho, NA))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = lab), size = 2.2) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0,
                       limits = c(-0.6, 0.6), na.value = "grey92", name = expression(rho),
                       guide = guide_cbar()) +
  labs(x = NULL, y = NULL, title = "TAS, age, stage, and six immune cell types (Spearman)") +
  theme_cns(base_size = 8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
save_pub(p, "FigS_tas_clinical_immune", OUT, w = 150, h = 125)

# ---- CellChat composite ---------------------------------------------------
tag <- function(img, t) {
  image_annotate(image_border(img, "white", "0x60"), t, size = 90, weight = 700,
                 gravity = "northwest", location = "+20+5")
}
rd <- function(f, w) image_scale(image_read(file.path(OUT, f)), as.character(w))
W <- 2100
a <- tag(rd("FigS_cellchat_outgoing.png", W / 2), "A")
b <- tag(rd("Fig12C_outgoing_circles.png", W), "B")
c_ <- tag(rd("Fig12E_outgoing_counts.png", W), "C")
dd <- tag(rd("Fig12D_pathway_circles.png", W), "D")
top <- image_append(c(a, image_blank(W / 2, image_info(a)$height, "white")))
comp <- image_background(image_append(c(top, b, c_, dd), stack = TRUE), "white")
image_write(comp, file.path(OUT, "FigS_cellchat.png"), density = 300)
image_write(comp, file.path(OUT, "FigS_cellchat.pdf"), format = "pdf", density = 300)
cat("saved: FigS_cellchat (png+pdf)\n")
