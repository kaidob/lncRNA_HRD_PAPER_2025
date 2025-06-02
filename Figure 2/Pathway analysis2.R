## Install required Bioconductor packages if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages
required_packages <- c("ReactomePA", "DOSE", "org.Hs.eg.db", "clusterProfiler", 
                       "enrichplot", "GOSemSim", "pathview")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

## Load necessary libraries
library(ReactomePA)
library(DOSE)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)

## Load gene cluster data
# Ersetzen Sie den Pfad mit Ihrem tatsächlichen Dateipfad
GenesClusterEnriched <- read.csv("/Volumes/GoogleDrive/My Drive/Work_MA/5 HDR_paper/HRD_prediction_2021/TCGAanalyze/cluster genes enrichment/GenesCLusterEnriched.csv")

## Preview data structure
print("Datenstruktur:")
head(GenesClusterEnriched)
print(paste("Anzahl Gene pro Cluster:"))
table(GenesClusterEnriched$cluster)

## Prepare gene lists for each cluster
# Gene IDs nach Cluster aufteilen
cluster_genes <- split(GenesClusterEnriched$number, GenesClusterEnriched$cluster)

# Konvertiere zu Character-Listen
cluster_genes <- lapply(cluster_genes, as.character)

# Namen für die Cluster setzen
names(cluster_genes) <- paste0("Cluster_", names(cluster_genes))

print("Cluster-Übersicht:")
sapply(cluster_genes, length)

## 1. KEGG Pathway Enrichment Analysis für alle Cluster
print("=== KEGG Pathway Enrichment Analysis ===")

# Einzelne KEGG-Analyse für jeden Cluster
kegg_results <- list()
for (i in names(cluster_genes)) {
  cat(paste("\nAnalysiere", i, "mit", length(cluster_genes[[i]]), "Genen...\n"))
  
  tryCatch({
    kegg_result <- enrichKEGG(
      gene = cluster_genes[[i]],
      organism = 'hsa',
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
    kegg_results[[i]] <- kegg_result
    
    if (nrow(kegg_result@result) > 0) {
      cat(paste("Gefunden:", nrow(kegg_result@result), "signifikante KEGG Pathways\n"))
    } else {
      cat("Keine signifikanten KEGG Pathways gefunden\n")
    }
  }, error = function(e) {
    cat(paste("Fehler bei", i, ":", e$message, "\n"))
    kegg_results[[i]] <- NULL
  })
}

# Vergleichende KEGG-Analyse
print("\n=== Vergleichende KEGG-Analyse ===")
compare_kegg <- compareCluster(
  geneCluster = cluster_genes,
  fun = "enrichKEGG",
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# Plots erstellen falls Ergebnisse vorhanden
if (nrow(compare_kegg@compareClusterResult) > 0) {
  # Barplot für Cluster-Vergleich
  print("Erstelle KEGG Barplot...")
  p1 <- enrichplot::barplot(compare_kegg, showCategory = 10) + 
    ggtitle("KEGG Pathway Enrichment - Cluster Vergleich") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p1)
  
  # Dotplot für Cluster-Vergleich
  print("Erstelle KEGG Dotplot...")
  p2 <- enrichplot::dotplot(compare_kegg, showCategory = 10) + 
    ggtitle("KEGG Pathway Enrichment - Dotplot")
  print(p2)
  
  # Heatplot falls verfügbar
  tryCatch({
    compare_kegg_sim <- pairwise_termsim(compare_kegg)
    p3 <- enrichplot::heatplot(compare_kegg_sim, showCategory = 10) + 
      ggtitle("KEGG Pathway Similarity Heatmap")
    print(p3)
  }, error = function(e) {
    cat("Heatplot konnte nicht erstellt werden:", e$message, "\n")
  })
} else {
  print("Keine signifikanten KEGG Pathways für Cluster-Vergleich gefunden")
}

## 2. Gene Ontology (GO) Enrichment Analysis
print("\n=== Gene Ontology Enrichment Analysis ===")

# GO Biological Process
compare_go_bp <- compareCluster(
  geneCluster = cluster_genes,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

if (nrow(compare_go_bp@compareClusterResult) > 0) {
  print("Erstelle GO Biological Process Plots...")
  p4 <- enrichplot::dotplot(compare_go_bp, showCategory = 8) + 
    ggtitle("GO Biological Process - Cluster Vergleich")
  print(p4)
  
  p5 <- enrichplot::barplot(compare_go_bp, showCategory = 8) + 
    ggtitle("GO Biological Process - Barplot") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p5)
}

# GO Molecular Function
compare_go_mf <- compareCluster(
  geneCluster = cluster_genes,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

if (nrow(compare_go_mf@compareClusterResult) > 0) {
  print("Erstelle GO Molecular Function Plots...")
  p6 <- enrichplot::dotplot(compare_go_mf, showCategory = 8) + 
    ggtitle("GO Molecular Function - Cluster Vergleich")
  print(p6)
}

## 3. Reactome Pathway Analysis
print("\n=== Reactome Pathway Analysis ===")

compare_reactome <- compareCluster(
  geneCluster = cluster_genes,
  fun = "enrichPathway",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

if (nrow(compare_reactome@compareClusterResult) > 0) {
  print("Erstelle Reactome Plots...")
  p7 <- enrichplot::dotplot(compare_reactome, showCategory = 8) + 
    ggtitle("Reactome Pathways - Cluster Vergleich")
  print(p7)
  
  p8 <- enrichplot::barplot(compare_reactome, showCategory = 8) + 
    ggtitle("Reactome Pathways - Barplot") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p8)
}

## 4. Disease Gene Network (DisGeNET) Analysis
print("\n=== Disease Gene Network Analysis ===")

compare_dgn <- compareCluster(
  geneCluster = cluster_genes,
  fun = "enrichDGN",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

if (nrow(compare_dgn@compareClusterResult) > 0) {
  print("Erstelle Disease Gene Network Plots...")
  p9 <- enrichplot::dotplot(compare_dgn, showCategory = 8) + 
    ggtitle("Disease Gene Network - Cluster Vergleich")
  print(p9)
}

