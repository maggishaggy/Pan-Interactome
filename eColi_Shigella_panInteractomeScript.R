library(WGCNA)
library(data.table)
library(splitstackshape)

ortho <- fread("~/Desktop/Documents/BNFO_620/interactomePipeline/pgat_ecoli_shigella_orthologyTable.txt")
colnames(ortho) <- as.character(ortho[1])
ortho <- ortho[-1]

o <- ortho[!`Gene name` %chin% "-"]
newTable <- o[,.(genbankNo = `Genbank accession`, geneName = `Gene name`, description = Description)]

for (col in 6:ncol(o)) {
  for (row in 1:nrow(o)) {
    if (o[[col]][row] == "") set(x = newTable, i = row, j = colnames(o)[col], value = 0)
    else if (o[[col]][row] != "") set(x = newTable, i = row, j = colnames(o)[col], value = 1)
  }
}

fwrite(newTable, "~/Desktop/Documents/BNFO_620/interactomePipeline/pgat_ecoli_shigella_orthologyTable_edited.csv")


lit <- fread("~/Desktop/Documents/BNFO_620/interactomePipeline/lit_curated_PPIs_gene_names.csv")
y2h <- fread("~/Desktop/Documents/BNFO_620/interactomePipeline/eColi_Y2H_PPIs_gene_names.csv")


proUnion <- funion(lit, y2h)
proUnion <- proUnion[!`gene symbol A` == `gene symbol B`]
colnames(proUnion) <- c("interactorA", "interactorB", "geneA", "geneB")
proUnion[,`:=` (interactorA = NULL, interactorB = NULL)]

truelength(proUnion)
alloc.col(proUnion, 1024)

fwrite(proUnion, "~/Desktop/Documents/BNFO_620/interactomePipeline/y2h_litCurated_eColiProteinInteractions.csv")


s1 <- cSplit(newTable, "geneName", sep = " ", direction = "long", stripWhite = TRUE)
s2 <- cSplit(s1, "geneName", sep = ",", direction = "long", stripWhite = TRUE)

s2[,geneName := gsub(pattern = "[^a-zA-Z0-9]", replacement = "", x = geneName)]

fwrite(s2, "~/Desktop/Documents/BNFO_620/interactomePipeline/pgat_ecoli_shigella_orthologyTable_geneNameStripEdited.csv")

#newTable <- copy(s2)

interactors <- 1:2
for (col in 4:ncol(newTable)) {
  colName <- colnames(newTable)[col]
  set(x = proUnion, i = NULL, j = colName,
      value = proUnion[,..interactors][,ifelse(Reduce(`&`, lapply(.SD, `%chin%`, newTable[get(colName) == 1, geneName])), 1, 0)])
  print(colName)
}

proUnion[,rowSums(.SD), .SDcols = 4:ncol(proUnion), keyby = .I][,.N, keyby = V1][V1 != 0][,sum(N)]

fwrite(proUnion, "~/Desktop/Documents/BNFO_620/interactomePipeline/y2h_litCurated_eColiShigella_interactomePresenceAbsence.csv")






superheat(proUnion[,..c], dist.method = "binary", pretty.order.rows = TRUE, pretty.order.cols = TRUE, bottom.label.text.angle = 90, bottom.label.text.size = 5, heat.col.scheme = "red", membership.cols = tCutColor)
superheat(proUnion[,..c], col.dendrogram = TRUE, dist.method = "binary", pretty.order.rows = TRUE, pretty.order.cols = TRUE, bottom.label.text.angle = 90, bottom.label.text.size = 2, heat.col.scheme = "red")
superheat(proUnion[counts != 0, ..c], dist.method = "binary", pretty.order.rows = TRUE, pretty.order.cols = TRUE, bottom.label.text.angle = 90, bottom.label.text.size = 4, heat.col.scheme = "red", membership.rows = cutColor, left.label.text.size = 3, membership.cols = tCutColor)

superheat(proUnion[counts != 0, ..c], dist.method = "binary", pretty.order.rows = TRUE, pretty.order.cols = TRUE, bottom.label.text.angle = 90, bottom.label.text.size = 4, membership.rows = cutColor, left.label.text.size = 3, membership.cols = tCutColor, column.title = "Strain Modules", row.title = "Interaction Modules", title = "Pan-Interactome Presence/Absence Heatmap with DTC Modules for E. Coli/Shigella Strains", heat.col.scheme = "red")




interaction_table <- fread("~/Desktop/PanInteractome/y2h_litCurated_eColiShigella_interactomePresenceAbsence.csv")
interactors <- 1:2
interaction_table[,!..interactors]
i_table <- interaction_table[,!..interactors]



