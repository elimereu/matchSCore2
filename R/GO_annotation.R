#' Compute the weighted GO enrichment of cluster-specific markers.
#'
#' This function computes the weighted GO enrichment across computational gene
#' markers by using their p-values.
#'
#' @param markers A data.frame of cluster specific gene markers like in the
#' Seurat output of the function FindAllMarkers. The used columns are
#' `c("p_val","cluster","gene")`
#' @param go.db Human (`org.Hs.eg.db`) or mouse (`org.Mm.eg.db`) GO database.
#' Corresponding R packages have to be installed.
#' @param species character indicating the species. Only 'human' or 'mouse'.
#' @param ontology.type Ontology family to be examined ("BP"= 'Biological Process',
#' "MF"='Molecular Function', "CC"='Cellular Component')
#' @param reformat.gene.names TODO
#' @param go.score.class TODO
#' @param p.val.threshold TODO
#' @param dag.file.prefix TODO
#' @param ngenes TODO
#'
#' @return A list with two elements:
#' - `GOenrich` contains data.frames with GO enrichments for each cluster.
#' - `genes_byGO` provides the set of genes from each resulting GO pathway that
#' are in overlap between your data and the pathway.
#'
#' @export
#'
#' @examples
#' # TODO
GOannotation <- function(markers,
                         go.db,
                         species = "mouse",
                         ontology.type = "BP",
                         reformat.gene.names = FALSE,
                         go.score.class = "weight01Score",
                         p.val.threshold = 0.05,
                         dag.file.prefix = FALSE,
                         ngenes = 20) {
  if (!ontology.type %in% c("BP", "MF", "CC")) {
    stop("Only 'BP', 'CC' and 'MF' are supported as ontology types")
  }
  if (!is.character(go.db)) {
    go.db <- go.db$packageName
  }

  if (reformat.gene.names) {
    message("Reformatting gene names (to lower case names starting with capital)")
    .simpleCap <- function(x) {
      s <- strsplit(x, " ")[[1]]
      paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " "
      )
    }
    gene.names <- sapply(tolower(gene.names), .simpleCap)
  }

  ks.test <- new(go.score.class,
    testStatistic = topGO::GOKSTest,
    name = "KS", scoreOrder = "decreasing"
  )
  results <- list()
  tables <- list()
  genes_byGO <- list()
  # selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
  allGO2genes <- annFUN.org(whichOnto = ontology.type, feasibleGenes = NULL, mapping = go.db, ID = "symbol")

  cluster <- factor(markers$cluster)

  if (species == "mouse") {
    db <- org.Mm.egALIAS2EG
    godb <- org.Mm.egGO2ALLEGS
  } else {
    db <- org.Hs.egALIAS2EG
    godb <- org.Hs.egGO2ALLEGS
  }
  db <- list2env(as.list(db))
  godb <- list2env(as.list(godb))

  ids <- unlist(lapply(mget(markers$gene, db, ifnotfound = NA), function(x) x[1]))
  rids <- names(ids)
  names(rids) <- ids
  # list all the ids per GO category
  go.env <- eapply(godb, function(x) as.character(na.omit(rids[x])))
  # library(GO.db)
  # desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
  # names(go.env) <- paste(names(go.env), desc)

  for (clus in levels(cluster)) {
    message(paste(
      "Computing GO enrichment for cluster:", clus,
      "\n"
    ))
    df <- markers[which(markers$cluster == clus), ]
    all.genes <- df$p_val
    names(all.genes) <- df$gene

    go.data <- new("topGOdata",
      ontology = ontology.type,
      allGenes = all.genes,
      annot = annFUN.GO2genes,
      GO2genes = allGO2genes,
      geneSel = function(x) {
        length(x)
      },
      nodeSize = 8
    )

    cor.p.val.threshold <- p.val.threshold
    ks.results <- topGO::getSigGroups(go.data, ks.test)
    if (is.character(dag.file.prefix)) {
      topic <- "" # TODOELI: do we need this?
      topGO::printGraph(go.data, ks.results,
        firstSigNodes = sum(score(ks.results) <
          cor.p.val.threshold), fn.prefix = paste(dag.file.prefix,
          "_T", topic, "_", ontology.type, "_DAG",
          sep = ""
        ),
        useInfo = "def", pdfSW = TRUE
      )
    }
    results <- append(results, list(score(ks.results)))
    t <- topGO::GenTable(go.data,
      weightKS = ks.results, numChar = 512,
      topNodes = sum(score(ks.results) < cor.p.val.threshold)
    )
    names(t)[3] <- "Total"
    names(t)[6] <- "p-Value"

    gt <- sapply(t[, "GO.ID"], function(x) intersect(go.env[[grep(x, names(go.env))]], df$gene))
    gt <- lapply(gt, function(x) x[1:ngenes])
    gt <- lapply(gt, function(x) x[!is.na(x)])
    gt2 <- unlist(lapply(gt, function(x) paste(x, collapse = ",")))
    t2 <- data.frame(t[c(1, 2, 3, 6)], Genes = gt2)
    tables <- append(tables, list(t2))
    genes_byGO <- append(genes_byGO, gt)
  }

  names(tables) <- levels(cluster)

  return(
    list(
      GOenrich = tables,
      genes_byGO = genes_byGO
    )
  )
}
