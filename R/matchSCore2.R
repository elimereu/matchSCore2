#' Compute the cellular phenotype matching across reference cell groups and
#' computational clusters
#'
#' This function computes the matchSCore across reference gene markers and
#' cluster-specific genes computationally obtained from the tool.
#'
#' @param gene_cl.ref A list of reference group markers. Each element of the
#' list represents the set of markers for each biological group.
#' @param gene_cl.obs A list of cluster-specific gene sets to match with `gene_cl.ref`.
#' @param ylab A character label (ylab) indicating the reference technology to plot
#' in the matchSCore table(e.g. "Smart-Seq2").
#' @param xlab A character label (xlab) indicating the tested technology to plot
#' in the matchSCore table (e.g. "Chromium").
#'
#' @return A list with the `matchSCore` matrix and the ggplot table at the
#' `ggplot` slot.
#'
#' @export
#'
#' @examples
#' # TODO
matchSCore2 <- function(gene_cl.ref,gene_cl.obs,ylab,xlab){

  score=0
  lab.ref=seq(1:length(gene_cl.ref))
  lab.obs=seq(1:length(gene_cl.obs))
  anno.lab=vector()
  max_ji=vector()
  ji_mat=vector()

  for(i in 1:length(lab.ref)){
    len1=length(gene_cl.ref[[i]])
    JI=vector()
    for(ind in 1:length(lab.obs)){
      len2=length(gene_cl.obs[[ind]])
      I=length(intersect(gene_cl.ref[[i]],gene_cl.obs[[ind]]))
      J=I/(len1+len2-I)
      JI=append(JI,J,length(JI))
    }
    score=sum(score,max(JI))
    anno.lab = append(anno.lab,which(JI==max(JI)),after = length(anno.lab))
    max_ji = append(max_ji,max(JI),after = length(max_ji))
    ji_mat =rbind(ji_mat,JI)
  }
  score=score/length(gene_cl.ref)
  colnames(ji_mat)=names(gene_cl.obs)
  rownames(ji_mat)=names(gene_cl.ref)
  gg<-summary_ggplot(data=ji_mat,ylab,xlab)

  return(list(matchScore=score,labels=anno.lab,max_JI=max_ji,JI.mat=ji_mat,ggplot=gg))

}

