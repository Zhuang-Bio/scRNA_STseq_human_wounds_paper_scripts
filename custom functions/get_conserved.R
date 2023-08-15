####----Create function to get conserved markers for any given cluster----####
####----This is very useful when two clusters are very similar----
get_conserved <- function(cluster){
  FindConservedMarkers(hswound.combined.sct,
                       ident.1 = cluster,
                       grouping.var = "sample", #need to change
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
# Iterate function across desired clusters, for example: Cluster 17, 20
#conserved_markers <- map_dfr(c(17,20), get_conserved) #