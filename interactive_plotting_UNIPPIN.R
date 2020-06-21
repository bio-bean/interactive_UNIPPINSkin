# required libraries
library(tidyverse)
library(visNetwork)
library(igraph)
# create node list with group size and node size (for node df)
# subppin is obtained from UNIPPIN_OmniPath_graph_generation.R
node_df <- get.vertex.attribute(subppin)
node_df2 <- as.data.frame(node_df)
node_df2$group <- 1
node_df2$node_size <- 1
# uniprot id queries - proteins that were searched in UNIPPIN_OmniPath_graph_generation.R
# this is only necessary if you want queries and neighbours to have different colours
queries <- c("P98066", "P00488", "P61587", "O60271")
# replace the node size with 2, will make these ones stand out
for (prot in queries) {
  r <- which(node_df2$name == prot)
  node_df2[r, "node_size"] <- 2
  node_df2[r, "group"] <- "Query"
}
node_df2$group[node_df2$group %in% "1"] <- "Neighbour"
node_df2 <- rename(node_df2, c("id" = "name", "label" = "prot_name"))
# add in uniprot ids where the name doesn't exist
for (r in 1:nrow(node_df2)) {
  uniprot_id <- node_df2[r, "id"]
  if (is.na(node_df2[r, "label"])) {
    node_df2[r, "label"] <- paste(uniprot_id)
  }
}
# create an edges_df using the igraph function
edge_df <- as_long_data_frame(subppin)
drops <- c("from", "to")
edge_df2 <- edge_df[, !(names(edge_df) %in% drops)]
# create new columns and populate them with the NA characters
edge_df2$arrows.from.enabled <- c(NA)
edge_df2$arrows.from.type <- c(NA)
edge_df2$arrows.to.enabled <- c(NA)
edge_df2$arrows.to.type <- c(NA)
# fill in the columns based on existing columns from the igraph object
for (r in 1:nrow(edge_df2)) {
  m <- which(edge_df2$omni_direction == "3")
  edge_df2[m, "arrows.from.enabled"] <- TRUE
  edge_df2[m, "arrows.from.type"] <- "arrow"
  edge_df2[m, "arrows.to.enabled"] <- TRUE
  edge_df2[m, "arrows.to.type"] <- "arrow"
  n <- which(edge_df2$omni_direction == "1")
  edge_df2[n, "arrows.from.enabled"] <- TRUE
  o <- which(edge_df2$omni_direction == "2")
  edge_df2[o, "arrows.to.enabled"] <- TRUE
  p <- which(edge_df2$omni_effect == "inhibition")
  edge_df2[p, "arrows.to.type"] <- "bar"
  edge_df2[p, "arrows.from.type"] <- "bar"
  q <- which(edge_df2$omni_effect == "stimulation")
  edge_df2[q, "arrows.to.type"] <- "arrow"
}
# rename columns - must be correct for visNetwork for automatic identification
edge_df2 <- rename(edge_df2, c("from" = "from_name", 
                               "to" = "to_name",
                               "label" = "omni_effect"))
edge_df2$title <- paste0("Interaction type : ", edge_df2$label, 
                         "<br> Database, PubMedID : ", edge_df2$omni_id)
# for reduced network containing only OmniPath edges
edge_df3 <- dplyr::filter(edge_df2, edge_df2$omni_direction > 0)
# plot - saves as a variable to write out later
network <- visNetwork(node_df2, # nodes data
                      edge_df3, # edges data
                      width = "100%", 
                      main = "UNIPPIN", # main title
                      submain = "", # write the proteins you were looking at here - subheading
                      footer = "Edges key: edges with arrows indicate a stimulation, 
           edges ending in a blunt end indicate an inhibition, 
           edges with no arrows describe an undefined interaction. PubMedIDs are available for interactions when edges are moused over") %>% # description, can be removed
  visPhysics(enabled = FALSE) %>% # turns off physics initially, can be toggled
  visOptions(highlightNearest = list(enabled = TRUE, # highlight the nearest nodes when a node is clicked
                                     hover = TRUE),
             selectedBy = list(variable = "group", # drop down option which allows highlighting of groups
                               multiple = TRUE)) %>% 
  visLegend(enabled = TRUE, main = "Key", zoom = TRUE) %>% # legend
  visConfigure(enabled = TRUE, filter = "physics") # toggle options on the html

# save - rename argument 2 to change name
visSave(network, "UNIPPIN_interactive.html", selfcontained = TRUE, background = "white")
