# load packages
library(igraph)
library(tidyverse)

# UNIPPIN --------
# Format: "Node1 ppi Node 2 DB1 DB2 ..."
network_file = "./UniPPIN2020_ACC_noSelfLoop_origin.sif"

# !!! The new file format has one line of header, so add header = TRUE !!!
unippin_db <- read.table(network_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Load the list of edges. Take 1st and 3rd column (UniProIDs).
edges_matrix <- as.data.frame(unippin_db[,c(1,3)])

# Merge other columns in a char string to get the origin (takes a few sec)
dbs <- colnames(unippin_db)[4:length(unippin_db[1,])]
unippin_db$origin <- apply( unippin_db[,dbs], 1, paste, collapse = " " )
unippin_db$origin <- lapply(unippin_db$origin, function(x) {gsub("\\. ", "", x)})
unippin_db$origin <- lapply(unippin_db$origin, function(x) {gsub(" \\.", "", x)})

# create a ppin from the dataframe and list the number of nodes and edges
ppin <- graph_from_data_frame(edges_matrix, directed = TRUE)
print(paste0('Number of nodes: ', vcount(ppin)))
print(paste0('Number of edges: ', ecount(ppin)))

# ADD edge attribute (databases of origin of each interaction)
ppin <- set_edge_attr(graph = ppin, name = "origin", value = unippin_db$origin)

# ADD protein names on nodes (take a few mins)
nodes_ID <- get.vertex.attribute(ppin)[[1]]
nodes_name <- array(data = NA, dim = length(nodes_ID))
names_db <- read.table("human_prot_names.txt", sep = '\t', stringsAsFactors = FALSE)
for (p in names_db[,1]) {
  pID <- strsplit(p,' ')[[1]][1]
  if (sum(nodes_ID == pID)) {
    nodes_name[nodes_ID == pID] <- strsplit(p,paste0(pID," "))[[1]][2]
  }
}
ppin <- set_vertex_attr(graph = ppin, name = "prot_name", value = nodes_name)

# Example: get neighborood of a node  (order 1 = first shell, ...) # examples
n1 <- neighborhood(graph = ppin, order = 1, nodes = "P98066")[[1]] # tsg6
n2 <- neighborhood(graph = ppin, order = 1, nodes = "P00488")[[1]] # f13a1
n3 <- neighborhood(graph = ppin, order = 1, nodes = "P61587")[[1]] # rnd3
n4 <- neighborhood(graph = ppin, order = 1, nodes = "O60271")[[1]] # spag9

# Create subgraph only with those nodes
all_nodes <- c(n1, n2, n3, n4)
subppin <- induced.subgraph(graph = ppin, v = all_nodes)

# display the number of nodes and edges in the subgraph
print(paste0('Number of nodes: ', vcount(subppin)))
print(paste0('Number of edges: ', ecount(subppin)))

# plot initial graph
plot.igraph(x = subppin, 
            layout= layout_with_lgl, 
            vertex.label.color="purple",
            vertex.label = get.vertex.attribute(subppin,'prot_name'),
            edge.color="pink", 
            edge.width=0.1, 
            edge.arrow.size=0, 
            vertex.size = 5, 
            vertex.border = NA)


# OmniPath ------
# required packages
library(OmnipathR)
library(tidyverse)

# create a list of edge ids to check against the database for annotated interactions
edge_id <- as_edgelist(subppin, names = TRUE)

# import databases - all - see OmniPath.R vignette for instructions on how to download specific databases
# this step requires a stable internet connection and takes a while
# reccomended: download once and save to a df to call back when needed
interactions <- import_AllInteractions(from_cache_file = NULL,
                                       filter_databases = get_interaction_databases(),
                                        select_organism = 9606)
# remove any self interacting nodes - UNIPPIN does not contain these so it will mess up the arrays
interactions_dd <- interactions[!(interactions$target == interactions$source), ]
write.csv(interactions_dd, "interactions_omnipath_dd.csv", row.names = FALSE)
interactions <- read.csv("interactions_omnipath_dd.csv", header = TRUE)

# create a dataframe of the interactions present within the subnetworks, this includes the references from omnipath
omni_reactions <- data.frame()
for (i in 1:nrow(edge_id)) {
  A <- edge_id[i, 1]
  B <- edge_id[i, 2]
  true_rows <- dplyr::filter(interactions, (A == interactions$source | B == interactions$source) & 
                               (A == interactions$target | B == interactions$target))
  omni_reactions <- bind_rows(omni_reactions, true_rows)
}
# save as a csv file for reference later
write.csv(omni_reactions, "omni_reactions_subppin_spag-9.csv")
# omni_dir stores 0, 1, 2, 3 representing the direction of the arrow in the pathway
# 3 is a two-way interaction
# 2 is a stimulation
# 1 is an inhibition
# 0 is no recorded interaction
omni_dir <- array(data = NA, dim = nrow(edge_id))
for (r in 1:nrow(edge_id)) {
  A <- edge_id[r, 1]
  B <- edge_id[r, 2]
  true_rows <- dplyr::filter(interactions, (A == interactions$source | B == interactions$source) & 
                               (A == interactions$target | B == interactions$target))
  if (nrow(true_rows) > 1) {
    omni_dir[r] <- as.integer(paste("3"))
  }
  else if (sum(interactions$source == A & interactions$target == B) == 1) {
    omni_dir[r] <- as.integer(paste("2"))
  }
  else if (sum(interactions$source == B & interactions$target == A) == 1) {
    omni_dir[r] <- as.integer(paste("1"))
  }
  else {
    omni_dir[r] <- as.integer(paste("0"))
  }
}
# create labels from the omnipath database
# two way interaction is for when there was more than one entry of a pair (pairs appear once in the db per direction)
# stimulation/inhibition are recorded in separate columns
# labels are stored in the omni_effects array
omni_effects <- array(data = NA, dim = nrow(edge_id))
for (r in 1:nrow(edge_id)) {
  A <- edge_id[r, 1]
  B <- edge_id[r, 2]
  true_rows <- dplyr::filter(interactions, (A == interactions$source | B == interactions$source) 
                             & (A == interactions$target | B == interactions$target))
  if (nrow(true_rows) > 1) {
    omni_effects[r] <- paste("two-way interaction")
  }
  else if (sum(true_rows$consensus_stimulation) == 1) {
    omni_effects[r] <- paste("stimulation")
  }
  else if (sum(true_rows$consensus_inhibition) == 1) {
    omni_effects[r] <- paste("inhibition")
  }
}
# omni_ids is an array that PubMedIDs and other references are stored in from OmniPath
omni_ids <- array(data = NA, dim = nrow(edge_id))
for (r in 1:nrow(edge_id)) {
  A <- edge_id[r, 1]
  B <- edge_id[r, 2]
  true_rows <- dplyr::filter(interactions, (A == interactions$source | B == interactions$source) 
                             & (A == interactions$target | B == interactions$target))
  if (nrow(true_rows) > 1) {
    comb <- paste(true_rows$references, sep = "", collapse = ";")
    omni_ids[r] <- paste(comb)
  }
  else if (nrow(true_rows) == 1) {
    omni_ids[r] <- paste(true_rows$references)
  }
  else if (nrow(true_rows) == 0) {
    #do nothing
  }
}
# set the edge attributes
subppin <- set_edge_attr(graph = subppin, name = "omni_direction", value = omni_dir)
subppin <- set_edge_attr(graph = subppin, name = "omni_effect", value = omni_effects)
subppin <- set_edge_attr(graph = subppin, name = "omni_ids", value = omni_ids)
# igraph for R doesn't have a reliable function for exporting a graph, 
# would recommend not exporting the subppin graph and instead add code from interactive_plotting_UNIPPIN.R onto the end of this script
