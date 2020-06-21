# Interactive UNIPPIN code and examples

Welcome to the GitHub repository for code written to create interactive graphs for UNIPPIN interactions annotated with literature citations from OmniPath.

UNIPPIN details can be found at: http://fraternalilab.kcl.ac.uk/wordpress/unippin/

OmniPath details can be found at: http://omnipathdb.org/

## Contents
There are two scripts and 10 examples here: the two scripts detail how to create the initial igraph object from UNIPPIN data and annotate the edges using OmniPath, and the second details how to create the interactive graph. Both scripts are written in R and require light editing to specify query proteins and databases to interrogate.

## Using the networks
The examples show the final network as a .html file, these can be downloaded and explored. The networks have physics options to toggle at the bottom of the screen, the recommended setting to sort the networks is 'force2AtlasBased'. Nodes can be moved and placed if physics settings are not enabled, otherwise nodes will return to their original position. PubMedIDs can be found for edged that have an annotation (inhibition/stimulation/two-way interaction) when the edge is moused over. Please note that some edges to not have PubMedIDs even though they contain an annotation, in these cases the interaction is documented through the Human Signalling Network: http://www.bri.nrc.ca/wang/.

## Required packages
The required packages to use the R scripts are (including dependencies):

igraph: https://igraph.org/r/

visNetwork: https://datastorm-open.github.io/visNetwork/

tidyverse: https://www.tidyverse.org/

The scripts were built in R version 4.0.0 (2020-04-24) -- "Arbor Day" 
