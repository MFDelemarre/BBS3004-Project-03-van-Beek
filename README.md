Scripts for Post-GWAS analysis.
===============================
R-Studio is needed to run this
CytoScape is used for the networks.
===============================
There are 2 main pathways.
  1. FUMA preperation
  2. Network preperation and integration
     - DisGeNET analysis
===============================
- FUMA preperation
Is pretty straightforward
Just look up what the columns-names supposed to be 
===============================
- Network preperation and integration
In 6 main steps
    1. READ FUMA OUTPUT
    2. Intergrating different databases
        - Biomart
        - STRING DB
        - MSigDB/Gene Set Enrichtment Analysis
    3. Enrichment Analysis
        - Gene Ontology
        - Disease Ontology
    4. Creation of Cytoscape-network
        - Includes colouring and clustering
    5. Enriching and plotting per cluster
        - EnrichGO & EnrichDGN 
    6.  DisGeNET intergration
   
