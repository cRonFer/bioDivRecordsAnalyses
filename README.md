# bioDivRecordsAnalyses
This repository contains R Scripts to perform and reproduce all the analyses and figures of 'Exploring the impact of data curation criteria on the observed geographical distribution of mosses'
It also includes the shapefiles and data used.

**Cite as:**
Ronquillo, C., Stropp, J., Medina, N.G. & Hortal, J. (2023). Exploring the impact of data curation criteria on the observed geographical distribution of mosses. _Ecology and Evolution_

# Workflow
The following **diagram** describes the workflow in which the scripts are organized to reproduce the methods sections of the paper.

```mermaid
flowchart TD;
  A[GBIF / CNABH Records] --> B(1. Records pre-processing);
  C(1b. iDigBio Records \n Preprocessing)--> B;
  E[iDigBio Records] --> C;
  D(1c. BIEN Records \n Preprocessing)--> B;
  F[BIEN Records] --> D;
  
  B --> G(2. Taxonomic Standardization);
  G --> H[Clean Dataset];
  H --> I(3. Records Quantification);
  H --> K(4. Inventory Completeness Analyses);
  K --> L[Inventory Completeness \n Estimators];
  L --> M(5. Quantification of Well-Sample Cells);
  
  I --> O(7. Figure / Maps Code);
  L --> O(7. Figure / Maps Code);
  M --> O(7. Figure / Maps Code);
  
  L --> N(6. Piecewise Regression Analyses - \n Latitudinal Gradient);
  
  style A fill:#d67da9, stroke:#333,stroke-width:2px,color:#fff;
  style E fill:#d67da9, stroke:#333,stroke-width:2px,color:#fff;
  style F fill:#d67da9, stroke:#333,stroke-width:2px,color:#fff;
  style H fill:#f9e4a6, stroke:#333,stroke-width:2px,color:#333;
  style L fill:#f9e4a6, stroke:#333,stroke-width:2px,color:#333;
  style C fill:#c4b1ee, stroke:#c4b1ee,stroke-width:2px,color:#fff;
  style D fill:#c4b1ee, stroke:#c4b1ee,stroke-width:2px,color:#fff;
    
  subgraph Raw datasets of records occurrences by biodiversity database
  A & E & F
  end;
```
