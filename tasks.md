# Our tasks in this problem are:

## Task 1: multimodal spectral clustering

## Inputs: 
1-	Two modalities, 1) images, and 2) texts (annotations or tags that describe images), with their labels.
2-	Very limited correspondences samples between two modalities.

## What is the expected:
We expect that spectral clustering using approximate common eigenbases must be more efficient compared to spectral clustering on every single modality independently (using eigenbases of each modality).

## What we do:
1- Getting very limited correspondences (10% of the data).
2-	Expanding very limited correspondences to the sufficient number of new matching and mismatching samples.
3-	Approximating common eigenbases considering all predetermined and new expanded matching and mismatching samples.
4-	Performing the spectral clustering on every single modality independently, first using eigenbases of each modality and then using common eigenbases approximate on each modality (the best accuracy of spectral clustering with approximate common eigenbases of each modality is reported as the final accuracy of multimodal data).

##Evaluation criteria:
The accuracy of these clustering models is computed using the micro-averaged accuracy and normalized mutual information.
