# Our tasks in this problem are:

## Task 1: multimodal spectral clustering

### Inputs: 
- 1- Two modalities including images and texts (annotations or tags that describe images), with their labels.
- 2- Very limited correspondences samples between two modalities.
 
### What is the expected:
We expect that spectral clustering using approximate common eigenbases must be more efficient compared to spectral clustering on every single modality independently (using eigenbases of each modality).

### What we do:
- 1- Getting very limited correspondences (10% of the data). 
- 2- Expanding very limited correspondences to the sufficient number of new matching and mismatching samples.
- 3- Approximating common eigenbases considering all predetermined and new expanded matching and mismatching samples.
- 4- Performing the spectral clustering on every single modality independently, first using eigenbases of each modality and then using common eigenbases approximate on each modality (the best accuracy of spectral clustering with approximate common eigenbases of each modality is reported as the final accuracy of multimodal data).
 
### Evaluation criteria:
The accuracies of these clustering models are computed using the micro-averaged accuracy and normalized mutual information.


## Task 2: Multimodal classification

### Inputs: 
- 1- Two modalities including images and texts (annotations or tags that describe images), with their labels.
- 2-	Very limited correspondences samples between two modalities.

### What is the expected:
We expect that classification using approximate common eigenbases must be more efficiently compared to classification on each single modality independently (using eigenbases of each modality).

### What we do:
- 1-	Getting very limited correspondences (10% of the data).
- 2-	Expanding very limited correspondences to the sufficient number of new matching and mismatching samples.
- 3-	Approximating common eigenbases considering all predetermined and new expanded matching and mismatching samples.
- 4-	Computing the diffusion distances using the sufficient number of the first eigenbases of 1) Laplacian of each modality, and 2) approximate common eigenbases of the joint Laplacian eigenspaces.
- 5-	Performing the classification by applying the k NN classifier based on the diffusion distances on every single modality independently, first using its eigenbases and then using its approximate common eigenbases (the best accuracy of classification with approximate common eigenbases of each modality is reported as the final classification accuracy of multimodal data).

### Evaluation criteria:
The standard criterion used to measure the classification performances is the mean accuracy and standard deviation.
