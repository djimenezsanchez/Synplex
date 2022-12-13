# Synplex: A synthetic simulator of highly multiplexed histological images.
***Summary:*** Multiplex tissue immunostaining is a technology of growing relevance as it can capture in situ the complex interactions existing between the elements of the tumor microenvironment. The existence and availability of large, annotated image datasets is key for the objective development and benchmarking of bioimage analysis algorithms. Manual annotation of multiplex images, is however, laborious, often impracticable. In this paper, we present Synplex, a simulation system able to generate multiplex immunostained in situ tissue images based on user-defined parameters. This includes the specification of structural attributes, such as the number of cell phenotypes, the number and level of expression of cellular markers, or the cell morphology. Synplex consists of three sequential modules, each being responsible for a separate task: modeling of cellular neighborhoods, modeling of cell phenotypes, and synthesis of realistic cell/tissue textures. Synplex flexibility and accuracy are demonstrated qualitatively and quantitatively by generating synthetic tissues that simulate disease paradigms found in the real scenarios. Synplex is publicly available for scientific purposes, and we believe it will become a valuable tool for the training and/or validation of multiplex image analysis algorithms. See our [*paper*](https://ieeexplore.ieee.org/document/9508562) for further description of Synplex.  

© [Daniel Jiménez Sánchez - CIMA University of Navarra](https://cima.cun.es/en/research/research-programs/solid-tumors-program/research-group-preclinical-models-preclinical-tools-analysis) - This code is made available under the GNU GPLv3 License and is available for non-commercial academic purposes. 

## Usage
Synplex consists of three independent modules. They can be tuned to create custom multiplex immunofluorescence image. For each module a guide is provided:

- [1. Modeling cellular neighborhoods](README_1.md#Modeling-cellular-neighborhoods)
- [2. Modeling cell phenotypes](README_2.md#Modeling-cell-phenotypes)
- [3. Cell expression and virtual microscopy](README_3.md#Cell-expression-and-virtual-microscopy)

## Example data
Run the main.m file using MATLAB to generate synthetic multiplexed histological images. With the default parameters that are set in Parameters_Neighborhoods.m, Parameters_Phenotypes.m, and Parameters_Texture.m the simulator will generate images like the following one: 

![image](https://user-images.githubusercontent.com/43730952/183055649-571b5bfc-83c6-4d77-a594-8415fe30af33.png)

where blue is DAPI (staining nuclear DNA), red is cytokeratin (staining tumor cells), and green is CD8 (staining T cells)

## Simulation of Endometrial Carcinoma from real data
We use real multiple immunofluorescence data (https://zenodo.org/record/4630664#.Y5iFFHbMJD8) to extract tumor features and use them to generate synthetic data.
[Modeling from real data](README_4.md#Modeling-from-real-data)

## Citation
Please cite this paper in case our method or parts of it were helpful in your work.
```diff
@article{jimenez2021synplex,
  title={Synplex: a synthetic simulator of highly multiplexed histological images},
  author={Jiménez-Sánchez, Daniel and Ariz, Mikel and Ortiz-de-Solórzano, Carlos},
  journal={2021 IEEE EMBS International Conference on Biomedical and Health Informatics (BHI)},
  year={2021}
}
```
