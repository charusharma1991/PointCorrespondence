# Simplicial Complex based Point Correspondence between Images warped onto Manifolds

### About
This is a source code for our [ECCV 2020](https://eccv2020.eu/) paper: [Simplicial Complex based Point Correspondence between Images warped onto Manifolds](https://arxiv.org/pdf/2007.02381.pdf). This solves assignment problem using higher order structures between images warped onto manifolds.

### Requirements
This code is developed in 
```
MATLAB R2017b
```
### Dataset
We used ten different datasets for Link Prediction Task and three datasets for Node Classification. The dataset names are:
- Power
- C.elegans
- USAir

### Usage
1. To run end to end matching model to generate results for "Desktop" dataset, do the following:
```
end2end.m
```
2. To reproduce results mentioned in Table 2, do the following:
```
main.m
```
This will generate results for our method (OurWarped). The results would be generated in the output files.

### Citation
Please cite the paper if you use this code.
```
@article{sharma2020simplicial,
  title={Simplicial Complex based Point Correspondence between Images warped onto Manifolds},
  author={Sharma, Charu and Kaul, Manohar},
  journal={arXiv preprint arXiv:2007.02381},
  year={2020}
}
```
