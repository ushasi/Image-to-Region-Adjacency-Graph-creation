# Image-to-Region-Adjacency-Graph-creation
Convertion of an RGB image to a Region Adjacency Graph (RAG) using SLIC super-pixel based segmentation technique.

1. segmentation.m - contains the type of segmentation algorithm. The code by default uses SLIC superpixel based segmentation. However, a graph-cut based segmentation implementation can also be found in the commented section.
2. Extractfeaturevec.m - contains the features to be extracted from each nodes. 
3. filterimage.m - This is the main file. We have used the single labelled UC-Merced dataset, containing 21 classes. Mention the path to these images in the "addpath" line. Totgra saves the final matrix containing the node features, wighted adjacency matrix information, etc.

Requirements-
Set up vlfeat library for using the SLIC super-pixel based segmentation.

Find the UC-Merced dtaset from http://bigearth.eu/datasets.html

If you are using this code, please cite the paper: https://arxiv.org/abs/1904.04794 (CMIR-NET : A Deep Learning Based Model For Cross-Modal Retrieval In Remote Sensing)
