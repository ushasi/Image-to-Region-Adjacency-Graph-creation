# Image-to-Region-Adjacency-Graph-creation
[Paper](https://reader.elsevier.com/reader/sd/pii/S1077314219300578?token=FF18DF6BD33340CB07179AE964A960F224B8A29AC597C8D4875C71AF366407364D877984BA4E1BA4FF97548B3C83FB2Ahttps://reader.elsevier.com/reader/sd/pii/S1077314219300578?token=FF18DF6BD33340CB07179AE964A960F224B8A29AC597C8D4875C71AF366407364D877984BA4E1BA4FF97548B3C83FB2A) | [MATLAB](https://www.mathworks.com/products/matlab.html)

Convertion of an RGB image to a Region Adjacency Graph (RAG) using SLIC super-pixel based segmentation technique.

1. segmentation.m - contains the type of segmentation algorithm. The code by default uses SLIC superpixel based segmentation. However, a graph-cut based segmentation implementation can also be found in the commented section.
2. Extractfeaturevec.m - contains the features to be extracted from each nodes. 
3. filterimage.m - This is the main file. We have used the single labelled UC-Merced dataset, containing 21 classes. Mention the path to these images in the "addpath" line. Totgra saves the final matrix containing the node features, wighted adjacency matrix information, etc.

#Requirements-
Set up vlfeat library for using the SLIC super-pixel based segmentation.

Find the UC-Merced dtaset from http://bigearth.eu/datasets.html


### Paper

*    The paper is also available at: [Siamese graph convolutional network for content based remote sensing image retrieval](https://reader.elsevier.com/reader/sd/pii/S1077314219300578?token=FF18DF6BD33340CB07179AE964A960F224B8A29AC597C8D4875C71AF366407364D877984BA4E1BA4FF97548B3C83FB2A)

*   If the work is any help to you, please feel free to cite the author:

```
@article{chaudhuri2019siamese,
  title={Siamese graph convolutional network for content based remote sensing image retrieval},
  author={Chaudhuri, Ushasi and Banerjee, Biplab and Bhattacharya, Avik},
  journal={Computer Vision and Image Understanding},
  volume={184},
  pages={22--30},
  year={2019},
  publisher={Elsevier}
}
