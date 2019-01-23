## Beyond Landscapes: An Exemplar-based Image colorization method ##
[Saulo Pereira](https://github.com/saulo-p)

![Results1](https://github.com/saulo-p/ImageColorization/blob/master/results/results1.png)
![Results2](https://github.com/saulo-p/ImageColorization/blob/master/results/results3.jpg)
![Results3](https://github.com/saulo-p/ImageColorization/blob/master/results/results4.jpg)

### Overview ###
This repository contains the source code for my proposed colorization method.
It consists of an exemplar-based method which is based on a classification framework and works at superpixel level for enhanced coherence and performance.

It was built upon the project I forked from, but at this stage it holds no resemblance.  

The method's pipeline includes third party code:
* The superpixel segmentation <Turbo>
* The saliency maps in <Saliency>
* Histogram matching functions 
* The scribble propagation algorithm <Levin>
* Parts from the colorization algorithm of <Gupta>


#### Usage ####
* The script ```single_colorization``` performs the colorization for the input pair and parameters defined in the file <./input/single.in>
* Features weights should be adjusted in order to achieve better results for each input pair.

The code was developed for experimental purposes, do not expect more than this.
