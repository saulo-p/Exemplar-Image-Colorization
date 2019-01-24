## Beyond Landscapes: An Exemplar-based Image colorization method ##
[Saulo Pereira](https://github.com/saulo-p)

![Results1](https://github.com/saulo-p/ImageColorization/blob/master/results/results1.png)
![Results2](https://github.com/saulo-p/ImageColorization/blob/master/results/results3.png)
![Results3](https://github.com/saulo-p/ImageColorization/blob/master/results/results4.png)

### Overview ###
This repository contains the source code for my proposed colorization method.
It consists of an exemplar-based method which is based on a classification framework and works at superpixel level for enhanced coherence and performance.

It was built upon the project I forked from, but at this stage it holds no resemblance.  

The method's pipeline includes third party code:
* The superpixel segmentation in:
> Levinshtein, A., Stere, A., Kutulakos, K. N., Fleet, D. J., Dickinson, S. J., & Siddiqi, K. (2009). Turbopixels: Fast superpixels using geometric flows

* The saliency maps in:
> Yang, C., Zhang, L., Lu, H., Ruan, X., & Yang, M. H. (2013). Saliency detection via graph-based manifold ranking.

* Histogram matching functions
> http://cvhci.anthropomatik.kit.edu/~bschauer/

* The scribble propagation algorithm:
> Levin, A., Lischinski, D., & Weiss, Y. (2004, August). Colorization using optimization.

* Parts from the colorization algorithm of:
> Gupta, R. K., Chia, A. Y. S., Rajan, D., Ng, E. S., & Zhiyong, H. (2012, October). Image colorization using similar images


#### Usage ####
* The script ```single_colorization``` performs the colorization for the input pair and parameters defined in the file <./input/single.in>
* Features weights should be adjusted in order to achieve better results for each input pair.

The code was developed for experimental purposes, do not expect more than this.
