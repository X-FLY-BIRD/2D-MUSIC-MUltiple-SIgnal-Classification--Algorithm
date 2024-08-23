# 2D-MUSIC-MUltiple-SIgnal-Classification--Algorithm
MUISC algorithm for two-dimensional direction-of-arrival (DOA) estimation.
## Introduction
MUSIC (**MU**ltiple **SI**gnal **C**lassification) is one of super-resolution algorithms for estimating DoA (Direction of Arrival) of multiple narrowband sources arriving at a antenna/sensor array. ULA (Uniform Linear Array) can only detect angle of incident wave with the line along ULA. 2D-MUSIC-Algorithm, however, can detect azimuth and pitch of incident wave simultaneously, which requires that the sensor array is a planar array (e.g., L-shaped, square, circular array) or a more complex spatial array. A common MUSIC algorithm includes the following steps: computing covariance matrix of receiving signals, eigenvalue decomposition and 2D-Spectral-Peak searching. 2D-Spectral-Peak searching typically requires a large amount of computational resources and storage capacity, especially in high-resolution situation. Based on L-shaped sensor array, this demo can significantly reduce computational and storage load by converting 2D-Searching into two 1D-Searching and adopting "local searching" strategy.  
## Algorithm Interpretation
Assuming a L-shaped sensor array is composed of two mutually perpendicular ULAs, $Mx$ sensors located on x-axis and $My$ sensors on y-axis. Computing the noise subspace of all these $(Mx + My - 1)$ signals and performing 2D-Searching will resolve DoA of sources. Now consider the noise subspace of $Mx$ signals received by ULA on x-axis and perform 1D-Searching will resolve the angle of sources with x-axis, namely $\alpha$. Angle of sources with y-axis, $\beta$, will also be resolved when signals from ULA on y-axis are used. $\alpha$ and $\beta$ can be expressed as 
<div align="center">
  <p> $\alpha = \langle sources, \hspace{4pt} \boldsymbol{x} \rangle$ <br> 
  $\beta = \langle sources, \hspace{4pt} \boldsymbol{y} \rangle$ </p>
</div>

The angle pair ($\alpha$, $\beta$) will uniquely determine azimuth $\phi$ and pitch $\theta$ of sources in the upper-half-space of 3D space, through the following transformation
<div align="center">
  <p> $\phi = arctan(\frac{cos(\beta)}{cos(\alpha)})$ <br>
  $\theta = arccos(\sqrt{1 - cos^2(\alpha) - cos^2(\beta)})$ </p>
</div>

By performing two 1D-Searching, a simple transformation and local searching, we resolve the DoA of sources.  
## Main Program Interpretation
In this main program, a L-shaped sensor array with sampling frequency of $10^7 \hspace{2pt} HZ$ is composed of $Mx$ sensors on x-axis and $My$ sensors on y-axis. It is assumed that there are $N$ independent sources that generate narrowband complex exponential signals with constant frequency of $10^6 \hspace{2pt} HZ$ and random amplitude following a standard Gaussian distribution.  
_Steps a ~ e_ in the main program provide an intuitive demonstration of 2D-MUSIC algorithm and plot 2D MUSIC pseudo-spectrum function with respect to $\phi$ and $\theta$. _Step f and g_ demonstrates how to estimate $\alpha$ and $\beta$ by using signals received from ULA on x-axis and y-axis, respectively. _Step h_ performs local fine-searching on 2D-spectrum function by adopting "four directions" searching strategy, of which the principle is also described in main program.  
## Acknowledgements
Thanks to the course project of repository [doa-estimation-music](https://github.com/msamsami/doa-estimation-music) for providing sources' data emulation and code framework. Contact _xkyshangjiaoda@gmail.com_ if there is any confusions or improvment suggestions. 123
