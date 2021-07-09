## Dependencies
Download the following and place them in folders adjacent to this repository:  
1. [slra](https://github.com/slra/slra "SLRA Git Repository")  
2. [Matlab version](https://www.mathworks.com/help/matlab/release-notes.html "Matlab Release Notes"): R2018a (MATLAB 9.4), February 23, 2018  

Example folder structure: 
>       folder/
>>      slra-opt/
>>      slra-master/
              
## Instructions
Clone this repository
```
git clone https://github.com/PericlesPet/slra-opt.git
cd slra-opt
```
And then:  

### Quickest way:

Run:
```
runEverything.m
```

### Otherwise

Run in the following order: 
```matlab
slraInit
GDescMain
panoc_lbfgs
panoc_fminlbfgs
custom_fminlbfgs
almTest
fminconTests
```


## Sources

Stella, Lorenzo, Andreas Themelis, Pantelis Sopasakis, and Panagiotis Patrinos. “A Simple and Efficient Algorithm for Nonlinear Model Predictive Control.” ArXiv:1709.06487 [Math], September 19, 2017. http://arxiv.org/abs/1709.06487 .

Kul-Optec/Nmpc-Codegen. C. 2017. Reprint, OPTEC, 2021. https://github.com/kul-optec/nmpc-codegen.

Nocedal, Jorge, and Stephen J. Wright. Numerical Optimization. 2nd ed. Springer Series in Operations Research. New York: Springer, 2006.

GitHub. “Jkaardal/Augmented-Lagrangian-Matlab-Octave.” (2017). https://github.com/jkaardal/augmented-lagrangian-matlab-octave.

Markovsky, Ivan, and Konstantin Usevich. “Software for Weighted Structured Low-Rank Approximation.” Journal of Computational and Applied Mathematics 256 (January 15, 2014): 278–92. <https://doi.org/10.1016/j.cam.2013.07.048>.

N. Boumal, B. Mishra, P.-A. Absil and R. Sepulchre, "Manopt a Matlab toolbox for optimization on manifolds", J. Mach. Learn. Res., vol. 15, pp. 1455-1459, Jan. 2014.

Nocedal, Jorge. “Updating Quasi-Newton Matrices with Limited Storage.” Mathematics of Computation 35, no. 151 (1980): 773–82. <https://doi.org/10.2307/2006193>.

Dirk-Jan Kroon (2021). FMINLBFGS: Fast Limited Memory Optimizer (<https://www.mathworks.com/matlabcentral/fileexchange/23245-fminlbfgs-fast-limited-memory-optimizer>), MATLAB Central File Exchange. Retrieved June 22, 2021.
