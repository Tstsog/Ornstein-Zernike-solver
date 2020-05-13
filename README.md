# Ornstein-Zernike-solver

## Publisher
Author: Tsogbayar Tsednee (PhD)
Contact: tsog215@gmail.com

## Introduction: 
  
A matlab code oz_lj_ts_run_opt_one_parm.m computes a free parameter in Verlet bridge function using dp_c = dp_v consistency criteria for the one-component Lennard-Jones fluid by solving the Ornstein-Zernike equation with a Picard iteration technique. 
     
The code OZ_solver_LJ_py.m solves the OZ equation for the one-component Lennard-Jones liquid in the Percus‐Yevick approximation. 

The code OZ_solver_one_component_plasma_hnc.m solves the OZ equation for one-component plasma in the hypernetted chain approximation. 
 
## Requirement: 
Any version of Matlab 

## Implementation details and running

The code oz_lj_ts_run_opt_one_parm.m uses the fminsearch.m function from Matlab software to find a optimal value of a free parameter. A code oz_lj_ts.m solves the Ornstein-Zernike integral equation for the Lennard-Jones potential together with a one-parameter Verlet-modified closure. The code oz_lj_ts.m  uses a code lsint.m which calculates the sine transform of a function.    
You may download them and run the oz_lj_ts_run_opt_one_parm.m directly. You can download OZ_solver_LJ_py.m and OZ_solver_one_component_plasma_hnc.m, and run them directly. 

## Copyright / License 

These codes are published as freeware. Basically, you have the right to use them and modify any of them. 

Find the GNU General Public License at:
https://www.gnu.org/licenses/gpl-3.0.en.html
    
### Acknowledgments
Author is very thankful for advices of Professor T. Luchko at California State University, Northridge regarding to an usage of a fminsearch function of Matlab software. 






