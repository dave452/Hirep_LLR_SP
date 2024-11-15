The is repository contains a modified version of HiRep (https://github.com/claudiopica/HiRep) containing the spn branch with modifications for the heatbath LLR algorthim.
<br/>
This repository contains the code for the LLR method used in the paper: Bennett, E., Lucini, B., Mason, D., Piai, M., Rinaldi, E., & Vadacchino, D. (2024). The density of states method for symplectic gauge theories at finite temperature. arXiv preprint arXiv:2409.19426. <br/>
<br/>
To replicate the results of the paper, download the data-release (10.5281/zenodo.13807598). <br/>
Build LLR_HB. <br/>
Run llr_hb on the relevant in the relevant folder containing 'input_file' and the folders 'Rep_*' where *=0,1,..,N_replicas-1 <br/>
The number of tasks = to the number of subdomains * the number of replicas <br/>
<br/>
To run more generally: <br/>
Build LLR_HB. <br/>
<br/>
Copy the files 'input_file', 'input_file_rep' and 'setup_replicas.sh' from the the folder 'LLR_HB/example/' to the desired location. <br/>
Create the file pre.dat containing: <br/>
E_0, a_0, dE <br/>
E_1, a_1, dE <br/>
... <br/>
E_N, a_N, dE <br/>
where E_n is the value of the centre of the interval n, a_n is the initial a_n value for the interval n and dE is the size of the interval.  <br/>
<br/>
In input_file:  <br/>
Set the lattice size with :GLB_* where *=T,X,Y,Z <br/>
The number of subdomains in each direction: NP_* where *=T,X,Y,Z <br/>
The number of replicas N_REP <br/>
 <br/>
In input_file_rep: <br/>
Set the number of initial thermalisations with: therm  <br/>
Set the number of RM iterations with: last conf <br/>
Set llr:make = true <br/>
Set llr:nmc to the number of configurations that should be measured each RM iteration. <br/>
Set llr:nth to the number of thermalisation steps each RM iteration. <br/>
Set llr:sfreq_fxa to be the number of fixed a_n iterations before an replica exchange is considered. <br/>
Set llr:nfxa to be the replpica exchange swaps in the fixed a_n iterations that are carried out. <br/>
Set llr:Smin to E_0 and llr:Smax to E_N. <br/>
<br/>
Run setup_replicas.sh <br/>
This will create a folder containing an input file for each replica Rep_0,..., Rep_N_replicas. <br/>
<br/>
Run llr_hb in this folder with the number of tasks = number of subdomains * the number of replicas. <br/>

