% LTE system simulator main simulator file. Check the LTE_sim_batch files
% to check how to launch the simulator.
% [] = LTE_sim_main()
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2008/08/11
% last changes: 2008/09/02  Colom Ikuno HARQ is implemented
%               2008/09/04  Bosanska    changed to user_biterror and cell_biterrors... 
%                                       for multi-user implementation
%               2008/09/08  Bosanska    added structure Results.UE_specific/CELL_specific
%               2008/09/10  Bosanska    HARQ adapted for multiple users
%               2008/09/15  Bosanska    changed structure of BS_output
%               2008/09/18  Bosanska    changed the name of chan_output -> ChanMod_output
%                                       added function [ChanMod_output] = LTE_channel_matrix(ChanMod);
%               2008/10/02  Bosanska    changed structure to multiple call of function LTE_RX 
%                                       for multi-user scenario (parallel receivers)
%                                       changed function LTE_RX, new inputs -> [1x1]struct BS_UE_specific
%                                                                              [1x1]double uu
%               2008/10/21  Bosanska    changed function LTE_channel_matrix, new inputs -> [1xnUE]struct UE_output (also output) 
%                                                                                          [1x1]double uu, [1x1]double SNR
%                                       prepared loop over number of users for LTE_channel_matrix 
%                                       (multi-user scenario with different channels for users) 
%               2008/12/04  Bosanska    commented loop over number of users for LTE_channel_matrix
%               2009/02/02  Simko       scheduler plot
%               2009/03/03  Simko       multi user channel- every user has different channel
%               2009/03/09  Jcolom      Added uplink delay
%               2009/16/09  Jcolom      Splitted tracing and plotting from the main simulation loop (better code readability)
%               2009/05/15  Jcolom      Public release of the simulator (r400).
%
%
%
% By using this simulator, you agree to the license terms stated in the license agreement included with this work.
% If you are using the simulator for your scientific work, please reference:
%
% BibTeX:
% @InProceedings{EUSIPCO2009,
%   author =        {Christian Mehlf\"uhrer and Martin Wrulich and Josep Colom Ikuno and Dagmar Bosanska and Markus Rupp},
%   title =         {Simulating the Long Term Evolution Physical Layer},
%   booktitle =     {Proc. of the 17th European Signal Processing Conference (EUSIPCO 2009)},
%   month =         aug,
%   year =          2009,
%   address =       {Glasgow, Scotland},
%   note =          {accepted for publication},
% }
% 
% ASCII
% C. Mehlführer, M. Wrulich, J. C. Ikuno, D. Bosanska and M. Rupp, "Simulating the Long Term Evolution Physical Layer,"
% in Proc. of the 17th European Signal Processing Conference (EUSIPCO 2009), Aug. 2008, Glasgow, Scotland

fprintf('Vienna LTE Link Level simulator\n');
fprintf('(c) 2008, INTHFT, TU Wien\n');
fprintf(' This work has been funded by A1 Telekom Austria AG and the Christian Doppler Laboratory for Design Methodology of Signal Processing Algorithms.\n\n');
fprintf('  By using this simulator, you agree to the license terms stated in the license agreement included with this work\n');
fprintf('  Contains code from:\n');
fprintf('    - pycrc (CRC checking)\n');
fprintf('    - The Coded Modulation Library (convolutional coding & SISO decoding)\n');
fprintf('  Convolutional coding & SISO decoding MEX files under the GNU lesser GPL license\n\n');

% Parallel toolbox
switch LTE_params.simulation_type
    case 'normal'
        LTE_sim_main_single;
    case 'parallel'
        LTE_sim_main_par;
end

clear tmp_results
