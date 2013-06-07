RequireVersion  ("2.11");

fprintf (stdout, 
"
============================= FUBAR v1.0 ============================

This file will perform a Fast Unbiased AppRoximate Bayesian (FUBAR)
analysis of a coding sequence alignment to determine whether some
sites have been subject to pervasive purifying or diversifying
selection. For details of the method and explanations of various
settings, please see http://www.hyphy.org/wiki/FUBAR

Please note that a FUBAR analysis generates many files in the same
directory as the original alignment. HyPhy needs to have write
privileges to this directory. For example if the original file is in
/home/sergei/FUBAR/data/pol.nex then at the end of a FUBAR run, there
will also exist files such as /home/sergei/FUBAR/data/pol.nex.samples,
/home/sergei/FUBAR/data/pol.nex.gridInfo etc. Many of these files can
be further examined for diagnostic and other purposes. They also 
provide checkpointing so that a partially completed analysis can be
restarted.

=====================================================================

");

/*
     Load the nucleotide file with multiple 
     trees/partitions if present 
*/

_cachingOK = 1;

ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{HYPHY_LIB_DIRECTORY[0][Abs(HYPHY_LIB_DIRECTORY)-2],"TemplateBatchFiles","TemplateModels","chooseGeneticCode.def"}}));
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");

_runAsFunctionLibrary = 0;
ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{HYPHY_LIB_DIRECTORY[0][Abs(HYPHY_LIB_DIRECTORY)-2],"TemplateBatchFiles","_MFReader_.ibf"}}));

filePaths = {"Base": LAST_FILE_PATH,
             "Nucleotide fit suffix": ".gtr_fit",
             "Codon fit suffix": ".codon_fit",
             "Grid information": ".grid_info",
             "Converged prior": ".thetaBar",
             "Output": ".fubar.csv",
             "SimGrid": ".sim_grid_info",
             "SimFitFile": ".sim_codon_fit"};
             
fprintf (stdout, "\n\nFUBAR will write intermediate and result files to\n", filePaths["Base"], ".extension\n\n");

//----------------------------------------------------------------------------
// PHASE 1: nucleotide fit
//----------------------------------------------------------------------------

_fubarNucFitLocation = filePaths["Base"] + filePaths["Nucleotide fit suffix"];

if (_cachingOK && !_fubarNucFitLocation) {
// file exists 
     fprintf (stdout, "[CACHED] FUBAR found the self-contained nucleotide fit file at ", _fubarNucFitLocation, "\n"); 
}
else
{
    _cachingOK = 0;
    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_1.bf"}}), {"0" : _fubarNucFitLocation});
    fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote the self-contained nucleotide fit file to ", _fubarNucFitLocation, "\n"); 
}

//----------------------------------------------------------------------------
// PHASE 2: branch length scaling and grid calculation
//----------------------------------------------------------------------------

fprintf (stdout, "\n\n");

_fubarCodonFitLocation = filePaths["Base"] + filePaths["Codon fit suffix"];
_fubarGridInfoLocation = filePaths["Base"] + filePaths["Grid information"];

if (_cachingOK && !_fubarGridInfoLocation && !_fubarCodonFitLocation) {
     fprintf (stdout, "[CACHED] FUBAR found the self-contained codon fit file at ", _fubarCodonFitLocation, "\n"); 
     fprintf (stdout, "[CACHED] FUBAR found the site likelihoods file at ", _fubarGridInfoLocation, "\n"); 
}
else
{
    _cachingOK = 0;
    _grid_points = prompt_for_a_value ("Number of grid points per dimension (total number is D^2)",20,5,50,1);
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will use a ", _grid_points , "X", _grid_points, " grid\n"); 
    
    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_2.bf"}}), {"0" : _fubarNucFitLocation,
                                                                                                                                 "1" : _fubarCodonFitLocation,
                                                                                                                                 "2" : _fubarGridInfoLocation,
                                                                                                                                 "3" : "" + _grid_points});
    fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote the self-contained codon fit file to ", _fubarCodonFitLocation, "\n"); 
    fprintf (stdout, "[DIAGNOSTIC] FUBAR wrote the the site likelihoods file to ", _fubarGridInfoLocation, "\n"); 
}

//----------------------------------------------------------------------------
// PHASE 3: MCMC
//----------------------------------------------------------------------------

fprintf (stdout, "\n\n");

_fubarWeightLocation = filePaths["Base"] + filePaths["Converged prior"];

if (_cachingOK && !_fubarWeightLocation) {
// file exists 
     fprintf (stdout, "[CACHED] FUBAR found the learned weights at ", _fubarWeightLocation, "\n"); 
}
else
{	
	_fubarMCMCBaseLocation = filePaths["Base"];
    _cachingOK = 0;
    _fubarPriorShape = prompt_for_a_value ("The concentration parameter of the Dirichlet prior",0.5,0.001,1,0);    
    fprintf (stdout, "[DIAGNOSTIC] FUBAR will use the Dirichlet prior concentration parameter of ", _fubarPriorShape, "\n"); 

    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_3.bf"}}), {"0" : _fubarMCMCBaseLocation,
                                                                                                                                 "1" : _fubarGridInfoLocation,
                                                                                                                                 "2" : "" + _fubarPriorShape
                                                                                                                                  });
}

//START FROM HERE
//----------------------------------------------------------------------------
// PHASE 4: PROCESSING & FDR Simulation
//----------------------------------------------------------------------------

fprintf (stdout, "\n\n");

_fubarResultLocation = filePaths["Base"] + filePaths["Output"];
_fubarSimGrid        = filePaths["Base"] + filePaths["SimGrid"];
_fubarSimFitFile     = filePaths["Base"] + filePaths["SimFitFile"];
_fubarConvergedPriorFile = filePaths["Base"] + filePaths["Converged prior"];


_fubar_do_simulations = 0;

if (_cachingOK && !_fubarResultLocation) {
     fprintf (stdout, "[CACHED] FUBAR found the results file at ",_fubarResultLocation  ,"\n"); 
}
else
{

    ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_PHASE_4.bf"}}), {"0" : _fubarNucFitLocation,
                                                                                                                                 "1" : _fubarGridInfoLocation,
                                                                                                                                 "2" : _fubarConvergedPriorFile,
                                                                                                                                 "3" : _fubarResultLocation
                                                                                                                                  });                                                                                                                                     
    fprintf (stdout, "\n[DIAGNOSTIC] FUBAR wrote the results of its analysis to ", _fubarResultLocation, "\n"); 
}


fubar_data = (ReadCSVTable (_fubarResultLocation, 1))[1]%4;

ExecuteAFile (Join(DIRECTORY_SEPARATOR,{{PATH_TO_CURRENT_BF[0][Abs(PATH_TO_CURRENT_BF)-2],"FUBAR_HBL","FUBAR_tools.ibf"}}));

fprintf (stdout, "\n[RESULTS] At posterior probability >= 0.9 ");

idx = Rows(fubar_data);
mean_pp = 0;

p_i = {};

while (fubar_data[idx-1][4] >= 0.9 && idx > 0) {
    mean_pp += (1-fubar_data[idx-1][4]);
    p_i + (1-fubar_data[idx-1][4]);
    idx += -1;
}

if (idx == Rows(fubar_data) ) {
    fprintf (stdout, "there were no sites under diversifying positive selection\n");
} else {
    detected = Rows(fubar_data)-idx;
    fprintf (stdout, "there were ", detected, " sites under diversifying positive selection.\n");
        fprintf (stdout, "\nCodon\tProb[dN/dS>1]\t\tEBF[dN/dS]>1\t\tE[Beta]/E[Alpha]");
        for (idx2 = Rows(fubar_data)-1; idx2 >= idx; idx2 += -1) {
            fprintf (stdout, "\n", fubar_data[idx2][0], "\t",  fubar_data[idx2][4], "\t",  fubar_data[idx2][5],"\t",fubar_data[idx2][2]/fubar_data[idx2][1]); 
        }    
    }
    fprintf (stdout, "\n");
}
