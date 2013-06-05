fscanf              (stdin, "String", nuc_fit_file);
fscanf              (stdin, "String", grid_file);
fscanf              (stdin, "String", weights_file);
fscanf              (stdin, "String", results_file);

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("WriteDelimitedFiles");

fscanf (grid_file, REWIND, "NMatrix,Raw", grid, site_probs);
fscanf (weights_file, REWIND, "NMatrix", learntWeights);

site_probs = Eval (site_probs);
sites   = Columns (site_probs["conditionals"]);

transWeights = Transpose(learntWeights);

positive_selection_stencil = {points,sites} ["grid[_MATRIX_ELEMENT_ROW_][0]<grid[_MATRIX_ELEMENT_ROW_][1]"];
negative_selection_stencil = {points,sites} ["grid[_MATRIX_ELEMENT_ROW_][0]>grid[_MATRIX_ELEMENT_ROW_][1]"];
diag_alpha = {points,points}["grid[_MATRIX_ELEMENT_ROW_][0]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"];
diag_beta  = {points,points}["grid[_MATRIX_ELEMENT_ROW_][1]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"];
    
norm_matrix         = (transWeights*site_probs["conditionals"]);
pos_sel_matrix      = (transWeights*(site_probs["conditionals"]$positive_selection_stencil) / norm_matrix);
neg_sel_matrix      = (transWeights*(site_probs["conditionals"]$negative_selection_stencil) / norm_matrix);
alpha_matrix        = ((transWeights*diag_alpha*site_probs["conditionals"])/norm_matrix);
beta_matrix         = ((transWeights*diag_beta*site_probs["conditionals"])/norm_matrix);

bySitePosSel = {sites,5};
for (s = 0; s < sites; s+=1) {
    	SetParameter (STATUS_BAR_STATUS_STRING, "Tabulating results for site "+ s + "/" + sites + " " + _formatTimeString(Time(1)-t0),0);
    	bySitePosSel [s][0] = alpha_matrix[s]; 
    	bySitePosSel [s][1] = beta_matrix[s];
    	bySitePosSel [s][2] = neg_sel_matrix[s];
    	bySitePosSel [s][3] = pos_sel_matrix[s];
    }

 
fubarRowCount     = Rows (bySitePosSel);
site_counter = {};
for (currentFubarIndex = 0; currentFubarIndex < fubarRowCount; currentFubarIndex += 1) {
    site_counter + (currentFubarIndex+1);
}

WriteSeparatedTable (results_file, {{"Codon","alpha","beta","Prob[alpha>beta]", "Prob[beta>alpha]"}}, bySitePosSel, site_counter, ",");
