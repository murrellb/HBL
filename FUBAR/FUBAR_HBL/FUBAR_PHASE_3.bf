fscanf              (stdin,"String", _sampleFile);
fscanf              (stdin,"String", _gridInfo);
fscanf              (stdin,"Number", _concentration);

ExecuteAFile        (PATH_TO_CURRENT_BF + "FUBAR_tools.ibf");
baseFilePath  		= PATH_TO_CURRENT_BF + "spool/"+_in_FilePath;

fscanf (_gridInfo, REWIND, "NMatrix,Raw", grid, gridInfo);
gridInfo = Eval(gridInfo);
    

points             = Rows(grid);
sites              = Columns(gridInfo["conditionals"]);
normalize_by_site  = ({1,points}["1"])*(gridInfo["conditionals"]);
normalized_weights = (gridInfo["conditionals"])*({sites,sites}["1/normalize_by_site[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"]);
sum_by_site        = normalized_weights * ({sites,1}["1"]);

//Trying out a Dirichlet without symmetry
priorvec		= {points,1}["_concentration"];
//priorvec		= {points,1}["3*_concentration*(grid[_MATRIX_ELEMENT_ROW_][0]==grid[_MATRIX_ELEMENT_ROW_][1])+_concentration"];


   
weights = {1,points}["1"];
weights = weights * (1/(+weights));
oldweights = Transpose(weights);

diffSum = 1;
iters=1;
t0 = Time (1);
while (diffSum > 0.0000001) {
	 phiUN = {points,sites}["normalized_weights[_MATRIX_ELEMENT_ROW_][_MATRIX_ELEMENT_COLUMN_]*weights[_MATRIX_ELEMENT_ROW_]"];
	 phiNormalizers  = ({1,points}["1"])*phiUN;
	 phi = phiUN*({sites,sites}["1/phiNormalizers[_MATRIX_ELEMENT_ROW_]*(_MATRIX_ELEMENT_ROW_==_MATRIX_ELEMENT_COLUMN_)"]);
	 //Would have thought this would be faster, but NUUUUUUU.
	 //phi = {points,sites}["phiUN[_MATRIX_ELEMENT_ROW_][_MATRIX_ELEMENT_COLUMN_]/phiNormalizers[_MATRIX_ELEMENT_COLUMN_]"];
	 weights = (phi * ({sites,1}["1"]))+priorvec;
	 weights = weights * (1/(+weights));
	 diffVec = weights - oldweights;
	 diffSum = +({points,1}["diffVec[_MATRIX_ELEMENT_ROW_]*diffVec[_MATRIX_ELEMENT_ROW_]"]);
	 SetParameter (STATUS_BAR_STATUS_STRING, "Iteration: "+ iters + " ------ delta:" + diffSum + " ----- Time:" + _formatTimeString(Time(1)-t0),0);
	 oldweights = weights;
	 iters = iters+1;
	 }


fprintf (stdout, "\nTime taken:" + _formatTimeString(Time(1)-t0), "\n");

convergefile = _sampleFile + ".thetaBar";
fprintf (convergefile,CLEAR_FILE, weights);
fprintf (_sampleFile +"_posterSurface.csv",CLEAR_FILE,vectorToMatrixCSVstring(oldweights,20),"\n");
fprintf (_sampleFile +"_priorSurface.csv",CLEAR_FILE,vectorToMatrixCSVstring(priorvec*(1/(+priorvec)),20),"\n");