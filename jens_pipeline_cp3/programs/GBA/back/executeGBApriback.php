<html>
<head>
<link type="text/css" rel="stylesheet" href="style.css">
</head>
<body>
<?php
	
	set_time_limit(0);
	//Get the time from the epoch
	$reqTime = $_SERVER['REQUEST_TIME'];
	// So $resultFile stores the name of the file that will contans the result. 
	
	//Variables
	$storeDir = 'filestore/';
	$swissDir = 'swissprotLearnedMatrices/wForgetRate/normalized/';

	if('text' == $_REQUEST['inputMode'])
	{
		$sequence =$_REQUEST['sequence'];
		
		//So $trimmedseq contains the formatted string
		// Now save this in a file
		//The name should be consistent: destFile
		$destFile = $storeDir.'destFile.txt';
		//open the File
		@$fp = fopen($destFile, 'w');
		//Handling error message
		if(!$fp)
		{
			echo '<p><strong> Your request cann\'t be processed at this time.<br />
			Please try again later.</p>';
		}
		//Write $sequence into the file
		fwrite($fp, $sequence, strlen($sequence));
		//Close it
		fclose($fp);
	}
	else if('file' == $_REQUEST['inputMode'])
	{
		if($_FILES['datafile']['error'])
		{
			echo 'Problem: ';
			switch($_FILES['datafile']['error'])
			{
				case 1: echo 'File exceeded upload_max_filesize'; break;
				case 2: echo 'File exceeded max_file_size'; break;
				case 3: echo 'File only partially uploaded'; break;
				case 4: echo 'No file uploaded'; break;
			}
			exit;
		}
		//Does the file have right MIME type?
		
		if($_FILES['datafile']['type'] != 'text/plain')
		{
			echo "Error type: ".$_FILES['datafile']['error'];
			echo 'Problem: file is not plain text';
			exit;
		}
		

		$destFile = $storeDir.$_FILES['datafile']['name'];
		if(is_uploaded_file($_FILES['datafile']['tmp_name']))
		{
			if(!move_uploaded_file($_FILES['datafile']['tmp_name'], $destFile))
			{
				echo 'Problem: could not move file to the destination directory';
				exit;
			}
			
		}
		else
		{
			echo 'Problem: Possible file upload attack. Filename: ';
			echo $_FILES['datafile']['name'];
			exit;
		}

		//File uploaded successfully
		//echo 'File uploaded successfully';
	}
	//Get the other variables
	$swiss =$_REQUEST['swiss'];
	$t1 = $_REQUEST['t1'];
	$t2 = $_REQUEST['t2'];
	$t3 = $_REQUEST['t3'];
	$eaddr = $_REQUEST['eaddr'];
	//Swissprot

	$swissProt = $swissDir.$swiss;
	
	//create the partpath
	
	$tenpath = $_SERVER['PHP_SELF'];
	$pathpieces = explode('/', $tenpath);
	$partpath="";
	for($iter=1; $iter < sizeof($pathpieces) -2; $iter=$iter+1)
	{
		$partpath=$partpath."/".$pathpieces[$iter];
	}
	$server=$_SERVER['SERVER_NAME'];
	$totalpath = 'http://'.$server.$partpath.'/resultDir/';
	//echo "totalpath: ".$totalpath."\n";

	//Now save five variables in a files.
	//Open the arguement fiile
	$argstr=$destFile." ".$swissProt." ".$t1." ".$t2." ".$t3." ".$eaddr." ".$totalpath;
	$commandstr="php executeGBAsec.php ".$argstr." 2>&1 > /dev/null &";
	
	@$fscript = fopen("script", "w");
	if(!$fscript)
		{
			echo '<p><strong> Your request cann\'t be processed at this time.<br />
			Please try again later.</p>';
		}
	fwrite($fscript, $commandstr, strlen($commandstr));
	fclose($fscript);
	exec("sh script 2>&1 >/dev/null &");
	//echo $commandstr;
	echo "<br />";
	//exec ($commandstr);
	echo '<h2 align="center">Your query is being processed. A link to the result will be mailed to your email address.</h2>';
	//We have to save all the variable names in a file in order and have to 
	//recover them in some other format
		
		
?>
</body>
</html>




