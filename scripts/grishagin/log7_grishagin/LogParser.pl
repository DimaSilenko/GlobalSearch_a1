@files = <*.txt>;


open(RF, ">comm.csv");
if(-e $fileOut)
{
	print "add to end\n";
	open(W, ">>comm.csv");
}
else
{
	print "create new\n";
	open(W, ">comm.csv");
	print W "TaskNumber;N;FuncNumber;eps;r;min;IsSolvePerebor;timePerebor;IterPerebor;IsSolveAGP;timeAGP;IterAGP;\n";
}

foreach $f(@files)
{
	open (R, "$f");
	
	$n = -1;
	$TaskNumber     = -1;
	$FuncNumber     = -1;
	$timePerebor   = -1.0;
	$timeAGP   = -1.0;
	$min = -1.0;
	$IterPerebor = -1;
	$IterAGP = -1;
	$r = -1.0;
	$eps = -1.0;
	$isSolvePerebor = 0;
	$isSolveAGP = 0;
	
	$err = "";
	
	$lib = "";
	while(<R>)
	{
		chomp;
		$tmp = $_;
		
		if (($tmp =~ /^libPath = [.]?[\/]?(\w+)/))
		{
			$lib = $1;
		}
		if (($tmp =~ /^srun: error: (\w+)/))
		{
			$err = $1;
		}		
		if (($tmp =~ /^Dimension = (\d+)/))
		{
			$n = $1;
		}
		
		if (($tmp =~ /^Iteration Perebor = (\d+)/))
		{
			$IterPerebor = $1;
		}
		
		if (($tmp =~ /^Iteration AGP = (\d+)/))
		{
			$IterAGP = $1;
		}
		
		if($tmp =~ /^Solve time Perebor = (\d+[.]?\d*)/)
		{
			$timePerebor        = $1;
		}
		if($tmp =~ /^Solve time AGP = (\d+[.]?\d*)/)
		{
			$timeAGP        = $1;
		}
		
		if($tmp =~ /^min = ([+-]?\d+[.]?\d*)/)
		{
			$min        = $1;
		}
		
		if (($tmp =~ /^Task_number = (\d+)/))
		{
			$TaskNumber = $1;
		}
		if (($tmp =~ /^Function_number = (\d+)/))
		{
			$FuncNumber = $1;
		}
		
		if($tmp =~ /^r = (\d+[.]?\d*)/)
		{
			$r        = $1;
		}
		if($tmp =~ /^Epsilon = (\d+[.]?\d*)/)
		{
			$eps        = $1;
		}

		if (($tmp =~ /^Perebor found Global optimum!/))
		{
			$isSolvePerebor = 1;
		}
		if (($tmp =~ /^AGP found Global optimum!/))
		{
			$isSolveAGP = 1;
		}
		
		#iterPointsSavePath = outLogFiles\1_80_54_0.70.log
		#libConfigPath = outLogFiles\grishaginc_conf1_80_54_0.7.xml
		#if($tmp =~ /^libConfigPath = outLogFiles\grishaginc_conf\\(\d+)_(\d+)_(\d+)_(\d+[.]?\d*)/)
		#if($tmp =~ /^libConfigPath = outLogFiles\\grishaginc_conf(\d+)_(\d+)_(\d+)_(\d+)_(\d+[,]?\d*).xml/)
		#{
		#	$cf        = $1;
		#	$f0 = $2;
		#	$f1 = $3;
		#	$f2 = $4;
		#	$delta = $5;
		#	#print "$cf $f0 $f1 $delta\n";
		#}
	}
	
	$tmp = $Reord + $VReordB;
	#if ($dif < $eps)
	#{
	#	$isSolv = 1;
	#}

	
	
	print W "$TaskNumber; $n; $FuncNumber; $eps; $r; $min; $isSolvePerebor; $timePerebor; $IterPerebor; $isSolveAGP; $timeAGP; $IterAGP; \n";
	close(R);
}
close(W);