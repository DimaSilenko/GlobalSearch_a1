#use List::Util qw(min max);
use XML::LibXML;
print("Hellow world!!!\n");

##@dists=( "0.9", "0.9", "0.66", "0.9", "0.66", "0.9", "0.66", "0.66", "0.9", "0.9", "0.9" );
##@radiuss=( "0.2", "0.1", "0.2", "0.2", "0.2", "0.2", "0.3", "0.2", "0.2", "0.2", "0.2" );
##@dims=( "2", "2", "3", "3", "4", "4", "5", "5", "6", "7", "8" );
##@epss=( "0.01", "0.01", "0.01", "0.01", "0.03", "0.03", "0.03", "0.03", "0.0464", "0.0517", "0.0562" );
##@shs=( "simple", "hard", "simple", "hard", "simple", "hard", "simple", "hard", "simple", "simple", "simple" );
##@rs=( "4.5", "5.6", "4.5", "5.6", "4.5", "5.6", "4.5", "5.6", "5.6", "5.6", "5.6" );
##@nls=(1, 2);
##@kdims=(0, 0, 1);
##$l=9;
##$npl=256;
$zad=100;
##$rm=6;
$N=2;
##$l1=1;
##$l2=1;

##$k=0;
##$numlev=2;
##$r0=2.5;
##$r1=4.0;
##$nt1 = 15;
##$NumThread=1;

for ( $i=4 ; $i<=4; $i++ )
{
	##$dist=$dists[$i];
	##$radius=$radiuss[$i];
	##$N=$dims[$i];
	##$eps=$epss[$i];
	##$sh0=$shs[$i];
	##$r=$rs[$i];

	##$m=10;
	
	for ( $func=1	; $func<=$zad; $func++ )
	#for ( $func=7 ; $func<=8; $func++ )
	{
		##$nt1 = 15;
		##$npl = 8;
		##for ( $NumThread=8 ; $NumThread<=$npl; $NumThread=$NumThread * 2 ) 
		##{
		my $doc    = XML::LibXML::Document->new('1.0', 'utf-8');
		my $create = $doc->createElement('config');
		$doc->setDocumentElement($create);

		my $name_doc = $doc->createElement('function_number');
		$name_doc->appendText($func);
		$create->appendChild($name_doc);

		open my $out, '>', 'file_grishagin.xml';
		print {$out} $doc->toString(2);		

			#$nt1 = $nt1 + 1;
			#for ( $r=$r0 ; $r<=$r1; $r= $r+0.7  )
			{
			
				#$argum1="examin.exe -lib gklsCGpu.dll  -StepPrintMessages 100000 -NumPoints $NumThread -NumThread $NumThread -TypeMethod 0 -m 10 -TypeCalculation 0 -function_number $func -global_dist $dist -global_radius $radius -Dimension $N -Epsilon $eps -r $r -Comment AGP -MaxNumOfPoints 1000000_2_1000_1000 -stopCondition 1 >log2\\AGP_${i}_${func}_${NumThread}_${r}.txt";
			
			    #print($argum1);
				#system $argum1;
				#print("\n");
			
				$argum1="..\\..\\_bin\\GlobalSearch.exe grishagin file_grishagin.xml >log7_grishagin\\AS1_${i}_${func}_${NumThread}_${r}.txt";
				
				print($argum1);
				system $argum1;
				print("\n");


			}	
		##}
		
		##$npl = 8;
		##$r=3.3;
		##for ( $NumThread=8 ; $NumThread<=$npl; $NumThread=$NumThread * 2 ) 
		##{

			#$nt1 = $nt1 + 1;
			#for ( $r=$r0 ; $r<=$r1; $r= $r+0.2  )
			##{
			
				##$argum1="examin.exe -lib gklsCGpu.dll  -StepPrintMessages 100000 -NumPoints $NumThread -NumThread $NumThread -TypeMethod 9 -m 10 -TypeCalculation 1 -function_number $func -global_dist $dist -global_radius $radius -Dimension $N -Epsilon $eps -r $r -Comment ASgpu1 -MaxNumOfPoints 10000_2_1000_1000 -stopCondition 1 -nl 2 -dl 2_2 -cl ${NumThread}_0 >log7\\ASgpu1_${i}_${func}_${NumThread}_${r}.txt";
				##print($argum1);
				##system $argum1;
				##print("\n");
				
				
				##$argum1="\"C:\\Program Files (x86)\\MPICH\\mpd\\bin\\MPIRun.exe\" -np 2 examin.exe -lib gklsCGpu.dll  -StepPrintMessages 100000 -NumPoints $NumThread -NumThread $NumThread -TypeMethod 6 -m 10 -TypeCalculation 4 -function_number $func -global_dist $dist -global_radius $radius -Dimension $N -Epsilon $eps -r $r -Comment ASgpu1 -MaxNumOfPoints 10000_2_1000_1000 -stopCondition 1 -nl 2 -dl 2_2 -cl ${NumThread}_0 -mt mpRotated -ml 2_1 -mpl 2_1 >log7\\MMASgpu1_${i}_${func}_${NumThread}_${r}.txt";
				##print($argum1);
				##system $argum1;
				##print("\n");
				


			##}	
		##}
	}

}