###############################################################################
# make sure the position is exactly matching, need to subtract 1 for starting at 0
# $Replacement[$j]=$Replacement[$j]-1; 
sub CreateRSD_ScaledSeq
{
	my ($scale,$geneID,$isoform,$length1,$ref_Synonymous,$ref_Replacement,$ref_damaging)=@_;

	my @new2;

	my $length=int($length1/$scale);
	for (my $i=0;$i<$length;$i++)
		{
		 push @new2,"*";
		}

	my @Replacement=@$ref_Replacement;
	my @Synonymous=@$ref_Synonymous;
	my @Damaging=@$ref_damaging;

	my $OrignalR= scalar(@Replacement);
	my $OrignalS= scalar(@Synonymous);
	#Silent first
	for (my $j=0;$j<scalar(@Synonymous);$j++)
		{
		 $Synonymous[$j]=$Synonymous[$j]-1;#make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Synonymous[$j]/$scale);
		 $new2[$t]="S";
		}

#Damaging second
	for (my $j=0;$j<scalar(@Damaging);$j++)
		{
		 $Damaging[$j]=$Damaging[$j]-1;#make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Damaging[$j]/$scale);
		 $new2[$t]="D";
		}
#Finally Replacement
	for (my $j=0;$j<scalar(@Replacement);$j++)
		{
		 $Replacement[$j]=$Replacement[$j]-1; #make sure the position is exactly matching, need to subtract 1 for starting at 0
		 my $t= int($Replacement[$j]/$scale);
		 $new2[$t]="R";
		}



	my $countR=0;
	my $countS=0;
	my $countD=0;


	for (my $ii=0;$ii<$length;$ii++)
		{
		 if ($new2[$ii] eq "R")
			 {
			  $countR=$countR+1;
			 }

		 if ($new2[$ii] eq "S")
			 {
			  $countS=$countS+1;
			 }

		 if ($new2[$ii] eq "D")
			 {
			  $countD=$countD+1;
			 }			 
		}

	printf "The total number of Replacement, Synonymous, Damaging after scaling is %d and %d.\n\n",$countR, $countS, $countD;
	printf ">%s_%s.%d_R_%d_S_%d_D_%d_len_%d_R_%d_S_%d\n%s\n",$geneID,$isoform,$length,$countR,$countS,$countD,$length1,$OrignalR,$OrignalS,join("",@new2);
	return @new2;
 
}
