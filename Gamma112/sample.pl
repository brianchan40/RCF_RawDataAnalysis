@flowCents = (
#"Data_160",
#"Data_161",
#"Data_162"
"Data_08",
"Data_09",
"Data_10",
"Data_11",
"Data_12",
"Data_13",
"Data_14",
"Data_15",
"Data_16"
);

#foreach $eachCent (@flowCents) {
  foreach(glob ("/star/u/brian40/PicoDst/PicoDst_NoMaker/Gamma112/condor_files/sched*_*_*.condor"))
  {
    my $eachScript = $_;
    open (prototypeMacro, "$eachScript");
    chomp $eachScript;
    @fields = split(/\//,$eachScript) ;
    $eachScriptNoPath= $fields[$#fields];
    ($eachJOBID = $eachScriptNoPath) =~ s/\.condor// ;
    chomp $eachJOBID;

    @names = split(/_/,$eachJOBID) ;
    $eachNumber= $names[$#names];
  	print "$eachNumber\n";
  	$eachName = $names[0];
  	print "$eachName\n";

    $count =0.;
    while ($eachLine = <prototypeMacro>) {

      if(($count % 17)==0) {open (macroFile,">Sample/$eachName\_$eachNumber.condor");}

      print macroFile $eachLine;
      $count = $count+1;

      if(($count % 17)==0) {close macroFile; $eachNumber = $eachNumber -1;}
    }
    
     close prototypeMacro;
  }
#}

