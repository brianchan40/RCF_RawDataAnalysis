#! /opt/star/bin/perl -w

use File::Basename;
use Getopt::Std;
use Cwd 'abs_path';     # aka realpath()
my $resubmit = "resub.sh";

@jobDirs = (
"/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/condor_files"
);

@outputDirs = (
"/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/output/lam_low"
);

@flowCents = (
#"Data_160",
#"Data_161"
#,
#"Data_162"
#"Data_07",
#"Data_08",
#"Data_09_",
#"Data_10_",
#"Data_11_"
#"Data_12_",
#"Data_13_",
"Data_14_"
#"Data_15_",
#"Data_16_"
);

@Centrality = (
#"0",
#"1",
#"2",
#"3",
"4"
#"5",
#"6",
#"7",
#"8"
);


foreach $eachCent (@flowCents) {
    foreach $centt (@Centrality) {
    
        $totalJobs=0.;
        $finishedJobs=0.;
        $finishedJobs1=0.;
    
        foreach $jobDir (@jobDirs) {
        
            foreach $outputDir (@outputDirs) {
            
                if (-e "$outputDir/${eachCent}${centt}/$resubmit") { print `rm $outputDir/${eachCent}${centt}/$resubmit \n`;}
            
                foreach (glob ("$jobDir/${eachCent}${centt}/sched*.csh") ){
                
                    my $eachScript = $_;
                    chomp $eachScript;
                    @fields = split(/\//,$eachScript) ;
                    $eachScriptNoPath= $fields[$#fields];
                    ($eachJOBID = $eachScriptNoPath) =~ s/\.csh// ;
                    chomp $eachJOBID;
                
                    #if ( (!(-s "$outputDir/${eachCent}${centt}/$eachJOBID\_plam.root")) || (!(-s "$outputDir/${eachCent}${centt}/${eachJOBID}cen${centt}.weight\_112\_module\_new.root"))) {
                    #if ( (!(-s "$outputDir/${eachCent}${centt}/$eachJOBID\_plam.root")) || (!(-s "$outputDir/${eachCent}${centt}/${eachJOBID}cen${centt}.gamma112_fullEP\_eff\_pT02\_module.root"))) {
                    if ( (!(-s "$outputDir/${eachCent}${centt}/$eachJOBID\_plam.root")) || (!(-s "$outputDir/${eachCent}${centt}/$eachJOBID\_tree.root"))) {
                        #print "$outputDir/${eachCent}${centt}/${eachJOBID}cen0.weight\_112\_module\_new.root";
                    
                        if (-e "$jobDir/${eachCent}${centt}/$eachJOBID.run.log") { print `rm $jobDir/${eachCent}${centt}/$eachJOBID.run.log \n`; }
                    
                        $submitCommand = "cd $jobDir\/${eachCent}${centt}; condor_submit $jobDir\/${eachCent}${centt}\/$eachJOBID.condor";
                    
                        open(macroFile,">>$outputDir/${eachCent}${centt}/$resubmit");
                        print macroFile "$submitCommand \n";
                    }
                }
            
                my @tempScripts = glob( "$jobDir/${eachCent}${centt}/sched*.csh" );
                my @tempFinished = glob( "$outputDir/${eachCent}${centt}/*plam.root" ) ;
                my @tempFinished1 = glob( "$outputDir/${eachCent}${centt}/*tree.root" ) ;
            
                #print scalar(@tempFinished1);
            
                $totalJobs +=scalar(@tempScripts);
                $finishedJobs +=scalar(@tempFinished);
                $finishedJobs1 +=scalar(@tempFinished1);
            }
        
            my $perCentFinished = 0.;
            my $perCentFinished1 = 0.;
            if ($totalJobs != 0) {$perCentFinished = ($finishedJobs/$totalJobs)*100.;}
            if ($totalJobs != 0) {$perCentFinished1 = ($finishedJobs1/$totalJobs)*100.;}
            print "for $jobDir/${eachCent}${centt}, finished $finishedJobs jobs out of $totalJobs, ".$perCentFinished."% completed \n";
            print "for $jobDir/${eachCent}${centt}, finished $finishedJobs1 jobs out of $totalJobs, ".$perCentFinished1."% completed \n";
        }
    }
}

exit;

