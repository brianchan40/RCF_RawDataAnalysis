#! /opt/star/bin/perl -w

use File::Basename;
use Getopt::Std;
use Cwd 'abs_path';     # aka realpath()
my $resubmit = "resub.sh";

@jobDirs = (
"/star/u/brian40/PicoDst/PicoDst_NoMaker/AMPT/condor_files/Gamma112_"
);

@outputDirs = (
"/star/u/brian40/PicoDst/PicoDst_NoMaker/AMPT/output/antilam/Gamma112_"
);

@flowCents = (
#"Data_160",
#"Data_161",
#"Data_162"
"0",
"1",
"2",
"3",
"4",
"5",
"6",
"7",
"8"
);


foreach $eachCent (@flowCents) {
    
    $totalJobs=0.;
    $finishedJobs=0.;
    
    foreach $jobDir (@jobDirs) {
        
        foreach $outputDir (@outputDirs) {
            
            if (-e "${outputDir}${eachCent}/$resubmit") { print `rm ${outputDir}${eachCent}/$resubmit \n`;}
            
            foreach (glob ("${jobDir}${eachCent}/sched*.csh") ){
                
                my $eachScript = $_;
                chomp $eachScript;
                @fields = split(/\//,$eachScript) ;
                $eachScriptNoPath= $fields[$#fields];
                ($eachJOBID = $eachScriptNoPath) =~ s/\.csh// ;
                chomp $eachJOBID;
                
                #if ( (!(-s "${outputDir}${eachCent}/${eachJOBID}cen${eachCent}.weight\_112\_module\_new.root")) ) {
                if ( (!(-s "${outputDir}${eachCent}/${eachJOBID}cen${eachCent}.gamma112\_fullEP\_eff\_pT02\_module.root")) ) {
                #if ( (!(-s "${outputDir}${eachCent}/${eachJOBID}.root")) ) {
                    
                    #print "${outputDir}${eachCent}/${eachJOBID}cen${eachCent}.weight\_112\_module\_new.root";
                    #print "${outputDir}${eachCent}/${eachJOBID}cen${eachCent}.gamma112\_fullEP\_eff\_pT02\_module.root";
		    #print "${outputDir}${eachCent}/${eachJOBID}.root";                    

                    if (-e "${jobDir}${eachCent}/$eachJOBID.run.log") { print `rm ${jobDir}${eachCent}\/$eachJOBID.run.log \n`; }
                    
                    $submitCommand = "cd ${jobDir}${eachCent}; condor_submit ${jobDir}${eachCent}\/$eachJOBID.condor";
                    #$submitCommand = "./run_gamma.csh 1 1 ./condor_files/$eachJOBID.list $eachJOBID";

                    
                    open(macroFile,">>${outputDir}${eachCent}/$resubmit");
                    print macroFile "$submitCommand \n";
                }
            }
            
            my @tempScripts = glob( "${jobDir}${eachCent}/sched*.csh" );
            #my @tempFinished = glob( "${outputDir}${eachCent}/*.root" ) ;
            my @tempFinished = glob( "${outputDir}${eachCent}/*gamma112\_fullEP\_eff\_pT02\_module.root" ) ;
            #my @tempFinished = glob( "${outputDir}${eachCent}/*.weight\_112\_module\_new.root" ) ;

            $totalJobs +=scalar(@tempScripts);
            $finishedJobs +=scalar(@tempFinished);
        }
        
        my $perCentFinished = 0.;
        if ($totalJobs != 0) {$perCentFinished = ($finishedJobs/$totalJobs)*100.;}
        print "for $jobDir/$eachCent, finished $finishedJobs jobs out of $totalJobs, ".$perCentFinished."% completed \n";
    }
}

exit;

