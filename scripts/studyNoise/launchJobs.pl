#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

use Env;

#PG lettura dei parametri da cfg file
#PG --------------------------------
print "reading ".$ARGV[0]."\n" ;

open (USERCONFIG,$ARGV[0]) ;

while (<USERCONFIG>)
{
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
#    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;
}

$username = getpwuid( $< );

$BASEDir      = $User_Preferences{"BASEDir"};
$CMSSWDir     = $User_Preferences{"CMSSWDir"};
$EXEName      = $User_Preferences{"EXEName"};
$CFGTemplate  = $User_Preferences{"CFGTemplate"};
$BJOBTemplate = $User_Preferences{"BJOBTemplate"};
$BJOBList     = $User_Preferences{"BJOBList"};
$QUEUEId      = $User_Preferences{"QUEUEId"};

print "BASEDir      = ".$BASEDir."\n" ;
print "CMSSWDir     = ".$CMSSWDir."\n" ;
print "EXEName      = ".$EXEName."\n" ;
print "CFGTemplate  = ".$CFGTemplate."\n" ;
print "BJOBTemplate = ".$BJOBTemplate."\n" ;
print "BJOBList     = ".$BJOBList."\n" ;
print "QUEUEId      = ".$QUEUEId."\n" ;
print("\n");






$sampleJobListFile = "./lancia.sh";
open(SAMPLEJOBLISTFILE, ">", $sampleJobListFile);



open(BJOBLIST,$BJOBList) ;
while(<BJOBLIST>)
{
  @chars = split("",$_);
  $firstChar = $chars[0];
  if( $firstChar eq "#" )
  {
    next;
  }
  
  ($label,$inputFileList_DA,$inputFileList_MC,$PUDistrFile,$PUReweighting) = split(" ");
  print("*** ".$label." ***\n");
  
  
  
  $jobDir = $BASEDir."/output/".$label;
  $command = "mkdir ".$jobDir;
  system($command);

  $jobFile = $jobDir."/bjob_".$label.".sh";
  open JOBFILE, ">", $jobFile;
  
  $newInputFileList_DA = $BASEDir."/data/tmpList_".$label."_DA.txt";
  open NEWINPUTFILELISTDA, ">", $newInputFileList_DA;
  $newInputFileList_MC = $BASEDir."/data/tmpList_".$label."_MC.txt";
  open NEWINPUTFILELISTMC, ">", $newInputFileList_MC;
  
  
  
  print JOBFILE "df -h\n";
  print JOBFILE "pwd\n";
  
  $eos = "/afs/cern.ch/project/eos/installation/0.2.5/bin/eos.select";
  
  print JOBFILE "mkdir ./DA_".$label."\n";
  open(INPUTFILELISTDA,$BASEDir."/data/".$inputFileList_DA);
  while(<INPUTFILELISTDA>)
  {
    chomp($_);
    print JOBFILE $eos." cp -r /eos/cms".$_." ./DA_".$label."/\n";
  }
  print NEWINPUTFILELISTDA "./DA_".$label."/*.root\n";
  
  print JOBFILE "mkdir ./MC_".$label."\n";
  open(INPUTFILELISTMC,$BASEDir."/data/".$inputFileList_MC);
  while(<INPUTFILELISTMC>)
  {
    chomp($_);
    print JOBFILE $eos." cp -r /eos/cms".$_." ./MC_".$label."/\n";
  }
  print NEWINPUTFILELISTMC "./MC_".$label."/*.root\n";
  
  print JOBFILE "cd ".$CMSSWDir."/src\n";
  print JOBFILE "pwd\n";
  print JOBFILE "eval \`scramv1 runtime -sh\`\n";
  print JOBFILE "cd -\n";
  print JOBFILE "pwd\n";
  print JOBFILE "cd ".$BASEDir."\n";
  print JOBFILE "pwd\n";
  print JOBFILE "source scripts/setup.sh\n";
  print JOBFILE "cd -\n";
  print JOBFILE "pwd\n";
  print JOBFILE "ls -lh ./*\n";
  print JOBFILE "unbuffer ".$BASEDir."/bin/".$EXEName." ".$BASEDir."/output/".$label."/studyNoise_".$label.".cfg > ".$jobDir."/out.txt\n";
  print JOBFILE "rm ".$newInputFileList_DA."\n";
  print JOBFILE "rm ".$newInputFileList_MC."\n";
  print JOBFILE "rm -rf ./DA_".$label."\n";
  print JOBFILE "rm -rf ./MC_".$label."\n";
  
  $command = "chmod 777 ".$jobFile;
  system($command);
  
  
  
  $cfgFile = $jobDir."/studyNoise_".$label.".cfg";
  system ("cat ".$CFGTemplate." | sed -e s%BASEDIR%".$BASEDir.
                            "%g | sed -e s%LABEL%".$label.
                            "%g | sed -e s%OUTPUTDIR%".$jobDir.
                            "%g | sed -e s%PUREWEIGHTING%".$PUReweighting.
                            "%g | sed -e s%PUDISTRFILE%".$PUDistrFile.
                            "%g > ".$cfgFile);
  
  $command = "chmod 777 ".$cfgFile;
  system($command);
  
  
  
  print SAMPLEJOBLISTFILE "echo \"bsub -cwd ".$jobDir." -q ".$QUEUEId." ".$jobFile."\"\n";  
  print SAMPLEJOBLISTFILE "bsub -cwd ".$jobDir." -q ".$QUEUEId." ".$jobFile."\n";  
}  
