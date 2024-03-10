######################################################################################################################################
##
##		TITLE:				n1stats.pl
##
##		PURPOSE:			This program is designed to read information from the outputs of 
##							n1parse.pl, match the scheduled images to the received images,
##							and compare that information, exporting to files for individual days
##							and one large file for cumulative stats, as well as a list of images
##							where at least one category had an issue.
##
##
##
##							Sub-procedures are as follows:
##						
##							<Stuff for that goes here>
##
##
##
##		AUTHOR: 			Ryan Krieger, Co-op Student, 23 Jan, 2019
##
##		
		my $VER =            '4.00'; # Initial Release
##
##
##
#######################################################################################################################################

#!usr/bin/perl
use strict;
use feature 'state';
use List::MoreUtils 'pairwise';
use Config::Tiny;
use File::Copy;
use lib 'C:\Perl\lib\Date\Calc';
use Date::Calc qw(:all);

# Store the home volume & directory
my ($homevolume, $homedir, $homefilename) = File::Spec->splitpath(File::Spec->rel2abs(__FILE__));
push @INC, File::Spec->catpath($homevolume, $homedir, "");

require 'n1parse.pl';

my $config = Config::Tiny->read('stats.config');

my @start_date; # $start_date[0] = year, $start_date[1] = day of year
my @end_date;	# $start_date[0] = year, $start_date[1] = day of year
my @ssc_files;
my @fits_files;
my @totals = ();
my %file_match = ();
my %fits_hash = (); #temporary hash containing all the fits information for a single day
my %ssc_hash = ();	#temporary hash containing all the ssc information for a single day
my $ssc_key;
my $fits_key;
my $year;
my $doy; #read as "day of year"
my $num_args_in = 0;
my $monthly = 'No';
my @month_requested = ();
my $use_daily = 0;
my $use_obsv = 0;
my $cat_req = '';
my $filename_base = '';
my $range_type = '';


## --- DIRECTORIES --- ##
my $read_dir = $config->{directories}->{mid_dir}; #mid_dir stands for middle directory, the files are stored there from parser program, then read into here
my $out_dir = $config->{directories}->{stats_out}; #directory that the files are output to 
my $debug_dir = $config->{directories}->{deb_dir};
my $compl_dir = $config->{directories}->{compl_dir}; #directory that files with no missing images will be copied to

my $ssc_dir = $config->{directories}->{ssc_in}; 
my $fits_dir = $config->{directories}->{fits_in}; 
my $prev_ssc_dir = $config->{directories}->{prev_ssc};
my $prev_fits_dir = $config->{directories}->{prev_fits};

unless (-d $out_dir) {
   File::Path::make_path($out_dir) or die "Fatal error: Cannot create $out_dir\n";
}
unless (-d $compl_dir) {
   File::Path::make_path($compl_dir) or die "Fatal error: Cannot create $compl_dir\n";
}
unless (-d $debug_dir) {
   File::Path::make_path($debug_dir) or die "Fatal error: Cannot create $debug_dir\n";
}
unless (-d $read_dir) {
   File::Path::make_path($read_dir) or die "Fatal error: Cannot create $read_dir\n";
}

## --- COMPLETION INFO FILE --- ##
my $compl_info_out = File::Spec->catfile($compl_dir, "compl_stats.info");
if (! -e $compl_info_out) {
   open COMPL_INFO, ">$compl_info_out" or die "Cannot write to $compl_info_out: $!\n";
   close COMPL_INFO;
}
my $compl_info = Config::Tiny->read($compl_info_out);


## --- RANGES --- ##
my $fp_dec_range = $config->{acceptable_ranges}->{fp_dec_range};
my $fp_ra_range = $config->{acceptable_ranges}->{fp_ra_range};
my $fp_roll_range = $config->{acceptable_ranges}->{fp_roll_range};
my $fp_rate_range = $config->{acceptable_ranges}->{fp_rate_range};

my $fs_dec_range = $config->{acceptable_ranges}->{fs_dec_range};
my $fs_ra_range = $config->{acceptable_ranges}->{fs_ra_range};
my $fs_roll_range = $config->{acceptable_ranges}->{fs_roll_range};
my $fs_rate_range = $config->{acceptable_ranges}->{fs_rate_range};

my $def_ra_range = $config->{acceptable_ranges}->{default_ra_range};
my $def_dec_range = $config->{acceptable_ranges}->{default_dec_range};
my $def_roll_range = $config->{acceptable_ranges}->{default_roll_range};
my $def_rate_range = 10 * $fs_rate_range;  # Default rate tolerance not defined in config file; setting to 10x larger than Fine Slew rate tolerance

my $sec_range = $config->{acceptable_ranges}->{start_range};
my $dur_range = $config->{acceptable_ranges}->{dur_range};


my %cat_idcs = ( 				#used for if a specific category_bad or if a specific Observer is requested, change if the output format is changed
	'ASTRO' => 5, 'HEOSS' => 5,												#Check: Whichever observer is specified
	'MISSING' => 6,															#Check: NO										(NO means image is not available)
	'COMPLETED' => 7,														#Check: BAD
	'START' => 9,															#Check: BAD
	'DURATION' => 13,														#Check: BAD
	'MODE' => 17,															#Check: BAD
	'FINE_POINT' => 18, 'FINE_SLEW' => 18, 'DESAT' => 18, 					#Check: Whichever mode is specified
	'RATE' => 20,															#Check: BAD
	'DEC' => 24,															#Check: BAD
	'RA' => 28,																#Check: BAD
	'ROLL' => 32, 															#Check: BAD
	'SIZE' => 38, 															#Check: BAD
);

my @cat_modifiers = ('','','');

foreach (@ARGV){	

	if ($_ =~ /-parse/i) {
                print "Calling n1parse.pl\n";	
		system('n1parse.pl');
	        &stats_log_parse($ssc_dir, $fits_dir, $read_dir, $prev_ssc_dir, $prev_fits_dir, $debug_dir);	
	
	} 

}

my $stats_debug = File::Spec->catfile($debug_dir, "stats_debug.txt");
open DEBUG, ">$stats_debug" or die "Cannot write to $stats_debug: $!";

opendir(FILES, "$read_dir");
my @filenames = sort (grep(/\.txt$/, readdir FILES)); #grabs all the filenames from the working directory and sorts
closedir FILES;

print "Files are: @filenames\n";

$filenames[0] =~ /^ (\d{4}) _ (\d{3}) .* /isx;
print "\nDefault start date: $1/$2\n";

$filenames[$#filenames] =~ /^ (\d{4}) _ (\d{3}) .* /isx;
print "Default end date: $1/$2\n\n";

{ 
	
	$filenames[0] =~ /^ (\d{4}) _ (\d{3}) .* /isx;
	@start_date = ($1, $2);

	$filenames[$#filenames] =~ /^ (\d{4}) _ (\d{3}) .* /isx;
	@end_date = ($1, $2);
	
	if (defined $ARGV[0]) { 
	
		my $month_requested = '';
		my @range = ();
	
		$num_args_in = scalar @ARGV;
	
		foreach (@ARGV) {
		
			if ($_ =~ /^ -range /ix) {
			
				if (($_ !~ /^ -range= \d{4} - \d{1,2} - \d{1,2} , \d{4} - \d{1,2} - \d{1,2} $/xi) && ($_ !~ /^ -range= \d{4} - \d{1,3} , \d{4} - \d{1,3} $/ix)) {
				
					print "USAGE: -range=[YYYY-MM-DD,YYYY-MM-DD or YYYY-DOY,YYYY-DOY] to specify a range\n";
					exit;
				
				}
				my @split = split /=/, $_;
				@range = split /,/, $split[1];

			
			} elsif ($_ =~ /-monthly/i) {
			
				$monthly = 'Yes';
				unless ($_ =~ /-monthly$/i) {
				
					if ($_ !~ / -monthly: \d{4} - \d{2} $/x) {
					
						print "USAGE: -monthly<:YYYY-MM> to specify a desired month, otherwise leave empty for previous month\n";
						exit;
					
					}
				
					my @split = split /:/, $_;
					$month_requested = $split[1];
					print $month_requested;
				
				} else {
				
					my $sys_time = localtime();
					print "$sys_time\n";
					
					my @sys_time = split ' ', $sys_time;
					my $month = Decode_Month("$sys_time[1]",1);
					
					$month = $month - 1;
					($sys_time[4] -= 1) if ($month == 0);
					($month = 12) if ($month == 0);
					$month_requested = $sys_time[4] . "-$month";
				
				}
			
			} elsif ($_ =~ /-daily/i) {
			
				unless ($_ =~ /-daily.+/i) {
				
					$use_daily = 1;
					last;
				
				} elsif ($_ =~ /-daily:.+/i) {
				
					my @split = split /:/, $_;
					
					if ($split[1] =~ /^ (\d{4}) - (\d{1,3}) $/x) {
					
						@start_date = ($1, $2);
						@end_date = ($1, $2);
					
					} elsif ($split[1] =~ /^ (\d{4}) - (\d{1,2}) - (\d{1,2}) $/x) {
					
						my $days = Day_of_Year($1, $2, $3);
						@start_date = ($1, $days);
						@end_date = ($1, $days);
					
					} else {
					
						print "USAGE: -daily{:YYYY-MM-DD or :YYYY-DDD} to specify a date, otherwise just -daily for all days\n";
						exit;
					
					}
					
					$use_daily = 1;
				
				}
			
			} elsif ($_ =~ /-cat_req:/x) {
			
				my @split = split /:/, $_;
				$cat_req = "$split[1]";
				if ($cat_req =~ /^MODE/) {
				
					unless ($cat_req =~ /^MODE$/) {
					
						@cat_modifiers = split /,/, $cat_req;
						$cat_req = $cat_modifiers[0];
						print "$cat_req => $cat_modifiers[1] , $cat_modifiers[2]\n";
					
					}
				
				}
				
				if (not exists $cat_idcs{$cat_req}) {
				
					print "USAGE: -cat_req:[CATEGORY] where category is:\n\tASTRO\n\tHEOSS\n\tMISSING\n\tCOMPLETED\n\tSTART\n\tDURATION\n\tMODE\n\tRATE\n\tDEC\n\tRA\n\tROLL\n\tSIZE\n\n";
					exit;
				
				} elsif (($cat_req eq 'ASTRO')||($cat_req eq 'HEOSS')) {
				
					$use_obsv = 1;
				
				}
			
			} elsif ($_ =~ /-help$/i) {
				
				print "INPUT HELP:\n\n";
				print "<> Indicates optional parameters\n{} Indicates optional parameters with optional formats\n[] Indicates required parameters but with optional formats\n\n";
				print "Commands:\n\n";
				print "-parse\n\tProgram will parse any log files within the specified directory\n\tLog files then moved to appropriate \"previous\" folder\n\n";
				print "-daily:{YYYY-MM-DD or YYYY-DDD}\n\t{} to specify a specific day\n\t-daily with nothing after for all days available\n\tCalling just fits_stats.pl has same effect as -daily with no parameters\n\n";
				print "-range=[YYYY-MM-DD,YYYY-MM-DD or YYYY-DDD,YYYY-DDD]\n\t[] to specify a desired range\n\tOutputs bad images to a range file\n\n";
				print "-monthly:<YYYY-MM>\n\t<> to specify a desired month\n\t-monthly with nothing after for prev month (based on system time)\n\tOutputs bad image stats to a monthly_bad file\n\tOutputs all image stats to monthly_all\n\n";
				print "-cat_req:[CATEGORY]\n\tCreates an extra third file specific to CATEGORY\n\tFor ASTRO/HEOSS, outputs all images for that observer\n\tFor all else, exports only when CATEGORY result is BAD\n\tMODE has modifiers:\n\t\tMODE_TYPE: FINE_SLEW, FINE_POINT, or DESAT\n\t\tWHICH: BOTH, NOT, or BOTH-NOT\n\t\tType it out as -cat_req:MODE,MODE_TYPE,WHICH\n\tCATEGORY is replaced with:\n\t\tASTRO\n\t\tHEOSS\n\t\tMISSING\n\t\tCOMPLETED\n\t\tSTART\n\t\tDURATION\n\t\tMODE\n\t\tRATE\n\t\tDEC\n\t\tRA\n\t\tROLL\n\t\tSIZE\n\n";
				
				exit;
			
			}
		
		}
		
		if (scalar @range != 0) {
		
			if (scalar @range != 2) {
			
				print "Too many input arguments for range, please try again\n";
				exit;
			
			}
			
			@start_date = split /-/, $range[0];
			@end_date = split /-/, $range [1];
			
			if (scalar @start_date == 3) {
			
				my $day_of_year = Day_of_Year($start_date[0], $start_date[1], $start_date[2]);
				pop @start_date; pop @start_date;
				push @start_date, $day_of_year;
				print "Start_Date = $start_date[0]/$start_date[1]\n";
			
			}
			
			if (scalar @end_date == 3) {
			
				my $day_of_year = Day_of_Year($end_date[0], $end_date[1], $end_date[2]);
				pop @end_date; pop @end_date;
				push @end_date, $day_of_year;
				print "End_Date = $end_date[0]/$end_date[1]\n";
			
			}
		
		} elsif (($month_requested ne '') && (scalar @range == 0)) {
		
			@month_requested = split /-/, $month_requested;
			$month_requested[1] =~ s/^0//;
			my $start_doy = Day_of_Year($month_requested[0], $month_requested[1], 1);
			my $days_in_month = Days_in_Month($month_requested[0], $month_requested[1]);
			my $end_doy = Day_of_Year($month_requested[0], $month_requested[1], $days_in_month);
			
			@start_date = ($month_requested[0], $start_doy);
			@end_date = ($month_requested[0], $end_doy);
			
			($start_date[1] = "00" . $start_date[1]) if ($start_date[1] < 10);
			($start_date[1] = "0" . $start_date[1]) if ($start_date[1] >= 10 and $start_date[1] < 100);
			($end_date[1] = "00" . $end_date[1]) if ($end_date[1] < 10);
			($end_date[1] = "0" . $end_date[1]) if ($end_date[1] >= 10 and $end_date[1] < 100);
			# print "Start_Date = $start_date[0]/$start_date[1]\n";
			# print "End_Date = $end_date[0]/$end_date[1]\n";
		
		}
		

	}

	print "Start_Date = $start_date[0]/$start_date[1]\nEnd_Date = $end_date[0]/$end_date[1]\n";

	#this foreach loop checks the dates on each of the files in the specified directory, then
	#if they are within the desired interval, pushes the name to either the SSC file name array, or the FITS file name array
	
	foreach (@filenames) {
		
		#these two states are used to compare the dates of the files, for instance 2016 001 (jan 1st 2016)
		#would be 736345 (2016*365.25 + 1) and 2016 010 (jan 10th 2016) would be 736354 (2016*365.25 + 10), 
		#an easy numeric start and end number to compare rather than trying to compare two sets of two different numbers
		
		state $abs_start = int($start_date[0]*365.25 + $start_date[1]);
		state $abs_end = int($end_date[0]*365.25 + $end_date[1]);
		
		
		if ($_ =~ /^ (\d{4})_(\d{1,3})_(\w{3,4})/isx) {
		
			#converts the date at the front of the file into the format used to define the absolute start and end date
			#then compares the new date to the start/end date, and sends to either a fits or ssc specific array
			
			my $abs_date = int($1*365.25 + $2);
			(push @ssc_files, $_) if (($abs_start <= $abs_date && $abs_date <= $abs_end) && $3 eq 'SSC');
			(push @fits_files, $_) if (($abs_start <= $abs_date && $abs_date <= $abs_end) && $3 eq 'FITS');
		
		}

	}

	#this next block of code checks if there are no files in the @ssc_files or @fits_files arrays
	#then gives the user the option to change the start/end dates if they would like, else exits the code
	
	if (!@ssc_files || !@fits_files){

		#backup in case there are no files at all in the range, as that would mean nothing is really happening
		undef @start_date;
		undef @end_date;
		print "ERROR: There seem to be no files of one type corresponding to that time interval\n";
		#chomp(my $choice = <STDIN>);
		#(redo) if ($choice =~ /y/i);
		#(exit) if ($choice =~ /n/i);
		exit

	}

}

#Matches FITS files to their associated SSC files
foreach my $ssc (@ssc_files) {

	#Finds the year and day from the SSC name
	$ssc =~ /^ (\d{4})_(\d{1,3})_\w{3,4} /isx;
	($year, $doy) = ($1, $2);
	
	
	foreach my $fits (@fits_files) {
	
		#matches SSC file names to FITS file names and adds the pair to the hash %file_match
		$fits =~ /^ (\d{4})_(\d{1,3})_\w{3,4} /isx;
		
		next unless ($year == $1 && $doy == $2);
		$file_match{$ssc} = $fits;
	
	}

}

#This block of code creates two different hashes, one for SSC and one for FITS, for each specific day, and then sends
#that to the subroutine &stats_creating to create/print the stats to files corresponding to their specific day
#certain totals are returned, added to previous totals (initially 0) and used to create final statistics over the whole time interval

foreach my $keys (sort keys %file_match) {

	#Opens the file corresponding to an individual date
	open SSCFILE,  File::Spec->catfile($read_dir,$keys) or die "could not open SSC: $!";
	open FITSFILE, File::Spec->catfile($read_dir,$file_match{$keys}) or die "could not open FITS: $!";
	
	$keys =~ /^ (\d{4})_(\d{1,3})_\w{3,4} /isx;
	($year, $doy) = ($1, $2);

	
	#the next two while loops parse individual day files, and export the data to either %ssc_hash or %fits_hash
	
	while (<SSCFILE>) {
		
		chomp;
		
		if ($_ =~ /^\d{4}-\d{3}-\d{2}:\d{2}:\d{2}.*/s) {
		
			#removes the _SSC from the end of the times
			$_ =~ s/_SSC$//;
			($ssc_key) = $_;

		} elsif ($_ =~ /^;/) {
			
			#removes the semicolon from the beginning of the category names
			$_ =~ s/^;//;
			chomp(my ($sub_key, $value) = split /=>/, $_);
	
			$ssc_hash{$ssc_key}{$sub_key} = $value;
			
		} elsif ($_ =~ /\s+/isg){
			
			#skips the whitespace lines
			next;
		
		}
		
	}

	while (<FITSFILE>) {
		
		chomp;
		
		if ($_ =~ /^\d{4}-\d{3}-\d{2}:\d{2}:\d{2}.*/s) {
		
			#removes the _FITS from the end of the times
			$_ =~ s/_FITS$//;
			($fits_key) = $_;

		} elsif ($_ =~ /^;/) {

			#removes the semicolon from the beginning of the category names
			$_ =~ s/^;//;
			chomp(my ($sub_key, $value) = split /=>/, $_);
	
			$fits_hash{$fits_key}{$sub_key} = $value;
			
		} elsif ($_ =~ /\s+/isg){
			
			#skips the whitespace lines
			next;
		
		}
		
	}
	
	close SSCFILE;
	close FITSFILE;
	
	
	#sends the two hashes as well as what specific day/year the stuff is from to stats_creating subroutine
	#which will return values of how much a specific thing failed, and total (total number of images, number of failures)
	
	my @returned_vals = stats_creating(\%ssc_hash, \%fits_hash, $year, $doy);
	
	#adds the values from @returned_vals to their corresponding element in @totals, which will be used for final statistics
	@totals = pairwise { $a + $b } @totals, @returned_vals;
	
	#undefines both hashes so there is no data leftover from previous files/days when loop begins again
	undef %ssc_hash;
	undef %fits_hash;
	
}


if ($monthly eq 'Yes') {

	my @months = ('element zero', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec');
	my $used_month = $month_requested[1];
	$used_month = $months[$used_month];
	print "$used_month\n";
	$range_type = 'MONTHLY';
	
	$filename_base = "NEOS-$year-$used_month\_";
	
	range_stats(\@totals, $filename_base);

}

#sends @totals to the final_stats subroutine for use in creating a file containing the end statistics for the entire interval
if (($num_args_in != 0) && ($monthly eq 'No') && ($use_daily != 1)) {
	
	$filename_base = "NEOS_$start_date[0]-$start_date[1]_$end_date[0]-$end_date[1]_";
	$range_type = 'RANGE';
	
	range_stats(\@totals, $filename_base);

}



## STATS_CREATING creates the statistics per-date as well as adds all the cases 
## where something went wrong to an end file for use in subroutine FINAL_STATS
sub stats_creating { 

	## --- VARIABLE DECLARATION --- ##
	my $ssc_ref = shift;
	my $fits_ref = shift;
	my $year = shift;
	my $doy = shift;
	my $STATSFILE;
	my $total = 0;
	my $total_mode_bad = 0;
	my $total_dec_bad = 0;
	my $total_ra_bad = 0;
	my $total_roll_bad = 0;
	my $total_dur_bad = 0;
	my $total_fine_point = 0;
	my $total_missing = 0;
	my $total_incomplete = 0;
	my $incorrect_time = 0;
	my $total_size_bad = 0;
	my $total_rate_bad = 0;
	my $total_astro = 0;
	my $total_heoss = 0;
	my @return_totals = ();
	my @fine_slew_stats = ();
	my @fine_point_stats = ();
	my $missing_check = 0;

	## --- CODE FOLLOWS --- ##
	
	my $stats_out_name = "NEOS_$year-$doy\_daily.stats";
	my $stats_out = File::Spec->catfile($out_dir, $stats_out_name);
	if ((exists $compl_info->{completed}->{"$year-$doy"}) && (-e $stats_out)) {
	
		my $values = $compl_info->{completed}->{"$year-$doy"};
		@return_totals = split /-/, $values;
		return @return_totals;
	
	}
	
	#opens a file named after the date of the file the information is coming from

	open ($STATSFILE, ">$stats_out") or die "cannot open $stats_out, $!";
	my $header = '#' x 150; #a header, just a bunch of # in a row
	
	#prints the header, can add information to it later, currently contains:
	#Filename, Input_files, and Acceptable Ranges (for Right Ascension, Declination, Roll, Duration, and Start_time)
	print_sub('HEADER', $STATSFILE, $stats_out_name, $year, $doy);
	print_sub('COLUMN_NAMES', $STATSFILE);
	
	#this loop does everything except print to the file, will go into more detain inside
	foreach my $ssc_times (sort keys %{$ssc_ref}) {
	
		my $next = undef;
		my $test_times = undef;
		my @test_times = ();
		my $fits_times = $ssc_times;
		my $diff_sec = 0;
		$missing_check = 0;
		my $ssc_dec;
		my $fits_dec;
		my $diff_dec = 0;
		my $ssc_ra;
		my $fits_ra;
		my $diff_ra = 0;
		my $ssc_roll;
		my $fits_roll;
		my $diff_roll = 0;
		my $fits_rate;
		my $diff_rate = 0;
		my ($year, $doy, $time) = split /-/, $ssc_times;
		my $expected;
		my $actual;
		my $fits_latitude = '~';
		my $fits_longitude = '~';
	
		printf $STATSFILE "%-6d %5d %10s ", $year, $doy, $time;
		
		$total += 1;
		my $ssc_source = $ssc_ref->{$ssc_times}{'SOURCE_FILE'};
		my $observer = $ssc_ref->{$ssc_times}{'OBSERVER'};
		my $object = $ssc_ref->{$ssc_times}{'OBJECT'};
		($total_astro += 1) if ($observer eq 'ASTRO');
		($total_heoss += 1) if ($observer eq 'HEOSS');
		
		$ssc_source =~ s/\.SSCLOG$//i;
		printf $STATSFILE "%26s %26s %9s ", $ssc_source, $object, $observer;
		
		#checks if the fits hash has a key corresponding to the same time as one from an SSC log
		if (not exists $fits_ref->{$ssc_times}) {
			
			#if it does not exist, it checks the time that the following scheduled image supposed to take place in the SSC_LOG
			#and sets it to $next
			$next = $ssc_ref->{$ssc_times}{'NEXT_TIME'};
			$test_times = $ssc_times; #the times to be tested will be set to this variable, so that $ssc_times doesn't get overwritten
			
			unless ($next eq 'final') { #next will be set to ~ at the last thing in a file
			
				#converts the times from one long string of numbers, to numbers seperated into YYYY-DDD-HH-MM-SS
				#this is then split by - and assigned to @test_times so that it's easier to change values
				$diff_sec = 0; #sets $diff_sec to 0 and iterates it up 1 every time the next while loop loops
				
				#sends to subroutine &time_iterator
				($diff_sec, $incorrect_time, $fits_times) = time_iterator($next, $test_times, $fits_ref, $incorrect_time, $sec_range, $ssc_ref);
				
				#This works because we've already checked that $ssc_times does not work as a key in the %fits_hash
				#Because $fits_times was not changed and is thus equal to $ssc_times, the image file is assumed missing
				if ($fits_times eq $ssc_times) {#With all the - and : in the time, it is moreso a string than a number, string comparater eq needed
				
					print_sub('MISSING', $STATSFILE, $ssc_ref, $ssc_times);
					
					$total_missing += 1;
					$missing_check = 1;
					next;
				
				} else {
				
					printf $STATSFILE "%10s ", 'YES';
				
				}
				
			
			} else { #runs if the value for $next is '~', as it means there is no next value to use as an end
			
				($diff_sec, $incorrect_time, $fits_times) = time_iterator($next, $test_times, $fits_ref, $incorrect_time, $sec_range);
				
				if ($fits_times eq $ssc_times) {#With all the - and : in the time, it is moreso a string than a number, string comparater eq needed
				
					print_sub('MISSING', $STATSFILE, $ssc_ref, $ssc_times);
					$total_missing += 1;
					$missing_check = 1;
					next;
				
				} else {
				
					printf $STATSFILE "%10s ", 'YES';
				
				}
			
			}
		
		} else {
		
			printf $STATSFILE "%10s ", 'YES';
		
		}

		#checks the value for completion for 100, prints GOOD if it matches or BAD if it doesn't and increases $total_incomplete by 1 if BAD
		my $compl_percent = $fits_ref->{$fits_times}{'%Compl'};

		(printf $STATSFILE "\t%8s %9s ", 'GOOD', "$compl_percent") if ($compl_percent =~ /100/);
		(printf $STATSFILE "\t%8s %9s ", 'BAD', "$compl_percent") if ($compl_percent !~ /100/);
		($total_incomplete += 1) if ($compl_percent !~ /100/);		
		
		#Checks whether the start time was on time, or late	by greater than the acceptable range
		$expected = join('',$year,$doy,$time);
		$expected =~ s/[:]//g;
		@test_times = split /[-:]/g, $fits_times;
		$actual = join('', @test_times);
		
		my $exp_start = $fits_ref->{$fits_times}{'Exposure_start'};
		$exp_start =~ / \d{6} (\d{3}) $/x;
		$diff_sec = "$diff_sec"."\.$1";
		
		(printf $STATSFILE "\t%12s %15s %15s %7.3f ", 'GOOD', $expected, $actual, $diff_sec, "+-$sec_range\s") if ($diff_sec <= $sec_range);
		(printf $STATSFILE "\t%12s %15s %15s %7.3f ", 'BAD', $expected, $actual, $diff_sec, "+-$sec_range\s") if ($diff_sec > $sec_range);
		
		#checks whether the duration is the same, and by how much it differs (0 if the same)
		my $ssc_dur = $ssc_ref->{$ssc_times}{'DUR'};
		my $fits_dur = $fits_ref->{$fits_times}{'DUR'};
		my $diff_dur = abs($ssc_dur - $fits_dur);
		if ($diff_dur > $dur_range) {
		
			$total_dur_bad += 1;
			
			printf $STATSFILE "\t%7s %9s %9s %7.2f ", 'BAD', $ssc_dur, $fits_dur, $diff_dur;
		
		} else {
		
			printf $STATSFILE "\t%7s %9s %9s %7.2f ", 'GOOD', $ssc_dur, $fits_dur, $diff_dur;
			
		}
		
		
		#checks the value for completion for 100, prints GOOD if it matches or BAD if it doesn't and increases $total_incomplete by 1 if BAD
		# my $compl_percent = $fits_ref->{$fits_times}{'%Compl'};
		
		#This next stuff used to check whether the mode is correct, and if the mode is FINE_POINT for $ssc_mode
		my $ssc_mode = $ssc_ref->{$ssc_times}{'MODE'};
		my $fits_mode = $fits_ref->{$fits_times}{'MODE'};
		
		unless ( ($ssc_mode eq '~') 
			 || ( ( index($ssc_mode, 'FINE_SETTLE') >= 0 )  && ( index($fits_mode, 'FINE_POINT') >= 0 )  ) ) {
		
			if (index($ssc_mode, $fits_mode) < 0) {  # Cannot find $fits_mode string in $ssc_mode string

				printf $STATSFILE "\t%6s %15s %15s ", 'BAD', "$ssc_mode", "$fits_mode";
				$total_mode_bad += 1;
				
			} else { # $fits_mode matches $ssc_mode 

				printf $STATSFILE "\t%6s %15s %15s ", 'GOOD', "$ssc_mode", "$fits_mode";
				
			}		
		
		} else {
		
			printf $STATSFILE "\t%6s %15s %15s ", 'GOOD', "$ssc_mode", "$fits_mode";
		
		}
		
		my $rate_bad = 0;
		my $ssc_rate = $ssc_ref->{$ssc_times}{'RATE'};
		my $fits_ra_vel = $fits_ref->{$fits_times}{'RA_VEL'};
		my $fits_dec_vel = $fits_ref->{$fits_times}{'DEC_VEL'};
		
		($fits_rate, $diff_rate, $rate_bad) = rate_conv($ssc_rate, $fits_ra_vel, $fits_dec_vel, $ssc_mode, $STATSFILE); #includes mode to determine range to use
		$total_rate_bad = $total_rate_bad + $rate_bad;
		
		if (index($ssc_mode, $fits_mode) >= 0) {  # use "index" rather than regular expression, because we might see FINE_SLEW-s (for steady-state FINE_SLEW)
		
			if (index($ssc_mode, 'FINE_POINT') >= 0) {
			
				$fine_point_stats[0] += 1;
				$fine_point_stats[7] += $diff_rate;
				($fine_point_stats[8] = abs($diff_rate)) if ($fine_point_stats[8] < abs($diff_rate));
			
			} elsif (index($ssc_mode, 'FINE_SLEW') >= 0) {
			
				$fine_slew_stats[0] += 1;
				$fine_slew_stats[7] += $diff_rate;
				($fine_slew_stats[8] = abs($diff_rate)) if ($fine_slew_stats[8] < abs($diff_rate));
			
			}
		
		}
		
		
		my $ssc_ra = $ssc_ref->{$ssc_times}{'RA'};
		my $ssc_dec = $ssc_ref->{$ssc_times}{'DEC'};
		my $ssc_roll = $ssc_ref->{$ssc_times}{'ROLL'};
		
		($total_fine_point += 1) 	unless ( (index($ssc_mode, 'DESAT') >= 0) 
									or ($ssc_ra =~ /^ .+ ~ .+ $/x) 
									or (($ssc_dec =~ /^ .+ ~ .+ $/x) 
									or ($ssc_roll =~ /^ .+ ~ .+ $/x) ) );
		unless (index($ssc_mode, 'DESAT') >= 0) { 
		
			my @add = ();
			#declaring variables and making calculations for declination, right ascension, and roll
			
			@add = ra_dec_roll($ssc_ref, $fits_ref, $ssc_times, $fits_times, $STATSFILE);
			$total_dec_bad = $total_dec_bad + $add[0];
			$total_ra_bad = $total_ra_bad + $add[2];
			$total_roll_bad = $total_roll_bad + $add[4];
			
			if ((index($ssc_mode, 'FINE_POINT') >= 0) && (index($fits_mode, 'FINE_POINT') >= 0 ) ) {
				
				$fine_point_stats[1] += $add[1]; #DEC
				($fine_point_stats[2] = abs($add[1])) if ($fine_point_stats[2] < abs($add[1]));
				
				$fine_point_stats[3] += $add[3]; #RA
				($fine_point_stats[4] = abs($add[3])) if ($fine_point_stats[4] < abs($add[3]));
				
				$fine_point_stats[5] += $add[5]; #ROLL
				($fine_point_stats[6] = abs($add[5])) if ($fine_point_stats[6] < abs($add[5]));
			
			} elsif ((index($ssc_mode, 'FINE_SLEW') >= 0) && (index($fits_mode, 'FINE_SLEW') >= 0))  {
				
				$fine_slew_stats[1] += $add[1]; #DEC
				($fine_slew_stats[2] = abs($add[1])) if ($fine_slew_stats[2] < abs($add[1]));
				
				$fine_slew_stats[3] += $add[3]; #RA
				($fine_slew_stats[4] = abs($add[3])) if ($fine_slew_stats[4] < abs($add[3]));
				
				$fine_slew_stats[5] += $add[5]; #ROLL
				($fine_slew_stats[6] = abs($add[5])) if ($fine_slew_stats[6] < abs($add[5]));
			
			}
			
		} else {
		
			print_sub('DEC', $STATSFILE, '~', '~', '~', '~');
			print_sub('RA', $STATSFILE, '~', '~', '~', '~');
			print_sub('ROLL', $STATSFILE, '~', '~', '~', '~');
		
		}
		
		if (exists $fits_ref->{$fits_times}{'LAT[N]'}) {
		
			$fits_latitude = $fits_ref->{$fits_times}{'LAT[N]'};
		
		}
		
		if (exists $fits_ref->{$fits_times}{'LONG[E]'}) {
		
			$fits_longitude = $fits_ref->{$fits_times}{'LONG[E]'};
		
		}
		
		
		printf $STATSFILE "\t%7s %8s", $fits_latitude, $fits_longitude;	
		
		#checks whether the sizes match, prints GOOD if they do, BAD if they don't, and increases $total_size_bad by 1 if BAD
		my $ssc_size = $ssc_ref->{$ssc_times}{'SIZE'};
		my $fits_size = $fits_ref->{$fits_times}{'Size'};
		(printf $STATSFILE "\t%11s %11s %11s ", 'GOOD', "$ssc_size", "$fits_size") if ($ssc_size eq $fits_size);
		(printf $STATSFILE "\t%11s %11s %11s ", 'BAD', "$ssc_size", "$fits_size") if ($ssc_size ne $fits_size);
		($total_size_bad += 1) if ($ssc_size ne $fits_size);
		
		print $STATSFILE "\n";
	}
	
	#used to add totals for all of them to the overall one, for final statistics
	push @return_totals, 
		($total, #0
		$total_fine_point, #1
		$total_mode_bad, #2
		$total_dec_bad, #3
		$total_ra_bad, #4
		$total_roll_bad, #5
		$total_dur_bad, #6
		$total_missing, #7
		$total_incomplete, #8
		$incorrect_time, #9
		$total_size_bad, #10
		$total_rate_bad, #11
		$total_astro, #12
		$total_heoss, #13
		);
	#Closing and returning values
	
	print_sub('FOOTER', $STATSFILE, \@return_totals, \@fine_point_stats, \@fine_slew_stats);
	
	my $info_totals = join('-', @return_totals); #used to add the totals to the compl_stats.info file
	if ($missing_check == 0) {
	
		$compl_info->{completed}->{"$year-$doy"} = $info_totals;
		my $compl_stats_out = File::Spec->catfile($compl_dir, $stats_out_name);
		$compl_info->write( $compl_info_out );
		copy($stats_out, $compl_stats_out);
		print "Writing $compl_info_out \n Copying $stats_out to $compl_stats_out\n";
	
	} else {
		print "End of data @ $year-$doy - Missing check ($missing_check) - Totals  ($info_totals)\n";
	}
	
	close $STATSFILE;
	
	$ssc_ref = undef;
	$fits_ref = undef;
	
	if ($monthly eq 'Yes') {
	
		my ($file_year, $file_month, $file_day) = Add_Delta_Days($year, 1, 1, $doy-1);
		
		if (($file_year == $month_requested[0]) && ($file_month == $month_requested[1])) {
		
			return @return_totals;
		
		}
	
	} else {
	
		return @return_totals;
	
	}
	
}




sub ra_dec_roll { 
	
	#a subroutine to calculate whether or not RA, DEC, or ROLL are good or bad, print it to the thing, and return certain values
	my $ssc_ref = shift;
	my $fits_ref = shift;
	my $ssc_times = shift;
	my $fits_times = shift;
	my $filehandle = shift;
	my $result = '';
	my $add = 0;
	my @categories = ('DEC', 'RA', 'ROLL');
	my @returned_vals = ();
	
	foreach (@categories) {

		my $ssc = $ssc_ref->{$ssc_times}{"$_"};
		my $fits = $fits_ref->{$fits_times}{"$_"};
		my $ssc_mode = $ssc_ref->{$ssc_times}{'MODE'};
		my $fits_mode = $fits_ref->{$fits_times}{'MODE'};
		my $range;
		
		if ($ssc =~ /^ .+ ~ .+ $/x) {
		
			print_sub("$_", $filehandle, 'N/A', $ssc, $fits, 'N/A');
			push @returned_vals, 0;

			if ($_ eq 'ROLL') {
			
				return @returned_vals;
			
			} else { next; }
			
		}
		
		
		
		if ((($ssc_mode, 'FINE_POINT') >= 0) or ( (index($ssc_mode, 'FINE_SETTLE') >= 0) and index($fits_mode, 'FINE_POINT') >= 0)) {
		
			($range = $fp_dec_range) if ($_ eq 'DEC');
			($range = $fp_ra_range) if ($_ eq 'RA');
			($range = $fp_roll_range) if ($_ eq 'ROLL');
		
		} elsif (($ssc_mode, 'FINE_SLEW') >= 0) {
		
			($range = $fs_dec_range) if ($_ eq 'DEC');
			($range = $fs_ra_range) if ($_ eq 'RA');
			($range = $fs_roll_range) if ($_ eq 'ROLL');	
		
		} else {
		
			($range = $def_ra_range) if ($_ eq 'RA');
			($range = $def_dec_range) if ($_ eq 'DEC');
			($range = $def_roll_range) if ($_ eq 'ROLL');
		
		}
		
		my $neg_range = $range * (-1);
		
		my $diff = ($ssc*3600) - ($fits*3600);
		unless ((($diff >= $neg_range) and ($diff <= $range)) or $diff == 0) {
		
			$result = 'BAD';
			push @returned_vals, 1, $diff;
		
		} else {
		
			$result = 'GOOD';
			push @returned_vals, 0, $diff;
		
		}
		
		print_sub("$_", $filehandle, $result, $ssc, $fits, $diff);
	
	}

	return @returned_vals;
	
}



sub time_iterator {

	#used to increase the time per second in an accurate way by splitting the test time initially into @test_times = (YYYY, DDD, HH, MM, SS)
	#and iterating the second element by 1, and changing itself and other values when it hits 60, this array then is joined and assigned to 
	#variable $test_times without changing the array, and $test_times is checked with the fits hash to see if it exists in it

	my $next = shift; #time of next image being taken
	my $test_times = shift; #used to check whether a time exists in the fits hash
	my $fits_ref = shift;
	my $incorrect_time = shift;
	my $sec_range = shift;
	my $ssc_ref = shift;
	
	#$next =~ s/[-:]//g;
	my @test_times = split /[-:]/, $test_times;
	my $fits_times = $test_times;
	my $diff_sec = 0;
	my @returned_vals = ();

	while (($next >= $test_times)||($diff_sec <= 10 && $next eq 'final')) {
		$diff_sec += 1;
		
		#next several if conditionals used to just make sure the numbers count up the right way for times (0-60 vs 0-100)
		$test_times[4] += 1;
		$test_times[4] = sprintf("%02d", $test_times[4]); #all sprintf are used to keep the correct # of zeroes in front
			if ($test_times[4] == 60) {
			
				$test_times[4] = '00'; #sets the second value back to 0, increases the minute value up 1
				$test_times[3] += 1;
				$test_times[3] = sprintf("%02d", $test_times[3]);
				
				if ($test_times[3] == 60) {
				
					$test_times[3] = '00'; #sets the minute value back to 0, increases the hour value up 1
					$test_times[2] += 1;
					$test_times[2] = sprintf("%02d", $test_times[2]);
					
					if ($test_times[2] == 24) {
					
						$test_times[2] = '00'; #sets the hour value back to 0, increases the day of year by 1
						$test_times[1] += 1;
						$test_times[1] = sprintf("%03d", $test_times[1]);
						
						if ($test_times[1] == 366) {
						
							$test_times[1] = 1; #sets the day of year to 1 (january 1st), and increases the year value by 1
							$test_times[0] += 1;
						
						}
					
					}
				
				}
			
			}## --- the time/date related iteration stuff ends here --- ##
		
		$test_times = join('-', @test_times); #sets $test_times to the joined values of @test_times, does not alter the array
		$test_times =~ s/^ (\d{4}) - (\d{3}) - (\d{2}) - (\d{2}) - (\d{2}) $/$1-$2-$3:$4:$5/x;
		#checks to see if there's a key that corresponds to the new value of $test_value
		if (exists $fits_ref->{$test_times} && $test_times ne $next) {
			
			if (exists $ssc_ref->{$test_times} && $diff_sec != 0) { #if another thing with the same time exists but the difference in seconds is 0
			
				last;
			
			}
			#changes $fits_times from being equal to $ssc_times to the now working $test_times
			$fits_times = $test_times;
			($incorrect_time += 1) if ($diff_sec > $sec_range);
			last;
		
		} elsif ($test_times == $next) {
		
			#ends the iterations if $test_times is equal to $next, as the next image is scheduled
			last;
		
		}
	
	}

	push @returned_vals, $diff_sec;
	push @returned_vals, $incorrect_time;
	push @returned_vals, $fits_times;
	return @returned_vals;
	
}



sub rate_conv {

	my $ssc_rate = shift;
	my $fits_ra_vel = shift;
	my $fits_dec_vel = shift;
	my $mode = shift;
	my $filehandle = shift;
	my $range = $def_rate_range;  # default rate tolerance
	my $fits_rate = 0;
	my $diff_rate = 0;
	my $rate_bad = 0;
	my @returned_vals = ();
	
	
	if (index($mode, 'FINE_POINT') >= 0) {
		$range = $fp_rate_range;
	}
	elsif (index($mode, 'FINE_SLEW') >= 0) {
		$range = $fs_rate_range;
	} 

	my $fits_rate = sqrt($fits_ra_vel**2 + $fits_dec_vel**2);
	my $diff_rate = abs($ssc_rate - $fits_rate);
	
	unless (($diff_rate <= $range) || ($diff_rate == 0) || ($mode eq 'DESAT') || ($mode eq 'RATE_SLEW')) {
	
		$rate_bad += 1;
		printf $filehandle "\t\t%6s %9.2f %7.2f %6.2f ", 'BAD', $ssc_rate, $fits_rate, $diff_rate;
	
	} elsif ($mode eq 'DESAT' || $mode eq 'RATE_SLEW') {
		
		printf $filehandle "\t\t%6s %9s %7s %6s ", '~', '~', '~', '~';
		
	} else {
	
		printf $filehandle "\t\t%6s %9.2f %7.2f %6.2f ", 'GOOD', $ssc_rate, $fits_rate, $diff_rate;
	
	}
	
	push @returned_vals, $fits_rate;
	push @returned_vals, $diff_rate;
	push @returned_vals, $rate_bad;
	return @returned_vals;

}



#This differentiates between the input (monthly or range) and any desired category requests
sub range_stats {

	my $aref = shift;
	my $filename_base = shift;
	my @returned_vals = @{$aref};
	my @matched_files = ();
	my $RANGEBAD;
	my $RANGEALL;
	my $RANGECAT;
	my $cat_idx = 0;
	my $check = '';
	my $counter = 0;

	
	opendir(FILES, "$out_dir");
	my @stats_files = grep(/daily\.stats$/, readdir FILES);
	closedir FILES;
	
	if ($range_type eq 'MONTHLY') {
	
		foreach (sort @stats_files) {
		
			$_ =~ /^NEOS_ (\d{4}) - (\d{3}) _daily\.stats $/ix;
			
			my ($year, $month, $day) = Add_Delta_Days($1, 1, 1, $2-1);
			#print "$year-$month-$day ignored $_\n";
			
			if (($year == $month_requested[0]) && ($month == $month_requested[1])) {
			
				push @matched_files, $_;

			}
		
		}
	
	} elsif ($range_type eq 'RANGE') {
	
		foreach (sort @stats_files) {
		
			my @stats_date = split /[-_]/, $_;
			my $stats_date = int($stats_date[1]*365.25 + $stats_date[2]);
			my $abs_start = int($start_date[0]*365.25 + $start_date[1]);
			my $abs_end = int($end_date[0]*365.25 + $end_date[1]);
			
			next if ($stats_date < $abs_start || $stats_date > $abs_end);
			
			push @matched_files, $_;
		
		}
	
	}

	# Define files 
        my $BAD_out = File::Spec->catfile($out_dir, $filename_base . $range_type . "_BAD.stats");	
        my $ALL_out = File::Spec->catfile($out_dir, $filename_base . $range_type . "_ALL.stats");	
        my $CATREQ_out = File::Spec->catfile($out_dir, $filename_base . $cat_req . ".stats");	
        my $CATMOD_out = File::Spec->catfile($out_dir, $filename_base  
		                                       . $cat_modifiers[0] . "," . $cat_modifiers[1] . "," . $cat_modifiers[2] 
						       . ".stats");	
        my $TEMP_out = File::Spec->catfile($out_dir, "temp.txt");

	if ($use_obsv == 0) {
		open ($RANGEBAD, ">$BAD_out") 
		or die "cannot open/create $BAD_out: $!";
		
		open ($RANGEALL, ">$ALL_out") 
		or die "cannot open/create $ALL_out: $!";
		
		print_sub('TOTAL_HEADER', $RANGEBAD, $range_type, $aref);
		print_sub('TOTAL_HEADER', $RANGEALL, $range_type, $aref);
		print_sub('COLUMN_NAMES', $RANGEBAD);
		print_sub('COLUMN_NAMES', $RANGEALL);
		
		if ($cat_req ne '') {
		
			$cat_idx = $cat_idcs{$cat_req};
			($cat_idx = $cat_idcs{$cat_modifiers[1]}) if ($cat_req eq 'MODE' and $cat_modifiers[1] ne '');
			
			if ($cat_idx == 6) {
		
				$check = 'NO';
			
			} elsif ($cat_idx == 18) {
			
				$check = $cat_modifiers[1];
			
			} else {
			
				$check = 'BAD';
			
			}
		
			unless ($cat_modifiers[1] ne '') {
			
				open ($RANGECAT, ">$CATREQ_out") 
				or die "cannot open/create $CATREQ_out: $!";
			
			} else {
			
				open ($RANGECAT, ">$CATMOD_out")
				or die "cannot open/create $CATMOD_out: $!";
			
			}
			
			print_sub('TOTAL_HEADER', $RANGECAT, 'RANGE', $aref);
			print_sub('COLUMN_NAMES', $RANGECAT);
		
		}
		
	} elsif ($use_obsv == 1) {
	
		open ($RANGECAT, ">$CATREQ_out") 
		or die "cannot open/create $CATREQ_out: $!";
		
		open TEMP, ">$TEMP_out" or die "Cannot open $TEMP_out: $!";
		
		$cat_idx = $cat_idcs{$cat_req};
		$check = $cat_req; 
		
		@returned_vals = ();
	
	}

	foreach my $stats_file (@matched_files) {
	        my $stats_in = File::Spec->catfile($out_dir, $stats_file);	
		open STATS_IN, "<$stats_in" or die "Could not open $stats_in: $!";
		
		while (<STATS_IN>) {
		
			chomp;
			
			if ($_ =~ /^ \d /x) { #prints lines to $RANGEBAD that contain the word "BAD" or "YES" which only appears when one is missing
				
				if ($use_obsv == 0) {
				
					if ($_ =~ /bad|no/ig) {
					
						print $RANGEBAD "$_\n";
						
					}			
				
					print $RANGEALL "$_\n";
					
					if ($cat_req ne '') {
					
						my @split = split /\s+/, $_;
						
						if ($cat_idx == 18 and $cat_modifiers[2] eq 'BOTH') {
						
							(print $RANGECAT "$_\n") if ($split[$cat_idx] eq $check and $split[$cat_idx + 1] eq $check);
						
						} elsif ($cat_idx == 18 and $cat_modifiers[2] eq 'NOT') {
						
							(print $RANGECAT "$_\n") if ($split[$cat_idx] ne $check);
						
						} elsif ($cat_idx == 18 and $cat_modifiers[2] eq 'BOTH-NOT') {
						
							(print $RANGECAT "$_\n") if ($split[$cat_idx] ne $check and $split[$cat_idx + 1] ne $check);
						
						} else {
						
							(print $RANGECAT "$_\n") if ($split[$cat_idx] eq $check);
						
						}
					}
				
				} elsif ($use_obsv == 1) {
				
					my @split = split /\s+/, $_;
					
					if ($split[$cat_idx] eq $check) {
					
						print TEMP "$_\n";
						$returned_vals[0] += 1;
						
						unless ($split[6] eq 'NO') {
							
							$returned_vals[1] += 0; #missing
							($returned_vals[2] += 1) if ($split[7] eq 'BAD'); #completed
							($returned_vals[3] += 1) if ($split[9] eq 'BAD'); #start time 
							($returned_vals[4] += 1) if ($split[13] eq 'BAD'); #Duration
							($returned_vals[5] += 1) if ($split[17] eq 'BAD'); #Mode
							($returned_vals[6] += 1) if ($split[20] eq 'BAD'); #Rate
							
							($returned_vals[7] += 1) if (($split[18] ne 'DESAT') #number of non-desat, except where the expected RA/DEC/ROLL are ranges not specific numbers (initial~final)
														and (($split[25] !~ /^ .+ ~ .+ $/x) 	#index 25 is expected Declination
														or ($split[29] !~ /^ .+ ~ .+ $/x) 	#index 29 is expected RA
														or ($split[33] !~ /^ .+ ~ .+ $/x) )); #index 33 is expected Roll
													
							($returned_vals[8] += 1) if ($split[24] eq 'BAD'); #Declination 
							($returned_vals[9] += 1) if ($split[28] eq 'BAD'); #Right Ascension 
							($returned_vals[10] += 1) if ($split[32] eq 'BAD'); #Roll 
							($returned_vals[11] += 1) if ($split[36] eq 'BAD'); #Size 
						
						} else {
						
							$returned_vals[1] += 1;
							$returned_vals[2] += 0;
							$returned_vals[3] += 0;
							$returned_vals[4] += 0;
							$returned_vals[5] += 0;
							$returned_vals[6] += 0;
							$returned_vals[7] += 0;
							$returned_vals[8] += 0;
							$returned_vals[9] += 0;
							$returned_vals[10] += 0;
							$returned_vals[11] += 0;
						
						}
						($returned_vals[12] += 1) if (($_ =~ /BAD/) || $split[6] eq 'NO');
					}
				}
			}
		}
	
		close STATS_IN;
		
	}
	
	close TEMP;
	
	if ($use_obsv == 0) {
	
		close $RANGEBAD;
		close $RANGEALL;

		(close $RANGECAT) if ($cat_req ne '');
	
	} elsif ($use_obsv == 1) {
	
		open TEMP, "<$TEMP_out" or die "Cannot open $TEMP_out: $!";
		print_sub('CATREQ HEADER', $RANGECAT, $cat_req, \@returned_vals);
		print_sub('COLUMN_NAMES', $RANGECAT);
		print "$cat_req\n";
		
		while (<TEMP>) {
		
			print $RANGECAT "$_";
		
		}
		
		close TEMP;
		close $RANGECAT;
		unlink $TEMP_out;
	
	}
	
}


sub print_sub { #can use this to easily organize the various columns in any order wanted, since it's all in one place, also prints really large blocks of prints

	my $usage = shift;
	my $filehandle = shift;
	
	if ($usage eq 'COLUMN_NAMES') { ##Prints the column names, two lines used to keep unnecessary length down and keep things close together
	
		printf $filehandle "%-6s %5s %10s %26s %26s %9s %10s ", 'YEAR', 'DAY', 'TIME', 'SSCLOG', 'SSCLOG', 'SSCLOG', 'IMAGE';	#VARIOUS
		printf $filehandle "\t%8s %9s ", 'RESLT', 'PRCNT';															#COMPLETION
		printf $filehandle "\t%12s %15s %15s %7s ", 'RESLT', 'EXPECTED', 'ACTUAL', 'DIFF';							#START TIME
		printf $filehandle "\t%7s %9s %9s %5s ", 'RESLT', 'EXPECTED', 'ACTUAL', 'DIFF';								#DURATION	
		printf $filehandle "\t%6s %15s %15s ", 'RESLT', 'EXPECTED', 'ACTUAL';										#MODE
		printf $filehandle "\t\t%6s %9s %7s %6s ", 'RESLT', 'EXPECTED', 'ACTUAL', 'DIFF';							#RATE
		printf $filehandle "\t%6s %15s %10s %13s ", 'RESLT', 'EXPECTED', 'ACTUAL', 'DIFF';							#DEC
		printf $filehandle "\t%6s %15s %10s %13s ", 'RESLT', 'EXPECTED', 'ACTUAL', 'DIFF';							#RA
		printf $filehandle "\t%6s %15s %10s %13s ", 'RESLT', 'EXPECTED', 'ACTUAL', 'DIFF';							#ROLL
		printf $filehandle "\t%7s %8s", 'LAT', 'LONG';																#LATITUDE AND LONGITUDE
		printf $filehandle "\t%11s %11s %11s ", 'RESLT', 'EXPECTED', 'ACTUAL';										#SIZE
		printf $filehandle "\n";
	
		printf $filehandle "%-6s %5s %10s %26s %26s %9s %10s ", 'YYYY', 'DDD', 'HH:MM:SS', 'SOURCE_FILE', 'OBJECT', 'OBSERVER', 'AVAILABLE';	#VARIOUS
		printf $filehandle "\t%8s %9s ", '%COMPL', 'COMPL';																						#COMPLETION
		printf $filehandle "\t%12s %15s %15s %7s ", 'START_TIME', 'TIME', 'TIME', 'TIME';														#START TIME
		printf $filehandle "\t%7s %9s %9s %5s ", 'DUR', 'DUR', 'DUR', 'DUR';																	#DURATION		
		printf $filehandle "\t%6s %15s %15s ", 'MODE', 'MODE', 'MODE';																			#MODE
		printf $filehandle "\t\t%6s %9s %7s %6s ", 'RATE', 'RATE', 'RATE', 'RATE';																#RATE
		printf $filehandle "\t%6s %15s %10s %13s ", 'DEC', 'DEC', 'DEC', 'DEC(arcsec)';															#DEC
		printf $filehandle "\t%6s %15s %10s %13s ", 'RA', 'RA', 'RA', 'RA(arcsec)';																#RA
		printf $filehandle "\t%6s %15s %10s %13s ", 'ROLL', 'ROLL', 'ROLL', 'ROLL(arcsec)';														#ROLL
		printf $filehandle "\t%7s %8s", 'North', 'East';																						#LAT MEASURED NORTH, LONG MEASURED EAST
		printf $filehandle "\t%11s %11s %11s ", 'RAST_SIZE', 'SIZE', 'SIZE';																	#SIZE
		printf $filehandle "\t\n;\n";
	
	} elsif ($usage eq 'MISSING') {
	
		my $ssc_ref = shift;
		my $ssc_start = shift;
		my $exp_start;
		my $exp_size;
		my $exp_mode;
		my $exp_dur;
		my $exp_rate;
		my $exp_dec;
		my $exp_ra;
		my $exp_roll;
		
		my ($year, $doy, $time) = split /-/, $ssc_start;
		$exp_start = join('',$year,$doy,$time);
		$exp_start =~ s/[:]//g;
		$exp_size = $ssc_ref->{$ssc_start}{'SIZE'};
		$exp_mode = $ssc_ref->{$ssc_start}{'MODE'};
		$exp_dur = $ssc_ref->{$ssc_start}{'DUR'};
		$exp_rate = $ssc_ref->{$ssc_start}{'RATE'};
		
		if (index($exp_mode, 'FINE_POINT') >= 0) {
		
			$exp_dec = $ssc_ref->{$ssc_start}{'DEC'};
			$exp_ra = $ssc_ref->{$ssc_start}{'RA'};
			$exp_roll = $ssc_ref->{$ssc_start}{'ROLL'};
		
		} else {
		
			$exp_dec = '~';
			$exp_ra = '~';
			$exp_roll = '~';
		
		}

		printf $filehandle "%10s ", 'NO';
		printf $filehandle "\t%8s %9s ", '~', '~'; 												#Percent Completion
		printf $filehandle "\t%12s %15s %15s %7s ", '~', $exp_start, '~', '~'; 					#Start_time/Expected/Actual/Diff
		printf $filehandle "\t%7s %9s %9s %5s ", '~', $exp_dur, '~', '~'; 						#duration/expected/actual/difference
		printf $filehandle "\t%6s %15s %15s ", '~', $exp_mode, '~'; 							#mode/Expected/actual
		printf $filehandle "\t\t%6s %9s %7s %6s ", '~', $exp_rate, '~', '~'; 					#rate/expected/actual/difference
		printf $filehandle "\t%6s %15s %10s %13s ", '~', $exp_dec, '~', '~'; 					#declination/expected/actual/difference
		printf $filehandle "\t%6s %15s %10s %13s ", '~', $exp_ra, '~', '~'; 					#right_ascension/expected/actual/difference
		printf $filehandle "\t%6s %15s %10s %13s ", '~', $exp_roll, '~', '~';					#roll/expected/actual/difference
		printf $filehandle "\t%7s %8s", '~', '~';												#LATITUDE/LONGITUDE
		printf $filehandle "\t%11s %11s %11s ", '~', $exp_size, '~'; 							#raster_size/expected/actual
		printf $filehandle "\n";
	
	} elsif ($usage eq 'RA' || $usage eq 'DEC' || $usage eq 'ROLL') {
		
		my $state = shift; #BAD or GOOD
		my $expected = shift;
		my $actual = shift;
		my $diff = shift;
	
		(printf $filehandle "\t%6s %15s %10s %13.2f ", $state, $expected, $actual, $diff) if ($state ne '~');
		(printf $filehandle "\t%6s %15s %10s %13s ", $state, $expected, $actual, $diff) if ($state eq '~');
	
	} elsif ($usage eq 'FOOTER') {
	
		my $rt_aref = shift;
		my $fp_aref = shift;
		my $fs_aref = shift;
		
		my @return_totals = @{$rt_aref};
		my @fine_point_stats = @{$fp_aref};
		my @fine_slew_stats = @{$fs_aref};
		
		my $header = '#' x 150; #a header, just a bunch of # in a row
	
		print $filehandle ";\n$header\n$header\n$header\n;\n";
		print $filehandle ";DAY STATS: TOTAL_IMAGES = $return_totals[0]\n";
		if ($return_totals[0] > 0) {
			printf $filehandle "%-15s %12s \t%6.2f% \n", "\t;TOTAL_ASTRO:", "$return_totals[12]/$return_totals[0] -> ", ($return_totals[12]/$return_totals[0])*100;
			printf $filehandle "%-15s %12s \t%6.2f% \n", "\t;TOTAL_HEOSS:", "$return_totals[13]/$return_totals[0] -> ", ($return_totals[13]/$return_totals[0])*100;
			printf $filehandle "%-15s %12s \t%6.2f% \t%-30s \n", "\t;START_TIME_BAD:", "$return_totals[9]/$return_totals[0] -> ", ($return_totals[9]/$return_totals[0])*100, "(Acceptable Variance: $sec_range seconds)";
			printf $filehandle "%-15s \t%12s \t%6.2f% \n", "\t;MISSING:", "$return_totals[7]/$return_totals[0] -> ", ($return_totals[7]/$return_totals[0])*100;
			printf $filehandle "%-15s \t%12s \t%6.2f% \n", "\t;INCOMPLETE:", "$return_totals[8]/$return_totals[0] -> ", ($return_totals[8]/$return_totals[0])*100;
			printf $filehandle "%-15s \t%12s \t%6.2f% \n", "\t;MODE_BAD:", "$return_totals[2]/$return_totals[0] -> ", ($return_totals[2]/$return_totals[0])*100;
			printf $filehandle "%-15s \t%12s \t%6.2f% \t%-30s \n", "\t;DUR_BAD:", "$return_totals[6]/$return_totals[0] -> ", ($return_totals[6]/$return_totals[0])*100, "Acceptable_Variance: $dur_range seconds";
			printf $filehandle "%-15s \t%12s \t%6.2f% \n", "\t;SIZE_BAD:", "$return_totals[10]/$return_totals[0] -> ", ($return_totals[10]/$return_totals[0])*100;
			printf $filehandle "%-15s \t%12s \t%6.2f% \t%-30s \n", "\t;RATE_BAD:", "$return_totals[11]/$return_totals[0] -> ", ($return_totals[11]/$return_totals[0])*100, 
			                   "Acceptable_Variance: (FINE_POINT: $fp_rate_range arcsec/sec, FINE_SLEW: $fs_rate_range arcsec/sec)";
		}

		if ($return_totals[1] > 0) {

			printf $filehandle "%-15s \t%12s \t%6.2f% \t%-30s \n", "\t;RA_BAD:", "$return_totals[4]/$return_totals[1] -> ", ($return_totals[4]/$return_totals[1])*100, 
			                    "Acceptable_Variance: (FINE_POINT: $fp_ra_range arcsec, FINE_SLEW: $fs_ra_range arcsec, DEFAULT: $def_ra_range arcsec)";
			printf $filehandle "%-15s \t%12s \t%6.2f% \t%-30s \n", "\t;DEC_BAD:", "$return_totals[3]/$return_totals[1] -> ", ($return_totals[3]/$return_totals[1])*100, 
				            "Acceptable_Variance: (FINE_POINT: $fp_dec_range arcsec, FINE_SLEW: $fs_dec_range arcsec, DEFAULT: $def_dec_range arcsec)";
			printf $filehandle "%-15s \t%12s \t%6.2f% \t%-30s \n", "\t;ROLL_BAD:", "$return_totals[5]/$return_totals[1] -> ", ($return_totals[5]/$return_totals[1])*100, 
			                    "Acceptable_Variance: (FINE_POINT: $fp_roll_range arcsec, FINE_SLEW: $fs_roll_range arcsec, DEFAULT: $def_roll_range arcsec)";
		}
		if ($fine_point_stats[0] > 0) {
		
			print $filehandle ";\n;FINE_POINT: $fine_point_stats[0] IMAGES\n";
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcseconds,\tMax_Deviation = %-12.2f Arcseconds\n", "Fine_Point_RA:", ($fine_point_stats[3]/$fine_point_stats[0]), $fine_point_stats[4];
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcseconds,\tMax_Deviation = %-12.2f Arcseconds\n", "Fine_Point_Dec:", ($fine_point_stats[1]/$fine_point_stats[0]), $fine_point_stats[2];
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcseconds,\tMax_Deviation = %-12.2f Arcseconds\n", "Fine_Point_Roll:", ($fine_point_stats[5]/$fine_point_stats[0]), $fine_point_stats[6];
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcsec/sec,\tMax_Deviation = %-12.2f Arcsec/sec\n", "Fine_Point_Rate:", ($fine_point_stats[7]/$fine_point_stats[0]), $fine_point_stats[8];
		
		} 
		
		if ($fine_slew_stats[0] > 0) {
		
			print $filehandle ";\n;FINE_SLEW: $fine_slew_stats[0] IMAGES\n";
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcseconds,\tMax_Deviation = %-12.2f Arcseconds\n", "Fine_Slew_RA:", ($fine_slew_stats[3]/$fine_slew_stats[0]), $fine_slew_stats[4];
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcseconds,\tMax_Deviation = %-12.2f Arcseconds\n", "Fine_Slew_Dec:", ($fine_slew_stats[1]/$fine_slew_stats[0]), $fine_slew_stats[2];
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcseconds,\tMax_Deviation = %-12.2f Arcseconds\n", "Fine_Slew_Roll:", ($fine_slew_stats[5]/$fine_slew_stats[0]), $fine_slew_stats[6];
			printf $filehandle "\t;%-16s \tAverage Deviation = %-8.2f Arcsec/sec,\tMax_Deviation = %-12.2f Arcsec/sec\n", "Fine_Slew_Rate:", ($fine_slew_stats[7]/$fine_slew_stats[0]), $fine_slew_stats[8];
		
		} 
		
	
	} elsif ($usage eq 'HEADER') {
	
		my $stats_out_name = shift;
		my $year = shift;
		my $doy = shift;
	
		my $header = '#' x 150; #a header, just a bunch of # in a row
		
		print $filehandle "$header\n#\n#\tFilename:\t\t\t$stats_out_name\n#\n#\tInput_Files:\t\t\t$year\_$doy\_SSC.txt and $year\_$doy\_FITS.txt\n#\n#\tVersion:\t\t\t$VER\n#";
		print $filehandle "\n$header\n";
		
		print $filehandle ";ACCEPTABLE_VARIANCES:\n;\n";
		print $filehandle ";START_TIME:\t\t$sec_range seconds\n;DURATION:\t\t$dur_range seconds\n";
		print $filehandle ";RIGHT_ASCENSION:\t(FINE_POINT: $fp_ra_range arcseconds, FINE_SLEW: $fs_ra_range arcseconds, Default: $def_ra_range arcseconds)\n";
		print $filehandle ";DECLINATION:\t\t(FINE_POINT: $fp_dec_range arcseconds, FINE_SLEW: $fs_dec_range arcseconds, Default: $def_dec_range arcseconds)\n";
		print $filehandle ";ROLL:\t\t\t(FINE_POINT:  $fp_roll_range arcseconds, FINE_SLEW: $fs_roll_range arcseconds, Default: $def_roll_range arcseconds)\n";
		print $filehandle ";RATE:\t\t\t(FINE_POINT: $fp_rate_range arcsec/sec, FINE_SLEW: $fs_rate_range arcsec/sec)\n;\n;\n";
	
	} elsif ($usage eq 'TOTAL_HEADER') {
	
		my $usage2 = shift;
		my $aref = shift;
		my @returned_vals = @{$aref};
		
		my $image_total += $returned_vals[0];
		my $image_fine_point += $returned_vals[1];
		my $image_mode_bad += $returned_vals[2];
		my $image_dec_bad += $returned_vals[3];
		my $image_ra_bad += $returned_vals[4];
		my $image_roll_bad += $returned_vals[5];
		my $image_dur_bad += $returned_vals[6];
		my $image_missing += $returned_vals[7];
		my $image_incomplete += $returned_vals[8];
		my $image_start_bad += $returned_vals[9];
		my $image_size_bad += $returned_vals[10];
		my $image_rate_bad += $returned_vals[11];
		my $image_astro += $returned_vals[12];
		my $image_heoss += $returned_vals[13];
		
		if ($usage2 eq 'MONTHLY') {
			
			my $used_month = shift;
		
			printf $filehandle "%s\n#\n#\tMonthly_Statistics_Log\n#\n#\tINTERVAL: %17s-%-3s\n#\n#\tTOTAL_IMAGES: %12s\n#\n#\tVERSION: %17s\n#\n%s\n", 
			'#' x 147, $start_date[0], $used_month, $image_total, $VER,'#' x 147;
		
		} elsif ($usage2 eq 'RANGE') {
		
			printf $filehandle "%s\n#\n#\tCumulative_Statistics_Log\n#\n#\tINTERVAL_START: %10s-%-3s\n#\n#\tINTERVAL_END: %12s-%-3s\n#\n#\tTOTAL_IMAGES: %12s\n#\n#\tVERSION: %17s\n#\n%s\n", 
			'#' x 147, $start_date[0], $start_date[1], $end_date[0], $end_date[1], $image_total, $VER,'#' x 147;			
		
		}
		
		
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%\n", ";OBSERVER:", "ASTRO:", "$image_astro/$image_total ->", ($image_astro/$image_total)*100, "HEOSS:", "$image_heoss/$image_total -> ", ($image_heoss/$image_total)*100;
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%\n", ";IMAGE:", "AVAILABLE:", ($image_total-$image_missing) . "/$image_total ->", (1-($image_missing/$image_total))*100, "MISSING:", "$image_missing/$image_total -> ", ($image_missing/$image_total)*100;
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%\n", ";COMPLETION:", "COMPLETE:", ($image_total-$image_incomplete) . "/$image_total ->", (1-($image_incomplete/$image_total))*100, "INCOMPLETE:", "$image_incomplete/$image_total -> ", ($image_incomplete/$image_total)*100;
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%,\t%-30s\n", ";START_TIME:", "GOOD:", ($image_total-$image_start_bad) . "/$image_total ->", (1-($image_start_bad/$image_total))*100, "BAD:", "$image_start_bad/$image_total -> ", ($image_start_bad/$image_total)*100, "Acceptable_Variance: $sec_range seconds";
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%,\t%-30s\n", ";DURATION:", "GOOD:", ($image_total-$image_dur_bad) . "/$image_total ->", (1-($image_dur_bad/$image_total))*100, "BAD:", "$image_dur_bad/$image_total -> ", ($image_dur_bad/$image_total)*100, "Acceptable_Variance: $dur_range seconds";
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%\n", ";MODE:", "CORRECT:", ($image_total-$image_mode_bad) . "/$image_total ->", (1-($image_mode_bad/$image_total))*100, "INCORRECT:", "$image_mode_bad/$image_total -> ", ($image_mode_bad/$image_total)*100;
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%,\t%-30s\n", ";RATE:", "GOOD:", ($image_total-$image_rate_bad) . "/$image_total ->", (1-($image_rate_bad/$image_total))*100, "BAD:",  "$image_rate_bad/$image_total -> ", ($image_rate_bad/$image_total)*100, "Acceptable_Variance: (FINE_POINT: $fp_rate_range arcsec/sec, FINE_SLEW: $fs_rate_range arcsec/sec)";
		(printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%,\t%-30s\n", ";DECLINATION:", "GOOD:", ($image_fine_point-$image_dec_bad) . "/$image_fine_point ->", (1-($image_dec_bad/$image_fine_point))*100, "BAD:", "$image_dec_bad/$image_fine_point -> ", ($image_dec_bad/$image_fine_point)*100, "Acceptable_Variance: (FINE_POINT: $fp_dec_range arcsec, FINE_SLEW: $fs_dec_range arcsec, DEFAULT: $def_dec_range arcsec)") if ($image_fine_point > 0);
		(printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%,\t%-30s\n", ";RIGHT_ASC:", "GOOD:", ($image_fine_point-$image_ra_bad) . "/$image_fine_point ->", (1-($image_ra_bad/$image_fine_point))*100, "BAD:", "$image_ra_bad/$image_fine_point -> ", ($image_ra_bad/$image_fine_point)*100, "Acceptable_Variance: (FINE_POINT: $fp_ra_range arcsec, FINE_SLEW: $fs_ra_range arcsec, DEFAULT: $def_ra_range arcsec)") if ($image_fine_point > 0);
		(printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%,\t%-30s\n", ";ROLL:", "GOOD:", ($image_fine_point-$image_roll_bad) . "/$image_fine_point ->", (1-($image_roll_bad/$image_fine_point))*100, "BAD:", "$image_roll_bad/$image_fine_point -> ", ($image_roll_bad/$image_fine_point)*100, "Acceptable_Variance: (FINE_POINT: $fp_roll_range arcsec, FINE_SLEW: $fs_roll_range arcsec, DEFAULT: $def_roll_range arcsec)") if ($image_fine_point > 0);
		printf $filehandle "%-16s %12s %14s %7.2f%, %12s %14s %7.2f%\n", ";RASTER_SIZE:", "CORRECT:", ($image_total-$image_size_bad) . "/$image_total ->", (1-($image_size_bad/$image_total))*100, "INCORRECT:", "$image_size_bad/$image_total -> ", ($image_size_bad/$image_total)*100;
		printf $filehandle "%s\n", '~' x 350;
	
	} elsif ($usage eq 'CATREQ HEADER') {
	
		my $usage2 = shift;
		my $aref = shift;
		my @returned_vals = @{$aref};
		
		my $image_total += $returned_vals[0];
		my $image_missing += $returned_vals[1];
		my $image_incomplete += $returned_vals[2];
		my $image_start_bad += $returned_vals[3];
		my $image_dur_bad += $returned_vals[4];
		my $image_mode_bad += $returned_vals[5];
		my $image_rate_bad += $returned_vals[6];
		my $image_fine_point += $returned_vals[7];
		my $image_dec_bad += $returned_vals[8];
		my $image_ra_bad += $returned_vals[9];
		my $image_roll_bad += $returned_vals[10];
		my $image_size_bad += $returned_vals[11];
		my $total_cat_bad += $returned_vals[12];
		
		printf $filehandle "%s\n#\n#\t$usage2\_Statistics_Log\n#\n#\tINTERVAL_START: %10s-%-3s\n#\n#\tINTERVAL_END: %12s-%-3s\n#\n#\tTOTAL_IMAGES: %12s\n#\n#\tVERSION: %17s\n#\n%s\n", 
		'#' x 147, $start_date[0], $start_date[1], $end_date[0], $end_date[1], $image_total, $VER,'#' x 147;
		
		printf $filehandle "%-16s %14s %7.2f%,\t%-30s\n", ";START_TIME_BAD:", "$image_start_bad/$image_total -> ", ($image_start_bad/$image_total)*100, "Acceptable_Variance: $sec_range seconds";
		printf $filehandle "%-16s %14s %7.2f% \n", ";MISSING:", "$image_missing/$image_total -> ", ($image_missing/$image_total)*100;
		printf $filehandle "%-16s %14s %7.2f% \n", ";INCOMPLETE:", "$image_incomplete/$image_total -> ", ($image_incomplete/$image_total)*100;
		printf $filehandle "%-16s %14s %7.2f% \n", ";MODE_BAD:", "$image_mode_bad/$image_total -> ", ($image_mode_bad/$image_total)*100;
		printf $filehandle "%-16s %14s %7.2f%,\t%-30s\n", ";DUR_BAD:", "$image_dur_bad/$image_total -> ", ($image_dur_bad/$image_total)*100, "Acceptable_Variance: $dur_range seconds";
		printf $filehandle "%-16s %14s %7.2f% \n", ";SIZE_BAD:", "$image_size_bad/$image_total -> ", ($image_size_bad/$image_total)*100;
		printf $filehandle "%-16s %14s %7.2f%,\t%-30s\n", ";RATE_BAD:", "$image_rate_bad/$image_total -> ", ($image_rate_bad/$image_total)*100, "Acceptable_Variance: (FINE_POINT: $fp_rate_range arcsec/sec, FINE_SLEW: $fs_rate_range arcsec/sec)";
		(printf $filehandle "%-16s %14s %7.2f%,\t%-30s\n", ";DEC_BAD:", "$image_dec_bad/$image_fine_point -> ", ($image_dec_bad/$image_fine_point)*100, "Acceptable_Variance: (FINE_POINT: $fp_dec_range arcsec, FINE_SLEW: $fs_dec_range arcsec, DEFAULT: $def_dec_range arcsec)") if ($image_fine_point > 0);
		(printf $filehandle "%-16s %14s %7.2f%,\t%-30s\n", ";RA_BAD:", "$image_ra_bad/$image_fine_point -> ", ($image_ra_bad/$image_fine_point)*100, "Acceptable_Variance: (FINE_POINT: $fp_ra_range arcsec, FINE_SLEW: $fs_ra_range arcsec, DEFAULT: $def_ra_range arcsec)") if ($image_fine_point > 0);
		(printf $filehandle "%-16s %14s %7.2f%,\t%-30s\n", ";ROLL_BAD:", "$image_roll_bad/$image_fine_point -> ", ($image_roll_bad/$image_fine_point)*100, "Acceptable_Variance: (FINE_POINT: $fp_roll_range arcsec, FINE_SLEW: $fs_roll_range arcsec, DEFAULT: $def_roll_range arcsec)") if ($image_fine_point > 0);
		printf $filehandle "%-16s %14s %7.2f% \n", ";$usage2\_BAD:", "$total_cat_bad/$image_total -> ", ($total_cat_bad/$image_total)*100;
		printf $filehandle "%s\n", '~' x 350;
	
	}
	
	return;

}

close DEBUG; close SSCFILE; close FITSFILE;
