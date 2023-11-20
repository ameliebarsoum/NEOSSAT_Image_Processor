#!usr/bin/perl


use strict;
use List::MoreUtils qw(first_index);
use IO::Handle;
use Config::Tiny;
use File::Copy;
use File::Path;
	
#chdir "C:\\Users\\RKrieger\\Desktop\\test_directory" or die "cannot change to that directory"; #using this as an export for debug purposes, easy to access

# Global variables
my %data_hash = ();
my $data_hash = 0;

# Directories, populated from config file (but passed into function by parameter)
my $ssc_dir    = undef;  #SSCLOG directory, used to read SSCLOG files from
my $fits_dir   = undef;  #FITSLOG directory, used to read FITSLOG files from
my $out_dir    = undef;  #Intermediate directory where "stats" files are
my $prev_fits  = undef;  #Directory where consumed FITSLOG files are copied to
my $prev_ssc   = undef;  #Directory where consumed SSCLOG files are copied to
my $debug_dir  = undef;  #Directory where debugging logs can be stored (not used presently)

# Main subroutine
sub stats_log_parse {

	## --- Populating global variables from input parameters (really from config file, but no need to read it twice
	#
        ($ssc_dir, $fits_dir, $out_dir, $prev_ssc, $prev_fits, $debug_dir) = @_;	

	## --- Local variables
	#
	my @ssc_filenames; #array of all filenames in both directories
	my @fits_filenames; #array of all filenames in both directories
	

	## --- Creating missing directories, if any
	unless (-d $ssc_dir) {
		File::Path::make_path($ssc_dir) or die "Cannot create $ssc_dir: $!\n";
	}
        unless (-d $prev_fits) {
		File::Path::make_path($prev_fits) or die "Cannot create $prev_fits: $!\n";
	}
        unless (-d $prev_ssc) {
		File::Path::make_path($prev_ssc) or die "Cannot create $prev_ssc: $!\n";
	}

	## --- Directory input and filename output section --- ##

	#my $parse_debug = File::Spec->catfile($debug_dir, "parse_debug.txt");
	#open DEBUG, ">$parse_debug" or die "Cannot write to $parse_debug: $!";

	opendir(SSCDIR, "$ssc_dir");
	@ssc_filenames = grep(/\.SSCLOG$/, readdir SSCDIR); #grabs the SSCLOG filenames from $ssc_dir
	closedir SSCDIR;

	opendir(FITSDIR, "$fits_dir");
	@fits_filenames = grep(/\.FITSLOG$/, readdir FITSDIR); #grabs the FITSLOG files from $fits_dir
	closedir FITSDIR;


	## --- File parsing section --- ##
	# Code reads the name of each of the files, changes to specified directory and calls the appropriate subroutine to parse the file #
	
	print "Parsing SSCLOGS...\n";
	(print "\tNo SSCLOGs to parse\n") if (scalar(@ssc_filenames) == 0);
	(printf "\t%d SSCLOGs being parsed\n", scalar(@ssc_filenames)) if (scalar(@ssc_filenames) > 0);
	
	chdir "$ssc_dir" or die "Cannot change to $ssc_dir"; #changes to the directory indicated by $ssc_dir
	foreach my $filename (@ssc_filenames) {
		
		ssc_parse($filename); #calls the ssc_parse subroutine and adds returned values to %data_hash
		
	}
	print "\tSSCLOG parsing completed\n\n";
	
	print "Parsing FITSLOGS...\n";
	(print "\tNo FITSLOGs to parse\n") if (scalar(@fits_filenames) == 0);
	(printf "\t%d FITSLOGs being parsed\n", scalar(@fits_filenames)) if (scalar(@fits_filenames) > 0);
	
	chdir "$fits_dir" or die "Cannot change to $fits_dir"; #changes to the directory indicated by $fits_dir
	foreach my $filename (@fits_filenames) {
	
		fits_parse($filename); #calls the fits_parse subroutine and adds returned values to %data_hash
		
	}
	print "\tFITSLOG parsing completed\n\n";
	
        return 1;	
}

sub ssc_parse {

	my @keys = ();
	my @values = ();
	my $keys_length;
	my $vals_length;
	my $size_index = 0;
	my $prev_event_time = undef;
	my $current_year = 0;
	my $current_doy = 0;
	my $start_span = 0;
	my $end_span = 0;
	my $start_line = 0;
	my $temp_indicator = 0;
	open FILE, "@_" or die "could not find/open @_";
	my $name = "@_";
	
	my $file_name_out = '';
	my $file_name_temp = '';
	my $temp_out = '';
        my $creation_date = 0;	
        my $start_time_index = 0;
        my $type_index = 0;
        my $rate_index = 0;
	my $event_times = 0;
	
	while (<FILE>) {
	
		#print "LINE: $_ \n";

		#Checks for a line that says file creation time for use later, otherwise skips lines that begin with whitespace or semicolons
		unless (/file.*?creation.*?time/isg) {
		
			unless (eof) {
			
				next if (/^;/ or /^\s+/);
			
			}

		
		} else {
			
			#used to extract the file creation date from one of the header lines, for use in overwriting older iterations of stuff on the same day/time
			chomp(my @c_date = split /\s+/, $_);
			$c_date[1] =~ s/(0{3}?)|([-:\.]?)//sig; #gets rid of the 3 trailing millisecond digits, as well as any punctuation
			$creation_date = $c_date[1];
			print "\t Created on $creation_date\n";
			
			next;
		
		}
		
		#if conditional checks for the line that contains the headers

		if ($_ =~ /observer/i) {
			
			#splits the line and assigns the values to @keys
			#Finds the index for the start time "START", used later, and then removes that element from @keys (to avoid duplication)
			#finds the indices for "TYPE" and "RATE" so that an adjustment can be made if needed later on
			#Pushes "CREATION_DATE" into @keys, and then checks the length of @keys
			
			chomp(@keys = split /\s+/, $_);
			$start_time_index = first_index {$_ eq 'START'} @keys;
			splice @keys, $start_time_index, 1;
			$type_index = first_index {$_ eq 'TYPE'} @keys; #both type and rate index used to insert a ~ in the index where rate is potentially empty
			$rate_index = first_index {$_ eq 'RATE'} @keys;
			$size_index = first_index {$_ =~ /^SIZE/i} @keys;
			push @keys, 'CREATION_DATE';
			push @keys, 'SOURCE_FILE';
			$keys_length = scalar @keys;
			$keys[$size_index] = 'SIZE';
	
		} elsif ($_ =~ /^;Estimated_Image_Volume:/i) {
		
			print "Last line of file @_ reached\n";
		
		} else {
	
			#Splits the line and assigns values to @values, uses previously found $start_time_index to find the start time
			#removes punctuation and the millisecond zeroes, then assigns that to $event_times, and appends _SSC to it 
			#the start time element is then removed from @values to avoid duplication ($event time is used as the keys to the hash)
			#$creation_date is pushed to @values, then the length of the array is taken
			
			chomp(@values = split /\s+/, $_);
			$values[$start_time_index] =~ s/\.\d{3}?//si;
		
			$event_times = $values[$start_time_index];		##Gets the value of the date/time from the correct index of @values
			$event_times = "$event_times" . "_SSC";
			splice @values, $start_time_index, 1;		##removes the date/time from @values cause it is now part of $event_times
			push @values, $creation_date;
			push @values, "@_";
			$vals_length = scalar @values;
			
			if ($start_line == 0) {
			
				$start_span = $event_times;
				$start_span =~ s/\D//g;
				$start_line = 1;
			
			} else {
			
				$end_span = $event_times;
				$end_span =~ s/\D//g;
			
			}
			
			#If conditional checks for a match to the mode DARK, which if true it then checks if length of @values
			#is not equal to @keys, splices a ~ in where the rate index is as it's likely missing a value in some files
			
			if ($values[$type_index] =~ /(DARK)/i) {
			
				(splice @values, $rate_index, 0, '~') if ($vals_length != $keys_length);
			
			}
			
			#making the size match the format used in the FITSLOGS
			
			if ($values[$size_index] =~ /^ \d+ x .+ x .+ $/xgi) {
			
				$values[$size_index] =~ /^ (\d+) x (\d+) .* \[(\d+) , \d+ \] .? $/xi;
				my $size_fix = $1 + $3;
				$values[$size_index] = "$size_fix" . "x$2";	
			
			} else {
			
				$values[$size_index] =~ /^ (\d+) x (\d+) .* \[(\d+) , \d+ \] .? $/xi;
				$values[$size_index] = "$1x$2";
			
			}

			
			$vals_length = scalar @values;

			unless (eof) {
			
				next if ($vals_length != $keys_length);
			
			}
			
			
			EMPTY: {	#naked loop used to cut down on repetitive code
			
				$event_times =~ / ^(\d{4}) - (\d{3}) - \d{2} : \d{2} : \d{2} _ \w+$ /x; 
				if ($current_year == 0 && $current_doy == 0) {
				
					undef %data_hash;
					$current_year = $1;
					$current_doy = $2;

					$file_name_out  = File::Spec->catfile($out_dir, $current_year . "_" . $current_doy . "_SSC.txt");
					$file_name_temp = File::Spec->catfile($out_dir, $current_year . "_" . $current_doy . "_SSC_temp.txt");
					$temp_out = File::Spec->catfile($out_dir, $current_doy . "-temp.txt");
					
					if (-e "$file_name_out") { #if the file already exists, just append
						
						print "$file_name_out already exists, creating $temp_out\n";
						$temp_indicator = 0;
						open OUT, '>', $temp_out or die "cannot append to temp.txt from file @_: $!\n";
						$temp_indicator = 1;
						
					} else {  # if file did not exist, create it
					
						open OUT, '>', $file_name_out or die "cannot create $file_name_out: $!\n";
					
					}	
					
				} elsif ($current_year != 0 && $current_doy != 0) {
				
					if ($current_year != $1 || ($current_year == $1 && $current_doy != $2) || eof) {
					
						if (eof FILE) {
						
							unless ($_ =~ /^;/){
							
								my $i = 0;
								while ($i < $vals_length) { 
								
									$data_hash{$event_times}{$keys[$i]} = $values[$i];
									$i += 1;
								
								}	
								
								$data_hash{$event_times}{'NEXT_TIME'} = 'final'; #used so that the last one for the day doesn't have an empty place
								
								#Adds the event time for this to the previous one, for use with the stats module finding matches to fits files with diff titles
								if (defined $prev_event_time && $prev_event_time ne $event_times) {
									
									my $next_time = $event_times;
									$next_time =~ s/_SSC$//i;
									$data_hash{"$prev_event_time"}{'NEXT_TIME'} = $next_time;

								
								}	
								
								$prev_event_time = $event_times;
								
							}
						
						}
					
						foreach my $start_times (sort keys %data_hash) {
			
							print OUT "$start_times\n"; #start times used as a main category to differentiate them
			
							foreach my $subheaders (sort keys %{$data_hash{$start_times}}) { #subheaders are categories
				
								print OUT ";$subheaders=>$data_hash{$start_times}{$subheaders}\n"; #prints categories and associated values

							}
		
							print OUT "\n";
		
						}
					
						close OUT;
						
						if ($temp_indicator == 1) {

							my $old_time = 0;
							my $old_line;
							my $old_year;
							my $old_doy;
							
							file_update( $start_span, $end_span, $file_name_out, $file_name_temp, $temp_out, $name );

							
							$temp_indicator = 0;
							
						}
						

						undef %data_hash;
						$current_year = 0;
						$current_doy = 0;
						undef $prev_event_time;
						(redo EMPTY) if (not eof (FILE)); #restarts the naked loop, which begins just before the previous $event_times matching
				
					
					}
				
				} 
				
			}
			
			#Initializes $i to 0, and runs the while loop until $i is equal to the length of vals_length (index of last element + 1)
			#the while loop uses $event_times as the key, where the value is a hash comprising of the previously determined keys and values
			#I'm still not entirely sure how array/hash references work, and since I haven't been able to get one to work yet, this is a 
			#much cleaner solution than just hardcoding in the hash, and allows the size of them to change more fluidly
			
			my $i = 0;
			while ($i < $vals_length) { 
			
				($values[$i] =~ s/^~//) if ((($keys[$i] eq 'ROLL')||($keys[$i] eq 'RA')||($keys[$i] eq 'DEC'))&&($values[$i] =~ /^~/));
				if ($keys[$i] eq 'ROLL') {
				
					($values[$i] = $values[$i] - 360) if ($values[$i] > 180);
					($values[$i] = $values[$i] + 360) if ($values[$i] < -180);
					redo if (($values[$i] < -180)||($values[$i] > 180));
				
				}
			
				$data_hash{$event_times}{$keys[$i]} = $values[$i];
				$i += 1;
			
			}	
			
			$data_hash{$event_times}{'NEXT_TIME'} = 'final'; #used so that the last one for the day doesn't have an empty place
			
			#Adds the event time for this to the previous one, for use with the stats module finding matches to fits files with diff titles
			if (defined $prev_event_time && $prev_event_time ne $event_times) {
				
				my $next_time = $event_times;
				$next_time =~ s/_SSC$//i;
				$data_hash{"$prev_event_time"}{'NEXT_TIME'} = $next_time;

			}	
			
			$prev_event_time = $event_times;
			
		}
	
	}
	
	close FILE;
	undef %data_hash;
	move(File::Spec->catfile($ssc_dir, @_), File::Spec->catfile($prev_ssc,@_));
	#unlink @_;
	
}



sub fits_parse {

	open FILE, "@_" or die "could not find/open @_";
	
	my @keys = ();
	my @values = ();
	my $keys_length = 0; 
	my $vals_length = 0;
	my $mode_index = 0;
	my $dur_index = 0;
	my $offset_index;
	my $stddev_index;
	my $roll_index;
	my $ra_vel_index;
	my $dec_vel_index;
	my $roll_vel_index;
	
        my $file_name_index;
        my $start_time_index;
        my $event_times; 

	while (<FILE>) {
		
		if ($_ =~ /observer/i) {
		
			#if line contains one of the keys, reads line, removes leading whitespace, splits line by whitespace
			#locates the indices corresponding to key values, stores them for later use, then removes the element
			#containing the filename, which has the date/time as part of it, as it's not useful to have twice
			
			$_ =~ s/^\s+//;
			chomp(@keys = split /\s+/isgx, $_);
			$start_time_index = first_index {$_ eq 'Exposure_start'} @keys;
			$file_name_index = first_index {$_ eq 'File_name'} @keys;
			splice @keys, $file_name_index, 1;
			push @keys, 'SOURCE_FILE';
			$keys_length = scalar @keys;
	
			$mode_index = first_index {$_ eq 'Pointing'} @keys; #using these next 4 lines to change the names of some keys, so they match with SSC
			$dur_index = first_index {$_ eq 'Expos.'} @keys;
			$keys[$mode_index] = 'MODE';
			$keys[$dur_index] = 'DUR';
	
	
		} else {
		
			#Skips all lines that don't start with the FITS image file name format
			#lines that do start with that are split into values
			#the element corresponding to the "Exposure_start" index in @keys is fixed to remove unnecessary values/punctuation
			#the element corresponding to the "File_name" index in @keys is assigned to variable $event_times, which is then cleaned update
			#The element corresponding to "File_name" is removed to conform with @keys
			
			next unless ($_ =~ /^ NEOS_SCI_ /isgx); 

			if ($start_time_index eq 'TBD') {
				next;
			}
			if ($mode_index eq 'XX-N/A') {
				next;
			}

			
			chomp(@values = split /\s+/s, $_);		
			$values[$start_time_index] =~ s/(^\d{4}-\d{2}-\d{2}T)|( [-:\.] ) //sigx;
			
			$event_times = $values[$file_name_index];		
			$event_times =~ s/NEOS_SCI_//igs;
			$event_times =~ s/^ (\d{4}) (\d{3}) (\d{2}) (\d{2}) (\d{2}) $/$1-$2-$3:$4:$5/x;
			$event_times = "$event_times" . "_FITS";
			splice @values, $file_name_index, 1;
			push @values, "@_";
			$vals_length = scalar @values;

			
			DIFFLENGTH: {
			
				my $offset_index;
				my $stddev_index;
				my $roll_index;
				my $ra_vel_index;
				my $dec_vel_index;
				my $roll_vel_index;
				my $index_count = 0;
			
				if ($vals_length != $keys_length) {
						
						my $obj_index = first_index {$_ eq 'OBJECT'} @keys;
						my $obsv_index = first_index {$_ eq 'Observer'} @keys;
						(print DEBUG "Object $obj_index\nObserver $obsv_index\n") if ($event_times =~ /2019-021-19:00:12/);
						
						
						(splice @values, $obj_index, 0, '~') if (($values[$obj_index] =~ /VIEW|DARK/i) && ($values[$obsv_index] =~ /^ \d+ x \d+ $/x));
						(splice @values, $obsv_index, 0, '~') if ($values[$obsv_index] =~ /VIEW|DARK/i);
						
						$offset_index = first_index {$_ eq 'Offset'} @keys;
						$stddev_index = first_index {$_ eq 'StdDev'} @keys;
						$roll_index = first_index {$_ eq 'ROLL'} @keys;
						$ra_vel_index = first_index {$_ eq 'RA_VEL'} @keys;
						$dec_vel_index = first_index {$_ eq 'DEC_VEL'} @keys;
						$roll_vel_index = first_index {$_ eq 'ROL_VEL'} @keys;
						
						
						$index_count = 0;
						foreach (@values) {
						
							if (($_ =~ /^ \d+ \. \d+ \. \d $/x) && ($index_count == $offset_index)) {
							
								$_ =~ /^(\d+ \. \d) (\d+ \. \d) $/x;
								print DEBUG "offset and stddev $1|$2 from file @_\n";
								$values[$index_count] = $1;
								splice @values, $stddev_index, 0, $2;
							
							} elsif (($_ =~ /^ -? \d+ \. \d{3} -? \d+ \. \d{3} $/x) && ($index_count == $roll_index)) {
							
								$_ =~ /^ (-? \d+ \. \d{3}) (-? \d+ \. \d{3}) $/x;
								print DEBUG "roll and ra $1|$2 from file @_\n";
								$values[$roll_index] = $1;
								splice @values, $ra_vel_index, 0, $2;
							
							} elsif (($_ =~ /^ -? \d+ \. \d{3} -? \d+ \. \d{3} $/x) && ($index_count == $ra_vel_index)) {
							
								$_ =~ /^ (-? \d+ \. \d{3}) (-? \d+ \. \d{3}) $/x;
								print DEBUG "ra and dec $1|$2 from file @_\n";
								$values[$ra_vel_index] = $1;
								splice @values, $dec_vel_index, 0, $2;
							
							} elsif (($_ =~ /^ -? \d+ \. \d{3} -? \d+ \. \d{3} $/x) && ($index_count == $dec_vel_index)) {
							
								$_ =~ /^ (-? \d+ \. \d{3}) (-? \d+ \. \d{3}) $/x;
								print DEBUG "dec and roll $1|$2 from file @_\n";
								$values[$dec_vel_index] = $1;
								splice @values, $roll_vel_index, 0, $2;
							
							}
						
							
							$index_count += 1;
						
						}
						
						$vals_length = scalar @values;
						(redo DIFFLENGTH) if ($vals_length < $keys_length);
						
						
						print DEBUG "Time: $event_times from file @_\n";
					
				
				}
				
			}
			#Initializes $i to 0, and runs the while loop until $i is equal to the length of vals_length (index of last element + 1)
			#the while loop uses $event_times as the key, where the value is a hash comprising of the previously determined keys and values
			#I'm still not entirely sure how array/hash references work, and since I haven't been able to get one to work yet, this is a 
			#much cleaner solution than just hardcoding in the hash, and allows the size of them to change more fluidly
			
			$values[$mode_index] =~ s/^ (.{3}) (\w+) $/$2/ix; #removes the stuff at the beginning of the mode name
			$values[$dur_index] =~ s/s$//i;
			
			my $i = 0; #used in the following loop
			while ($i < $vals_length) { 
				
				if ($keys[$i] eq 'RA') {
				
					my ($hour, $min, $sec) = split /:/, $values[$i];
					$values[$i] = (($sec/3600) + ($min/60) + $hour);
					$values[$i] = sprintf("%.3f", $values[$i]);
				
				} elsif ($keys[$i] eq 'DEC') {
				
					my ($deg, $min, $sec) = split /:/, $values[$i];
					($values[$i] = ($sec/3600) + ($min/60) + ($deg) ) if ($deg >= 0);
					($values[$i] = ($deg) - ($min/60) - ($sec/3600) ) if ($deg < 0);
					$values[$i] = sprintf("%.3f", $values[$i]);
				
				} elsif ($keys[$i] eq 'ROLL') {
				
					($values[$i] = $values[$i] - 360) if ($values[$i] > 180);
					($values[$i] = $values[$i] + 360) if ($values[$i] < -180);
					redo if (($values[$i] < -180)||($values[$i] > 180));
				
				}
				
				$data_hash{$event_times}{$keys[$i]} = $values[$i];
				$i += 1;
			
			}
			
		}
		
	}
	
	close FILE;
	move(File::Spec->catfile($fits_dir, @_), File::Spec->catfile($prev_fits,@_));
	#unlink @_;
	
	
	my @keys_array = sort keys %data_hash;
	
	my $counter = 1;
	foreach (@keys_array) {
	
		$counter += 1;
	
	}
	$keys_array[0] =~ /^(\d{4})-(\d{3})-\d{2}:\d{2}:\d{2}_\w+$/;
	my ($year, $doy) = ($1,$2); #priming the first doy/year combo to start a new file

	#the outer foreach loop deals with each of the $event_times keys

	foreach my $start_times (sort keys %data_hash) {
		
			
			#checks $start_times for values, assigns them to memory, then checks if they are not the same as one or all of the previous values
			#for $year, $doy, and $file_type, or if the filehandle OUT is not open. If any of these are true, the code block runs and opens a
			#filehandle that exports to a file that matches the day of year, year, or filetype (in this case, )
			
			$start_times =~ / ^(\d{4}) - (\d{3}) - \d{2} : \d{2} : \d{2} _ (\w+)$ /x; 
				
			## this if block used to discern if it needs to close the OUT filehandle, change $file_name_out, and re-open the filehandle to a different output filename ##
			if ($doy != $2 || $year != $1 || not(OUT->opened())) { 
			
				(close OUT) if (OUT->opened()); #closes the filehandle OUT if it is still open
				
				($year, $doy) = ($1, $2); #changes the values to the new doy/year

				#make $file_name_out the new combination of year and doy for use in the new file
				my $file_name_out = File::Spec->catfile($out_dir, $year . "_" . $doy . "_FITS.txt"); 
				
					
				if (-f $file_name_out) { 
					
					#checks if the filename is already in use, and if it is, the new filehandle is an append one as opposed to overwriting it
					open OUT, ">>$file_name_out" or die "cannot append to $file_name_out: $!\n";
					
				} else {
					
					#if the filename is not already in use, just opens a new filehandle by that name
					open OUT, ">$file_name_out" or die "cannot create $file_name_out: $!\n";
					
				}
					
			}
			
			print OUT "$start_times\n"; #start times used as a main category to differentiate them
			
			foreach my $subheaders (sort keys %{$data_hash{$start_times}}) { #subheaders are categories
				
				print OUT ";$subheaders=>$data_hash{$start_times}{$subheaders}\n"; #prints categories and associated values

			}
		
		print OUT "\n";
		
	}
	
	close OUT;
	undef %data_hash;
	
	

}

sub file_update {

	my $old_time = 0;
	my $old_line;
	my $start_span = shift;
	my $end_span = shift;
	my $file_name_out = shift;
	my $file_name_temp = shift;
	my $temp_out = shift;
	my $name = shift;
	
	open (TEMP, ">", $file_name_temp) or die "cannot create $file_name_temp for file $name: $!";
	open (OLD_IN, "<", $file_name_out) or die "cannot open $file_name_temp for file $name: $!";
	open (NEW_IN, "<", $temp_out) or die "cannot open $temp_out for file $name: $!";
	
	print DEBUG "Reading from file $name\n\tFile $file_name_temp:\n";
	
	while (<OLD_IN>) {
	
		if ($_ =~ /^\d/) {
		
			$old_time = $_;
			$old_time =~ s/\D//g;
		
		}
		
		if ($old_time < $start_span) {
		
			print TEMP "$_";
		
		} elsif ($old_time >= $start_span) {
		
			last;
		
		}
	
	
	}
	#print DEBUG "\tNew info starts at $start_span\n";
	
	while (<NEW_IN>) {
	
		print TEMP "$_";
	
	}
	
	#print DEBUG "\tNew info ends at $end_span\n";
	while (<OLD_IN>) {
	
		if ($_ =~ /^\d/) {
		
			$old_time = $_;
			$old_time =~ s/\D//g;
		
		}
		
		if ($old_time > $end_span) {
		
			print TEMP "$_";
		
		} elsif ($old_time <= $end_span) {
		
			next;
		
		}
	
	}
	
	close OLD_IN; close NEW_IN; close TEMP;
	
	my $unlinked = unlink $file_name_out, $temp_out;
	rename $file_name_temp, $file_name_out;

}

1;
