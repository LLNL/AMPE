#!/usr/bin/perl

#--------------------------------------------------------------------------
# Define configuration constants
#--------------------------------------------------------------------------

# @REFERERS contains the host names and IP addresses of the domains which
# are allowed to use this copy of FormHandler to parse their forms.

@REFERERS = ('localhost:8090');
#@REFERERS = ('www.llnl.gov');

# Host name of the web server.

$WEB_SERVER = 'localhost:8090';
#$WEB_SERVER = 'www.llnl.gov';

$sunmapURL = '../sunmap/sunmap.html';

# %CONFIG defines which form fields should be considered configuration
# fields rather than standard data fields. Each of the default variables
# defined in the array below have special meaning to FormHandler and are
# usually set using hidden fields. Default values used in the array will
# be overridden by form fields with the same name. Any variable that should
# be considered a configuration variable must be defined in this array.

%CONFIG = ('realname',                '', 
           'city',                    '',
           'country',                 '',
           'lat',                     '',
           'lng',                     '',
           'photo',                   '',
           'new_id',                  '',
           'old_id',                  '',
           'action',                  '');

# Directory where all lock files will be placed.

$LOCK_DIR = 'lock/';

# Max number of seconds the lock script will wait before overiding the lock file.

$MAX_WAIT = 5;

# Directory in which all of your required files are placed.

$REQUIRE_DIR = './';

# Push the $REQUIRE_DIR onto the @INC array for include file directories.

push(@INC, $REQUIRE_DIR) if $REQUIRE_DIR;

#--------------------------------------------------------------------------
# Check that the form is coming from a web site that is included in the
# @REFERERS array.  If not, sent out an error. Otherwise, continue.
#--------------------------------------------------------------------------

# Set the flag to 0 so the referer will fail by default.

$check_referer = "0";

# Get the hostname out of the referring URL.  If there is no 
# referring URL, set flag to 1, since some browsers don't pass the 
# HTTP_REFERER environment variable.

if ($ENV{'HTTP_REFERER'} =~ m#http://([^/]+)#) {
    $hostname = $1;
}
else {
    $check_referer = 1;
}

# If there is a hostname found, check it against the hostnames in the 
# @REFERERS array until a match is found.  Set flag to 1 if match is 
# found.

if ($hostname) {
    foreach $referer (@REFERERS) {
        if ($hostname =~ /$referer/i) {
            $check_referer = 1;
            last;
        }
    }
}

# If flag is not set to 1, throw an error to the &error subroutine.

if (!$check_referer) {
    &error("$ENV{'HTTP_REFERER'} is not allowed access to this program.");
}


#--------------------------------------------------------------------------
# Parse the form contents and put configuration fields into %CONFIG.
#--------------------------------------------------------------------------

if (!&parse_form) {
    &error($Error_Message);
}

#--------------------------------------------------------------------------
# Update XML file depending on action type
#--------------------------------------------------------------------------

$XML_file = 'sunmarkers.xml';
$NEW_file = 'sunmarkers_tmp.xml';

# Lock XML file

if (&lock($XML_file, $LOCK_DIR, $MAX_WAIT)) {
    &error($Error_Message);
}


if ($CONFIG{'action'} eq 'add') {

    open(LOG_IN, "< $XML_file");
    open(LOG_OUT, "> $NEW_file");

    while ($line = <LOG_IN>) {
        $id_start = index($line, "</markers>");
        if($id_start != -1) {
            print LOG_OUT <<EOL;
 <marker lat="$CONFIG{'lat'}" lng="$CONFIG{'lng'}" name="$CONFIG{'realname'}" city="$CONFIG{'city'}" country="$CONFIG{'country'}" photo="$CONFIG{'photo'}" id="$CONFIG{'new_id'}" />
</markers>
EOL
        } else {
            print LOG_OUT $line;
        }
    }

    close(LOG_IN);
    close(LOG_OUT);    
    rename($NEW_file, $XML_file);

}

if ( ($CONFIG{'action'} eq 'replace') || ($CONFIG{'action'} eq 'delete') ) {

    open(LOG_IN, "< $XML_file");
    open(LOG_OUT, "> $NEW_file");

    $MSG_FILE = 'msg.txt';
    open(DUMP, "> $MSG_FILE");

    print DUMP "action is $CONFIG{'action'}";

    while ($line = <LOG_IN>) {

        $new_line = $line;

        $id_start = index($line, "id=");
        if($id_start != -1) {
           $id_start = $id_start + 4;
           $id_end = rindex($line, "\"");
           $length = $id_end-$id_start; 
           $this_id = substr($line, $id_start, $length);            
           if($this_id eq "$CONFIG{'old_id'}") {
               print DUMP "this is the line...  $CONFIG{'action'} \n";
               if($CONFIG{'action'} eq 'replace') {
                   $new_line = qq( <marker lat="$CONFIG{'lat'}" lng="$CONFIG{'lng'}" name="$CONFIG{'realname'}" city="$CONFIG{'city'}" country="$CONFIG{'country'}" photo="$CONFIG{'photo'}" id="$CONFIG{'new_id'}" />\n);
               }
               if($CONFIG{'action'} eq 'delete') {
                   $new_line = "";
               }
           }
        }

        print LOG_OUT $new_line;

        print DUMP $new_line;

    }


    close(DUMP);


    close(LOG_IN);
    close(LOG_OUT);    
    rename($NEW_file, $XML_file);
}


# Unlock the log file
    
&unlock($XML_file, $LOCK_DIR);

#--------------------------------------------------------------------------
# Redirect output
#--------------------------------------------------------------------------

print "Location: $sunmapURL \n\n";

exit;

###########################################################################
#
#  LOCAL SUBROUTINES
#
###########################################################################

#--------------------------------------------------------------------------
# Print out the Error Messages and exit.
#--------------------------------------------------------------------------

sub error {

    # Localize the error variable.
    
    local($error) = $_[0];
    print "Content-type: text/html\n\n";
    
    # If the error is because of missing_required_fields, use the
    # error_template that was specified or print out a generic response.
    
    if ($error eq 'missing_required_fields') {
    
        # Prepare the error_fields config field so that users can use it 
        # in their templates.
        
        $CONFIG{'error_fields'} = "<UL>\n";
        foreach $missing_required_field (@missing_required_fields) {
            if ($ALT_NAME{$missing_required_field}) {
                $CONFIG{'error_fields'} .= "<LI>$ALT_NAME{$missing_required_field}</LI>\n";
            }
            else {
                $CONFIG{'error_fields'} .= "<LI>$missing_required_field</li>\n";
            }
        }
        $CONFIG{'error_fields'} .= "</UL>";
        
        # Print out formatted template to user.
        
        if (!&parse_template($CONFIG{'error_template'}, *STDOUT)) {
            $error = "Can't open $CONFIG{'error_template'} ($!).";
        }
        else { exit }
        
    }
    
    # For any other errors, just print a title and heading which supplies 
    # the error.
    
    print <<HTML_END;
<HTML>
   <HEAD>
      <TITLE>$error</TITLE>
   </HEAD>
   <BODY BGCOLOR=#FFFFFF TEXT=#000000>
      <CENTER>
      <H4>$error</H4>
      </CENTER>
   </BODY>
</HTML>
HTML_END
    exit;
}

#--------------------------------------------------------------------------
# parse_form
#--------------------------------------------------------------------------
#
# Function:      Takes form field data from a POST or GET request and      #
#                converts it to name/value pairs in the %FORM array or,    #
#                if a corresponding entry in the %CONFIG array is defined, #
#                in %CONFIG.                                               #
#                                                                          #
# Usage:         &parse_form;                                              #
#                                                                          #
# Variables:     None                                                      #
#                                                                          #
# Returns:       0 if invalid request method                               #
#                1 if successful                                           #
#                                                                          #
# Uses Globals:  Sets %CONFIG with name/value pairs if corresponding entry #
#                  in %CONFIG is defined                                   #
#                Otherwise sets entries in %FORM                           # 
#                $Error_Message for descriptive error messages
#

sub parse_form {
    local($name, $value, $pair, $buffer, @pairs);
    
    # Check for request method and handle appropriately
    
    if ($ENV{'REQUEST_METHOD'} eq 'GET') {
        @pairs = split(/&/, $ENV{'QUERY_STRING'});
    }
    elsif ($ENV{'REQUEST_METHOD'} eq 'POST') {
        read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'});
        @pairs = split(/&/, $buffer);
    }
    else {
        $Error_Message = "Bad request method ($ENV{'REQUEST_METHOD'}).  Use POST or GET";
       return(0);
    }

    # Convert the data to its original format
    
    foreach $pair (@pairs) {
        ($name, $value) = split(/=/, $pair);

        $name =~ tr/+/ /;
        $name =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
        $name =~ s/\n//g;
        $value =~ tr/+/ /;
        $value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
        $value =~ s/\n//g;

        # If they try to include server side includes, erase them, so they
        # arent a security risk if the HTML gets returned.  Another
        # security hole plugged up.
        
        $value =~ s/<!--(.|\n)*-->//g;

        # Store name/value pair in %CONFIG if the corresponding entry is
        # defined
        
        if ($CONFIG{$name}) {
            $CONFIG{$name} .= ",$value";
        }
        elsif (defined($CONFIG{$name})) {
            $CONFIG{$name} = $value;
        }
        
        # Otherwise store in %FORM
        
        elsif ($FORM{$name}) {
            $FORM{$name} .= ",$value";
        }
        else {
            $FORM{$name} = $value;
        }
    }
    return(1);
}

#--------------------------------------------------------------------------
# lock
#--------------------------------------------------------------------------
#
# Function:      Creates an exclusive lock for a file. The lock will       #
#                only work if other programs accessing the file are also   #
#                using this subroutine.                                    #
#                                                                          #
# Usage:         &lock($filename, $LOCK_DIR[, $MAX_WAIT]);                 #
#                                                                          #
# Variables:     $filename --   Name of file being locked.                 #
#                               Example "filename.html"                    #
#                $LOCK_DIR --   Path of directory to store lock files      #
#                               Should be "/tmp/" on UNIX sytems           #
#                               Example "/home/lockdir/"                   #
#                $MAX_WAIT --   Maximum seconds to wait if the file is     #
#                               already locked                             #
#                                                                          #
# Returns:       0 if successful                                           #
#                1 if $LOCK_DIR/$filename.tmp could not be created         #
#                2 if $filename is currently in use                        #
#                3 if lock file could not be created or opened             #
#                                                                          #
# Uses Globals:  $Error_Message for descriptive error messages             #
#                $NAME_LEN for maximum filename length                     #
#                                                                          #
# Files Created: Creates $LOCK_DIR/$filename.tmp                           #
#                Creates $LOCK_DIR/$filename.lok (exists only while file   #
#                  is locked)                                              #
#

sub lock {
    
    # Initialize variables
    
    local($filename, $LOCK_DIR, $MAX_WAIT) = @_; 
    local($wait, $lock_pid);
    local($temp_file) = "$LOCK_DIR$$.tmp";
    $Error_Message = '';
    
    local($lock_file) = $filename;
    $lock_file =~ tr/\/\\:.//d;         # Remove file separators/periods
    if ($NAME_LEN && ($NAME_LEN < length($lock_file))) {
        $lock_file = substr($lock_file, -$NAME_LEN);
    }
    $lock_file = "$LOCK_DIR$lock_file.lok";
    
    # Create temp file with PID
    
    if (!open(TEMP, ">$temp_file")) {
        $Error_Message = "Could not create $temp_file ($!).";
        return(1);
    }           
    print TEMP $$;
    close(TEMP);
    
    # Test for lock file
    
    if (-e $lock_file) {

        # Wait for unlock if lock file exists
        
        for ($wait = $MAX_WAIT; $wait; --$wait) {
            sleep(1);
            last unless -e $lock_file;
        }
    }
    
    # Check to see if there's still a valid lock
    
    if ((-e $lock_file) && (-M $lock_file < 0)) {
        
        # The file is still locked but has been modified since we started
                    
        unlink($temp_file);
        $Error_Message = "The file \"$filename\" is currently in use. Please try again later.";
        return(2);
    }
    else {

        # There is either no lock or the lock has expired
        
        if (!rename($temp_file, $lock_file)) { 

            # Problem creating the lock file
        
            unlink($temp_file);
            $Error_Message = "Could not lock file \"$filename\" ($!).";
            return(3);
        }
        
        # Check to make sure the lock is ours

        if (!open(LOCK, "<$lock_file")) {
            $Error_Message = "Could not verify lock for file \"$filename\" ($!).";
            return(3);
        }
        $lock_pid = <LOCK>;
        close(LOCK);        
        if ($lock_pid ne $$) { 
            $Error_Message = "The file \"$filename\" is currently in use. Please try again later.";
            return(2);
        }
        else { return(0) }
    }
}


#--------------------------------------------------------------------------
# unlock
#--------------------------------------------------------------------------
#
# Function:      Unlocks a file that has been locked using lock().         #
#                                                                          #
# Usage:         &unlock($filename, $LOCK_DIR);                            #
#                                                                          #
# Variables:     $filename --   Name of file being locked.                 #
#                               Example "filename.html"                    #
#                $LOCK_DIR --   Path of directory to store lock files      #
#                               Should be "/tmp/" on UNIX sytems           #
#                               Example "/home/lockdir/"                   #
#                                                                          #
# Returns:       0 if successful                                           #
#                1 if the lock file could not be deleted                   #
#                                                                          #
# Uses Globals:  $Error_Message for descriptive error messages             #
#                $NAME_LEN for maximum filename length                     #
#                                                                          #
# Files Created: Removes $LOCK_DIR/$filename.lok                           #
#

sub unlock {
    
    # Initialize variables
    
    local($filename, $LOCK_DIR) = @_;
    local($lock_file) = $filename;
    $Error_Message = '';
    
    $lock_file =~ tr/\/\\:.//d;         # Remove file separators/periods
    if ($NAME_LEN < length($lock_file)) {
        $lock_file = substr($lock_file, -$NAME_LEN);
    }
    $lock_file = "$LOCK_DIR$lock_file.lok";
    
    # Check to make sure the lock is ours

    if (!open(LOCK, "<$lock_file")) {
        $Error_Message = "Could not access the lock file for \"$filename\" ($!).";
        return(1);
    }
    $lock_pid = <LOCK>;
    close(LOCK);        
    if ($lock_pid ne $$) { 
        $Error_Message = "The file \"$filename\" is locked by another process.";
       return(2);
    }

    # Release the lock by unlinking the lock file
    
    if (!unlink($lock_file)) {
        $Error_Message = "Could not unlock file \"$filename\" ($!).";
        return(3);
    }
    return(0);
}

