#!/usr/bin/perl

#--------------------------------------------------------------------------
# Define configuration constants
#--------------------------------------------------------------------------


# %CONFIG defines which form fields should be considered configuration
# fields rather than standard data fields. Each of the default variables
# defined in the array below have special meaning to FormHandler and are
# usually set using hidden fields. Default values used in the array will
# be overridden by form fields with the same name. Any variable that should
# be considered a configuration variable must be defined in this array.


%CONFIG = ('password', '',
           'action',   '');

#--------------------------------------------------------------------------
# Parse the form contents and put configuration fields into %CONFIG.
#--------------------------------------------------------------------------

&parse_form;

#--------------------------------------------------------------------------
# Check password
#--------------------------------------------------------------------------

if ($FORM{'CheckPassword'}) {
    
    $pass = "porsche";
    $validurl = "../sunmap/editmap.html";
    $invalidurl = "../sunmap/sunmap.html";
    
    if ($CONFIG{'password'} eq $pass) {
        
        # Password OK
        
        open(XML_IN, "< sunmarkers.xml"); 
        @fileContent= <XML_IN>;
        close(XML_IN); 
        
        print <<HTML_END;
<HTML>
<HEAD><TITLE>Edit sunmarkers.xml</TITLE></HEAD>
<BODY BGCOLOR=#FFFFFF TEXT=#000000>
$CONFIG{'action'}
  <FORM NAME="EditXML" action="passwords.cgi">
  <textarea rows="20" cols="140" style="overflow:scroll; line-height:18pt;" wrap="off">@fileContent</textarea>
  <br><br>
  <INPUT TYPE="submit" NAME="submit" VALUE="Submit" CLASS="button" onClick="submitAdd(); return true;"></INPUT>
  <INPUT TYPE="reset" VALUE="Cancel" CLASS="button"></INPUT>
  <INPUT TYPE="hidden" ID="action" NAME="action" VALUE="Edit"></INPUT>
  </FORM>
</BODY>
</HTML>
HTML_END

   }
    else {
        # Password not OK - goto the invalid password URL
        print "Location: $invalidurl\n\n";
    }

}

#--------------------------------------------------------------------------
# Save modified file
#--------------------------------------------------------------------------

if ($FORM{'EditXML'}) {
    print "Location: www.cnn.com";
}


exit;

#--------------------------------------------------------------------------
# parse_form
#--------------------------------------------------------------------------
#
# Function:      Takes form field data from a POST or GET request and
#                converts it to name/value pairs in the %FORM array or,
#                if a corresponding entry in the %CONFIG array is defined,
#                in %CONFIG.
#
# Usage:         &parse_form;
#
# Variables:     None
#
# Returns:       0 if invalid request method
#                1 if successful
#
# Uses Globals:  Sets %CONFIG with name/value pairs if corresponding entry
#                  in %CONFIG is defined
#                Otherwise sets entries in %FORM
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

