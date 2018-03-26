#!/usr/bin/perl
use strict;
use CGI::Cookie;
use CGI qw(:standard);

# Author: John Gardner
# Date:   5 January 2005

# ******************************************************************************

#  THIS SECTION IS USER MODIFIABLE
#  ===============================

# Set up the password strings and associated URLs. Note that the elements of all
# but the last of this hash list are separated by commas.

my %urlList = ("xyzzy" => "http://www.braemoor.co.uk/software/valid.html",
               "abcdef" => "http://www.braemoor.co.uk/software/valid.html",
               "123456" => "http://www.braemoor.co.uk/software/valid.html"
              );

# Set up the invalid URL - you also need to set this up in the JavaScript  
   
my $invalidurl = "http://www.braemoor.co.uk/software/invalid.html";

# ******************************************************************************

# Pick up the name of the password from the query string
my $password = param ('password');
 
my $q = new CGI;
if (exists($urlList{$password})) {

  # Password OK - goto the required URL with a cookie to show that the password
  # has been given. 
  my $validurl = $urlList{$password};
  my $cookie = $q->cookie(-name => "validpassword", -value => "0", -path => "/");
  print $q->redirect (-url =>$validurl, -cookie => $cookie);
}
else {

  # Password not OK - goto the invalid password URL
  print $q->redirect (-url =>$invalidurl);
}
