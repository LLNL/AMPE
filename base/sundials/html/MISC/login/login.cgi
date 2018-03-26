#!/usr/bin/perl
# (Put the address to the location of PERL on your system.  Find
#  it with 'which perl')
######################################################################
#                                                                    #
#           Login      Version 1.5.7                                 #
#                                                                    #
#               By Andrew Cantino (cantino@tectonicdesigns.com)      #
#               (c) 1998 by Andrew Cantino.    Do not edit without   #
#                permission.  Thank you.  If you like this program,  #
#                please e-mail me!                                   #
#                                                                    #
#                Please send me bug reports as well.                 #
#								     #
#		 Get updates at:				     #
#			http://tectonicdesigns.com/freecgi           #
#                                                                    #
#                                                                    #
#     What it does:                                                  #
#       Login is a simple password protection script that I have     #
#       written on request by many different people.  This script    #
#       is only intended to be a simple password protection system.  #
#       This script is not intended to be used to protect anything   #
#       that needs high level security.  I take no responsibility    #
#       for stolen or lost data because of this script.              #
#       Generally, just use this at your own risk.  Now enjoy!       #
#                                                                    #
#     How to set up:                                                 #
#      1) Put login.cgi where it can be run from the web.  This may  #
#         just be in your directory or it my be in your server's     #
#         cgi-bin area.                                              #
#      2) chmod 755 login.cgi     (make login.cgi executable)        #
#      3) configure login.cgi     (edit it)                          #
#      4) chmod 777 log           (make the log file editable)       #
#      5) If you use this script, please tell me your URL and I will #
#         make a link to you from my pages!                          #
#                                                                    #
######################################################################

# Address to this script.
$thisscript = "http://cantino.frognet.net/login/login.cgi";

# Page that the users should be sent to when they put in the correct
# password.
$forward = "http://cantino.frognet.net/design/";

#Name of the page that you are logging into.
$pagename = "Tectonic Designs";

# Password required to login.
$pass = "login";

#Send mail to YOU when someone logs in?
$send_mail = 0;

# UNIX path to the mail program on your system.
# This program must be able to take the -s option to send mail.  Such as:
# elm, Mail, etc.  If you run into problems, turn mail sending off.
$mail = "/usr/bin/Mail";

#Email address to send mail to (your personal e-mail address.)
#You MUST put a backslash (\) in front of the 'at' (@) sign in the e-mail
# address.
$to_email = "cantino\@tectonicdesigns.com";

# Do you wish to log logins?  (1/0)
# LOG file is NOT auto cleared.  You will have to edit it by hand.  If you
# delete it, remember to chmod the new file 777 when you re-make it.
$log = 1;

#Ask for name?  (Will be logged.)
$name = 1;

#Ask for an e-mail address?  (Will be logged.)
$email = 1;

# What is the address to the log file?  (Remember to create the file and
#                                         to chmod it 777)
$log_file = "/home/cantino/html/login/log";

# Path to your system's date program for logging.
$date_prog = "/bin/date";

# Settings for login page colors.
$text = "black";
$link = "blue";
$vlink = "blue";
$bgcolor = "white";

#  DO NOT EDIT PAST THIS POINT WITHOUT PERMISSION.

$date = `$date_prog '+%D %H:%M:%S'`;
%in = &getcgi;
print "Content-type: text/html\n\n";
if ($in{'pass'} ne $pass) {
	&print_login;
	exit;
}

if ($log == 1) {
	if (($name == 1) && ($in{'name'} eq "")) {
		&error("name");
	}
	if (($email == 1) && ($in{'email'} eq "")) {
		&error("E-Mail address");
	}
	unless ($in{'email'} =~ m/\@/i) {
		if ($email == 1) {
			&error("correct e-mail address");
		}
	}

	if (! -e "$log_file") {
		open(FILE, ">$log_file");
		print FILE "File START $date\n";
		close(FILE);
	}
	open(FILE, ">>$log_file");
	print FILE "Login: $ENV{'REMOTE_ADDR'} (with $ENV{'HTTP_USER_AGENT'}) $date";
	if ($name == 1) {
		print FILE "  Name: $in{'name'}\n";
	}
	if ($email == 1) {
		print FILE "  E-mail: $in{'email'}\n";
	}
	close(FILE);
}

if ($send_mail == 1) {
        if (-x $mail) {
                open(MAIL, "|$mail -s 'Login Detected' $to_email");
                print MAIL "Someone has used Login.\n"; 
		print MAIL "$ENV{'REMOTE_ADDR'} (with $ENV{'HTTP_USER_AGENT'})\n";
		print MAIL "$date\n";
	        if ($name == 1) {
	                print MAIL "  Name: $in{'name'}\n";
	        }
	        if ($email == 1) {
	                print MAIL "  E-mail: $in{'email'}\n";
	        }
		close(MAIL);
	}
}


print "<html><head><META HTTP-EQUIV=\"REFRESH\" CONTENT=\"0;URL=$forward\"></head></html>\n";
print "<a href=\"$forward\">Click here to contine!</a>";

sub print_login {
print <<"html";
<html><head><title>Login</title></head>
<body bgcolor=$bgcolor text=$text link=$link vlink=$vlink>
<center>
<font size=5>Please login to <b>$pagename</b>:
</center>
</font>
<form method=post action=\"$thisscript\">
html
if ($name == 1) {
	print "Name: <input type=text name=name><br>\n";
}
if ($email == 1) {
	print "E-mail Address: <input type=text name=email><br>\n";
}
print <<"html";
Password: <input type=password name=pass><br>
<input type=submit value=Submit name=submit>
</form>
<hr>
Login Script \&copy; 1998, by <a
href=\"http://tectonicdesigns.com/\" target=\"_top\">Andrew Cantino</a>
</body></html>
html
return 1;
}

sub getcgi {
    local($in, %in);
    local($name, $value);

        # If REQUEST_METHOD is POST, use CONTENT_LENGTH. Else, use QUERY_STRING.
    if ($ENV{'REQUEST_METHOD'} eq 'POST') {
        if ($ENV{'CONTENT_TYPE'}=~ m#^application/x-www-form-urlencoded#i) {
            read(STDIN, $in, $ENV{'CONTENT_LENGTH'});
                }
    }
    else {
                $in= $ENV{'QUERY_STRING'};
        }
        # Resolve and unencode name/value pairs into %in
    foreach (split('&', $in)) {
        s/\+/ /g;
        ($name, $value)= split('=', $_, 2);
        $name=~ s/%(..)/sprintf("%c",hex($1))/ge;
        $value=~ s/%(..)/sprintf("%c",hex($1))/ge;
        $in{$name}.= $value;
    }
    return %in;
}

sub error {
	local($error) = $_[0];
print <<"html";
<html><head<title>Error</title></head>
<body bgcolor=white>
<center><font size=5>Error</font><p><b>
You must enter a $error!
</body></html>
html
	exit;
	return 1;
}
