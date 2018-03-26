<?php

// PHP is too old on www.llnl.gov server for file_get_contents :(
//$theFile = file_get_contents('../cgi/sunmarkers.xml');

$filename = '../cgi/sunmarkers.xml';

$handle = fopen($filename, "r");
$theFileSize = filesize($filename);
$theFile = fread($handle, $theFileSize);
fclose($handle);

header("Content-Type: text/xml; charset=ISO-8859-1");

// Force browser to always reload the file, instead of caching.
header("Expires: Mon, 26 Jul 1997 05:00:00 GMT");
header("Last-Modified: ". gmdate("D, d M Y H:i:s") ." GMT");
header("Cache-Control: no-store, no-cache, must-revalidate");
header("Cache-Control: post-check=0, pre-check=0", false);
header("Pragma: no-cache");

header('Content-length: '.$theFileSize);

print $theFile;

?>
