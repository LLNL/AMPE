// Routine to get a nominated cookie value
function getCookieValue (cookieName) {
  var exp = new RegExp (cookieName + "=([^;]+)");
  if (exp.test (document.cookie + ";")) {
    exp.exec (document.cookie + ";");
    return unescape(RegExp.$1);
  }
  else return false
}

// Get the "password" cookie value, otherwise goto the default invalid page

var invalidpassword  = "http://braemoor.co.uk/software/invalid.html";
if (!getCookieValue ("validpassword")) {
  location.replace (invalidpassword);
}
else {
  // Get the cookie value
  var myCookie = getCookieValue ("password");

  // If it isn't 0, then it is the URL of the invalid password URL
  if (myCookie != "0") {location.replace (myCookie);}
}