// JavaScript April 2006 -- Radu Serban, LLNL

var map;
var pointName = [];
var pointLat = [];
var pointLng = [];
var pointCity = [];
var pointCountry = [];
var pointDetail = [];
var pointPhoto = [];
var pointID = [];
var npdm = 0;

var mode = 0;        // 0: navigate, 1: add

var new_marker;

var show_list = 0;   // 0: hide, 1: show

var add_location = 0;

var have_selection = 0;
var last_marker;
var last_row;


// -------------------------------------------------
// Functions for the SunMapCookie
// -------------------------------------------------

// Create a cookie with the realname + creation date
// and store its value in a hidden form field

function createSunmapCookie() {

  var date = new Date();
  date.setTime(date.getTime()+(365*24*60*60*1000));
  var realname = document.getElementById('realname').value;  
  var cookiename = realname+" "+String(date.getTime());
  document.cookie = "SunMapCookie=" + cookiename + "; expires=" + date.toGMTString() + "; path=/";

  document.getElementById('new_id').value = cookiename;

}

// Erase the SunMapCookie

function deleteSunmapCookie() {

    var date = new Date();
    date.setTime(date.getTime()-10000000);
    document.cookie = "SunMapCookie=; expires="+date.toGMTString()+"; path=/";

}

// Attempt to read the SunMapCookie

function readSunmapCookie() {

  var nameEQ = "SunMapCookie=";
  var ca = document.cookie.split(';');  
  for(var i=0;i < ca.length;i++) {
    var c = ca[i];
    while (c.charAt(0)==' ') c = c.substring(1,c.length);
    if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
  }

  return null;

}

// -------------------------------------------------
// Functions to toggle display of the list
// -------------------------------------------------

function toggleList() {
  show_list = 1-show_list;
  displayList(show_list);
}


function displayList(what) {

  var lst = document.getElementById('listing');

  if(what) { 

    document.getElementById('list').value="Hide List";
    lst.style.display='block';

  } else {

    map.closeInfoWindow();
    document.getElementById('list').value="Show List";
    lst.style.display='none';

  }

}

// -------------------------------------------------
// Functions to switch mode NAVIGATE <-> ADD
// -------------------------------------------------

function toggleAdd() {

  if(mode==0) {

    // We are currently in NAVIGATE mode

    var proceed = proceed2add();    

    if(proceed==1) {
      document.getElementById('add').value="Cancel add";
      switch2add();        
    }

  } else {

    document.getElementById('add').value="Add";
    switch2nav();
  }

}

// Test if we can proceed with adding a new marker

function proceed2add() {

  // By default, action is ADD

  document.getElementById('action').value = "add";

  // Try to read the SunMapCookie

  var cval = readSunmapCookie();

  // If no cookie found, proceed to ADD

  if(cval == null) return 1;

  // A cookie was found, seach for a marker with ID = cval

  var myidx = -1;
  for(var i=0; i<npdm; i++) {
    if(pointID[i] == cval) myidx = i;
  }

  // If no such marker was found, proceed to ADD

  if(myidx == -1) return 1;

  // A matching marker was found, ask the user what to do

  var action = confirm('A previous marker was found. Click OK to replace it.');

  // If the user cancels, do not proceed to ADD

  if(action==false) return 0;

  // The user wants to replace the marker. Change action to REPLACE.

  document.getElementById('action').value = "replace";
  return 1;

}

// Change mode to ADD

function switch2add() {

  mode = 1;

  // Parse all PdMarkers and set showDetailOnClick on FALSE
  var marker = map.getFirstMarker();
  while (marker != null) {
    marker.closeDetailWin();
    marker.showDetailOnClick = false;
    marker = map.getNextMarker();
  }

  // Close map blowup
  map.closeInfoWindow();

  // Show the 'invisible' marker
  new_marker.setImage("img/markergreen.png");

  // Set new message
  document.getElementById("message").innerHTML = "Click on the map to specify the position.";

}   

// Change mode to NAVIGATE

function switch2nav() {

  mode = 0;

  // Parse all PdMarkers and set showDetailOnClick on TRUE
  var marker = map.getFirstMarker();
  while (marker != null) {
    marker.showDetailOnClick = true;
    marker = map.getNextMarker();
  }

  // Hide add_form table
  tbl = document.getElementById('add_form');
  tbl.style.display='none';

  // Hide the 'invisible' marker
  new_marker.setImage("img/markerinvisible.png");
  new_marker.setPoint(new GLatLng(90,0));

  // Clear message
  document.getElementById("message").innerHTML = "";

} 

// -------------------------------------------------
// Functions to add new marker
// -------------------------------------------------

function submitAdd() {

  // Store location of new_marker in hidden form fields

  point = new_marker.getPoint();

  document.getElementById('lat').value = String(point.lat());
  document.getElementById('lng').value = String(point.lng());

  // Create a cookie that expires in one year with the name and creation date

  createSunmapCookie();

}

// -------------------------------------------------
// Functions to create the markers database
// -------------------------------------------------

// Read the sunmarkers.xml and create markers

function readXML() {

  var request = GXmlHttp.create();

  // Open XML file through PHP so that we can force reloading
  request.open("GET", "loadxml.php", true);

  request.onreadystatechange = function() {
    if(request.readyState == 4) {

      var xmlDoc = request.responseXML;

      // obtain the array of markers and loop through it
      var markers = xmlDoc.documentElement.getElementsByTagName("marker");
          
      for(var i = 0; i < markers.length; i++) {
        // obtain the attribues of each marker
        var lat = parseFloat(markers[i].getAttribute("lat"));
        var lng = parseFloat(markers[i].getAttribute("lng"));
        var point = new GLatLng(lat,lng);
        var name = markers[i].getAttribute("name");
        var city = markers[i].getAttribute("city");
        var country = markers[i].getAttribute("country");
        var photo = markers[i].getAttribute("photo");
        var id = markers[i].getAttribute("id");
        pointName[i] = name;
        pointPhoto[i] = photo;
        pointLat[i] = point.lat();
        pointLng[i] = point.lng();
        pointCity[i] = city;
        pointCountry[i] = country;
        pointDetail[i] = "<b>"+name+"</b><br>"+city+", "+country;
        if(pointPhoto[i].length > 0) {
          pointDetail[i] = pointDetail[i]+"<br><img src=\""+pointPhoto[i]+"\" height=80 border=1>";
        }
        pointID[i] = id;

        // create the marker
        var marker = createMarker(i, point);

        map.addOverlay(marker);

        // update listing table
        addEntry(i, marker);

      }

      npdm = markers.length;

    } // end if readyState

  } // end request.onreadystatechange function

  request.send(null);

}

// Function to create and set-up a new marker

function createMarker(i, point) {
  var marker = new PdMarker(point);
  marker.setTitle(pointName[i]);
  marker.setTooltip(pointName[i]);
  marker.setHoverImage("img/markeryellow.png");
  marker.setDetailWinHTML(pointDetail[i]);
  marker.setCursor("help");
  GEvent.addListener(marker, "mouseover", function() {
     if(mode==0) {
       var center = marker.getPoint();
       document.getElementById("message").innerHTML = pointName[i]+"'s location: "+center.toString();
     }
  }); 
  GEvent.addListener(marker, "mouseout", function() {
     var center = marker.getPoint();
     document.getElementById("message").innerHTML = "";
  }); 
  return marker;
}

// Add listing entry

function addEntry(i, marker) {

  // See if the user owns this marker
  var mine = 0;
  var cval = readSunmapCookie();
  if(cval == pointID[i]) {
     mine = 1;  
     document.getElementById('old_id').value = pointID[i];
  }

  var table = document.getElementById('listing_table');

  var tr          = document.createElement('TR');

  var td_i        = document.createElement('TD');
  var td_name     = document.createElement('TD');
  var td_city     = document.createElement('TD');
  var td_country  = document.createElement('TD');
  var td_lat      = document.createElement('TD');
  var td_lng      = document.createElement('TD');
  var td_del      = document.createElement('TD');

  td_lat.align="right";
  td_lng.align="right";
  td_del.align="center";

  td_name.onclick = "zoom2point("+pointLat[i]+","+pointLng[i]+")";

  var txt_i       = document.createTextNode(String(i+1));
  var txt_name    = document.createTextNode(" "+pointName[i]);
  var txt_city    = document.createTextNode(pointCity[i]);
  var txt_country = document.createTextNode(pointCountry[i]);
  var txt_lat     = document.createTextNode(Math.round(10000*pointLat[i])/10000);
  var txt_lng     = document.createTextNode(Math.round(10000*pointLng[i])/10000);

  var img_zoom = document.createElement('IMG');
  img_zoom.src = "img/zoom.png";
  img_zoom.style.cursor = "pointer";
  img_zoom.onclick = function() {
     zoomToMarker(marker,i+1);
     have_selection = 1;
     last_marker = marker;
     last_row = i+1;
  }

  var img_del = document.createElement('IMG');
  if (mine == 1) {
    img_del.src = "img/cross.png";
    img_del.style.cursor = "pointer";
    img_del.onclick = function() {
      var action = confirm('Are you sure you want to delete this marker?');
      if(action==false) return 1;
      document.getElementById('action').value = "delete";
      document.getElementById('old_id').value = pointID[i];
      document.forms.sunMapAdd.submit.click();
     }
  } else {
    img_del.src = "img/blank.png";
  }  

  table.appendChild(tr);

  tr.appendChild(td_i);
  tr.appendChild(td_name);
  tr.appendChild(td_city);
  tr.appendChild(td_country);
  tr.appendChild(td_lat);
  tr.appendChild(td_lng);
  tr.appendChild(td_del);

  td_i.appendChild(txt_i);
  td_name.appendChild(img_zoom);
  td_name.appendChild(txt_name);
  td_city.appendChild(txt_city);
  td_country.appendChild(txt_country);
  td_lat.appendChild(txt_lat);
  td_lng.appendChild(txt_lng);
  td_del.appendChild(img_del);
  
}

// Zoom to selected entry

function zoomToMarker(marker,row) {

  // Scroll to top of page

  document.documentElement.scrollTop=0;

  // Get rows of list table

  var table = document.getElementById('listing_table');
  var rows = table.getElementsByTagName("tr");

  // Reset previous selection

  if(have_selection) {
    // Restore the Z index of the previous marker
    last_marker.restoreMarkerZIndex();  
    // Restore table row color of previous selection
    var cols = rows[last_row].getElementsByTagName("td");
    cols[0].className="default";
  }

  // Zoom to marker for selected entry

    var point = marker.getPoint();
    map.setCenter(point, 6);  
    marker.topMarkerZIndex();
    marker.showMapBlowup();

  // Highlight selected row

  var cols = rows[row].getElementsByTagName("td");
  cols[0].className="selected";

}

// -------------------------------------------------
// Map listener functions
// -------------------------------------------------

function addMapListeners() {

  // Click event when adding new location

  GEvent.addListener(map, "click", function(marker,point){

  if(mode) {

      tbl = document.getElementById('add_form');
      tbl.style.display='block';

      // Move last marker
      new_marker.setPoint(point);

      // Set new message
      document.getElementById("message").innerHTML = "You can continue changing the location on the map until you click 'Submit' or 'Cancel add'.";

  }
  });


  // Zoomend event - hide any opened DetailWindows

  GEvent.addListener(map, 'zoomend', function (oldzoomlevel,newzoomlevel) {

    var marker = map.getFirstMarker();
    for(var i = 0; i < npdm; i++) {
      marker.closeDetailWin();
      marker = map.getNextMarker();
    }  

  }); 
}

// -------------------------------------------------
// -------------------------------------------------


function onPageLoad() {

  // Create the map

  map = new GMap2(document.getElementById("map"));

  map.addControl(new GLargeMapControl());
  map.addControl(new GMapTypeControl());
  map.addControl(new GScaleControl());

//  var Tsize = new GSize(200, 100);
//  map.addControl(new GOverviewMapControl(Tsize));

  map.setCenter(new GLatLng(30,0), 2, G_HYBRID_MAP);

  // Read XML file and create markers database

  readXML();

  // Create the 'invisible' marker at the North Pole

  var point = new GLatLng(90,0);
  new_marker = new PdMarker(point);
  map.addOverlay(new_marker);
  new_marker.setImage("img/markerinvisible.png");
  new_marker.setDetailWinHTML("huh");

  // Add map listener functions

  addMapListeners();

}

