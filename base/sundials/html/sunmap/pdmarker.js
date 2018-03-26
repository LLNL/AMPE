/*
 PdMarker

 Purpose: extends Google Map API GMap and GMarker (hover effects, image swapping, moving)
 Details: http://wwww.pixeldevelopment.com/pdmarker.asp
 Updated: [see getPdMarkerRevisionInfo]
 Author:  Peter Jones
 Notes:   Relies on undocumented features of the Google Map API which may change.
	    Based on my own PJToolTip and ideas from GxMarker, TLabel and the Google Maps API forum.

 Contact www.pixeldevelopment.com for your custom Google Map needs
*/

function getPdMarkerRevisionInfo() {
var cr = "<br/>";
var s =
"1.99a 04/18/06 - revised for Google Maps API Version 2, GMap2 required." + cr +
"0.99c 01/30/06 - added setDetailWinClass and resetDetailWinClass." + cr +
"0.99a 10/12/05 - now handles maps in containers with undefined widths" + cr +
"define a div with id 'pdmarkerwork' to reduce flicker" + cr +
"0.99  10/03/05 - added setImageEnabled, allowLeftTooltips (global)" + cr +
"0.98  09/30/05 - fixed zoomToMarkers" + cr +
"0.97  09/24/05 - added setHoverImage, setShowDetailOnClick, setDetailWinHTML, showDetailWin, closeDetailWin" + cr +
"0.96  09/22/05 - added setTooltipHiding, getTooltipHiding" + cr +
"0.95  09/20/05 - handle zoom for lingering tooltips mouseOutEnabled(false) " +
		   "disables setImage and restoreImage" + cr +
"0.94  09/20/05 - added setTooltipClass and resetTooltipClass" + cr +
"0.93  09/19/05 - added slopPercentage [optional] parameter to zoomToMarkers" + cr +
"0.92  09/18/05 - added getMouseOutEnabled, setMouseOutEnabled" + cr +
"0.91  09/17/05 - fixed setOpacity";
return s;
}

function getPdMarkerVersion() {
	return getPdMarkerRevisionInfo().substring(0,15);
}

function getPdMarkerShortVersion() {
	return getPdMarkerRevisionInfo().substring(0,5);
}

function getGoogleMapsVersion() {
	var i, a, b, c;
	var v = "unknown";

	if (document.getElementsByTagName)
		for(i=0; (a = document.getElementsByTagName("script")[i]); i++)
			if(a.getAttribute("src"))
			{
				b = a.getAttribute("src");
				c = b.indexOf("/mapfiles/maps");
				if (c > 0)
					v = parseFloat(b.substring(c+14));
			}
	return v;
}

function latLongToPixel(map,coord,zoom) {
	var topLeft = map.getCurrentMapType().getProjection().fromLatLngToPixel(map.fromDivPixelToLatLng(new GLatLng(0,0),true),map.getZoom());
	var point = map.getCurrentMapType().getProjection().fromLatLngToPixel(coord,map.getZoom());
	return new GPoint(point.x - topLeft.x, point.y - topLeft.y);
}

var pdMarkerOpenList = [];

function PdMarkerAddDetailWinOpen(marker) {
	pdMarkerOpenList.push(marker);
}

function PdMarkerClose(id) {
	for (var i=0; i<pdMarkerOpenList.length; i++)
		if (pdMarkerOpenList[i].internalId == id)
			{
				pdMarkerOpenList[i].closeDetailWin();
				pdMarkerOpenList.splice(i,1);
			}
}

// GMap extension for walking through PdMarker list
// Note: some overlays are not markers, some may not be PdMarkers

function isPdMarker(a) {
	if (a.isMarker)
		return true;		
	return false;
}

function getPdMarkerCount(a) {
	if (a.pdMarkers)
		return a.pdMarkers.length;
	return 0;
}

GMap2.prototype.getMarkerById = function(id) {
	var count = getPdMarkerCount(this);
	for (var i = 0; i < count; i++)
		if (isPdMarker(this.pdMarkers[i]))
			if (this.pdMarkers[i].internalId == id)
			{
				this.cursor = i;
				return this.pdMarkers[i];
			}
	return null;
}

GMap2.prototype.getFirstMarker = function() {
	var count = getPdMarkerCount(this);
	for (var i = 0; i < count; i++)
		if (isPdMarker(this.pdMarkers[i]))
		{
			this.cursor = i;
			return this.pdMarkers[i];
		}
	return null;
}

GMap2.prototype.getNextMarker = function() {
	var count = getPdMarkerCount(this);
	if (count > 0)
		if (this.cursor >= 0)
			for (var i = this.cursor+1; i < count; i++)
				if (isPdMarker(this.pdMarkers[i]))
				{
					this.cursor = i;
					return this.pdMarkers[i];
				}
	return null;
}

GMap2.prototype.getNthMarker = function(nTh) {
	var count = getPdMarkerCount(this);
	for (var i = 0; i < count; i++)
		if (isPdMarker(this.pdMarkers[i]))
		{
			nTh--;
			if (nTh == 0)
			{
				this.cursor = i;
				return this.pdMarkers[i];
			}
		}
	return null;
}

GMap2.prototype.getMarkerCount = function() {
	return getPdMarkerCount(this);
}

GMap2.prototype.boxMap = function(center, span) {
	var spec = this.spec;
	var zoom = spec.getLowestZoomLevel(center, span, this.viewSize);
	this.centerAndZoom(new GPoint(center.x, center.y), zoom);
}

GMap2.prototype.zoomToMarkers = function(slopPercentage, heightOffsetPct) {
	var count = 0;
	var thePoint, x, y, minX, maxX, minY, maxY, span;
	var marker = this.getFirstMarker();
	while (marker != null)
	{
		thePoint = marker.getPoint();
		// x = thePoint.x; y = thePoint.y;
		x = thePoint.lat(); y = thePoint.lng();
		if (count == 0)
		{
			minX = x; maxX = x; minY = y; maxY = y;
		}
		else
		{
			if (x < minX) minX = x;
			if (x > maxX) maxX = x;
			if (y < minY) minY = y;
			if (y > maxY) maxY = y;
		}
		marker = this.getNextMarker();
		count++;
	}
	if (count == 1)
		this.setCenter(new GLatLng(x,y), this.getZoom());
	else if (count > 1)
	{
		var center = new GLatLng((minX + maxX) / 2, (minY + maxY) / 2)
		span = new GSize(Math.abs(maxX - minX), Math.abs(maxY - minY));
		slopWid = 0;
		slopHgt = 0;
		if (typeof slopPercentage != "undefined")
		{
			slopWid = span.width * slopPercentage / 200;
			slopHgt = span.height * slopPercentage / 200;
			span.width  *= 1 + slopPercentage / 100;
			span.height *= 1 + slopPercentage / 100;
		}
		deltaHgt = 0;
		if (typeof heightOffsetPct != "undefined")
		{
			deltaHgt = span.height * heightOffsetPct / 100;
			center = new GLatLng(center.lat() + deltaHgt, center.lng());
		}
		// needs slop
		var bounds = new GLatLngBounds(new GLatLng(minX-slopHgt, minY-slopWid), new GLatLng(maxX+slopHgt, maxY+slopWid)); // sw, ne
		var zoom = this.getBoundsZoomLevel(bounds);
		this.setCenter(center, zoom);
	}
}

function getPoweredBy() {
	var i, a, b, c;

	if (document.getElementsByTagName)
	{
		var imgs = document.getElementsByTagName("img");
		for(i=0; i < imgs.length; i++)
		{
			var a = imgs[i];
			if(a.getAttribute("src")) // only returned for Firefox (png uses filter in IE)
			{
				b = a.getAttribute("src");
				c = b.indexOf("/poweredby.png");
				if (c > 0)
					return a;
			}
		}
	}
	return null;
}

function setPoweredBy(map) {
	if (!map.poweredByObj)
		map.poweredByObj = getPoweredBy();
	if (map.poweredByObj)
	{
		var newTitle = "GMap " + getGoogleMapsVersion() + " & PdMarker " + getPdMarkerShortVersion();
		map.poweredByObj.setAttribute("title", newTitle);
		map.poweredByObj.setAttribute("alt", newTitle);
	}
}


// PdMarker code


function PdMarkerNamespace() {

var userAgent = navigator.userAgent.toLowerCase();
var n4=(document.layers);
var n6=(document.getElementById&&!document.all);
var ie=(document.all);
var o6=(userAgent.indexOf("opera") != -1);
var safari=(userAgent.indexOf("safari") != -1);
var msie  = (userAgent.indexOf("msie") != -1) && (userAgent.indexOf("opera") == -1);

var nextMarkerId = 10;
var permitLeft = true;

icon = new GIcon();
icon.shadow = "http://www.google.com/mapfiles/shadow50.png";
icon.iconSize = new GSize(20, 34);
icon.shadowSize = new GSize(37, 34);
icon.iconAnchor = new GPoint(9, 34);
icon.infoWindowAnchor = new GPoint(9, 2);
icon.infoShadowAnchor = new GPoint(18, 25);
icon.image = "http://www.google.com/mapfiles/marker.png";

// Globals - careful of multiple maps

function PdMarker(a, b, tooltip) {
	this.inheritFrom = GMarker;
	if (typeof b == "undefined") // pmj oct 23, 2005
		b = icon;
	this.inheritFrom(a,b);
	if (typeof tooltip != "undefined")
		this.pendingTitle = tooltip;
	else
		this.pendingTitle = "";
	this.internalId = nextMarkerId;
	nextMarkerId += 1;
	this.zIndexSaved = false;
	this.oldImagePath = "";
	this.pendingCursor = "";
	this.percentOpacity = 70;
	this.mouseOutEnabled = true;
	this.setImageOn = true;
	this.hidingEnabled = true;
	this.showDetailOnClick = true;
	this.detailOpen = false;
	this.userData = "";
}

// PdMarker.prototype = new GMarker;
PdMarker.prototype = new GMarker(new GLatLng(1, 1));


function addMarkerToMapList(map,marker) {
	try {
		if (map.pdMarkers.length) ;
	}
	catch(e) {
		map.pdMarkers = new Array();
	}
	// add to list
	map.pdMarkers.push(marker);
}

function removeMarkerFromMapList(map,marker) {
	var id = marker.internalId;
	for (var i=0; i<map.pdMarkers.length; i++)
		if (map.pdMarkers[i].internalId == id)
		{
			map.pdMarkers.splice(i,1);
			return;
		}
}

PdMarker.prototype.initialize = function(a) {
	if (typeof a == "GMap")
	{
		GLog.write("PdMarker requires GMap2");
		return;
	}
	addMarkerToMapList(a,this);
	try
	{
		GMarker.prototype.initialize.call(this, a);
		this.isMarker = true;
		if (this.pendingTitle.length > 0)
			this.setTitle(this.pendingTitle);
		if (this.pendingCursor.length > 0)
			this.setCursor(this.pendingCursor);

		this.map = a;
		setPoweredBy(a);

		GEvent.bindDom(this, "mouseover", this, this.onMouseOver);
		GEvent.bindDom(this, "mouseout",  this, this.onMouseOut);
		GEvent.bindDom(this, "click",  this, this.onClick);
		GEvent.bind(this.map, "zoom", this, this.reZoom);
	}
	catch(e) {
		alert("PdMarker initialize error: " + e);
	}
}

PdMarker.prototype.allowLeftTooltips = function(a){
	permitLeft = a;
}

PdMarker.prototype.reZoom = function(){
	var didSet = false;
	if (this.tooltipObject)
		if (this.tooltipObject.style.display == "block")
		{
			setTTPosition(this);
			didSet = true;
		}
	if (this.detailObject)
	{
		if (!didSet)
			setTTPosition(this);
		this.detailObject.style.top  = this.ttTop + "px";
		this.detailObject.style.left  = this.ttLeft + "px";
	}
}

PdMarker.prototype.setId = function(id) {
	this.internalId = id;
}

PdMarker.prototype.getId = function() {
	return this.internalId;
}

PdMarker.prototype.setName = function(a) {
	this.name = a;
}

PdMarker.prototype.getName = function() {
	if (this.name)
		return this.name;
	else
		return null;
}

PdMarker.prototype.setUserData = function(a) {
	this.userData = a;
}

PdMarker.prototype.getUserData = function() {
	if (this.userData)
		return this.userData;
	else
//		return null;
		return "";
}

PdMarker.prototype.setUserData2 = function(a) {
	this.userData2 = a;
}

PdMarker.prototype.getUserData2 = function() {
	if (this.userData2)
		return this.userData2;
	else
		return "";
}

PdMarker.prototype.setImageEnabled = function(a) {
	this.setImageOn = a;
}

var PdMIN = "";

function PdCompPdMIN(marker) {
	if (PdMIN.length == 0)
		for (var i in marker)
			if (eval("typeof marker." + i) == "object")
				try {
					if (eval("typeof marker." + i + "[0].src") != "undefined")
						PdMIN = "this." + i + "[0]";
				}
				catch (e) {}
}

PdMarker.prototype.setImage = function(a) {
	var msFilter = 'progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + a + '")';
	if (this.mouseOutEnabled && this.setImageOn)
	{
		PdCompPdMIN(this);
		if (this.oldImagePath.length == 0)
			eval("this.oldImagePath = " + PdMIN + ".src");
		if (msie)
			eval(PdMIN + ".style.filter = msFilter");
		else
			eval(PdMIN + ".src = a");
	}
}

PdMarker.prototype.restoreImage = function() {
	if (this.mouseOutEnabled && this.setImageOn && this.oldImagePath.length > 0)
	{
		var a = this.oldImagePath;
		this.setImage(a);
		this.oldImagePath = "";
	}
}

PdMarker.prototype.setIcon = function(a) {
	this.remove();
	this.icon = a;
	this.initialize(this.map);
	this.redraw(true); 
}

PdMarker.prototype.setMarkerZIndex = function(a) {
	if (!this.zIndexSaved)
	{
		this.zIndexSaved = true;
		this.oldZIndex = eval(PdMIN + ".style.zIndex");
	}
	eval(PdMIN + ".style.zIndex = a")
	this.redraw(true);
}

PdMarker.prototype.topMarkerZIndex = function() {
	this.setMarkerZIndex (600000);
}

PdMarker.prototype.restoreMarkerZIndex = function() {
	if (this.zIndexSaved)
	{
		this.zIndexSaved = false;
		eval(PdMIN + ".style.zIndex = this.oldZIndex")
		this.redraw(true);
	}
}

PdMarker.prototype.onInfoWindowOpen = function() {
	this.hideTooltip();
	GMarker.prototype.onInfoWindowOpen.call(this);
}

PdMarker.prototype.setHoverImage = function(a) {
	this.hoverImage = a;
}

var inMouseOver = false;

// problem.  if we use new method of changing image (re-init, etc), then we keep getting mouse-over events.
// which truly blows.
// maybe force things to 2.45 and use images array, assuming a name of 'eb'

PdMarker.prototype.onMouseOver = function() {
	if (inMouseOver)
		return;
	inMouseOver = true;
	if (this.hoverImage)
		this.setImage(this.hoverImage);
	if (!this.detailOpen)
		this.showTooltip();
	inMouseOver = false;
}

PdMarker.prototype.onMouseOut = function() {
	if (this.hoverImage)
		this.restoreImage();
	if (!this.detailOpen)
		if (this.mouseOutEnabled)
			this.hideTooltip();
}

PdMarker.prototype.setMouseOutEnabled = function(a) {
	this.mouseOutEnabled = a;
}

PdMarker.prototype.getMouseOutEnabled = function() {
	return this.mouseOutEnabled;
}

PdMarker.prototype.setTooltipHiding = function(a) {
	this.hidingEnabled = a;
}

PdMarker.prototype.getTooltipHiding = function() {
	return this.hidingEnabled;
}

PdMarker.prototype.setTitle = function(a) {
	this.tooltipText = "";
	PdCompPdMIN(this);
	try {
		eval(PdMIN + ".title = a");
	}
	catch (e) {
		this.pendingTitle = a;
	}
}

PdMarker.prototype.setCursor = function(a) {
	PdCompPdMIN(this);
	try {
		eval(PdMIN + ".style.cursor = a");
	}
	catch (e) {
		this.pendingCursor = a;
	}
}

PdMarker.prototype.setTooltipClass = function(a) {
	this.pendingClassName = a;
	if (this.tooltipObject)
	{
		var showing = (this.tooltipObject.style.display != "none");
		this.deleteObjects();
		if (this.tooltipRaw)
			this.setTooltipNoResize(this.tooltipRaw);
		if (showing)
			this.showTooltip();

	}
	else
		if (this.tooltipRaw)
			this.setTooltipNoResize(this.tooltipRaw);
}

PdMarker.prototype.resetTooltipClass = function() {
	this.setTooltipClass("markerTooltip");
}

PdMarker.prototype.getTooltip = function() {
	try {
		return this.tooltipRaw;
	}
	catch (e)
	{
		return "";
	}
}

PdMarker.prototype.setTooltipNoResize = function(a) {
	this.setTitle("");
	var ttClass = "markerTooltip";
	if (this.pendingClassName)
		ttClass = this.pendingClassName;
	this.tooltipRaw = a;
	this.tooltipText = "<div class='" + ttClass + "'>" + a + "</div>";
	if (this.tooltipObject)
		this.tooltipObject.innerHTML = this.tooltipText;
}

PdMarker.prototype.setTooltip = function(a) {
	this.setTooltipNoResize(a);
	this.deleteObjects();
}

PdMarker.prototype.showTooltip = function() {
	if (this.tooltipText)
	{
		if (!this.tooltipObject)
			initTooltip(this);
		setTTPosition(this);
		this.tooltipObject.style.display = "block";
	}
}

PdMarker.prototype.hideTooltip = function() {
	if (this.tooltipObject)
		if (this.hidingEnabled)
			this.tooltipObject.style.display = "none";
}

PdMarker.prototype.onClick = function(a) {
	if (this.showDetailOnClick && this.detailWinHTML)
		this.showDetailWin();
}

PdMarker.prototype.setShowDetailOnClick = function(a) {
	this.showDetailOnClick = a;
}

PdMarker.prototype.setDetailWinHTML = function(a) {
	this.detailWinHTML = a;
}




PdMarker.prototype.setDetailWinClass = function(a) {
	this.pendingDetailClassName = a;
}

PdMarker.prototype.resetDetailWinClass = function() {
	this.setDetailWinClass("markerDetail");
}



PdMarker.prototype.showDetailWin = function() {
	if (this.detailOpen)
	{
		this.closeDetailWin();
		return;
	}
	this.hideTooltip();
	this.setMouseOutEnabled(false);

	var winClass = "markerDetail";
	if (this.pendingWinClassName)
		winClass = this.pendingWinClassName;

	var html = "<table><tr><td>" + this.detailWinHTML + "<\/td><td valign='top'><a class='markerDetailClose' href='javascript:PdMarkerClose(" + this.internalId + ")'><img src='http://www.google.com/mapfiles/close.gif' width='14' height='13'><\/a><\/td><\/tr><\/table>";
	html = "<div class='" + winClass + "'>" + html + "</div>";
	this.detailOpen = true;
	if (!this.tooltipText)
	{
		this.ttWidth = 150;
		this.ttHeight = 30;
		setTTPosition(this); // compute ttTop, ttLeft
	}
	initDetailWin(this, this.ttTop, this.ttLeft, html);
	PdMarkerAddDetailWinOpen(this);
}


PdMarker.prototype.closeDetailWin = function() {
	this.detailOpen = false;
	if (this.detailObject)
	{
		this.setMouseOutEnabled(true);
		this.onMouseOut();
		// GEvent.trigger(this, "mouseout");
	      this.map.getPane(G_MAP_FLOAT_PANE).removeChild(this.detailObject);
		this.detailObject = null;
	}
}

PdMarker.prototype.deleteObjects = function() {
	if (this.tooltipObject)
	{
	      this.map.getPane(G_MAP_FLOAT_PANE).removeChild(this.tooltipObject);
		this.tooltipObject = null;
	}
	if (this.detailObject)
	{
		this.map.getPane(G_MAP_FLOAT_PANE).removeChild(this.detailObject);
		this.detailObject = null;
	}
}

PdMarker.prototype.remove = function(a) {
	removeMarkerFromMapList(this.map, this);
	GMarker.prototype.remove.call(this);
	this.deleteObjects();
}

PdMarker.prototype.setOpacity = function(b) {
	if (b < 0)
		b=0;
	if (b >= 100)
		b=100;
	var c = b / 100;
	this.percentOpacity = b;
	var d = document.getElementById(this.objId);
	if (d)
	{
		if(typeof(d.style.filter)=='string'){d.style.filter='alpha(opacity:'+b+')';}
		if(typeof(d.style.KHTMLOpacity)=='string'){d.style.KHTMLOpacity=c;}
		if(typeof(d.style.MozOpacity)=='string'){d.style.MozOpacity=c;}
		if(typeof(d.style.opacity)=='string'){d.style.opacity=c;}
	}
}

PdMarker.prototype.setOpacityNew = function(b) {
	setObjOpacity(this.objId);
	this.percentOpacity = b;
}

// ***** Private routines *****

function setObjOpacity(objId, b) {
	if (b < 0)
		b=0;
	if (b >= 100)
		b=100;
	var c = b / 100;
	var d = document.getElementById(objId);
	if (d)
	{
		if(typeof(d.style.filter)=='string'){d.style.filter='alpha(opacity:'+b+')';}
		if(typeof(d.style.KHTMLOpacity)=='string'){d.style.KHTMLOpacity=c;}
		if(typeof(d.style.MozOpacity)=='string'){d.style.MozOpacity=c;}
		if(typeof(d.style.opacity)=='string'){d.style.opacity=c;}
	}
}

function idToElemId(id) {
	return "ttobj" + id;
}

function initTooltip(theObj) {
	theObj.objId = idToElemId(theObj.internalId);
	theObj.anchorLatLng = theObj.point;

	var b = document.createElement('span');
	theObj.tooltipObject = b;
	b.setAttribute('id',theObj.objId);
	b.innerHTML = theObj.tooltipText;

	// append to body for size calculations
	var c = document.body;
	var d = document.getElementById("pdmarkerwork");
	if (d)
		c = d;
	c.appendChild(b);
	b.style.position = "absolute";
	b.style.bottom = "5px";
	b.style.left = "5px";
	b.style.zIndex = 1;
	if (theObj.percentOpacity)
		theObj.setOpacity(theObj.percentOpacity);
	var tempObj = document.getElementById(theObj.objId);
	theObj.ttWidth  = tempObj.offsetWidth;
	theObj.ttHeight = tempObj.offsetHeight;
	c.removeChild(b);

	b.style.zIndex = 600000;
	b.style.bottom = "";
	b.style.left = "";
	// theObj.map.div.appendChild(b);
      theObj.map.getPane(G_MAP_FLOAT_PANE).appendChild(b);
}

function initDetailWin(theObj, top, left, html) {
	theObj.detailId = "detail" + theObj.internalId;
	var b = document.createElement('span');
	theObj.detailObject = b;
	b.setAttribute('id',theObj.detailId);
	b.innerHTML = html;
	b.style.display = "block";
	b.style.position = "absolute";
	b.style.top  = top + "px";
	b.style.left = left + "px";
	b.style.zIndex = 600001;
	map = theObj.map;
      map.getPane(G_MAP_FLOAT_PANE).appendChild(b);
}

function setTTPosition(theObj) {
	var map = theObj.map;
	var pt  = theObj.getPoint();
	var ttPos = latLongToPixel(map, pt, map.getZoom());
	var theIcon = theObj.getIcon();
	ttPos.x += Math.floor(theIcon.iconAnchor.x/2);
	ttPos.y -= Math.floor(theIcon.iconAnchor.y/2);

	var rightSide = true;
	var bounds = map.getBounds();
	var boundsSpan	= bounds.toSpan();
	var longSpan = boundsSpan.lng();
	var mapWidth = map.getSize().width;

	var tooltipWidthInDeg = (theObj.ttWidth + icon.iconSize.width + 6) / mapWidth * longSpan;
	if (pt.lng() + tooltipWidthInDeg > bounds.getNorthEast().lng() && permitLeft)
		rightSide = false;
	if (rightSide)
	{
		ttPos.y -= Math.floor(theObj.ttHeight/2);
		ttPos.x += icon.iconSize.width;
	}
	else
	{
		ttPos.y -= Math.floor(theObj.ttHeight/2);
		ttPos.x -= (icon.iconSize.width + theObj.ttWidth);
	}
	theObj.ttLeft = ttPos.x;
	theObj.ttTop  = ttPos.y;
	if (theObj.tooltipObject)
	{
		theObj.tooltipObject.style.left = ttPos.x + "px";
		theObj.tooltipObject.style.top  = ttPos.y + "px";
	}
}

function makeInterface(a) {
	var b = a || window;
	b.PdMarker = PdMarker;
}

makeInterface();
}


PdMarkerNamespace();
