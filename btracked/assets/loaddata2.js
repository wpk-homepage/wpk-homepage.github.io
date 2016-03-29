'use strict';

function json2data(json, name) {
	function bandpass(i,min,max){
		return function(x){
			return !((min!=undefined&&x[i]<min) || (max!=undefined&&x[i]>=max));
		}
	}
	
	function datapoint(item) {
		// morning: green
		// afternoon: orange
		// evening: purple
		var x = item[2];
		var y = item[3];
		var a = item[4];
		var h = item[5];
		//var c;
		//a = a < 1 ? 1 : (a < 25 ? (25-a)/24 : 0);
		//a = ('' + a).substr(0,5);
		//c = h < 12 ? '80,156,106' : (h < 15.5 ? '211,132,37' : '139,101,165');
		//c = 'rgba(' + c + ',' + a + ')';
		return {
			x: x * 100,	// %
			y: y * 3.6,	// km/h
			//color: c,
			//size: 6,
		};
	}
	
	var today     = new Date();
	var onejan    = new Date(today.getFullYear(),0,1);
	var day       = Math.ceil((today - onejan) / 86400000);
	
	var accuracy  = bandpass(4, 0, 13);
	var longago   = bandpass(5, undefined, day-1);
	var yesterday = bandpass(5, day-1, day);
	var today     = bandpass(5, day);
	
	var values    = json.filter(accuracy);
	
	var data      = [
		{
			key: name || 'Frank',
			values: values.map(datapoint),
		},
		//{
		//	key: name || 'Hele reis',
		//	values: values.filter(longago).map(datapoint),
		//},
		//{
		//	key: 'Gisteren',
		//	values: values.filter(yesterday).map(datapoint),
		//},
		//{
		//	key: 'Vandaag',
		//	values: values.filter(today).map(datapoint),
		//},
	];
	return data;
}


function loaddata(filename, cb) {
	var xhr = new XMLHttpRequest();
	xhr.onload = function() {
		var json = JSON.parse(xhr.responseText);
		var data = json2data(json);
		cb(data);
	};
	xhr.open('get', filename, !!cb);
	xhr.send();
}

