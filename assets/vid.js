//javascript:(function(){document.body.appendChild(document.createElement('script')).src='http://paulkoppen.com/assets/vid.js';})();
(function(){

/** dobutton
 * Press the "Human verification" button on the web page to submit the form.
 * Then analyze the HTML and extract the video link. This link is returned
 * using a callback function.
 *
 * @param  fcn  Callback function will be invoked with the video link
 */
function dobutton(fcn) {
	var btn  = $(':submit'),
		frm  = btn.parents('form'),
		data = frm.serialize() + '&' +
					btn.attr('name') + '=' + btn.val().replace(/ /g, '+');
	$.ajax({
		url:      location.href,
		type:     'POST',
		data:     data,
		dataType: 'text',
		error:    function(){alert('Failed');},
		success:  function(html){
			var p   = html.match(/eval\((function.+?)\)\s*<\/script/)[1],
				js  = eval('(' + p.replace(/\}/, '})'));
			var u   = js.match(/src="(.+?)"/) || js.match(/'file','(http.+?)'/);
			var ddl = u[1];
			fcn(ddl);
		},
	});
}

/** doembed
 * Find the video object in the dom and extract the link to the video file.
 *
 * @param  fcn  Callback function will be invoked with the video link
 */
function doembed(fcn) {
	fcn($('embed').attr('flashvars').split(/\bfile=/)[1].split('&')[0]);
}

/** dovideolink
 * Get link to video file based on videolink.php.
 *
 * @param  fcn  Callback function will be invoked with the video link
 */
function dovideolink(fcn) {
	var v  = /\bv=(\w+)/.exec(location.search)[1],
		vl = '/xml/videolink.php?v=' + v + '&width=1080&id=' +
					(new Date()).getTime() + '&u=undefined';
	$.ajax({
		url:      vl,
		dataType: 'text',
		error:    function(){alert('Failed');},
		success:  function(xml) {
			var item, maxw=0;
			$(xml).find('row').each(function(){
				var row  = {},
					attr = this.attributes,
					$row = $(this),
					emb  = $row.attr('embed'),
					w    = emb.replace(/^.*width%3D%22(\d+).*$/, '$1');
				if (w > maxw)
					for (var i=attr.length; i--;)
						if (attr[i].name != 'embed')
							row[attr[i].name] = attr[i].value;
				item = $.extend({}, row);
			});
			// submit item, permanently redirects to download
			var s   = codestring(item.un, item.k1, item.k2),
				ddl = 'http://www' + item.s + '.megavideo.com/files/' + s + '/';
			fcn(ddl);
		}
	});
}

/** codestring
 * Make new string from un using k1 and k2.
 *
 * @param  un  Hexadecimal string
 * @param  k1  Number
 * @param  k2  Number
 * @return s   Hexadecimal string
 */
function codestring(un, k1, k2) {
	var hex2bin = {}, bin2hex = {};
	for (var i=0; i<16; i++) {
		hex2bin[i.toString(16)]=('000'+i.toString(2)).substr(-4);
		bin2hex[hex2bin[i.toString(16)]]=i.toString(16);}
	var binary = [];
	for (var i=0,len=un.length; i<len; i++)
		binary.push(hex2bin[un[i]]);
	binary = binary.join('').split(/\B/);
	var keys = [];
	for (var i=0; i<384; i++) {
		k1 = (k1*11+77213) % 81371;
		k2 = (k2*17+92717) % 192811;
		keys.push((k1+k2) % 128);
	}
	for (i=257; i--;) {
		var idx1 = keys[i],
			idx2 = i % 128,
			tmp  = binary[idx1];
		binary[idx1] = binary[idx2];
		binary[idx2] = tmp;
	}
	for (i=0; i<128; i++)
		binary[i] ^= keys[i+256] & 1;
	binary = binary.join('');
	var bytes = [];
	for (var i=0,len=binary.length; i<len; i+=4)
		bytes.push(bin2hex[binary.substr(i,4)]);
	return bytes.join('');
}

function downloadVideo($) {
	var lr = function(u){location.replace(u);};
	switch (location.hostname) {
	case 'www.megavideo.com':
	case 'wwwstatic.megavideo.com':
		dovideolink(lr);
		break;
	case 'www.vidbux.com':
	case 'www.vidxden.com':
		if ($('embed').length)
			doembed(lr);
		else dobutton(lr);
		break;
	default:
		alert('yes');
	}
}


/* Load jQuery, then download video */
if (typeof jQuery != 'undefined') downloadVideo(jQuery);
else {
	document.body.appendChild(document.createElement('script')).src =
				'http://code.jquery.com/jquery-1.6.min.js';
	(function init(){if(typeof jQuery=='undefined')
		return setTimeout(init,500);downloadVideo(jQuery)})();
}

})();

