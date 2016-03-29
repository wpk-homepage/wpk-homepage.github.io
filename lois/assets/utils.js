

/** preload(src)
 * Preloads an array of images.
 * @param {Array | string} src  List of image URL's or a single image URL.
 */
function preload(src) {
	if (typeof src == 'string')
		(new Image()).src = src;
	else
		for (var i=src.length; i--; )
			(new Image()).src = src[i];
}

