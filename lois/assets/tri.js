// Library of functions to deal with triangle meshes
// Author: Paul Koppen
//
// The following javascript libraries are required:
//

/**
 * Parse VRML to a trimesh object
 *
 * @param  vrml  String holding the contents of a VRML file (.wrl)
 * @param  path  String specifying the relative path of the VRML file to resolve
 *               the referenced texture image.
 * @return tri   Trimesh. Object with properties .texturefile, .vertices,
 *               .faces, .texturecoords, and .textureindices. Each (except
 *               texturefile obviously) holds an array of vectors, either vec3
 *               or vec2.
 *
 * @note This function only parses the file and does not compute things like
 *       normals and such. You will need to compute them separately.
 */
function triload(vrml, path) {
	var i1, i2, i3, i4, i5, txt, matches;
	var data;
	
	// texture file name
	//
	i1      = vrml.indexOf('appearance Appearance');
	i2      = vrml.indexOf('texture ImageTexture', i1);
	i3      = vrml.indexOf('{', i2);
	i4      = vrml.indexOf('}', i3);
	txt     = vrml.substring(i3+1, i4);
	matches = txt.match(/url\s+['"](.+?)['"]/);
	
	var texturefile = matches ? path + matches[1] : '';
	
	// vertices
	//
	i1      = vrml.indexOf('geometry IndexedFaceSet');
	i2      = vrml.indexOf('coord Coordinate', i1);
	i3      = vrml.indexOf('point', i2);
	i4      = vrml.indexOf('[', i3);
	i5      = vrml.indexOf(']', i4);
	txt     = vrml.substring(i4+1, i5);
	matches = txt.match(/[\d\.e-]+/g);	// nasty e, nasty e!!
	
	var vertices = Two.Arr([3,matches.length/3]);
	data = vertices.data;
	for (var i=matches.length; i--; )
		data[i] = parseFloat(matches[i]);
	
	// faces
	//
	i1      = vrml.indexOf('coordIndex', i5);
	i2      = vrml.indexOf('[', i1);
	i3      = vrml.indexOf(']', i2);
	txt     = vrml.substring(i2+1, i3);
	matches = txt.match(/[\d-]+/g);
	
	var faces = Two.Arr([3,matches.length/4], null, Uint16Array);	// N.B. four
	data = faces.data;
	for (var i=0,cnt=matches.length,j=0; i<cnt; i+=4,j+=3) {
		data[j+0] = parseInt(matches[i+0]);
		data[j+1] = parseInt(matches[i+1]);
		data[j+2] = parseInt(matches[i+2]);
	}
	
	// texture coordinates (UV)
	//
	i1      = vrml.indexOf('texCoord TextureCoordinate', i3);
	i2      = vrml.indexOf('point', i1);
	i3      = vrml.indexOf('[', i2);
	i4      = vrml.indexOf(']', i3);
	txt     = vrml.substring(i3+1, i4);
	matches = txt.match(/[\d\.e-]+/g);
	
	var texturecoords = Two.Arr([2,matches.length/2]);
	data = texturecoords.data;
	for (var i=matches.length; i--; )
		data[i] = parseFloat(matches[i]);
	
	// texture indices
	//
	i1      = vrml.indexOf('texCoordIndex', i4);
	i2      = vrml.indexOf('[', i1);
	i3      = vrml.indexOf(']', i2);
	txt     = vrml.substring(i2+1, i3);
	matches = txt.match(/[\d-]+/g);
	
	var textureindices = Two.Arr([3,matches.length/4], null, Uint16Array);
	data = textureindices.data;
	for (var i=0,cnt=matches.length,j=0; i<cnt; i+=4,j+=3) {
		data[j+0] = parseInt(matches[i+0]);
		data[j+1] = parseInt(matches[i+1]);
		data[j+2] = parseInt(matches[i+2]);
	}
	
	// trimesh
	//
	return {
		texturefile:    texturefile,
		vertices:       vertices,
		faces:          faces,
		texturecoords:  texturecoords,
		textureindices: textureindices,
	};
}


/**
 * Make vertices and texture coords referenced by same index
 *
 * @param  {Tri} tri  Trimesh with .vertices, .normals, .faces, .texturecoords,
 *                    and .textureindices. Will be modified by reference.
 * @return {Tri} tri  Input argument with .vertices and .texturecoords
 *                    organised such that .faces equals .textureindices.
 *                    This is needed for OpenGL.
 *
 * @note The input value will be modified.
 */
function trisingleindex(tri) {
	var key2newidx       = {},
		newvertices      = [],
		newtexturecoords = [],
		newfaces         = [],
		idxmax           = 0;
	
	var vdata  = tri.vertices.data,
		fdata  = tri.faces.data,
		tcdata = tri.texturecoords.data,
		tidata = tri.textureindices.data;
	for (var i=0,cnt=fdata.length; i<cnt; i++) {
		var vi  = fdata[i],
			ti  = tidata[i],
			key = vi + ',' + ti,
			idx = key2newidx[key];
		if (typeof idx == 'number') {
			newfaces.push(idx);
		} else {
			//alert(key + ': ' + idxmax);
			key2newidx[key] = idxmax;
			newfaces.push(idxmax++);
			var ovi = 3 * vi,
				oti = 2 * ti;
			newvertices.push(vdata[ovi+0]);
			newvertices.push(vdata[ovi+1]);
			newvertices.push(vdata[ovi+2]);
			newtexturecoords.push(tcdata[oti+0]);
			newtexturecoords.push(tcdata[oti+1]);
			//newtexturecoords.push(tcdata[oti+2]);
		}
	}
	var nfaces         = newfaces.length / tri.faces.size[0];
	tri.vertices       = Two.Arr([3,idxmax], newvertices);
	tri.texturecoords  = Two.Arr([2,idxmax], newtexturecoords);
	tri.faces          = Two.Arr([tri.faces.size[0],nfaces], newfaces, Uint16Array);
	tri.textureindices = Two.Arr([tri.faces.size[0],nfaces], newfaces, Uint16Array);
	return tri;
}


/**
 * Compute vertex normals
 *
 * @param  tri      Object with attributes .vertices and .faces, each of which
 *                  is an array of vec3 objects.
 * @return normals  Array of unit length vec3 objects, one per vertex.
 *
 * @note Algorithm described in Max, 1999.
 */
function trinormals(tri) {
	var fdata   = tri.faces.data,
		flen    = tri.faces.data.length,
		fdim    = tri.faces.size[0],		// fdim == 3
		fnum    = tri.faces.size[1],
		vdata   = tri.vertices.data,
		vlen    = tri.vertices.data.length,
		vdim    = tri.vertices.size[0],
		vnum    = tri.vertices.size[1],
		normals = Two.Arr([vdim,vnum]),
		ndata   = normals.data;
	
	var Arr = tri.vertices.type;
	
	for (var fi=0; fi<flen; fi+=fdim) {
		var vi0 = fdata[fi] * vdim,
			v0  = vdata.subarray(vi0, vi0+vdim),
			vi1 = fdata[fi+1] * vdim,
			v1  = vdata.subarray(vi1, vi1+vdim),
			vi2 = fdata[fi+2] * vdim,
			v2  = vdata.subarray(vi2, vi2+vdim),
			e01 = [],
			e12 = [],
			e20 = [],
			l01 = 0,
			l12 = 0,
			l20 = 0,
			ans = vec3.create(),
			scale = 0;
		
		// edges and lengths (squared)
		//
		for (var vi=0; vi<vdim; vi++) {
			e01[vi] = v1[vi] - v0[vi];
			l01    += e01[vi] * e01[vi];
			e12[vi] = v2[vi] - v1[vi];
			l12    += e12[vi] * e12[vi];
			e20[vi] = v0[vi] - v2[vi];
			l20    += e20[vi] * e20[vi];
		}
		
		vec3.cross(e01, e20, ans);
		vec3.negate(ans);
		vec3.scale(ans, 1 / (l01 * l20));
		
		Two.times(ans, scale, ans);	// can do this with matrix .* rowvector
		for (var i=nfaces; i--; ) {
			Two.add_(ncols[fcols[i][0]], acols[i], ncols[fcols[i][0]],
						vdim, 1, vdim, 1, vdim, 1);
		}
	}
}

function trinormals_(tri) {
	var vertices = tri.vertices.data,
		vcols    = tri.vertices.columns,
		vdim     = tri.vertices.size[0],
		nverts   = tri.vertices.size[1],
		faces    = tri.faces.data,
		fcols    = tri.faces.columns,
		fdim     = tri.faces.size[0],
		nfaces   = tri.faces.size[1],
		// vertices
		A        = Two.Arr([vdim,nfaces]),
		B        = Two.Arr([vdim,nfaces]),
		C        = Two.Arr([vdim,nfaces]),
		Adata    = A.data,
		Bdata    = B.data,
		Cdata    = C.data,
		// edges & lengths
		AB       = Two.Arr([vdim,nfaces]),
		BC       = Two.Arr([vdim,nfaces]),
		CA       = Two.Arr([vdim,nfaces]),
		lAB      = Two.Arr([1,nfaces]),
		lBC      = Two.Arr([1,nfaces]),
		lCA      = Two.Arr([1,nfaces]),
		// intermediate answer
		ans      = Two.Arr([vdim,nfaces]),
		acols    = ans.columns,
		scale    = Two.Arr([1,nfaces]),
		// output
		normals  = Two.Arr(tri.vertices.size),
		ncols    = normals.columns;
	
	// vertices in face-order
	//
	for (var i=0,cnt=faces.length; i<cnt; i+=3) {
		for (var j=3; j--; ) {
			Adata[i+j] = vcols[faces[i+0]][j];
			Bdata[i+j] = vcols[faces[i+1]][j];
			Cdata[i+j] = vcols[faces[i+2]][j];
		}
	}
	
	// edges
	//
	Two.sub(B, A, AB);
	Two.sub(C, B, BC);
	Two.sub(A, C, CA);
	
	// edge lengths (squared)
	//
	Two.dot(AB, AB, 0, lAB);
	Two.dot(BC, BC, 0, lBC);
	Two.dot(CA, CA, 0, lCA);
	
	// face normals
	// and add to vertex normals
	//
	Two.cross(AB, CA, 0, ans);			// A -> B & C
	Two.neg(ans, ans);
	Two.times(lAB, lCA, scale);
	Two.div(Two.ONE, scale, scale);
	Two.times(ans, scale, ans);	// can do this with matrix .* rowvector
	for (var i=nfaces; i--; ) {
		Two.add_(ncols[fcols[i][0]], acols[i], ncols[fcols[i][0]],
					vdim, 1, vdim, 1, vdim, 1);
	}
	
	Two.cross(BC, AB, 0, ans);			// B -> C & A
	Two.neg(ans, ans);
	Two.times(lBC, lAB, scale);
	Two.div(Two.ONE, scale, scale);
	Two.times(ans, scale, ans);
	for (var i=nfaces; i--; ) {
		Two.add_(ncols[fcols[i][1]], acols[i], ncols[fcols[i][1]],
					vdim, 1, vdim, 1, vdim, 1);
	}
	
	Two.cross(CA, BC, 0, ans);			// C -> A & B
	Two.neg(ans, ans);
	Two.times(lCA, lBC, scale);
	Two.div(Two.ONE, scale, scale);
	Two.times(ans, scale, ans);
	for (var i=nfaces; i--; ) {
		Two.add_(ncols[fcols[i][2]], acols[i], ncols[fcols[i][2]],
					vdim, 1, vdim, 1, vdim, 1);
	}
	
	// normalise all vectors
	// zero-length vectors will be normalised to (0,0,0)
	//
	Two.normalise(normals, 0, normals);
	
	return normals;
}


/**
 * N-th order neighbourhood vertices around each vertex
 *
 * @param {trimesh} tri  Trimesh. Object with .vertices and .faces, both arrays
 *                       vec3 objects.
 * @param {int} order    Optional neighbourhood order, default is 1 which gets
 *                       all directly connected vertices.
 * @return {array{array}} neighbours  An array holding for each vertex an array
 *                                    with the indices of all vertices in its
 *                                    neighbourhood.
 *
 * @note Although the implementation is as fast as I could get it, it is still a
 *       very expensive operation for reasonably large meshes.
 */
function trineighbourhood(tri, order) {
	if (!order)
		order = 1;
	
	var nverts = tri.vertices.length;
	var nfaces = tri.faces.length;
	
	// pre-allocate some variables
	//
	var firstorder    = [];		// first order (direct) neighbours
	var firstorderset = [];		// alternative representation because searching
								// elements in arrays is expensive
	
	for (var i=0; i<nverts; i++) {
		firstorder[i]       = [i];	// zero-th order is self
		firstorderset[i]    = {};
		firstorderset[i][i] = true;
	}
	
	// first order neighbourhood
	//
	for (var i=0; i<nfaces; i++) {
		var face = tri.faces[i],
			a    = face[0],
			b    = face[1],
			c    = face[2];
		if (!firstorderset[a][b]) {
			firstorderset[a][b] = true;
			firstorder[a].push(b);
			firstorderset[b][a] = true;
			firstorder[b].push(a);
		}
		if (!firstorderset[a][c]) {
			firstorderset[a][b] = true;
			firstorder[a].push(c);
			firstorderset[c][a] = true;
			firstorder[c].push(a);
		}
		if (!firstorderset[b][c]) {
			firstorderset[b][c] = true;
			firstorder[b].push(c);
			firstorderset[c][b] = true;
			firstorder[c].push(b);
		}
	}
	
	if (order == 1)
		return firstorder;
	
	// next order neighbourhoods
	//
	var nthorder    = [];				// complete neighbourhood
	var nthorderset = firstorderset;	// we can safely overwrite values here
	var expansion   = [];				// incremental neighbourhood
	
	for (var i=0; i<nverts; i++) {
		nthorder[i]  = firstorder[i].slice();
		expansion[i] = firstorder[i].slice(1);
	}
	
	for (var i=1; i<order; i++) {
		for (var j=0; j<nverts; j++) {
			var nb      = [],
				nth     = nthorder[j],
				nthset  = nthorderset[j],
				exp     = expansion[j],
				nextexp = [];
			
			for (var k=0,cnt=exp.length; k<cnt; k++) {
				nb.push.apply(nb, firstorder[exp[k]]);
			}
			for (var k=0,cnt=nb.length; k<cnt; k++) {
				if (!nthset[nb[k]]) {
					nthset[nb[k]] = true;
					nth.push(nb[k]);
					nextexp.push(nb[k]);
				}
			}
			expansion[j] = nextexp;
		}
	}
	
	return nthorder;
}


/**
 * Principal, Gaussian, and mean curvatures at each vertex
 *
 */
function tricurvature(tri, nb, normals) {
	if (!nb)
		nb      = trineighbourhood(tri, 2);
	if (!normals)
		normals = trinormals(tri);
	
	var vertices = tri.vertices,
		nverts   = vertices.length,
		yaxis    = vec3.createFrom(0, 1, 0),
		zaxis    = vec3.createFrom(0, 0, 1),
		localx   = vec3.create(),
		localy   = vec3.create(),
		basis    = mat3.create(),
		K        = Array(nverts);
	
	for (var i=0; i<nverts; i++) {
		// construct a local coordinate system with its xy plane tangent to the
		// mesh surface (so the vertex normal is z axis)
		var localz = normals[i];
		vec3.cross(normals[i], yaxis, localx);
		if (vec3.length(localx) < 0.1)
			vec3.cross(normals[i], zaxis, localx);
		vec3.normalize(localx);
		vec3.cross(normals[i], localx, localy);
		basis[0] = localx[0];
		basis[1] = localy[0];
		basis[2] = localz[0];
		basis[3] = localx[1];
		basis[4] = localy[1];
		basis[5] = localz[1];
		basis[6] = localx[2];
		basis[7] = localy[2];
		basis[8] = localz[2];
		
		// project all neighbouring vertices into the local basis
		// construct equation matrix
		var nbi  = nb[i],
			cnt  = nbi.length,
			ans  = vec3.create(),
			eqn  = Array(cnt),
			sol  = Array(cnt);
		for (var j=0; j<cnt; j++) {
			mat3.multiplyVec3(basis, vertices[nbi[j]], ans);
			// [x^2  xy  y^2  x  y  1]
			eqn[j] = [ans[0]*ans[0], ans[0]*ans[1], ans[1]*ans[1], ans[0], ans[1], 1];
			sol[j] = [ans[2]];
		}
		
		// find the parameters to the solution of the equation
		var eqnT = numeric.transpose(eqn);
		sol = numeric.dotMMsmall(eqnT, sol);
		eqn = numeric.inv(numeric.dotMMsmall(eqnT, eqn));
		var p = numeric.dotMMsmall(eqn, sol);
		
		var denom = (1 + p[3]*p[3] + p[4]*p[4]);
		K[i] = (4 * p[0] * p[2] - p[1]*p[1]) / (denom*denom);
	}
	
	return K;
}

