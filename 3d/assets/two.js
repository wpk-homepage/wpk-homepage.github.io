// Two: numerical operations on matrices, a javascript library
//
// @author Paul Koppen
// (c) Copyright Paul Koppen, 2012
//
// fix: normalise + write tests
//
// conceptualise matrix & vector operations similar to element-wise operations?
//
var Two = {
/**
 * Object with ".data" Float32Array and ".size" [h,w]
 */
Arr: function (size, data, type) {
	if (size.length == 1)
		size.push(1);
	else if (size.length != 2)
		throw new this.InvalidArgumentsError;
	if (!data) // no data? 'data' will be array length
		data = size[0] * size[1];
	if (!type)
		type = Float32Array;
	var typeddata,
		columns = [];
	if (type === Array && typeof(data) == 'object') {
		typeddata = data;
		typeddata.subarray = typeddata.slice;
		columns = null;	// unsupported feature for JS Array object
	} else {
		typeddata = new type(data);
		for (var i=0,cnt=typeddata.length; i<cnt; i+=size[0])
			columns.push(typeddata.subarray(i, i+size[0]));
	}
	return {
		columns:  columns,
		data:     typeddata,
		size:     size,
		toString: function(){return Two.toString(this);},
		type:     type,
	};
},
isArr: function(arr) {
	return arr.data && arr.size;
},
toString: function (arr) {
	var adata = arr.data,
		jsarr = [];
	//for (var i=0,cnt=adata.length; i<cnt; i++)
	//	jsarr.push(adata[i]);
	jsarr.push.apply(jsarr, adata);
	return 'Arr('+arr.size.join('x') + ') ['+jsarr.join(',')+']';
},
/**
 * Element-wise addition
 *
 * @param {Arr} a  First matrix / vector.
 * @param {Arr} b  Second matrix / vector. In each dimension, the size of "a"
 *                 and "b" must match or either must be 1, in which case the
 *                 data is repeated (singleton dimension expansion) to fit the
 *                 other.
 * @param {Arr} res  Output result by reference.
 * @return {Arr} A copy of "a" incremented by "b". Matrix / vector, a + b.
 */
add: function (arr, b, res) {
	return this.__element(arr, b, res, this.add_);
},
add_: function (adata, bdata, rdata, asize0, asize1, bsize0, bsize1, rsize0, rsize1) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length;
	
	// same size
	if (asize0 == bsize0 && asize1 == bsize1) {
		for (var i=rnumel; i--; )
			rdata[i] = adata[i] + bdata[i];
	}
	// at least one scalar
	else if (anumel == 1 || bnumel == 1) {
		var matrix = anumel == 1 ? bdata : adata,
			scalar = anumel == 1 ? adata[0] : bdata[0];
		for (var i=rnumel; i--; )
			rdata[i] = matrix[i] + scalar;
	}
	// matrix and column vector
	else if (asize0 != 1 && bsize0 != 1) {
		var matrix = asize1 == 1 ? bdata : adata,
			vector = asize1 == 1 ? adata : bdata,
			i      = rnumel;
		for (var col=rsize1; col--; )
			for (var row=rsize0; row--; ) {
				i--;
				rdata[i] = matrix[i] + vector[row];
			}
	}
	// matrix and row vector
	else if (asize1 != 1 && bsize1 != 1) {
		var matrix = asize0 == 1 ? bdata : adata,
			vector = asize0 == 1 ? adata : bdata,
			i      = rnumel,
			vcol;
		for (var col=rsize1; col--; ) {
			vcol = vector[col];
			for (var row=rsize0; row--; ) {
				i--;
				rdata[i] = matrix[i] + vcol;
			}
		}
	}
	// column vector and row vector
	else {
		var cvector = asize0 == 1 ? bdata : adata,
			rvector = asize0 == 1 ? adata : bdata,
			i       = rnumel;
		for (var col=rsize1; col--; )
			for (var row=rsize0; row--; )
				rdata[--i] = cvector[row] + rvector[col];
	}
},
/**
 * Vector cross product
 *
 * @param {Arr}    a    First matrix / vector. Size in "dim" must be 3. In the
 *                      other dimension "arr" and "b" must match size or be of
 *                      size 1, in which case it is replicated (singleton
 *                      dimension expansion) to match the other.
 * @param {Arr}    b    Second matrix / vector.
 * @param {number} dim  Optional dimension. Pass 1 for row vectors, default 0.
 * @param {Arr}    res  Pass Arr instance to get results by reference.
 * @return {Arr} Cross product of each vector in "a" with the corresponding
 *               vector in "b", a x b.
 */
cross: function (arr, b, dim, res) {
	var asize = arr.size,
		bsize = b.size,
		rsize = [];
	
	// find the right dimension to operate along
	if (typeof dim == 'number') {
		if (dim < 0 || dim > 1)
			throw new this.InvalidDimensionError;
		else if (asize[dim] != 3 || bsize[dim] != 3)
			throw new this.InvalidDimensionError;
	}
	else {
		if (asize[0] == 3 && bsize[0] == 3)
			dim = 0;
		else if (asize[1] == 3 && bsize[1] == 3)
			dim = 1;
		else
			throw new this.DimensionMismatchError;
	}
	
	if (asize[1-dim] != bsize[1-dim] && asize[1-dim] != 1 && bsize[1-dim] != 1)
		throw new this.DimensionMismatchError;
	
	rsize[dim]   = 3;
	rsize[1-dim] = Math.max(asize[1-dim], bsize[1-dim]);
	
	if (!res)
		res = this.Arr(rsize, rsize[0]*rsize[1]);
	else if (res.data.length == rsize[0]*rsize[1])
		res.size = rsize;
	else
		throw new this.DimensionMismatchError;
	
	this.cross_(arr.data, b.data, res.data, asize[0], bsize[0], rsize[0], dim);
	return res;
},
cross_: function(adata, bdata, rdata, asize0, bsize0, rsize0, dim) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length,
		ai     = 0,
		bi     = 0,
		ainc   = anumel > 3 ? 1 : 0,
		binc   = bnumel > 3 ? 1 : 0;
	
	if (dim == 0) {
		ainc *= 3;
		binc *= 3;
		for (var ri=0; ri<rnumel; ri+=3) {
			acol        = adata.subarray(ai, ai+3);
			bcol        = bdata.subarray(bi, bi+3);
			rdata[ri+0] = acol[1]*bcol[2] - acol[2]*bcol[1];
			rdata[ri+1] = acol[2]*bcol[0] - acol[0]*bcol[2];
			rdata[ri+2] = acol[0]*bcol[1] - acol[1]*bcol[0];
			ai         += ainc;
			bi         += binc;
		}
	}
	else {
		var a1 = adata.subarray(0, asize0),
			a2 = adata.subarray(asize0, 2*asize0),
			a3 = adata.subarray(2*asize0),
			b1 = bdata.subarray(0, bsize0),
			b2 = bdata.subarray(bsize0, 2*bsize0),
			b3 = bdata.subarray(2*bsize0),
			r1 = rdata.subarray(0, rsize0),
			r2 = rdata.subarray(rsize0, 2*rsize0),
			r3 = rdata.subarray(2*rsize0);
		for (var ri=0; ri<rsize0; ri++) {
			r1[ri] = a2[ai]*b3[bi] - a3[ai]*b2[bi];
			r2[ri] = a3[ai]*b1[bi] - a1[ai]*b3[bi];
			r3[ri] = a1[ai]*b2[bi] - a2[ai]*b1[bi];
			ai    += ainc;
			bi    += binc;
		}
	}
},
/**
 * Element-wise division
 *
 * @param {Arr}          a    Numerator.
 * @param {Arr | number} b    Denominator.
 * @param {Arr}          res  Result matrix / array by reference.
 * @return {Arr} A copy of "a" with each element divided by "b", a ./ b.
 */
div: function (arr, b, res) {
	return this.__element(arr, b, res, this.div_);
},
div_: function(adata, bdata, rdata, asize0, asize1, bsize0, bsize1, rsize0, rsize1) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length;
	
	// same size
	if (asize0 == bsize0 && asize1 == bsize1) {
		for (var i=rnumel; i--; )
			rdata[i] = adata[i] / bdata[i];
	}
	// scalar
	else if (anumel == 1 || bnumel == 1) {
		var matrix = anumel == 1 ? bdata : adata,
			scalar = anumel == 1 ? adata[0] : bdata[0];
		for (var i=rnumel; i--; )
			rdata[i] = matrix[i] / scalar;
	}
	// matrix and column vector
	else if (asize0 != 1 && bsize0 != 1) {
		var matrix = asize1 == 1 ? bdata : adata,
			vector = asize1 == 1 ? adata : bdata,
			i      = rnumel;
		for (var col=rsize1; col--; )
			for (var row=rsize0; row--; ) {
				i--;
				rdata[i] = matrix[i] / vector[row];
			}
	}
	// matrix and row vector
	else if (asize1 != 1 && bsize1 != 1) {
		var matrix = asize0 == 1 ? bdata : adata,
			vector = asize0 == 1 ? adata : bdata,
			i      = rnumel,
			vcol;
		for (var col=rsize1; col--; ) {
			vcol = vector[col];
			for (var row=rsize0; row--; ) {
				i--;
				rdata[i] = matrix[i] / vcol;
			}
		}
	}
	// column vector and row vector
	else {
		var cvector = asize0 == 1 ? bdata : adata,
			rvector = asize0 == 1 ? adata : bdata,
			i       = rnumel;
		for (var col=rsize1; col--; )
			for (var row=rsize0; row--; )
				rdata[--i] = cvector[row] / rvector[col];
	}
},
/**
 * Vector dot product
 *
 * @param {Arr}    a    First matrix / vector.
 * @param {Arr}    b    Second matrix / vector. "a" and "b" must match in size
 *                      on dimension "dim". In the other dimension either may be
 *                      of size 1 in which case it is replicated (singleton
 *                      dimension expansion) to match the other.
 * @param {number} dim  Optional dimension. Pass 1 for row vectors, default 0.
 * @param {Arr}    res  Pass Arr instance to get results by reference.
 * @return {Arr} Dot product of each vector in "a" with the corresponding vector
 *               in "b", a . b.
 */
dot: function(arr, b, dim, res) {
	var asize  = arr.size,
		bsize  = b.size,
		rsize  = [],
		rnumel = 0;
	
	// find the right dimension to operate along
	if (typeof dim == 'number') {
		if (dim < 0 || dim > 1)
			throw new this.InvalidDimensionError;
		else if (asize[dim] != bsize[dim])
			throw new this.InvalidDimensionError;
	}
	else {
		if (asize[0] == bsize[0])
			dim = 0;
		else if (asize[1] == bsize[1])
			dim = 1;
		else
			throw new this.DimensionMismatchError;
	}
	
	if (asize[1-dim] != bsize[1-dim] && asize[1-dim] != 1 && bsize[1-dim] != 1)
		throw new this.DimensionMismatchError;
	
	rsize[dim]   = 1;
	rsize[1-dim] = Math.max(asize[1-dim], bsize[1-dim]);
	rnumel       = rsize[0] * rsize[1];
	
	if (!res)
		res = this.Arr(rsize, rnumel);
	else if (res.data.length == rnumel)
		res.size = rsize;
	else
		throw new this.DimensionMismatchError;
	
	this.dot_(arr.data, b.data, res.data, asize[0], asize[1], bsize[0], dim);
	return res;
},
dot_: function(adata, bdata, rdata, asize0, asize1, bsize0, dim) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length,
		ai, bi, ri, vdot;

	if (dim == 0) {
		for (var ri=rnumel; ri--; ) {
			ai   = ai || anumel;
			bi   = bi || bnumel;
			vdot = 0;
			for (var j=asize0; j--; ) {
				vdot += adata[--ai] * bdata[--bi];
			}
			rdata[ri] = vdot;
		}
	}
	else {
		var ainc = asize0 == 1 ? anumel : (anumel - 1),
			binc = bsize0 == 1 ? bnumel : (bnumel - 1);
		ai       = anumel - 1;
		bi       = bnumel - 1;
		for (var ri=rnumel; ri--; ) {
			vdot = 0;
			for (var j=asize1; j--; ) {
				vdot += adata[ai] * bdata[bi];
				ai   -= asize0;
				bi   -= bsize0;
			}
			rdata[ri] = vdot;
			ai       += ainc;
			bi       += binc;
		}
	}
},
/**
 * Identity matrix
 *
 * @param {number} n  Size of the returned matrix along both dimensions.
 * @return {Arr} The n x n identity matrix (square matrix of ones along the
 *               diagonal and zeros elsewhere).
 */
eye: function (n) {
	var arr   = this.Arr([n,n]),
		adata = arr.data,
		step  = n + 1;
	for (var i=0,cnt=adata.length; i<cnt; i+=step)
		adata[i] = 1;
	return arr;
},
/**
 * Overwrite all values with a single new value
 *
 * @param {Arr}    a      Matrix / vector.
 * @param {number} value  Value written to each element of "a".
 * @return {Arr} The original "a" is returned for testing purposes.
 *
 * @note This function writes to the input Arr directly. If you desire a copy
 *       filled with "value", first run "var copy = Two.Arr(a.size);" and then
 *       fill the copy.
 */
fill: function (arr, value) {
	var adata = arr.data;
	for (var i=adata.length; i--; )
		adata[i] = value;
	return arr;
},
/**
 * Indices into the matrix running over the rows first
 */
indexrowsfirst: function (arr, type) {
	var adata  = arr.data,
		anumel = adata.length,
		nrows  = arr.size[0],
		ri     = 0,
		rdata;
	
	if (!type)
		type = anumel < 256 ? Uint8Array : (
				anumel < 65536 ? Uint16Array : Uint32Array);
	rdata = new type(anumel);
	
	for (var row=0; row<nrows; row++) {
		for (var ai=row; ai<anumel; ai+=nrows)
			rdata[ri++] = ai;
	}
	
	return rdata;
},
/**
 * Vector mean
 */
mean: function (arr, dim, res) {
	var asize  = arr.size,
		rsize  = [],
		rnumel = 0;
	
	// find the right dimension to operate along
	if (typeof dim == 'number') {
		if (dim < 0 || dim > 1)
			throw new this.InvalidDimensionError;
	}
	else {
		if (asize[0] == 1)
			dim = 1;
		else
			dim = 0;
	}
	
	rsize[dim]   = 1;
	rsize[1-dim] = asize[1-dim];
	rnumel       = asize[1-dim];
	
	if (!res)
		res = this.Arr(rsize, rnumel);
	else if (res.data.length == rnumel)
		res.size = rsize;
	else
		throw new this.DimensionMismatchError;
	
	this.mean_(arr.data, res.data, asize[0], asize[1], dim);
	return res;
},
mean_: function(adata, rdata, asize0, asize1, dim) {
	var anumel = adata.length,
		rnumel = rdata.length,
		vec;
	
	if (dim == 0) {
		var ri = 0;
		for (var col=0; col<anumel; col+=asize0) {
			vec = adata.subarray(col, col+asize0);
			for (var row=0; row<asize0; row++)
				rdata[ri] += vec[row];
			rdata[ri] /= asize0;
			ri++;
		}
	}
	else {
		for (var col=0; col<anumel; col+=asize0) {
			vec = adata.subarray(col, col+asize0);
			for (var row=0; row<asize0; row++)
				rdata[row] += vec[row];
		}
		for (var row=asize0; row--; )
			rdata[row] /= asize1;
	}
},
/**
 * Matrix multiplication
 *
 * @param {Arr}  a    NxM matrix.
 * @param {Arr}  b    MxP matrix or 1x1.
 * @param {Arr}  res  Pass Arr instance to get results by reference.
 * @return {Arr} "a" multiplied by "b". NxP matrix, a * b.
 *               If "b" is a scalar then a NxM matrix is returned.
 */
mul: function(arr, b, res) {
	var adata   = arr.data,
		anumel  = adata.length,
		bdata   = b.data,
		bnumel  = bdata.length,
		szleft  = arr.size[0],
		szinner = arr.size[1],
		szright = b.size[1],
		rnumel  = szleft * szright,
		rdata;
	
	// check dimensions
	if (szinner != b.size[0] && anumel != 1 && bnumel != 1)
		throw new this.DimensionMismatchError;
	
	// initialise 'res'
	if (!res)
		res = this.Arr([szleft,szright], rnumel);
	else if (res.data.length == rnumel)
		res.size = [szleft, szright];
	else
		throw new this.DimensionMismatchError;
	rdata  = res.data;
	
	// some trivial cases
	// zero data
	if (rnumel == 0)
		return res;
	else if (szinner == 0) {
		for (var i=rnumel; i--; )
			rdata[i] = 0;	// copy Matlab's behaviour
		return res;
	}
	
	// perform matrix multiplication
	this.mul_(adata, bdata, rdata, szleft, szinner, szright);
	return res;
},
mul_: function(adata, bdata, rdata, szleft, szinner, szright) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length,
		ri, bcol, vdot;
	
	// multiply matrix by matrix / column vector
	if (szleft > 1 && szinner > 1) {
		ri = 0;
		for (var col=0; col<bnumel; col+=szinner) {
			bcol = bdata.subarray(col, col+szinner);
			for (var row=0; row<szleft; row++) {
				vdot = 0;
				for (var ai=row,bi=0; ai<anumel; ai+=szleft,bi++)
					vdot += bcol[bi] * adata[ai];
				rdata[ri++] = vdot;
			}
		}
	}
	// multiply by scalar
	else if (anumel == 1 || bnumel == 1) {
		this.times_(anumel==1?bdata:adata, anumel==1?adata:bdata, rdata);
	}
	// multiply row vector by matrix / column vector
	else if (szleft == 1) {
		ri = 0;
		for (var bi=0; bi<bnumel; bi+=szinner) {
			bcol = bdata.subarray(bi, bi+szinner);
			vdot = 0;
			for (var j=szinner; j--; )
				vdot += bcol[j] * adata[j];
			rdata[ri++] = vdot;
		}
	}
	// multiply column vector by row vector (outer product)
	else {
		ri = rnumel;
		for (var bi=bnumel; bi--; )
			for (var ai=anumel; ai--; )
				rdata[--ri] = bdata[bi] * adata[ai];
	}
},
/**
 * Element-wise negation
 *
 * @param {Arr}    a    Matrix / vector /array.
 * @param {Arr}    res  Pass Arr instance to get results by reference.
 * @return {Arr} Copy of "a" with each element "a[i]" replaced by "-a[i]".
 */
neg: function (arr, res) {
	return this.__element(arr, 0, res, this.neg_);
},
neg_: function (adata, bdata, rdata, asize0, asize1, bsize0, bsize1, rsize0, rsize1) {
	for (var i=rdata.length; i--; )
		rdata[i] = -adata[i];
},
/**
 * Matrix / vector norms
 *
 * @param {Arr}    a  Matrix / vector / array.
 * @param {number} p  Type of norm. Choices are:
 *                    1:        1-norm (sum of absolute values)
 *                    2:        2-norm / Frobenius / Euclidean norm (square root
 *                              of sum of squares)
 *                    Infinity: Infinity / spectral norm (maximum of absolute
 *                              values)
 *                    p:        p-norm = sum_i(abs(a_i)^p) ^ (1/p)
 * @return {number} Norm of the input argument, |a|_p.
 */
norm: function(arr, p) {
	var adata = arr.data,
		asize = arr.size,
		n = 0;
	
	// trivial case
	if (adata.length == 0 || asize[0] == 0 || asize[1] == 0)
		return 0;
	
	// default is 2-norm
	if (!p)
		p = 2;
	
	if (p == 1)					// l1-norm
		return this.norm_1(adata);
	else if (p == 2)			// l2-norm / Frobenius norm / Euclidean norm
		return this.norm_2(adata);
	else if (p == Infinity)		// spectral norm
		return this.norm_inf(adata);
	else						// p-norm
		return this.norm_p(adata, p);
},
norm_1: function (adata) {
	var n = 0;
	for (var i=adata.length; i--; )
		n += Math.abs(adata[i]);
	return n;
},
norm_2: function (adata) {
	var n = 0;
	for (var i=adata.length; i--; )
		n += adata[i] * adata[i];
	return Math.sqrt(n);
},
norm_inf: function (adata) {
	var n = 0;
	for (var i=adata.length; i--; )
		n = Math.max(n, Math.abs(adata[i]));
	return n;
},
norm_p: function (adata, p) {
	var n = 0;
	for (var i=adata.length; i--; )
		n += Math.pow(Math.abs(adata[i]), p);
	return Math.pow(n, 1/p);
},
/**
 * Normalise vector lengths
 */
normalise: function(arr, dim, res) {
	var adata  = arr.data,
		anumel = adata.length,
		asize0 = arr.size[0],
		asize1 = arr.size[1],
		acol   = null,
		vnorm  = 0;
	
	if (dim == 1) {	// row vectors
	}
	else {			// column vectors
		for (var col=0; col<anumel; col+=asize0) {
			acol    = adata.subarray(col, col+asize0);
			invnorm = 1 / this.norm_2(acol);
			for (var row=0; row<asize0; row++)
				acol[row] *= invnorm;
		}
	}
},
/**
 * Create a 3x3 rotation matrix from specified axis and angle
 *
 * @param {Arr}    axis   Vector of size 3x1.
 * @param {number} theta  Rotation in radians.
 * @param {Arr}    res    Optional 3x3 matrix taking the result by reference.
 * @return {Arr} 3x3 matrix.
 */
rotationmatrix: function(axis, theta, res) {
	// http://en.wikipedia.org/wiki/Rotation_matrix
	var adata = axis.data || axis,
		cost  = Math.cos(theta),
		cost1 = 1 - cost,
		sint  = Math.sin(theta),
		// I * cos theta + sin theta * cross product matrix
		r     = this.Arr([3,3], [cost,sint*adata[2],sint*-adata[1],
							sint*-adata[2],cost,sint*adata[0],
							sint*adata[1],sint*-adata[0],cost]);
		// (1 - cos theta) * tensor product
		tp    = this.Arr([3,3], [
				cost1*adata[0]*adata[0],cost1*adata[0]*adata[1],cost1*adata[0]*adata[2],
				cost1*adata[1]*adata[0],cost1*adata[1]*adata[1],cost1*adata[1]*adata[2],
				cost1*adata[2]*adata[0],cost1*adata[2]*adata[1],cost1*adata[2]*adata[2]]);
	
	return this.add(r, tp, res);
},
/**
 * Element-wise subtraction
 *
 * @param {Arr} a  First matrix / vector.
 * @param {Arr} b  Second matrix / vector. In each dimension, the size of "a"
 *                 and "b" must match or either must be 1, in which case the
 *                 data is repeated (singleton dimension expansion) to fit the
 *                 other.
 * @param {Arr} res  Output result by reference.
 * @return {Arr} A copy of "a" incremented by "b". Matrix / vector, a - b.
 */
sub: function (arr, b, res) {
	return this.__element(arr, b, res, this.sub_);
},
sub_: function (adata, bdata, rdata, asize0, asize1, bsize0, bsize1, rsize0, rsize1) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length;
	
	// same size
	if (asize0 == bsize0 && asize1 == bsize1) {
		for (var i=rnumel; i--; )
			rdata[i] = adata[i] - bdata[i];
	}
	// at least one scalar
	else if (anumel == 1 || bnumel == 1) {
		if (anumel == 1)
			for (var i=rnumel; i--; )
				rdata[i] = adata[0] - bdata[i];
		else//if (bnumel == 1)
			for (var i=rnumel; i--; )
				rdata[i] = adata[i] - bdata[0];
	}
	// matrix and column vector
	else if (asize0 != 1 && bsize0 != 1) {
		var i = rnumel - 1;
		if (asize1 == 1)
			for (var col=rsize1; col--; )
				for (var row=rsize0; row--; i--)
					rdata[i] = adata[row] - bdata[i];
		else//if (bsize1 == 1)
			for (var col=rsize1; col--; )
				for (var row=rsize0; row--; i--)
					rdata[i] = adata[i] - bdata[row];
	}
	// matrix and row vector
	else if (asize1 != 1 && bsize1 != 1) {
		var i = rnumel - 1,
			vcol;
		if (asize0 == 1)
			for (var col=rsize1; col--; ) {
				vcol = adata[col];
				for (var row=rsize0; row--; i--)
					rdata[i] = vcol - bdata[i];
			}
		else//if (bsize0 == 1)
			for (var col=rsize1; col--; ) {
				vcol = bdata[col];
				for (var row=rsize0; row--; i--)
					rdata[i] = adata[i] - vcol;
			}
	}
	// column vector and row vector
	else {
		var i = rnumel;
		if (asize0 == 1)
			for (var col=rsize1; col--; )
				for (var row=rsize0; row--; )
					rdata[--i] = adata[col] - bdata[row];
		else//if (bsize0 == 1)
			for (var col=rsize1; col--; )
				for (var row=rsize0; row--; )
					rdata[--i] = adata[row] - bdata[col];
	}
},
/**
 * Element-wise multiplication
 *
 * @param {Arr}          a    Matrix / array.
 * @param {Arr | number} b    Multiplier.
 * @param {Arr}          res  Result matrix / array by reference.
 * @return {Arr} A copy of "a" with each element multiplied by "b", a .* b.
 */
times: function (arr, b, res) {
	return this.__element(arr, b, res, this.times_);
},
times_: function(adata, bdata, rdata, asize0, asize1, bsize0, bsize1, rsize0, rsize1) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length;
	
	// same size
	if (asize0 == bsize0 && asize1 == bsize1) {
		for (var i=rnumel; i--; )
			rdata[i] = adata[i] * bdata[i];
	}
	// scalar
	else if (anumel == 1 || bnumel == 1) {
		var matrix = anumel == 1 ? bdata : adata,
			scalar = anumel == 1 ? adata[0] : bdata[0];
		for (var i=rnumel; i--; )
			rdata[i] = matrix[i] * scalar;
	}
	// matrix and column vector
	else if (asize0 != 1 && bsize0 != 1) {
		var matrix = asize1 == 1 ? bdata : adata,
			vector = asize1 == 1 ? adata : bdata,
			i      = rnumel;
		for (var col=rsize1; col--; )
			for (var row=rsize0; row--; ) {
				i--;
				rdata[i] = matrix[i] * vector[row];
			}
	}
	// matrix and row vector
	else if (asize1 != 1 && bsize1 != 1) {
		var matrix = asize0 == 1 ? bdata : adata,
			vector = asize0 == 1 ? adata : bdata,
			i      = rnumel,
			vcol;
		for (var col=rsize1; col--; ) {
			vcol = vector[col];
			for (var row=rsize0; row--; ) {
				i--;
				rdata[i] = matrix[i] * vcol;
			}
		}
	}
	// column vector and row vector
	else {
		var cvector = asize0 == 1 ? bdata : adata,
			rvector = asize0 == 1 ? adata : bdata,
			i       = rnumel;
		for (var col=rsize1; col--; )
			for (var row=rsize0; row--; )
				rdata[--i] = cvector[row] * rvector[col];
	}
},
/**
 * Matrix transpose multiplication
 *
 * @param {Arr}  a    MxN matrix.
 * @param {Arr}  b    MxP matrix or 1x1.
 * @param {Arr}  res  Pass Arr instance to get results by reference.
 * @return {Arr} "a" multiplied by "b". NxP matrix, a^T * b.
 *               If "b" is a scalar then a NxM matrix is returned.
 */
tmul: function(arr, b, res) {
	var adata   = arr.data,
		anumel  = adata.length,
		bdata   = b.data,
		bnumel  = bdata.length,
		szleft  = arr.size[1],
		szinner = arr.size[0],
		szright = b.size[1],
		rnumel  = szleft * szright,
		rdata;
	
	if (szinner != b.size[0] && anumel != 1 && bnumel != 1)
		throw new this.DimensionMismatchError;
	
	if (!res)
		res = this.Arr([szleft,szright], rnumel);
	else if (res.data.length == rnumel)
		res.size = [szleft, szright];
	else
		throw new this.DimensionMismatchError;
	rdata   = res.data;
	
	// some trivial cases
	// zero data
	if (rnumel == 0)
		return res;
	else if (szinner == 0) {
		for (var i=rnumel; i--; )
			rdata[i] = 0;	// copy Matlab's behaviour
		return res;
	}
	
	// perform transpose multiplication
	this.tmul_(adata, bdata, rdata, szleft, szinner, szright);
	return res;
},
tmul_: function(adata, bdata, rdata, szleft, szinner, szright) {
	var anumel = adata.length,
		bnumel = bdata.length,
		rnumel = rdata.length,
		ri, acol, bcol, vdot;
	
	// multiply matrix by matrix / column vector
	if (szleft > 1 && szinner > 1) {
		ri = 0;
		for (var bi=0; bi<bnumel; bi+=szinner) {
			bcol = bdata.subarray(bi, bi+szinner);
			for (var ai=0; ai<anumel; ai+=szinner) {
				acol = adata.subarray(ai, ai+szinner);
				vdot = 0;
				for (var j=szleft; j--; )
					vdot += bcol[j] * acol[j];
				rdata[ri++] = vdot;
			}
		}
	}
	// multiply by scalar
	else if (anumel == 1 || bnumel == 1) {
		this.times_(anumel==1?bdata:adata, anumel==1?adata:bdata, rdata);
	}
	// multiply column vector by matrix / column vector
	else if (szleft == 1) {
		ri = 0;
		for (var bi=0; bi<bnumel; bi+=szinner) {
			bcol = bdata.subarray(bi, bi+szinner);
			vdot = 0;
			for (var j=szinner; j--; )
				vdot += bcol[j] * adata[j];
			rdata[ri++] = vdot;
		}
	}
	// multiply row vector by row vector / scalar (outer product)
	else {
		ri = rnumel;
		for (var bi=bnumel; bi--; )
			for (var ai=anumel; ai--; )
				rdata[--ri] = bdata[bi] * adata[ai];
	}
},
/**
 * Transpose
 *
 * @param {Arr}  a    MxN matrix.
 * @param {Arr}  res  Pass Arr instance to get results by reference.
 * @return {Arr} A transposed copy of "a". NxM matrix, a^T.
 */
transp: function (arr, res) {
	if (!res)
		res = this.Arr(arr.size);
	
	var adata  = arr.data,
		asize  = arr.size,
		rdata  = res.data,
		rnumel = rdata.length;
	
	res.size = [asize[1], asize[0]];
	
	// trivial case
	if (asize[0] == 0 || asize[1] == 0)
		return res;
	// transpose vector / scalar
	else if (asize[0] == 1 || asize[1] == 1) {
		for (var i=rnumel; i--; )
			rdata[i] = adata[i];
	}
	// full matrix transposition
	else {
		var nrows = asize[0],
			ri    = 0;
		for (var row=0; row<nrows; row++) {
			for (var ai=row; ai<rnumel; ai+=nrows)
				rdata[ri++] = adata[ai];
		}
	}
	
	return res;
},


/*******************************************************************************
 *    Wrappers for classes of functions
 ******************************************************************************/

/**
 * Element-wise operations like addition and subtraction
 */
__element: function (arr, b, res, fn) {
	if (typeof arr == 'number')
		arr = this.Arr([1,1], [arr]);
	if (typeof b == 'number')
		b = this.Arr([1,1], [b]);
	
	var asize  = arr.size,
		bsize  = b.size,
		rsize  = asize.slice(),
		rnumel = 1;
	
	// can only expand singleton dimensions
	if (asize[0] != bsize[0] && asize[0] != 1 && bsize[0] != 1)
		throw new this.DimensionMismatchError;
	else if (asize[1] != bsize[1] && asize[1] != 1 && bsize[1] != 1)
		throw new this.DimensionMismatchError;
	
	// expand to dimensions of b where applicable
	if (rsize[0] == 1) rsize[0] = bsize[0];
	if (rsize[1] == 1) rsize[1] = bsize[1];
	rnumel = rsize[0] * rsize[1];
	
	// initialise 'res'
	if (!res)
		res = this.Arr(rsize, rnumel);
	else if (res.data.length == rnumel)
		res.size = rsize;
	else
		throw new this.DimensionMismatchError;
	
	// trivial case
	if (rnumel == 0)
		return res;
	
	fn(arr.data, b.data, res.data, asize[0], asize[1], bsize[0], bsize[1], rsize[0], rsize[1]);
	return res;
},


/*******************************************************************************
 *    Errors
 ******************************************************************************/

DimensionMismatchError: function(msg) {
	this.name     = 'DimensionMismatchError';
	this.message  = msg || 'Array dimensions must match.';
},
InvalidDimensionError: function(msg) {
	this.name     = 'InvalidDimensionError';
	this.message  = msg || 'Cannot operate along that dimension.';
},
InvalidArgumentsError: function() {
	this.name     = 'InvalidArgumentsError';
	this.message  = 'Each Two.Arr must have exactly two dimensions.';
},
NotImplementedError: function() {
	this.name     = 'NotImplementedError';
	this.message  = 'That feature has not yet been implemented.';
},
};

Two.DimensionMismatchError.prototype = new Error();
Two.InvalidDimensionError.prototype  = new Error();
Two.InvalidArgumentsError.prototype  = new Error();
Two.NotImplementedError.prototype    = new Error();


/*******************************************************************************
 *    Convenience
 ******************************************************************************/

Two.ONE = Two.Arr([1,1], [1]);


/*******************************************************************************
 *    Aliases
 ******************************************************************************/

Two.transpose = Two.T      = Two.transp;
Two.scale                  = Two.times;
Two.multiply  = Two.mtimes = Two.mul;
Two.negate                 = Two.neg;
Two.normalize              = Two.normalise;
Two.subtract               = Two.sub;
Two.identity               = Two.eye;


